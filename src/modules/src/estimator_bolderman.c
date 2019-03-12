/*
CRAZYFLIE-FIRMWARE
this model implements the UKF with position measurements to perform a state estimation
for the quadcopter.
The model implemented is the following:
- INPUT
    - Thurst
    - Gyroscopic measurements
- OUTPUT
    - Position
    - Acceleration (in body coordinates + gravity)
*/

// Include all libraries / files used here
#include "estimator_bolderman.h"    // header file

#include "stm32f4xx.h"

#include "FreeRTOS.h"
#include "queue.h"          // Queueing
#include "task.h"           // Input, probably?
#include "sensors.h"        // Using sensors

#include "log.h"            // Using the logblocks
#include "param.h"          // Using the parameters

#include "math.h"           // Using the math library
#include "cf_math.h"        // Crazyflie math

// Define some parameters / variables useful later
#define DEG_TO_RAD (PI/180.0f)            // Degrees to radians
#define RAD_TO_DEG (180.0f/PI)            // Radians to degrees
#define GRAVITY_MAGNITUDE (9.81f)         // Gravitational constant
#define PREDICT_RATE RATE_100_HZ          // Frequency of which this filter is applied
#define CRAZYFLIE_WEIGHT_grams (27.0f)    // Weight of crazyflie
#define CONTROL_TO_ACC (GRAVITY_MAGNITUDE*60.0f/CRAZYFLIE_WEIGHT_grams/65536.0f)
                // variable defining the factor from control input to acceleration
#define N 12                              // Dimension of system
#define NOUT 6                            // Dimension of measurement
#define NSIGMA (2*N+1)                    // Number of sigma points
#define NCALC 12.0f                       // Dimension as a float
// Output = [ax, ay, az] --> acceleration in body coordinates
// Output = [x, y, z, ax, ay, az] --> position in global coordinates and acceleration in body coordinates
#define UPPERBOUND 10000.0f                  // Upperbound after which inversion is adjusted
#define UPPERBOUND_TRACE 100.0f
#define TS (1.0f/100.0f)
#define UPPERBOUND_DT 0.1f                // Upperbound of th timestep

// QUEUEING for position measurement
#define TWR_MEASUREMENT_QUEUE_LENGTH (10) // Number of measurements maximally in the queue
static xQueueHandle twrDataQueueu;

// Prototype of function to add measurement to queue
bool estimatorBoldermanEnqueueDistance(distanceMeasurement_t *measurement);
// Prototype of Boolean to receive measurement from queue
static inline bool estimatorHasTWRMeasurement(distanceMeasurement_t *dist) {
  return (pdTRUE == xQueueReceive(twrDataQueue, dist, 0));
}

/*******************************************************************************
*   Internal variables for the estimator

*******************************************************************************
*   Variables for debugging and logging
*******************************************************************************/
static bool isInit = false;                 // Boolean variable, whether system is initialized
static int32_t lastPrediction;              // Integer value, containing time sample when last prediction was done

// Accumulator variables
static Axis3f accAccumulator;               // Acceleration, summed over n samples (for average taking purposes)
static Axis3f gyroAccumulator;              // Gyroscopic, summed over n sample (for average taking purposes)
static Axis3f magAccumulator;               // Angular acceleration, summed over n samples (for average taking purposes)
static float thrustAccumulator;             // Thrust / input, summed over n samples (for average taking pursposes)
static uint32_t accAccumulatorCount;        // Count variable, number of samples of which acceleration is summed
static uint32_t gyroAccumulatorCount;       // Count variable, number of samples of which gyroscopic measurement is summed
static uint32_t magAccumulatorCount;        // Count variable, number of samples of which angular acceleration measurement is summed
static uint32_t thrustAccumulatorCount;     // Count variable, number of samples of which thrust is summed
// Extra variables used in the Unscented Kalman Filter
static float alpha = 0.1f;                // Defines spread of sigma points
static float kappa = 0.0f;                  // Secondary scaling factor
static float beta = 2.0f;                   // Distributation parameter, gaussion := 2
static float lambda;                        // Scaling parameter used for determining: Weights and Sigma points
static float wm[NSIGMA];                    // Weights for calculating the mean (initialized in estimaterBoldermanInit())
static float wc[NSIGMA];                    // Weights for calculating the covariance (initialized in estimaterBoldermanInit())
static float q[N][N];                       // Covariance of process noise (zero-mean Gaussian distributed)
static float r[NOUT][NOUT];                 // Covariance of measurement noise (zero-mean Gaussian distributed)
// State momentary (x, y, z, xdot, ydot, zdot, xddot, yddot, zddot, pitch, yaw, roll)
static float x[N];                          // State column, containing values of the states
static float sigmaX[N][NSIGMA];             // Matrix containing all sigmapoints
static float sigmaXplus[N][NSIGMA];         // Matrix containing updated sigmapoints
static float sigmaYplus[NOUT][NSIGMA];      // Matrix containing output of sigma points
static float xpred[N];                      // Predicted state vector
static float ypred[NOUT];                   // Predicted output of outputted sigmapoints
static float Pxx[N][N];                     // Covariance of state Estimation
static float Pxy[N][NOUT];                  // Covariance of state - outputs
static float Pyy[NOUT][NOUT];               // Covariance of ouput
static float Pyyinv[NOUT][NOUT];            // Inverse of Pyy
static float k[N][NOUT];                    // Kalman gain (Pxy * inv(Pyy))
static uint8_t fly = 0;                    // Boolean stating that the crazyflie flying or not
static float y[NOUT];                       // Measurement vectors (just to show how to use static variables)
static uint8_t ycounter = 0;                // Counter, to see how much postionmeasurements are available
static float anchor[3][NOUT];               // Anchor positions from where distance to drone is measured
static float omega[3];                      // Gyroscopic measurements
static float tracePxx;                      // Trace of the covariance matrix Pxx
// The following arrays help with initializing Pxx, and defining q and r
static float Pxxdiag[N]  = {1.0f,1.0f,1.0f, 1.0f,1.0f,1.0f, 1.0f,1.0f,1.0f, 1.0f,1.0f,1.0f};
static float qdiag[N]    = {TS*0.001f,TS*0.001f,TS*0.001f, TS*0.01f,TS*0.01f,TS*0.01f, TS*0.001f,TS*0.001f,TS*0.001f, TS*0.01f,TS*0.01f,TS*0.01f};
static float rdiag[NOUT] = {0.001f,0.001f,0.001f, 0.001f,0.001f,0.001f};

// Functions used for multiplication
static inline void mat_trans(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_trans_f32(pSrc, pDst)); }
static inline void mat_inv(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_inverse_f32(pSrc, pDst)); }
static inline void mat_mult(const arm_matrix_instance_f32 * pSrcA, const arm_matrix_instance_f32 * pSrcB, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_mult_f32(pSrcA, pSrcB, pDst)); }
static inline float arm_sqrt(float32_t in)
{ float pOut = 0; arm_status result = arm_sqrt_f32(in, &pOut); configASSERT(ARM_MATH_SUCCESS == result); return pOut; }



// Prototype of update function, which is used in estimatorBolderman function
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt, uint32_t osTick);
// Prototype of function combining prediction and measurement = actual update
static void estimatorBoldermanDynMeas(void);
// Prototype of function setting prediction to estimation (when no position measurements available / or not enough)
static void estimatorBoldermanPredEst(void);
// Prototype of function saving the updated state
static void estimatorBoldermanStateSave(state_t *state, uint32_t osTick);
// Prototype of function defining the dynamics
static void estimatorBoldermanPredict(float dt, float thrust);
// Prototype of function for resetting Pxx (when cholesky cannot be computed)
static void resetPxx(void);

// Prototypes of function By Marcus Greiff (cholesky decomposition)
int  assert_element(float val);
int  cholesky_decomposition(float (*A)[N], float (*R)[N], int n);



/*
HERE THE MODEL STARTS
---------------------------------------------------------------
*/
// Function receiving sensordata and control input
void estimatorBolderman(state_t *state, sensorData_t *sensors, control_t *control, const uint32_t tick)
{

  // Asks for tick at moment of start calling this function (estimatorBolderman() )
  uint32_t osTick = xTaskGetTickCount(); // would be nice if this had a precision higher than 1ms...

  // Average IMU measurements (in case the observer is run at a lower
  // rate than the controller). Convert accelerometer measurements from Gs to
  // ms^-2 and gyroscopic measurements from deg/sec to rad/sec
  if (sensorsReadAcc(&sensors->acc)) {
    accAccumulator.x += GRAVITY_MAGNITUDE*sensors->acc.x; // Adds acceleration to corresponding direction
    accAccumulator.y += GRAVITY_MAGNITUDE*sensors->acc.y;
    accAccumulator.z += GRAVITY_MAGNITUDE*sensors->acc.z;
    accAccumulatorCount++;                                // Increments counter variable
  }

  if (sensorsReadGyro(&sensors->gyro)) {
    gyroAccumulator.x += sensors->gyro.x * DEG_TO_RAD;    // Adds gyroscopic measurement to corresponding direction
    gyroAccumulator.y += sensors->gyro.y * DEG_TO_RAD;
    gyroAccumulator.z += sensors->gyro.z * DEG_TO_RAD;
    gyroAccumulatorCount++;                               // Increments counter variable
  }

  if (sensorsReadMag(&sensors->mag)) {
      magAccumulator.x += sensors->mag.x;                 // Adds angular acceleration measurement to corresponding direction
      magAccumulator.y += sensors->mag.y;
      magAccumulator.z += sensors->mag.z;
      magAccumulatorCount++;                              // Increments counter variable
  }

  // Average the thrust command from the last timestep, generated externally by the controller
  thrustAccumulator += control->thrust * CONTROL_TO_ACC; // thrust is in grams, we need ms^-2
  thrustAccumulatorCount++;

  // Run the system dynamics to predict the state forward.
  // This is specified for a certain frequency, lower than the tick frequency !!
  // Furthermore it can only be performed when there were sensor and thrust inputs
  if ((osTick-lastPrediction) >= configTICK_RATE_HZ/PREDICT_RATE // update at the PREDICT_RATE
      && gyroAccumulatorCount > 0
      && accAccumulatorCount > 0
      && thrustAccumulatorCount > 0
      && thrustAccumulatorCount > 0)                      // Conditions, note that there must be input
  {
    // Compute the average gyroscopic velocity
    gyroAccumulator.x /= gyroAccumulatorCount;
    gyroAccumulator.y /= gyroAccumulatorCount;
    gyroAccumulator.z /= gyroAccumulatorCount;

    // Compute average acceleration
    accAccumulator.x /= accAccumulatorCount;
    accAccumulator.y /= accAccumulatorCount;
    accAccumulator.z /= accAccumulatorCount;

    // Compute average angular acceleration
    magAccumulator.x /= magAccumulatorCount;
    magAccumulator.y /= magAccumulatorCount;
    magAccumulator.z /= magAccumulatorCount;

    // Compute average thrust
    thrustAccumulator /= thrustAccumulatorCount;

    // Computes the time step used
    float dt = (float)(osTick-lastPrediction)/configTICK_RATE_HZ;
    if (dt > UPPERBOUND_DT) {
      consolePrintf("dt exceeds upperbound \n");
    }


    // This is where we perform the time-integration of the state estimate
    estimatorBoldermanUpdate(state, thrustAccumulator, &accAccumulator, &gyroAccumulator, &magAccumulator, dt, osTick);

    // Reset the accumulators and counters to be used next sample
    lastPrediction = osTick;                  // Last moment in time we made a prediction == now
    accAccumulator = (Axis3f){.axis={0}};     // Reset, so accumulator can start over
    accAccumulatorCount = 0;                  // Set to zero, so accumalation can start over
    gyroAccumulator = (Axis3f){.axis={0}};
    gyroAccumulatorCount = 0;
    magAccumulator = (Axis3f){.axis={0}};
    magAccumulatorCount = 0;
    thrustAccumulator = 0;
    thrustAccumulatorCount = 0;
  }
}

// Here the integration and the state estimation is performed
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt, uint32_t osTick)
{
  // Store the measuremets y1 y2 y3 as rows in a matrix
  y[3] = acc->x;
  y[4] = acc->y;
  y[5] = acc->z;
  omega[0] = gyro->x;   // ACTUALLY MEASUREMENT, USED AS INPUT!
  omega[1] = gyro->y;
  omega[2] = gyro->z;

  // Get the distance to quadcopter and anchor positions
  // while loop -> multiple position measurements can be obtained during one timestep
  // When there are more than three distance measurements, we only take the three most recent
  distanceMeasurement_t dist;
  while (estimatorHasTWRMeasurement(&dist))
  {
    static uint8_t yindex;
    yindex = ycounter % 3;
    y[yindex] = dist.distance;
    anchor[0][yindex] = dist.x;
    anchor[1][yindex] = dist.y;
    anchor[2][yindex] = dist.z;
    ycounter++;
    //consolePrintf("Queued measurementÖ %f", (double)distance);
  }


  /* The following things are performed here:
      - a prediction of the state at the current moment is made, using the previous state (note that you need dt)
          * compute sigma points
          * use sigma points in dynamical and output equation
          * compute prediction using the weighted averages of sigma points
          * use prediction of sigmapoints to compute Pxy and Pyy
      - Update prediction
          * determine kalman gain
          * update estimation x, using measurement and predicted value
          * update the covariance matrix Pxx
      - estimation is saved in the state variable, so it can be used by controller
  */
  // Use dynamics to predict state at next interval
  estimatorBoldermanPredict(dt, thrust);
  // The update can only be performed with at least three position measurements
  if (ycounter >= 3) {
    // Use predicted state and measurement to improve Estimation
    estimatorBoldermanDynMeas();
  } else {
    estimatorBoldermanPredEst();
  }
  // Save new estimation in the state
  estimatorBoldermanStateSave(state, osTick);
}


// Dynamics to make a prediction of the next state
static void estimatorBoldermanPredict(float dt, float thrust)
{
  // Define the square root of Pxx (cholesky decomosition)
  float delta[N][N] = {{0}};
  int status = 0;
  status = cholesky_decomposition(Pxx,delta,N);
  if (status == 0) {
    // Then no good solution delta is found, so reinitialize Pxx
    resetPxx();
    status = cholesky_decomposition(Pxx,delta,N);
  }
  // Define the sigma points
  for (int ii = 0; ii<NSIGMA; ii++) {
    for (int jj = 0; jj<N; jj++) {
      // define sigmaX for each entry
      if (ii == 0) {
        // First sigmapoint is mean
        sigmaX[jj][ii] = x[jj];
      } else if ((ii>0) && (ii<N+1)) {
        // Sigmapoints is mean + difference
        sigmaX[jj][ii] = x[jj] + sqrtf(lambda+NCALC)*delta[ii-1][jj];
      } else if ((ii>N) && (ii<NSIGMA)) {
        sigmaX[jj][ii] = x[jj] - sqrtf(lambda+NCALC)*delta[ii-N-1][jj];
      }
    }
  }


  // integrate one timestep
  for (int ii = 0; ii<NSIGMA; ii++) {   // Integrate one step (somewhat different when position measurement is available)
    if ((thrust > GRAVITY_MAGNITUDE) || (sigmaX[2][ii]>0.0f)) {
      fly = 1;
    }

    // Adjust value of pitch IF == PI/2 + n*pi
    if (1.0f/cosf(sigmaX[10][ii]) > UPPERBOUND) {
      if (sinf(sigmaX[10][ii]) > 0.0f) {
        sigmaX[10][ii] = acosf(1.0f/UPPERBOUND);
      } else {
        sigmaX[10][ii] = -acosf(1.0f/UPPERBOUND);
      }
    } else if (1.0f/cosf(sigmaX[10][ii]) < -UPPERBOUND) {
      if (sinf(sigmaX[10][ii]) > 0.0f) {
        sigmaX[10][ii] = PI - acosf(1.0f/UPPERBOUND);
      } else {
        sigmaX[10][ii] = -PI + acosf(1.0f/UPPERBOUND);
      }
    }

    // Update the sigmapoints (discretized using forward euler)
    // Furthermore, note that the gyroscopic measurements are directly send through the dynamics as an input
    if (fly==1) {
      sigmaXplus[0][ii] = sigmaX[0][ii] + dt*sigmaX[3][ii];
      sigmaXplus[1][ii] = sigmaX[1][ii] + dt*sigmaX[4][ii];
      sigmaXplus[2][ii] = sigmaX[2][ii] + dt*sigmaX[5][ii];
      // Velocity
      sigmaXplus[3][ii] = sigmaX[3][ii] + dt*thrust*(cosf(sigmaX[11][ii])*sinf(sigmaX[10][ii])*cosf(sigmaX[9][ii]) + sinf(sigmaX[9][ii])*sinf(sigmaX[11][ii]));
      sigmaXplus[4][ii] = sigmaX[4][ii] + dt*thrust*(sinf(sigmaX[11][ii])*sinf(sigmaX[10][ii])*cosf(sigmaX[9][ii]) - cosf(sigmaX[11][ii])*sinf(sigmaX[9][ii]));
      sigmaXplus[5][ii] = sigmaX[5][ii] + dt*(-GRAVITY_MAGNITUDE + thrust*(cosf(sigmaX[9][ii])*cosf(sigmaX[10][ii])));
      // Acceleration (is the same as the additionstep above divided by dt)
      // This one might be tricky -> actually average acceleration from [k, k+1] -> IS THE SAME AS THE MEASUREMENT
      sigmaXplus[6][ii] = (sigmaXplus[3][ii]-sigmaX[3][ii])/dt;
      sigmaXplus[7][ii] = (sigmaXplus[4][ii]-sigmaX[4][ii])/dt;
      sigmaXplus[8][ii] = (sigmaXplus[5][ii]-sigmaX[5][ii])/dt;
      // Attitude
      sigmaXplus[9][ii]  = sigmaX[9][ii]  + (dt/cosf(sigmaX[10][ii])) * (cosf(sigmaX[10][ii])*omega[0] + sinf(sigmaX[10][ii])*sinf(sigmaX[9][ii])*omega[1] + sinf(sigmaX[10][ii])*cosf(sigmaX[9][ii])*omega[2]);
      sigmaXplus[10][ii] = sigmaX[10][ii] + (dt/cosf(sigmaX[10][ii])) * (cosf(sigmaX[10][ii])*cosf(sigmaX[9][ii])*omega[1] - cosf(sigmaX[10][ii])*sinf(sigmaX[9][ii])*omega[2]);
      sigmaXplus[11][ii] = sigmaX[11][ii] + (dt/cosf(sigmaX[10][ii])) * (sinf(sigmaX[9][ii])*omega[1] + cosf(sigmaX[9][ii])*omega[2]);
    }
    if (sigmaXplus[2][ii] < 0) {   // When landend (can only be done when position measurement is available, otherwise positions will drift)
      fly = 0;
    }
    if (fly == 0) {
      sigmaXplus[0][ii] = sigmaX[0][ii];
      sigmaXplus[1][ii] = sigmaX[1][ii];
      sigmaXplus[2][ii] = 0.0f;                // Quadcopter does not lift !
      // Velocity
      sigmaXplus[3][ii] = 0.0f;
      sigmaXplus[4][ii] = 0.0f;
      sigmaXplus[5][ii] = 0.0f;
      // Acceleration
      sigmaXplus[6][ii] = 0.0f;
      sigmaXplus[7][ii] = 0.0f;
      sigmaXplus[8][ii] = 0.0f;                // Hence no accelleration to (ground defines height)
      // Attitude
      sigmaXplus[9][ii]  = 0.0f;
      sigmaXplus[10][ii] = 0.0f;
      sigmaXplus[11][ii] = sigmaX[11][ii];  // Yaw is not influenced by ground
    }
    // Now compute the output vector
    sigmaYplus[0][ii] = sqrtf((sigmaXplus[0][ii]-anchor[0][0])*(sigmaXplus[0][ii]-anchor[0][0]) + (sigmaXplus[1][ii]-anchor[1][0])*(sigmaXplus[1][ii]-anchor[1][0]) + (sigmaXplus[2][ii]-anchor[2][0])*(sigmaXplus[2][ii]-anchor[2][0]));
    sigmaYplus[1][ii] = sqrtf((sigmaXplus[0][ii]-anchor[0][1])*(sigmaXplus[0][ii]-anchor[0][1]) + (sigmaXplus[1][ii]-anchor[1][1])*(sigmaXplus[1][ii]-anchor[1][1]) + (sigmaXplus[2][ii]-anchor[2][1])*(sigmaXplus[2][ii]-anchor[2][1]));
    sigmaYplus[2][ii] = sqrtf((sigmaXplus[0][ii]-anchor[0][2])*(sigmaXplus[0][ii]-anchor[0][2]) + (sigmaXplus[1][ii]-anchor[1][2])*(sigmaXplus[1][ii]-anchor[1][2]) + (sigmaXplus[2][ii]-anchor[2][2])*(sigmaXplus[2][ii]-anchor[2][2]));
    sigmaYplus[3][ii] = cosf(sigmaXplus[9][ii])*cosf(sigmaXplus[11][ii])*sigmaXplus[6][ii] + cosf(sigmaXplus[9][ii])*sinf(sigmaXplus[11][ii])*sigmaXplus[7][ii] - sinf(sigmaXplus[10][ii])*(sigmaXplus[8][ii]+GRAVITY_MAGNITUDE);
    sigmaYplus[4][ii] = (cosf(sigmaXplus[11][ii])*sinf(sigmaXplus[10][ii])*sinf(sigmaXplus[9][ii])-cosf(sigmaXplus[9][ii])*sinf(sigmaXplus[11][ii]))*sigmaXplus[6][ii] + (sinf(sigmaXplus[11][ii])*sinf(sigmaXplus[10][ii])*sinf(sigmaXplus[9][ii])+cosf(sigmaXplus[9][ii])*cosf(sigmaXplus[11][ii]))*sigmaXplus[7][ii] + (cosf(sigmaXplus[10][ii])*sinf(sigmaXplus[9][ii]))*(sigmaXplus[8][ii]+GRAVITY_MAGNITUDE);
    sigmaYplus[5][ii] = (cosf(sigmaXplus[11][ii])*sinf(sigmaXplus[10][ii])*cosf(sigmaXplus[9][ii])+sinf(sigmaXplus[9][ii])*sinf(sigmaXplus[11][ii]))*sigmaXplus[6][ii] + (sinf(sigmaXplus[11][ii])*sinf(sigmaXplus[10][ii])*cosf(sigmaXplus[9][ii])-sinf(sigmaXplus[9][ii])*cosf(sigmaXplus[11][ii]))*sigmaXplus[7][ii] + (cosf(sigmaXplus[10][ii])*cosf(sigmaXplus[9][ii]))*(sigmaXplus[8][ii]+GRAVITY_MAGNITUDE);
  }

  // Calculate the predicted state and the predicted output
  for (int ii=0; ii<N; ii++) {      // Predicted state
    xpred[ii] = 0.0f;
    for (int jj=0; jj<NSIGMA; jj++) {
      xpred[ii] += wm[jj]*sigmaXplus[ii][jj];
    }
  }
  // Exerting the bounds
  while ((xpred[9] > PI) || (xpred[9] < -PI)) {
    if (xpred[9] > PI) {
      xpred[9] -= 2*PI;
    } else if (xpred[9] < -PI) {
      xpred[9] += 2*PI;
    }
  }
  while ((xpred[10] > PI) || (xpred[10] < -PI)) {
    if (xpred[10] > PI) {
      xpred[10] -= 2*PI;
    } else if (xpred[10] < -PI) {
      xpred[10] += 2*PI;
    }
  }
  while ((xpred[11] > PI) || (xpred[11] < -PI)) {
    if (xpred[11] > PI) {
      xpred[11] -= 2*PI;
    } else if (xpred[11] < -PI) {
      xpred[11] += 2*PI;
    }
  }
  for (int ii=0; ii<NOUT; ii++) {   // Predicted output
    ypred[ii] = 0.0f;
    for (int jj=0; jj<NSIGMA; jj++) {
      ypred[ii] += wm[jj]*sigmaYplus[ii][jj];
    }
  }

  // Using the predicted state and the predicted output, calculate the covariances Pxx Pxy and Pyy
  /* Following formulas used:
      Pxx = sum_k( w_k(sigmaXplus_k - xpred)(sigmaXplus_k - xpred)^T ) + Q
      Pxy = sum_k( w_k(sigmaXplus_k - xpred)(sigmaYplus_k - ypred)^T )
      Pyy = sum_k( w_k(sigmaYplus_k - ypred)(sigmaYplus_k - ypred)^T ) + R
  */
  for (int ii = 0; ii<N; ii++) {  // Pxx, Pxy
    for (int jj = 0; jj<N; jj++) {  // Pxx
      Pxx[ii][jj] = q[ii][jj];
      for (int kk = 0; kk<NSIGMA; kk++) {
        Pxx[ii][jj] += wc[kk] * (sigmaXplus[ii][kk] - xpred[ii]) * (sigmaXplus[jj][kk] - xpred[jj]);
      }
    }
  }
}

// Perform the update rule, using the prediction and measurement
static void estimatorBoldermanDynMeas(void)
{
  // Calculate Pxy, Pyy and K. Performed here, because only necessary when position measurements are available
  for (int ii = 0; ii<N; ii++) {
    for (int jj = 0; jj<NOUT; jj++) {   // Pxy
      Pxy[ii][jj] = 0.0f;
      for (int kk = 0; kk<NSIGMA; kk++) {
        Pxy[ii][jj] += wc[kk] * (sigmaXplus[ii][kk] - xpred[ii]) * (sigmaYplus[jj][kk] - ypred[jj]);
      }
    }
  }
  for (int ii = 0; ii<NOUT; ii++) {   // Pyy
    for (int jj = 0; jj<NOUT; jj++) {   // Pyy
      Pyy[ii][jj] = r[ii][jj];
      for (int kk = 0; kk<NSIGMA; kk++) {
        Pyy[ii][jj] += wc[kk] * (sigmaYplus[ii][kk] - ypred[ii]) * (sigmaYplus[jj][kk] - ypred[jj]);
      }
    }
  }

  // Pyy needs to be inversed, to do so we need an extra variable. Due too the method being Gaussian elimination
  static float Pyyhelp[NOUT][NOUT];
  for (int ii=0; ii<NOUT; ii++) {
    for (int jj=0; jj<NOUT; jj++) {
      Pyyhelp[ii][jj] = Pyy[ii][jj];
    }
  }
  // Define the kalman gain ( K = Pxy * inv(Pyy) )
  static arm_matrix_instance_f32 Pyymhelp = {NOUT, NOUT, (float *)Pyyhelp};
  static arm_matrix_instance_f32 Pyyinvm = {NOUT, NOUT, (float *)Pyyinv};
  static arm_matrix_instance_f32 Pxym = {N, NOUT, (float *)Pxy};
  static arm_matrix_instance_f32 Km = {N, NOUT, (float *)k};
  // Do multiplications and inverse
  mat_inv(&Pyymhelp, &Pyyinvm);       // Inverse of Pyy
  mat_mult(&Pxym, &Pyyinvm, &Km);     // K = Pxy inv(Pyy)


  // UKF UPDATE OF THE STATE
  for (int ii=0; ii<N; ii++) {    // Difference with/without posiotion measurement not needed (K depends on measurement)
    x[ii] = xpred[ii];
    for (int jj=0; jj<NOUT; jj++) {
      x[ii] += k[ii][jj] * (y[jj]-ypred[jj]);
    }
  }

  // UKF UPDATE OF THE STATE COVARIANCE
  // Pxxnew = Pxx - K Pyy K'
  static float KPyy[N][NOUT];
  static arm_matrix_instance_f32 KPyym = {N, NOUT, (float *)KPyy};
  static float Ktrans[NOUT][N];
  static arm_matrix_instance_f32 Ktransm = {NOUT, N, (float *)Ktrans};
  static float KPyyK[N][N];
  static arm_matrix_instance_f32 KPyyKm = {N, N, (float *)KPyyK};
  static arm_matrix_instance_f32 Pyym = {NOUT, NOUT, (float *)Pyy};

  // Do multiplations and transpose
  mat_mult(&Km, &Pyym, &KPyym);
  mat_trans(&Km, &Ktransm);
  mat_mult(&KPyym, &Ktransm, &KPyyKm);
  // Determine the new Pxx
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      Pxx[ii][jj] -= KPyyK[ii][jj];
    }
  }
  tracePxx = 0.0f;
  for (int ii=0; ii<N; ii++) {
    tracePxx += Pxx[ii][ii];
  }
  if (tracePxx > UPPERBOUND_TRACE) {
    consolePrintf("Trace too high, resetting Pxx. \n");
    resetPxx();
    tracePxx = 0.0f;
    for (int ii=0; ii<N; ii++) {
      tracePxx += Pxx[ii][ii];
    }
  }
}
// Function used to set prediction to estimation, when no measurements available
static void estimatorBoldermanPredEst(void) {
  for (int ii=0; ii<N; ii++) {
    // Set predicted value to estimated, because no update can be performed, due too the lack of measuremets
    x[ii] = xpred[ii];
  }
}

// Save the estimated state in the variable state
static void estimatorBoldermanStateSave(state_t *state, uint32_t osTick)
{
  // Input here the method to save the improved estimation in the variable state...
  // Postion -> GLOBAL FRAME
  state->position = (point_t) {
    .timestamp = osTick,
    .x = x[0],
    .y = x[1],
    .z = x[2]
  };
  // Velocity -> GLOBAL FRAME
  state->velocity = (velocity_t) {
    .timestamp = osTick,
    .x = x[3],
    .y = x[4],
    .z = x[5]
  };
  // Acceleration -> GLOBAL FRAME (without the gravity and in unit G)
  state->acc = (acc_t) {
    .timestamp = osTick,
    .x = ((x[6]) / GRAVITY_MAGNITUDE),
    .y = ((x[7]) / GRAVITY_MAGNITUDE),
    .z = ((x[8]) / GRAVITY_MAGNITUDE)
  };
  // Attitude
  // STATED IN stabilizer_types.h AS LEGACY CF2 BODY COORDINATES, WHERE PITCH IS INVERTED????????
  state->attitude = (attitude_t) {
    .timestamp = osTick,
    .roll = (x[9] * RAD_TO_DEG),
    .pitch = -1.0f*(x[10] * RAD_TO_DEG),
    .yaw = (x[11] * RAD_TO_DEG)
  };
}






/* INITIALIZATION AND TESTING
*/
// Test function, can test whether the estimatorBolderman is readily initialized
bool estimatorBoldermanTest(void)
{
  return isInit;
}

// Sometimes Pxx needs to be reinitilized
static void resetPxx(void) {
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      if (ii == jj) {
        Pxx[ii][jj] = Pxxdiag[ii];
      } else {
        Pxx[ii][jj] = 0.0f;
      }
    }
  }
}


// Initializing function; is called first, necessary to perform the estimation
void estimatorBoldermanInit(void) {
  // Provide user with information about the initialization
  consolePrintf("Initializing the estimator");
  // Initializing is actually getting the first prediction, defined here
  lastPrediction = xTaskGetTickCount();

  // Define the measurement Queue
  twrDataQueue = xQueueCreate(TWR_MEASUREMENT_QUEUE_LENGTH, sizeof(distanceMeasurement_t));


  // Reset accumulators and counters
  // Often theseinverse 3x3 matrix are readily zero, but when starting over, these values need to be resetted
  accAccumulator = (Axis3f){.axis={0}};
  accAccumulatorCount = 0;
  gyroAccumulator = (Axis3f){.axis={0}};
  gyroAccumulatorCount = 0;
  magAccumulator = (Axis3f){.axis={0}};
  magAccumulatorCount = 0;
  thrustAccumulator = 0;
  thrustAccumulatorCount = 0;

  // LAMBDA scaling parameter
  lambda = alpha*alpha*(NCALC+kappa)-NCALC;

  // WEIGHTS for computing mean and covariance
  for (int ii = 0; ii<NSIGMA; ii++) {
    if (ii == 0) {
      wc[ii] = lambda/(NCALC+lambda) + (1-alpha*alpha+beta);
      wm[ii] = lambda/(NCALC+lambda);
    } else {
      wc[ii] = 1.0f/(2.0f*(NCALC+lambda));
      wm[ii] = 1.0f/(2.0f*(NCALC+lambda));
    }
  }

  // INITIAL state estimate
  /*
  x = {0,0,0, 0,0,0, 0,0,0, 0,0,0};
  */

  // COVARIANCES and INITIAL STATE
  // Estimation covariance
  for (int ii=0; ii<N; ii++) {
    x[ii] = 0.0f;
    for (int jj=0; jj<N; jj++) {
      if (ii == jj) {
        Pxx[ii][jj] = Pxxdiag[ii];
      } else {
        Pxx[ii][jj] = 0.0f;
      }
    }
  }
  // Process noise covariance
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      if (ii == jj) {
        q[ii][jj] = qdiag[ii];
      } else {
        q[ii][jj] = 0.0f;
      }
    }
  }
  // Measurement noise covariance
  for (int ii=0; ii<NOUT; ii++) {
    for (int jj=0; jj<NOUT; jj++) {
      if (ii == jj) {
        r[ii][jj] = rdiag[ii];
      } else {
        r[ii][jj] = 0.0f;
      }
    }
  }


  // Provide the boolean stating the initialization is performed
  isInit = true;
  // Provide user with information about the initialization
  consolePrintf(" -> Initializing the estimator finished. ");
}




/*  CHOLESKY DECOMPOSITION
  Written by Marcus Greiff
*/
int assert_element(float val){
  /***********************************************************************
   * Here we check that the diagonal element of R is greater than or
   * or equal to some limit, ensuring that the decomposition is done
   * correctly and allowing detection of non-PSD matrices A
   **********************************************************************/
  float limit = 0.000001f;
  if (val < limit){
    consolePrintf("Cholesky failed, lower than limit. \n");
    return 1;
  }
  if (isnan(val)){
    consolePrintf("Cholesky failed, not a number. \n");
    return 1;
  }
  if (isinf(val)){
    consolePrintf("Cholesky failed, infinite number. \n");
    return 1;
  }
  return 0;
}
int cholesky_decomposition(float (*A)[N], float (*R)[N], int n)
{
  /***********************************************************************
   * CHOLESKY DECOMPOSITION
   * Written by Marcus Greiff 17/02/2019 based on Numerical Linear
   * Algebra, Lloyd N. Trefethen.
   ***********************************************************************/
  int ii, jj, i, j, q;
  for (ii = 0; ii < n; ii++){
    for (jj = ii; jj < n; jj++){
      R[jj][ii] = 0;
      R[ii][jj] = A[ii][jj];
    }
  }
  for (q = 0; q < n; q++){
    for (j = q+1; j<n; j++){
      for (i = j; i < n; i++){
        if (assert_element(R[q][q])){return 0;}
        R[j][i] -= R[q][i]*R[q][j]/R[q][q];
      }
    }
    for (i = n-1; i>=q; i--){
      if (assert_element(R[q][q])){return 0;}
      R[q][i] /= sqrtf(R[q][q]);
    }
  }
  return 1;
}

// Function called from outside to add measurement to Queue
bool estimatorBoldermanEnqueueDistance(distanceMeasurement_t *measurement)
{
  ASSERT(isInit);
  return xQueueSend(twrDataQueue, measurement, 0);
}

// State all wanted outputs, which can be plotted on the crazyflie client
// Estimation group
LOG_GROUP_START(BOLDERMAN_est)
  LOG_ADD(LOG_FLOAT, tracePxx,   &tracePxx)
  LOG_ADD(LOG_FLOAT, pos_x_est,  &x[0])
  LOG_ADD(LOG_FLOAT, pos_y_est,  &x[1])
  LOG_ADD(LOG_FLOAT, pos_z_est,  &x[2])
  LOG_ADD(LOG_FLOAT, vel_x_est,  &x[3])
  LOG_ADD(LOG_FLOAT, vel_y_est,  &x[4])
  LOG_ADD(LOG_FLOAT, vel_z_est,  &x[5])
  LOG_ADD(LOG_FLOAT, acc_x_est,  &x[6])
  LOG_ADD(LOG_FLOAT, acc_y_est,  &x[7])
  LOG_ADD(LOG_FLOAT, acc_z_est,  &x[8])
  LOG_ADD(LOG_FLOAT, roll_est,   &x[9])
  LOG_ADD(LOG_FLOAT, pitch_est,  &x[10])
  LOG_ADD(LOG_FLOAT, yaw_est,    &x[11])
LOG_GROUP_STOP(BOLDERMAN_est)

// Measurement group
LOG_GROUP_START(BOLDERMAN_meas)
  LOG_ADD(LOG_FLOAT, acc_x_meas, &y[3])
  LOG_ADD(LOG_FLOAT, acc_y_meas, &y[4])
  LOG_ADD(LOG_FLOAT, acc_z_meas, &y[5])
  LOG_ADD(LOG_FLOAT, gyr_x_meas, &omega[0])
  LOG_ADD(LOG_FLOAT, gyr_y_meas, &omega[1])
  LOG_ADD(LOG_FLOAT, gyr_z_meas, &omega[2])
LOG_GROUP_STOP(BOLDERMAN_meas)

PARAM_GROUP_START(BOLDERMAN_param)
PARAM_ADD(PARAM_UINT8, force_fly, &fly)
PARAM_GROUP_STOP(BOLDERMAN_param)
