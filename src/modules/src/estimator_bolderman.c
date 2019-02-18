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
#define NOUT 3                            // Dimension of measurement
#define NSIGMA (2*N+1)                    // Number of sigma points
// Output = [ax, ay, az] --> acceleration in body coordinates
// Output = [x, y, z, ax, ay, az] --> position in global coordinates and acceleration in body coordinates
#define UPPERBOUND 10000                  // Upperbound after which inversion is adjusted

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
static float alpha = 0.001f;                // Defines spread of sigma points
static float kappa = 0.0f;                  // Secondary scaling factor
static float beta = 0.0f;                   // Distributation parameter, gaussion := 2
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
static bool fly = false;                    // Boolean stating that the crazyflie flying or not
static float y[3][4];                       // Measurement vectors (just to show how to use static variables)
// Functions used for multiplication
static inline void mat_inv(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_inverse_f32(pSrc, pDst)); }
static inline void mat_mult(const arm_matrix_instance_f32 * pSrcA, const arm_matrix_instance_f32 * pSrcB, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_mult_f32(pSrcA, pSrcB, pDst)); }



// Prototype of update function, which is used in estimatorBolderman function
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt, uint32_t osTick);
// Prototype of function combining prediction and measurement = actual update
static void estimatorBoldermanDynMeas(void);
// Prototype of function saving the updated state
static void estimatorBoldermanStateSave(state_t *state, uint32_t osTick);
// Prototype of function defining the dynamics
static void estimatorBoldermanPredict(float dt, float thrust);

// Prototypes of function By Marcus Greiff (cholesky decomposition)
int  assert_element(float val);
int  cholesky_decomposition(float (*A)[N], float (*R)[N], int n);



/*
HERE THE MODEL STARTS
---------------------------------------------------------------
*/
// Function receiving sensordata, control input and orchestrates update
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
    // When not flying, define when quadcopter lifts
    // Afterwards landing can only be detected WHEN USING POSITION MEASUREMENTS
    if ((fly == false) && (thrustAccumulator > GRAVITY_MAGNITUDE)) {
      fly = true;
    }

    // Computes the time step used
    float dt = (float)(osTick-lastPrediction)/configTICK_RATE_HZ;

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

// Here the integration is performed and the state estimation is performed
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt, uint32_t osTick)
{
  // Store the measuremets y1 y2 y3 as rows in a matrix
  // Acceleration
  y[0][0] = acc->x;
  y[0][1] = acc->y;
  y[0][2] = acc->z;
  // Magnetometer (not used in state estimation)
  y[1][0] = mag->x;
  y[1][1] = mag->y;
  y[1][2] = mag->z;
  // Gyroscope
  y[2][0] = gyro->x;
  y[2][1] = gyro->y;
  y[2][2] = gyro->z;
  // Position (if applicable) (not yet included)
  /*
  if (positionmeasurement == true) {
    y[3][0] = ...;
    y[3][1] = ...;
    y[3][2] = ...;
  }
  */

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
  // Use predicted state and measurement to improve Estimation
  estimatorBoldermanDynMeas();
  // Save new estimation in the state
  estimatorBoldermanStateSave(state, osTick);
}


// Dynamics to make a prediction of the next state
static void estimatorBoldermanPredict(float dt, float thrust)
{
  // Define the square root of Pxx (cholesky decomosition)
  float delta[N][N] = {0};
  int status = 0;
  status = cholesky_decomposition(Pxx,delta,N);
  if (status == 0) {
    // Then no good solution delta is found, so reinitialize Pxx
    reinitializePxx();
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
        sigmaX[jj][ii] = x[jj] + sqrt(lambda+n)*delta[ii-1][jj];
      } else if ((ii>N) && (ii<NSIGMA)) {
        sigmaX[jj][ii] = x[jj] - sqrt(lambda+n)*delta[ii-N-1][jj];
      }
    }
  }


  // integrate one timestep
  for (int ii = 0; ii<NSIGMA; ii++) {
    // Adjust value of pitch IF == PI/2 + n*pi
    if (1/sigmaX[10][ii] > UPPERBOUND) {
      if (sin(sigmaX[10][ii]) > 0.0) {
        sigmaX[10][ii] = acos(1/UPPERBOUND);
      } else {
        sigmaX[10][ii] = -acos(1/UPPERBOUND);
      }
    } else if (1/sigmaX[10][ii] < -UPPERBOUND) {
      if (sin(sigmaX[10][ii]) > 0.0) {
        sigmaX[10][ii] = PI - acos(1/UPPERBOUND);
      } else {
        sigmaX[10][ii] = -PI + acos(1/UPPERBOUND);
      }
    }
    // Update the sigmapoints (discretized using forward euler)
    // Furthermore, note that the gyroscopic measurements are directly send through the dynamics as an input
    if (positionmeasurement == false) {
      if (fly) {    // When flying
        // Position
        sigmaXplus[0][ii] = sigmaX[0][ii] + dt*sigmaX[3][ii];
        sigmaXplus[1][ii] = sigmaX[1][ii] + dt*sigmaX[4][ii];
        sigmaXplus[2][ii] = sigmaX[2][ii] + dt*sigmaX[5][ii];
        // Velocity
        sigmaXplus[3][ii] = sigmaX[3][ii] + dt*thrust*(cos(sigmaX[11][ii])*sin(sigmaX[10][ii])*cos(sigmaX[9][ii]) + sin(sigmaX[9][ii])*sin(sigmaX[11][ii]));
        sigmaXplus[4][ii] = sigmaX[4][ii] + dt*thrust*(sin(sigmaX[11][ii])*sin(sigmaX[10][ii])*cos(sigmaX[9][ii]) - sin(sigmaX[9][ii])*cos(sigmaX[11][ii]));
        sigmaXplus[5][ii] = sigmaX[5][ii] + dt*(thrust*cos(sigmaX[10][ii])*cos(sigmaX[9][ii]) - g);
        // Acceleration (is the same as the additionstep above divided by dt)
        sigmaXplus[6][ii] = (sigmaXplus[3][ii]-sigmaX[3][ii])/dt;
        sigmaXplus[7][ii] = (sigmaXplus[4][ii]-sigmaX[4][ii])/dt;
        sigmaXplus[8][ii] = (sigmaXplus[5][ii]-sigmaX[5][ii])/dt;
        // Attitude
        sigmaXplus[9][ii] = sigmaX[9][ii]   + (dt/cos(sigmaX[10][ii]))*(cos(sigmaX[10][ii])*y[2][0] + sin(sigmaX[10][ii])*sin(sigmaX[9][ii])*y[2][1] + sin(sigmaX[10][ii])*cos(sigmaX[9][ii])*y[2][2]);
        sigmaXplus[10][ii] = sigmaX[10][ii] + (dt/cos(sigmaX[10][ii]))*(cos(sigmaX[10][ii])*cos(sigmaX[9][ii])*y[2][1] - cos(sigmaX[10][ii])*sin(sigmaX[9][ii])*y[2][2]);
        sigmaXplus[11][ii] = sigmaX[11][ii] + (dt/cos(sigmaX[10][ii]))*(sin(sigmaX[9][ii])*y[2][1] + cos(sigmaX[9][ii])*y[2][2]);
      } else {      // When not flying (NO SLIP ASSUMED, however position is already unrealiable)
        // Position
        sigmaXplus[0][ii] = sigmaX[0][ii];
        sigmaXplus[1][ii] = sigmaX[1][ii];
        sigmaXplus[2][ii] = 0;                // Quadcopter does not lift !
        // Velocity
        sigmaXplus[3][ii] = 0;
        sigmaXplus[4][ii] = 0;
        sigmaXplus[5][ii] = 0;
        // Acceleration
        sigmaXplus[6][ii] = 0;
        sigmaXplus[7][ii] = 0;
        sigmaXplus[8][ii] = 0;                // Hence no accelleration to (ground defines height)
        // Attitude
        sigmaXplus[9][ii] = 0;
        sigmaXplus[10][ii] = 0;
        sigmaXplus[11][ii] = sigmaX[11][ii];  // Yaw is not influenced by ground
      }
      // Determine the predicted output (using the state)
      // Here only acceleration is measured (IN BODY COORDINATES)
      sigmaYplus[0][ii] = cos(sigmaXplus[9][ii])*cos(sigmaXplus[11][ii])*sigmaXplus[6][ii] + sin(sigmaXplus[11][ii])*cos(sigmaXplus[9][ii])*sigmaXplus[7][ii] - sin(sigmaXplus[10][ii])*sigmaXplus[8][ii];
      sigmaYplus[1][ii] = (cos(sigmaXplus[11][ii])*sin(sigmaXplus[10][ii])*sin(sigmaXplus[9][ii])-sin(sigmaXplus[11][ii])*cos(sigmaXplus[9][ii]))*sigmaXplus[6][ii] + (sin(sigmaXplus[11][ii])*sin(sigmaXplus[10][ii])*sin(sigmaXplus[9][ii])+cos(sigmaXplus[9][ii])*cos(sigmaXplus[11][ii]))*sigmaXplus[7][ii] + (cos(sigmaXplus[10][ii])*sin(sigmaXplus[9][ii]))*sigmaXplus[8][ii];
      sigmaYplus[2][ii] = (cos(sigmaXplus[11][ii])*sin(sigmaXplus[10][ii])*cos(sigmaXplus[9][ii])+sin(sigmaXplus[11][ii])*sin(sigmaXplus[9][ii]))*sigmaXplus[6][ii] + (sin(sigmaXplus[11][ii])*sin(sigmaXplus[10][ii])*cos(sigmaXplus[9][ii])-sin(sigmaXplus[9][ii])*cos(sigmaXplus[11][ii]))*sigmaXplus[7][ii] + (cos(sigmaXplus[10][ii])*cos(sigmaXplus[9][ii]))*sigmaXplus[8][ii];
    } else {    // Position measurement IS available
      if (fly) {    // When flying
        // Position
        sigmaXplus[0][ii] = sigmaX[0][ii] + dt*sigmaX[3][ii];
        sigmaXplus[1][ii] = sigmaX[1][ii] + dt*sigmaX[4][ii];
        sigmaXplus[2][ii] = sigmaX[2][ii] + dt*sigmaX[5][ii];
        // Velocity
        sigmaXplus[3][ii] = sigmaX[3][ii] + dt*thrust*(cos(sigmaX[11][ii])*sin(sigmaX[10][ii])*cos(sigmaX[9][ii]) + sin(sigmaX[9][ii])*sin(sigmaX[11][ii]));
        sigmaXplus[4][ii] = sigmaX[4][ii] + dt*thrust*(sin(sigmaX[11][ii])*sin(sigmaX[10][ii])*cos(sigmaX[9][ii]) - sin(sigmaX[9][ii])*cos(sigmaX[11][ii]));
        sigmaXplus[5][ii] = sigmaX[5][ii] + dt*(thrust*cos(sigmaX[10][ii])*cos(sigmaX[9][ii]) - g);
        // Acceleration (is the same as the additionstep above divided by dt)
        sigmaXplus[6][ii] = (sigmaXplus[3][ii]-sigmaX[3][ii])/dt;
        sigmaXplus[7][ii] = (sigmaXplus[4][ii]-sigmaX[4][ii])/dt;
        sigmaXplus[8][ii] = (sigmaXplus[5][ii]-sigmaX[5][ii])/dt;
        // Attitude
        sigmaXplus[9][ii] = sigmaX[9][ii]   + (dt/cos(sigmaX[10][ii]))*(cos(sigmaX[10][ii])*y[2][0] + sin(sigmaX[10][ii])*sin(sigmaX[9][ii])*y[2][1] + sin(sigmaX[10][ii])*cos(sigmaX[9][ii])*y[2][2]);
        sigmaXplus[10][ii] = sigmaX[10][ii] + (dt/cos(sigmaX[10][ii]))*(cos(sigmaX[10][ii])*cos(sigmaX[9][ii])*y[2][1] - cos(sigmaX[10][ii])*sin(sigmaX[9][ii])*y[2][2]);
        sigmaXplus[11][ii] = sigmaX[11][ii] + (dt/cos(sigmaX[10][ii]))*(sin(sigmaX[9][ii])*y[2][1] + cos(sigmaX[9][ii])*y[2][2]);
      } else if (fly==false || sigmaXplus[2][ii]<0.0) {      // When not flying (NO SLIP ASSUMED, however position is already unrealiable)
        fly = false;      // Position does not diverge, so quadcopter landed
        // Position
        sigmaXplus[0][ii] = sigmaX[0][ii];
        sigmaXplus[1][ii] = sigmaX[1][ii];
        sigmaXplus[2][ii] = 0;                // Quadcopter does not lift !
        // Velocity
        sigmaXplus[3][ii] = 0;
        sigmaXplus[4][ii] = 0;
        sigmaXplus[5][ii] = 0;
        // Acceleration
        sigmaXplus[6][ii] = 0;
        sigmaXplus[7][ii] = 0;
        sigmaXplus[8][ii] = 0;                // Hence no accelleration to (ground defines height)
        // Attitude
        sigmaXplus[9][ii] = 0;
        sigmaXplus[10][ii] = 0;
        sigmaXplus[11][ii] = sigmaX[11][ii];  // Yaw is not influenced by ground
      }
      // Determine the predicted output (using the state)
      // Here only acceleration is measured (IN BODY COORDINATES)
      sigmaYplus[0][ii] = sigmaXplus[0][ii];
      sigmaYplus[1][ii] = sigmaXplus[1][ii];
      sigmaYplus[2][ii] = sigmaXplus[2][ii];
      sigmaYplus[3][ii] = cos(sigmaXplus[9][ii])*cos(sigmaXplus[11][ii])*sigmaXplus[6][ii] + sin(sigmaXplus[11][ii])*cos(sigmaXplus[9][ii])*sigmaXplus[7][ii] - sin(sigmaXplus[10][ii])*sigmaXplus[8][ii];
      sigmaYplus[4][ii] = (cos(sigmaXplus[11][ii])*sin(sigmaXplus[10][ii])*sin(sigmaXplus[9][ii])-sin(sigmaXplus[11][ii])*cos(sigmaXplus[9][ii]))*sigmaXplus[6][ii] + (sin(sigmaXplus[11][ii])*sin(sigmaXplus[10][ii])*sin(sigmaXplus[9][ii])+cos(sigmaXplus[9][ii])*cos(sigmaXplus[11][ii]))*sigmaXplus[7][ii] + (cos(sigmaXplus[10][ii])*sin(sigmaXplus[9][ii]))*sigmaXplus[8][ii];
      sigmaYplus[5][ii] = (cos(sigmaXplus[11][ii])*sin(sigmaXplus[10][ii])*cos(sigmaXplus[9][ii])+sin(sigmaXplus[11][ii])*sin(sigmaXplus[9][ii]))*sigmaXplus[6][ii] + (sin(sigmaXplus[11][ii])*sin(sigmaXplus[10][ii])*cos(sigmaXplus[9][ii])-sin(sigmaXplus[9][ii])*cos(sigmaXplus[11][ii]))*sigmaXplus[7][ii] + (cos(sigmaXplus[10][ii])*cos(sigmaXplus[9][ii]))*sigmaXplus[8][ii];
    }
  }
  // Calculate the predicted state and the predicted output
  for (int ii=0; ii<N; ii++) {      // Predicted state
    xpred[ii] = 0;
    for (int jj=0; jj<NSIGMA; jj++) {
      xpred[ii] += wm[jj]*sigmaXplus[ii][jj];
    }
  }
  for (int ii=0; ii<NOUT; ii++) {   // Predicted output
    ypred[ii] = 0;
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
  for (int ii = 0; ii<N; ii++) {
    for (int jj = 0; jj<N; jj++) {
      Pxx[ii][jj] = q[ii][jj];
      for (int kk = 0; kk<NSIGMA; kk++) {
        Pxx[ii][jj] += wc[kk] * (sigmaXplus[ii][kk] - xpred[ii]) * (sigmaXplus[jj][kk] - xpred[jj]);
      }
    }
    for (int jj = 0; jj<NOUT; jj++) {
      Pxy[ii][jj] = 0;
      for (int kk = 0; kk<NSIGMA; kk++) {
        Pxy[ii][jj] += wc[kk] * (sigmaXplus[ii][kk] - xpred[ii]) * (sigmaYplus[jj][kk] - ypred[jj]);
      }
    }
  }
  for (int ii = 0; ii<NOUT; ii++) {
    for (int jj = 0; jj<NOUT; jj++) {
      Pyy[ii][jj] = r[ii][jj];
      for (int kk = 0; kk<NSIGMA; kk++) {
        Pyy[ii][jj] += wc[kk] * (sigmaYplus[ii][kk] - ypred[ii]) * (sigmaYplus[jj][kk] - ypred[jj]);
      }
    }
  }
  // Define the kalman gain ( K = Pxy * inv(Pyy) )
  static arm_matrix_instance_f32 Pyym = {NOUT, NOUT, (float *)Pyy};
  static arm_matrix_instance_f32 Pyyinvm = {NOUT, NOUT, (float *)Pyyinv};
  static arm_matrix_instance_f32 Pxym = {N, NOUT, (float *)Pxy};
  static arm_matrix_instance_f32 Km = {N, NOUT, (float *)k};
  mat_inv(&Pyym, &Pyyinvm);           // Inverse of Pyy
  mat_mult(&Pxym, &Pyyinvm, &Km);     // K = Pxy inv(Pyy)
}

// Perform the update rule, using the prediction and measurement
static void estimatorBoldermanDynMeas(void)
{
  // UKF UPDATE OF THE STATE
  // Slightly different depending on whether position information available
  if (positionmeasurement == false) {
    // Update using only acceleration measurements
    for (int ii = 0; ii<n; ii++) {
      x[ii] = xpred[ii] + (k[ii][0]*(y[0][0]-ypred[0]) + k[ii][1]*(y[0][1]-ypred[1]) + k[ii][2]*(y[0][2]-ypred[2]));
    }
  } else {
    // Update using acceleration and position measurements
    for (int ii = 0; ii<n; ii++) {
      x[ii] = xpred[ii] + (k[ii][0]*(y[3][0]-ypred[0]) + k[ii][1]*(y[3][1]-ypred[1]) + k[ii][2]*(y[3][2]-ypred[2]) + k[ii][3]*(y[0][0]-ypred[3]) + k[ii][4]*(y[0][1]-ypred[4]) + k[ii][5]*(y[0][2]-ypred[5]) );
    }
  }
  // UKF UPDATE OF THE STATE COVARIANCE
  // Pxxnew = Pxx - K Pyy K'
  static float KPyy[N][NOUT];
  static arm_matrix_instance_f32 KPyym = {N, NOUT, (float *)KPyy};
  static float Ktrans[NOUT][N];
  static arm_matrix_instance_f32 Ktransm = {NOUT, N, (float *)Ktrans};
  static float KPyyK[N][N];
  static arm_matrix_instance_f32 KPyyKm = {N, N, (float *)KPyyK}
  mat_mult(&Km, &Pyym, &KPyym);
  mat_trans(&Km, &Ktransm);
  mat_mult(&KPyym, &Ktransm, &KPyyKm);
  // Determine the new Pxx
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      Pxx[ii][jj] -= KPyyK[ii][jj];
    }
  }
}

// Save the estimated state in the variable state
static void estimatorBoldermanStateSave(state_t *state, uint32_t osTick)
{
  // Input here the method to save the improved estimation in the variable state...
  // Postion
  state->position = (point_t) {
    .timestamp = osTick,
    .x = x[0],
    .y = x[1],
    .z = x[2]
  };
  // Velocity
  state->velocity = (velocity_t) {
    .timestamp = osTick,
    .x = x[3],
    .y = x[4],
    .z = x[5]
  };
  // Acceleration (without the gravity and in unit G)
  state->acc = (acc_t) {
    .timestamp = osTick,
    .x = ((x[6]) / CONTROL_TO_ACC),
    .y = ((x[7]) / CONTROL_TO_ACC),
    .z = ((x[8] + GRAVITY_MAGNITUDE) / CONTROL_TO_ACC);
  };
  // Attitude
  // STATED IN stabilizer_types.h AS LEGACY CF2 BODY COORDINATES, WHERE PITCH IS INVERTED????????
  state->attitude = (attitude_t) {
    .timestamp = osTick,
    .roll = (x[9] * RAD_TO_DEG),
    .pitch = 1/(x[10] * RAD_TO_DEG),
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
static void reinitializePxx(void) {
  for (int ii = 0; ii<N; ii++) {
    for (int jj = 0; jj<N; jj++) {
      if (ii == jj) {
        Pxx[ii][jj] = 0.1;
      } else {
        Pxx[ii][jj] = 0.0;
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

  // Initialize the variables necessary for UKF (declared above)
  lambda = alpha*alpha*(N+kappa)-N;
  // Covariance of initial Estimation
  for (int ii = 0; ii<N; ii++) {
    Pxx[ii][ii] = 0.1;
    for (int jj = 0; jj<N; jj++) {
      Pxy[ii][jj] = 0;
    }
  }
  // Weights for computing mean and covariance
  for (int ii = 0; ii<NSIGMA; ii++) {
    if (ii == 0) {
      wc[ii] = lambda/(n+lambda) + (1-alpha*alpha+beta);
      wm[ii] = lambda/(n+lambda);
    } else {
      wc[ii] = 1/(2*(N+lambda));
      wm[ii] = 1/(2*(N+lambda));
    }
  }
  // Noise definitions
  // covariance process noise; q
  // Use same for loop to immediately initialize x and xpred
  // Also initialize the sigmapoints to zero
  for (int ii = 0; ii<N; ii++) {
    x[ii] = 0;
    xpred[ii] = 0;
    q[ii][ii] = (1/100)*0.01;
    for (int jj = 0; jj<NSIGMA; jj++) {
      sigmaX[ii][jj] = 0;
      sigmaXplus[ii][jj] = 0;
    }
  }
  // covariance measurement noise; r
  // Use same for loop to immediately initialize ypred, which should be done using output equations
  // Also initialize the output of the sigmapoints to zero
  for (int ii = 0; ii<NOUT; ii++) {
    r[ii][ii] = 0.001;
    ypred[ii] = 0;
    for (int jj = 0; jj<NOUT; jj++) {
      Pyy[ii][jj] = 0;
    }
    for (int jj = 0; jj<NSIGMA; jj++) {
      sigmaYplus = 0;
    }
  }


  // Provide the boolean stating the initialization is performed
  isInit = true;
  // Provide user with information about the initialization
  consolePrintf("Initializing the estimator finished");
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
  float limit = 0.0001;
  if (val < limit){ return 1; }
  if (isnan(val)){ return 1; }
  if (isinf(val)){ return 1; }
  return 0;
}
int cholesky_decomposition(float (*A)[N], float (*R)[N], int n)
{
  /***********************************************************************
   * CHOLESKY DECOMPOSITION
   * Written by Marcus Greiff 17/02/2019 based on Numerical Linear
   * Algebra, Lloyd N. Trefethen.
   ***********************************************************************/
  int ii, jj, i, j, k;
  for (ii = 0; ii < n; ii++){
    for (jj = ii; jj < n; jj++){
      R[jj][ii] = 0;
      R[ii][jj] = A[ii][jj];
    }
  }
  for (k = 0; k < n; k++){
    for (j = k+1; j<n; j++){
      for (i = j; i < n; i++){
        if (assert_element(R[k][k])){return 0;}
        R[j][i] -= R[k][i]*R[k][j]/R[k][k];
      }
    }
    for (i = n-1; i>=k; i--){
      if (assert_element(R[k][k])){return 0;}
      R[k][i] /= sqrtf(R[k][k]);
    }
  }
  return 1;
}

// State all wanted outputs, which can be plotted on the crazyflie client
// SAVE ALL STATE VARIABLES
LOG_GROUP_START(BOLDERMAN_ESTIMATION)
  LOG_ADD(LOG_FLOAT, pos_x, &x[0])
  LOG_ADD(LOG_FLOAT, pos_y, &x[1])
  LOG_ADD(LOG_FLOAT, pos_z, &x[2])
  LOG_ADD(LOG_FLOAT, vel_x, &x[3])
  LOG_ADD(LOG_FLOAT, vel_y, &x[4])
  LOG_ADD(LOG_FLOAT, vel_z, &x[5])
  LOG_ADD(LOG_FLOAT, acc_x, &x[6])
  LOG_ADD(LOG_FLOAT, acc_y, &x[7])
  LOG_ADD(LOG_FLOAT, acc_z, &x[8])
  LOG_ADD(LOG_FLOAT, roll, &x[9])
  LOG_ADD(LOG_FLOAT, pitch, &x[10])
  LOG_ADD(LOG_FLOAT, yaw, &x[11])
LOG_GROUP_STOP(BOLDERMAN_ESTIMATION)

// SAVE ALL MEASUREMENTS
LOG_GROUP_START(BOLDERMAN_MEASUREMENTS)
  LOG_ADD(LOG_FLOAT, acc_x, &y[0][0])
  LOG_ADD(LOG_FLOAT, acc_y, &y[0][1])
  LOG_ADD(LOG_FLOAT, acc_z, &y[0][2])
  LOG_ADD(LOG_FLOAT, mag_x, &y[1][0])
  LOG_ADD(LOG_FLOAT, mag_y, &y[1][1])
  LOG_ADD(LOG_FLOAT, mag_z, &y[1][2])
  LOG_ADD(LOG_FLOAT, gyr_x, &y[2][0])
  LOG_ADD(LOG_FLOAT, gyr_y, &y[2][1])
  LOG_ADD(LOG_FLOAT, gyr_z, &y[2][2])
LOG_GROUP_STOP(BOLDERMAN_MEASUREMENTS)
