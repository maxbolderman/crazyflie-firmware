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

/*******************************************************************************
*   Internal variables for the estimator

/*******************************************************************************
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
static const uint32_t n = 12;               // Dimension of dynamical model (used by Bolderman)
static bool positionmeasurement = false;    // Boolean variable stating whether there is observation of position (0 = no position measurement, 1 = position measurement)
// Output = [ax, ay, az] --> acceleration in body coordinates
// Output = [x, y, z, ax, ay, az] --> position in global coordinates and acceleration in body coordinates
static const uint32_t nout = 3;             // Dimension of output measurements (used by Bolderman); (acc and gyro)
static const uint32_t nsigma = 2*n+1;       // Number of sigma points used for UKF
static float alpha = 0.001;                 // Defines spread of sigma points
static float kappa = 0;                     // Secondary scaling factor
static float beta = 0;                      // Distributation parameter, gaussion := 2
static float lambda = alpha*alpha*(n+kappa)-n; // Scaling parameter used for determining: Weights and Sigma points
static float[nsigma] wm;                    // Weights for calculating the mean (initialized in estimaterBoldermanInit())
static float[nsigma] wc;                    // Weights for calculating the covariance (initialized in estimaterBoldermanInit())
static float[n][n] q;                       // Covariance of process noise (zero-mean Gaussian distributed)
static float[nout][nout] r;                 // Covariance of measurement noise (zero-mean Gaussian distributed)
// State momentary (x, y, z, xdot, ydot, zdot, xddot, yddot, zddot, pitch, yaw, roll)
static float[n] x;                          // State column, containing values of the states
static float[n][nsigma] sigmaX;             // Matrix containing all sigmapoints
static float[n][nsigma] sigmaXplus;         // Matrix containing updated sigmapoints
static float[nout][nsigma] sigmaYplus;      // Matrix containing output of sigma points
static float[n] xpred;                      // Predicted state vector
static float[nout] ypred;                   // Predicted output of outputted sigmapoints
static float[n][n] Pxx;                     // Covariance of state Estimation
static float[n][nout] Pxy;                  // Covariance of state - outputs
static float[nout][nout] Pyy;               // Covariance of ouput
static float[n][nout] k;                    // Kalman gain (Pxy * inv(Pyy))
static uint32_t tick;                       // Tick variable used for updating state
static bool fly = false;                    // Boolean stating that the crazyflie flying or not
static float y[3][4];                       // Measurement vectors (just to show how to use static variables)


// Prototype of update function, which is used in estimatorBolderman function
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt);
// Prototype of function combining prediction and measurement = actual update
static void estimatorBoldermanDynMeas();
// Prototype of function saving the updated state
static void estimatorBoldermanStateSave();
// Prototype of function defining the dynamics
static void estimatorBoldermanPredict();


// Function receiving sensordata, control input and orchestrates update
void estimatorBolderman(state_t *state, sensorData_t *sensors, control_t *control, const uint32_t tick)
{

  // Asks for tick at moment of start calling this function (estimatorBolderman() )
  uint32_t osTick = xTaskGetTickCount(); // would be nice if this had a precision higher than 1ms...
  tick = osTick;                         // Used when defining the updated state

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
    if (fly == false) && (thrustAccumulator > GRAVITY_MAGNITUDE) {
      fly = true;
    }

    // Computes the time step used
    float dt = (float)(osTick-lastPrediction)/configTICK_RATE_HZ;

    // This is where we perform the time-integration of the state estimate
    estimatorBoldermanUpdate(state, thrustAccumulator, &accAccumulator, &gyroAccumulator, &magAccumulator, dt);

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
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt)
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
  estimatorBoldermanStateSave(state);
}


// Dynamics to make a prediction of the next state
static void estimatorBoldermanPredict(float dt, float thrust)
{
  // Input here dynamics, to compute a new
  // Define sigmapoints
  // sigmaX = [xhat, xhat+sqrtm((n+lambda)*Pxx), xhat-sqrtm((n+lambda)*Pxx)]
  sigmaX = ....;
  // integrate one timestep
  for (int ii = 0; ii<nsigma; ii++) {
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
    } else {
    /* Define here the dynamics when position measurement is known (difference between not flying)
       Not yet necessary, since position measurement is not yet included.
       Here it is doable to detect landing
       COULD BE DONE DIFFERENTLY (NOT ALL DYNAMICS SHOULD BE SPLITTED I THINK, COMBINE WITH ABOVE)*/
    }
  }
  // Calculate the predicted state and the predicted output
  for (int ii=0; ii<n; ii++) {      // Predicted state
    xpred[ii] = 0;
    for (int jj=0; jj<nsigma; jj++) {
      xpred[ii] += wm[jj]*sigmaXplus[ii][jj];
    }
  }
  for (int ii=0; ii<nout; ii++) {   // Predicted output
    ypred[ii] = 0;
    for (int jj=0; jj<nsigma; jj++) {
      ypred[ii] += wm[jj]*sigmaYplus[ii][jj];
    }
  }

  // Using the predicted state and the predicted output, calculate the covariances Pxx Pxy and Pyy
  Pxx = .......;
  Pxy = .......;
  Pyy = .......;
  K = .......;

}

// Perform the update rule, using the prediction and measurement
static void estimatorBoldermanDynMeas()
{
  // Input here the UKF update...
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
}

// Save the estimated state in the variable state
static void estimatorBoldermanStateSave(state_t *state)
{
  // Input here the method to save the improved estimation in the variable state...
  // Postion
  state->position = (point_t) {
    .timestamp = tick,
    .x = x[0],
    .y = x[1],
    .z = x[2]
  };
  // Velocity
  state->velocity = (velocity_t) {
    .timestamp = tick,
    .x = x[3],
    .y = x[4],
    .z = x[5]
  };
  // Acceleration (without the gravity and in unit G)
  state->acc = (acc_t) {
    .timestamp = tick,
    .x = ((x[6]) / CONTROL_TO_ACC),
    .y = ((x[7]) / CONTROL_TO_ACC),
    .z = ((x[8] + GRAVITY_MAGNITUDE) / CONTROL_TO_ACC);
  };
  // Attitude
  // STATED IN stabilizer_types.h AS LEGACY CF2 BODY COORDINATES, WHERE PITCH IS INVERTED????????
  state->attitude = (attitude_t) {
    .timestamp = tick,
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
  // Covariance of initial Estimation
  for (int ii = 0; ii<n; ii++) {
    Pxx[ii][ii] = 0.1;
    for (int jj = 0; jj<n; jj++) {
      Pxy[ii][jj] = 0;
    }
  }
  // Weights for computing mean and covariance
  for (int ii = 0; ii<nsigma; ii++) {
    if (ii == 0) {
      wc[ii] = lambda/(n+lambda) + (1-alpha*alpha+beta);
      wm[ii] = lambda/(n+lambda);
    } else {
      wc[ii] = 1/(2*(n+lambda));
      wm[ii] = 1/(2*(n+lambda));
    }
  }
  // Noise definitions
  // covariance process noise; q
  // Use same for loop to immediately initialize x and xpred
  // Also initialize the sigmapoints to zero
  for (int ii = 0; ii<n; ii++) {
    x[ii] = 0;
    xpred[ii] = 0;
    q[ii][ii] = (1/100)*0.01;
    for (int jj = 0; jj<nsigma; jj++) {
      sigmaX[ii][jj] = 0;
      sigmaXplus[ii][jj] = 0;
    }
  }
  // covariance measurement noise; r
  // Use same for loop to immediately initialize ypred, which should be done using output equations
  // Also initialize the output of the sigmapoints to zero
  for (int ii = 0; ii<nout; ii++) {
    r[ii][ii] = 0.001;
    ypred[ii] = 0;
    for (int jj = 0; jj<nout; jj++) {
      Pyy[ii][jj] = 0;
    }
    for (int jj = 0; jj<nsigma; jj++) {
      sigmaYplus = 0;
    }
  }


  // Provide the boolean stating the initialization is performed
  isInit = true;
  // Provide user with information about the initialization
  consolePrintf("Initializing the estimator finished");
}

// State all wanted outputs, which can be plotted on the crazyflie client
LOG_GROUP_START(BOLDERMAN)
  LOG_ADD(LOG_FLOAT, acc_x, &y[0][0])
  LOG_ADD(LOG_FLOAT, acc_y, &y[0][1])
  LOG_ADD(LOG_FLOAT, acc_z, &y[0][2])
  LOG_ADD(LOG_FLOAT, mag_x, &y[1][0])
  LOG_ADD(LOG_FLOAT, mag_y, &y[1][1])
  LOG_ADD(LOG_FLOAT, mag_z, &y[1][2])
  LOG_ADD(LOG_FLOAT, gyr_x, &y[2][0])
  LOG_ADD(LOG_FLOAT, gyr_y, &y[2][1])
  LOG_ADD(LOG_FLOAT, gyr_z, &y[2][2])
LOG_GROUP_STOP(BOLDERMAN)
