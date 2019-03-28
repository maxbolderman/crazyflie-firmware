/*
--------------------------------------------------------------------------------
ESTIMATOR_BOLDERMAN.C
--------------------------------------------------------------------------------
This model uses the UKF (Unscented Kalman Filter) for state estimation of the Quadcopter.
Important aspects of this model:
- STATE and OUTPUT
    * x = [x,y,z, xdot,ydot,zdot, roll,pitch,yaw] -> position and velocity in GLOBAL COORDINATES
    * y = [d1,d2,d3] -> distance of quadcopter to a given ANCHOR POSITION (position of UWB sensor)
- INPUT
    * utrans = [ax,ay,az]eb + [0,0,g]eg -> Acceleration measurements from IMU -> acceleration in BODY COORDINATES, with gravity attribution
    * urot   = [wx,wy,wz]eb             -> Gyroscopic measurements from IMU   -> angular velocity in BODY COORDINATES
--------------------------------------------------------------------------------
WRITTEN BY; Max Bolderman
--------------------------------------------------------------------------------
*/

// INCLUDE all libraries necessary for the state estimation
#include "estimator_bolderman.h"    // header file of this model
#include "outlierFilter.h"
#include "stm32f4xx.h"
#include "FreeRTOS.h"
#include "queue.h"
#include "task.h"
#include "sensors.h"
#include "log.h"                    // Logging variables
#include "param.h"                  // Logging parameters (and changing)
#include "math.h"                   // Math library
#include "cf_math.h"                // Extended library for more involved computations

// DEFINE helpfull parameters
#define PREDICT_RATE (10)           // Update rate
#define DEG_TO_RAD (PI/180.0f)      // Degrees to Radians
#define RAD_TO_DEG (180.0f/PI)      // Radians to Degrees
#define GRAVITY_MAGNITUDE (9.81f)   // Gravitational acceleration
#define CRAZYFLIE_WEIGHT_grams (27.0f) // Weight of the Crazyflie
#define CONTROL_TO_ACC (GRAVITY_MAGNITUDE*60.0f/CRAZYFLIE_WEIGHT_grams/65536.0f)
                                    // Control input to acceleration
#define N (9)                        // Dimension of the state
#define NCALC (9.0f)                // Dimension of the state used for float calculations
#define NOUT (3)                    // Dimension of the measurement space
#define NSIGMA (2*N+1)              // Number of sigmapoints used for prediction
#define TS (1.0f/PREDICT_RATE)      // Time Step
#define MAX_DISTANCE (10.0f)        // Maximum distance for measurement which will be used
#define MIN_DISTANCE (0.01f)        // Minimum distance for measurement which will be used
#define LIMIT_CHOLESKY (0.00001f)
#define LIMIT_INVERT (0.0001f)      // Limit after which the numberical inversion is not feasable
#define LIMIT_TRACE (1000.0f)       // Limit of the covariance matrix trace for reset
#define MAX_ACCELERATION (20.0f)

// VARIABLES standard
static bool busy = false;           // Statemenet if the quadcopter is busy
static bool isInit = false;         // Boolean static if initialization occured yet
static int32_t lastPrediction;      // Time of the last prediction (not in seconds but in ticks)
static Axis3f accAccumulator;       // Acceleration accumulator from measurements
static Axis3f gyroAccumulator;      // Gyroscopic accumulator from measurements
static Axis3f magAccumulator;       // Magnetometer accumulator from measurements
static float thrustAccumulator;     // Thrust accumulator from control
static uint32_t accAccumulatorCount;// Counter of the acceleration measurements
static uint32_t gyroAccumulatorCount;// Counter of the gyroscopic measurements
static uint32_t magAccumulatorCount;// Counter of the magnetometer measurements
static uint32_t thrustAccumulatorCount;// Counter of the thrust inputs
// MEASUREMENT variables
static float acceleration[3];       // Contains acceleration measurements
static float omega[3];              // Contains gyroscopic measurements
static float y[NOUT] = {0.0f};      // Contains distance measurements
static uint32_t yCount = 0;         // Count the number of measurements received
static float anchor[3][NOUT] = {{0.5f}};// Will be filled with anchor positions -> initialization should just be not a real anchor position
// UKF variables
static float alpha = 0.1f;
static float kappa = 0.0f;
static float beta = 2.0f;
static float lambda;
//static float Rgb[3][3] = {{0.0f}};             // Rotation matrix
//static float Winv[3][3] = {{0.0f}};            // Angular rate matrix
static float Wm[NSIGMA];            // Weights for calculating the mean
static float Wc[NSIGMA];            // Weights for calculating the covariance
static float x[N] = {1.0f,4.0f,0.5f, 0.0f,0.0f,0.0f, 0.0f,0.0f,0.0f};
static float xpred[N];                // Predicted state, where measurements are not used to update yet
static float ypred[NOUT];
static float sigmaX[N][NSIGMA] = {{0.0f}};     // Sigmapoints
static float sigmaXplus[N][NSIGMA] = {{0.0f}}; // Sigmapoints propagated through dynamics
static float sigmaYplus[NOUT][NSIGMA] = {{0.0f}};
static float tracePxx = 0.0f;                  // Trace of state covariance matrix
static float Pxx[N][N] = {{0.0f}};             // Covariance state estimate
static float Pxxdiag[N] = {1.0f,1.0f,1.0f, 0.01f,0.01f,0.01f, 0.01f,0.01f,0.01f};
static float q[N][N] = {{0.0f}};
static float qdiag[N] = {0.5f*TS*TS,0.5f*TS*TS,0.5f*TS*TS, 0.5f*TS,0.5f*TS,0.5f*TS, 0.1f*TS,0.1f*TS,0.1f*TS};
static float Pyy[NOUT][NOUT] = {{0.0f}};
static float Pyyhelp[NOUT][NOUT] = {{0.0f}};
static float r[NOUT][NOUT] = {{0.0f}};
static float rdiag[NOUT] = {0.25f,0.25f,0.25f};
static float Pxy[N][NOUT] = {{0.0f}};
//static float determinantPyy;
static float Pyyinv[NOUT][NOUT] = {{0.0f}};
static float deltaPxx[N][N] = {{0.0f}};
static float K[N][NOUT] = {{0.0f}};
static float accglob[3];
static float omegaglob[3];
static float accext[3] = {0.0f,0.0f,0.0f};                        // Acceleration to externalize


// PROTOTYPES of functions used in estimator
// Orchestrates kalman algorithm
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt, uint32_t osTick);
// Makes the prediction
static void estimatorBoldermanPredict(float dt);
// Update with using distance measurements
static void estimatorBoldermanDynMeas(void);
// Update without using distance measurements
static void estimatorBoldermanPredEst(void);
// Externalize the state in the state variable used in control
static void estimatorBoldermanStateSave(state_t *state, uint32_t osTick);

// Functions used for computations
static void resetPxx(void);
//static void calculateRgb(float phi, float theta, float psi);
//static void calculateWinv(float phi, float theta);
static void safetyAngleBounds(void);
// Cholesky decomposition
int assert_element(float val);
int cholesky_decomposition(float (*A)[N], float (*R)[N], int n);

// Math functions
static inline void mat_inv(const arm_matrix_instance_f32 * pSrc, arm_matrix_instance_f32 * pDst)
{ configASSERT(ARM_MATH_SUCCESS == arm_mat_inverse_f32(pSrc, pDst)); }



/*
--------------------------------------------------------------------------------
START MODEL
--------------------------------------------------------------------------------
*/

// Called when a measurement is available
void estimatorBolderman(state_t *state, sensorData_t *sensors, control_t *control, const uint32_t tick)
{
  // Ask for moment this function is called
  uint32_t osTick = xTaskGetTickCount(); // Precision only 1ms

  // Collect all IMU measurements -> sum them -> average when update is called
  if (sensorsReadAcc(&sensors->acc)) {
    accAccumulator.x += GRAVITY_MAGNITUDE*sensors->acc.x;
    accAccumulator.y += GRAVITY_MAGNITUDE*sensors->acc.y;
    accAccumulator.z += GRAVITY_MAGNITUDE*sensors->acc.z;
    accAccumulatorCount++;
  }
  if (sensorsReadGyro(&sensors->gyro)) {
    gyroAccumulator.x += DEG_TO_RAD*sensors->gyro.x;
    gyroAccumulator.y += DEG_TO_RAD*sensors->gyro.y;
    gyroAccumulator.z += DEG_TO_RAD*sensors->gyro.z;
    gyroAccumulatorCount++;
  }
  // MAGNETOMETER MEASUREMENTS ARE NOT USED !!
  if (sensorsReadMag(&sensors->mag)) {
    magAccumulator.x += sensors->mag.x;
    magAccumulator.y += sensors->mag.y;
    magAccumulator.z += sensors->mag.z;
    magAccumulatorCount++;
  }
  thrustAccumulator += CONTROL_TO_ACC*control->thrust;
  thrustAccumulatorCount++;

  // UPDATE only when certain criterias are met
  if ((osTick-lastPrediction) >= configTICK_RATE_HZ/PREDICT_RATE
        && gyroAccumulatorCount > 0
        && accAccumulatorCount > 0
        && thrustAccumulatorCount > 0)
  {
    // AVERAGE the measurements
    thrustAccumulator /= thrustAccumulatorCount;
    accAccumulator.x /= accAccumulatorCount;
    accAccumulator.y /= accAccumulatorCount;
    accAccumulator.z /= accAccumulatorCount;
    gyroAccumulator.x /= gyroAccumulatorCount;
    gyroAccumulator.y /= gyroAccumulatorCount;
    gyroAccumulator.z /= gyroAccumulatorCount;
    magAccumulator.x /= magAccumulatorCount;
    magAccumulator.y /= magAccumulatorCount;
    magAccumulator.z /= magAccumulatorCount;

    /*
    if (accAccumulator.x > MAX_ACCELERATION || accAccumulator.y > MAX_ACCELERATION || accAccumulator.z > MAX_ACCELERATION
         || accAccumulator.x < -MAX_ACCELERATION || accAccumulator.y < -MAX_ACCELERATION || accAccumulator.z < -MAX_ACCELERATION) {
      consolePrintf("Max acceleration exceeded \n");
      if (accAccumulator.x > MAX_ACCELERATION) {
        accAccumulator.x = MAX_ACCELERATION;
      } else if (accAccumulator.x < -MAX_ACCELERATION) {
        accAccumulator.x = -MAX_ACCELERATION;
      }
      if (accAccumulator.y > MAX_ACCELERATION) {
        accAccumulator.y = MAX_ACCELERATION;
      } else if (accAccumulator.y < -MAX_ACCELERATION) {
        accAccumulator.y = -MAX_ACCELERATION;
      }
      if (accAccumulator.z > MAX_ACCELERATION) {
        accAccumulator.z = MAX_ACCELERATION;
      } else if (accAccumulator.z < -MAX_ACCELERATION) {
        accAccumulator.z = -MAX_ACCELERATION;
      }
    }
    */

    // CALL update
    float dt = (float)(osTick-lastPrediction)/configTICK_RATE_HZ;
    if (dt < (2.0f/PREDICT_RATE)) {
      busy = true;
      estimatorBoldermanUpdate(state, thrustAccumulator, &accAccumulator, &gyroAccumulator, &magAccumulator, dt, osTick);
      busy = false;
    } else {
      consolePrintf("Timestep to large \n");
    }

    // RESET all accumulator variables
    lastPrediction = osTick;
    thrustAccumulator = 0.0f;
    accAccumulator = (Axis3f){.axis={0}};
    gyroAccumulator = (Axis3f){.axis={0}};
    magAccumulator = (Axis3f){.axis={0}};
    thrustAccumulatorCount = 0;
    accAccumulatorCount = 0;
    gyroAccumulatorCount = 0;
    magAccumulatorCount = 0;
  }
}

/*
--------------------------------------------------------------------------------
*/
// ESTIMATORBOLDERMANUPDATE
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt, uint32_t osTick)
{
  // ENTER CODE HERE
  acceleration[0] = acc->x;
  acceleration[1] = acc->y;
  acceleration[2] = acc->z;
  omega[0] = gyro->x;
  omega[1] = gyro->y;
  omega[2] = gyro->z;

  // KALMAN ALGORITHM
  // PREDICT
  estimatorBoldermanPredict(dt);

  // UPDATE
  if (yCount >= 3) {
    // Combine dynamical model with measurements
    estimatorBoldermanDynMeas();
    // Reset variables
    //consolePrintf("Full rank distance measurement \n");
    yCount = 0;         // Reset counter
    y[0] = 0.0f;        // Reset distance measurements
    y[1] = 0.0f;
    y[2] = 0.0f;
    anchor[0][0] = 0.5f;    // Reset anchor positions DOES NOT MATTER WHAT THEY ARE, AS LONG AS THIS DIFFERS FROM THE REAL ANCHOR POSITIONS
    anchor[0][1] = 0.5f;
    anchor[0][2] = 0.5f;
    anchor[1][0] = 0.5f;
    anchor[1][1] = 0.5f;
    anchor[1][2] = 0.5f;
    anchor[2][0] = 0.5f;
    anchor[2][1] = 0.5f;
    anchor[2][2] = 0.5f;
  } else {
    if (isInit) {
      consolePrintf("No full rank distance measurement \n");
    }
    // Only use dynamical model, no full space measurements available
    estimatorBoldermanPredEst();
  }
  if (tracePxx > LIMIT_TRACE) {
    resetPxx();
  }

  // EXTERNALIZE
  estimatorBoldermanStateSave(state, osTick);
}

// PREDICT in X-variables
static void estimatorBoldermanPredict(float dt) {
  // DEFINE SIGMAPOINTS
  // Cholesky decomposition
  float delta[N][N] = {{0}};
  int status = 0;
  status = cholesky_decomposition(Pxx,delta,N);
  if (status == 0) {
    consolePrintf("Cholesky failed");
    // Then no good solution delta is found, so reinitialize Pxx
    resetPxx();
    memset(delta, 0, sizeof(delta));
    status = cholesky_decomposition(Pxx,delta,N);
    if (status == 0) {
      consolePrintf("Cholesky failes after resetting Pxx");
    }
  }
  // Fill in in sigmaX
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

  // Reset predicted x-value
  for (int ii=0; ii<N; ii++) {
    xpred[ii] = 0.0f;
  }
  // PROPOGATE THROUGH DYNAMICS
  for (int ii=0; ii<NSIGMA; ii++) {
    // Rotation
    //calculateRgb(sigmaX[6][ii], sigmaX[7][ii], sigmaX[8][ii]);
    //calculateWinv(sigmaX[6][ii], sigmaX[7][ii]);
    //accglob[0] = Rgb[0][0]*acceleration[0] + Rgb[1][0]*acceleration[1] + Rgb[2][0]*acceleration[2];
    //accglob[1] = Rgb[0][1]*acceleration[0] + Rgb[1][1]*acceleration[1] + Rgb[2][1]*acceleration[2];
    //accglob[2] = Rgb[0][2]*acceleration[0] + Rgb[1][2]*acceleration[1] + Rgb[2][2]*acceleration[2] - GRAVITY_MAGNITUDE;
    //omegaglob[0] = Winv[0][0]*omega[0] + Winv[0][1]*omega[1] + Winv[0][2]*omega[2];
    //omegaglob[1] = Winv[1][0]*omega[0] + Winv[1][1]*omega[1] + Winv[1][2]*omega[2];
    //omegaglob[2] = Winv[2][0]*omega[0] + Winv[2][1]*omega[1] + Winv[2][2]*omega[2];

    accglob[0] = cosf(sigmaX[6][ii])*cosf(sigmaX[8][ii])*acceleration[0] + (cosf(sigmaX[8][ii])*sinf(sigmaX[7][ii])*sinf(sigmaX[6][ii])-cosf(sigmaX[6][ii])*sinf(sigmaX[8][ii]))*acceleration[1] + (cosf(sigmaX[8][ii])*sinf(sigmaX[7][ii])*sinf(sigmaX[6][ii])+sinf(sigmaX[6][ii])*sinf(sigmaX[8][ii]))*acceleration[2];
    accglob[1] = cosf(sigmaX[6][ii])*sinf(sigmaX[8][ii])*acceleration[0] + (sinf(sigmaX[8][ii])*sinf(sigmaX[7][ii])*sinf(sigmaX[6][ii])+cosf(sigmaX[6][ii])*cosf(sigmaX[8][ii]))*acceleration[1] + (sinf(sigmaX[8][ii])*sinf(sigmaX[7][ii])*cosf(sigmaX[6][ii])-sinf(sigmaX[6][ii])*cosf(sigmaX[8][ii]))*acceleration[2];
    accglob[2] = -sinf(sigmaX[7][ii])*acceleration[0] + cosf(sigmaX[7][ii])*sinf(sigmaX[6][ii])*acceleration[1] + cosf(sigmaX[7][ii])*cosf(sigmaX[6][ii])*acceleration[2] - GRAVITY_MAGNITUDE;
    if (cosf(sigmaX[7][ii])>LIMIT_INVERT || cosf(sigmaX[7][ii])<LIMIT_INVERT) {
      omegaglob[0] = omega[0] + (sinf(sigmaX[7][ii])*sinf(sigmaX[6][ii])/cosf(sigmaX[7][ii]))*omega[1] + (sinf(sigmaX[7][ii])*cosf(sigmaX[6][ii])/cosf(sigmaX[7][ii]))*omega[2];
      omegaglob[1] = (cosf(sigmaX[7][ii])*cosf(sigmaX[6][ii])/cosf(sigmaX[7][ii]))*omega[1] - (cosf(sigmaX[7][ii])*sinf(sigmaX[6][ii])/cosf(sigmaX[7][ii]))*omega[2];
      omegaglob[2] = (sinf(sigmaX[6][ii])/cosf(sigmaX[7][ii]))*omega[1] + (cosf(sigmaX[6][ii])/cosf(sigmaX[7][ii]))*omega[2];
    } else {
      consolePrintf("Winv not calculated, under limit");
      omegaglob[0] = 0.0f;
      omegaglob[1] = 0.0f;
      omegaglob[2] = 0.0f;
    }
    // Dynamics
    sigmaXplus[0][ii] = sigmaX[0][ii] + dt*sigmaX[3][ii] + 0.5f*dt*dt*accglob[0];
    sigmaXplus[1][ii] = sigmaX[1][ii] + dt*sigmaX[4][ii] + 0.5f*dt*dt*accglob[1];
    sigmaXplus[2][ii] = sigmaX[2][ii] + dt*sigmaX[5][ii] + 0.5f*dt*dt*accglob[2];
    sigmaXplus[3][ii] = sigmaX[3][ii] + dt*accglob[0];
    sigmaXplus[4][ii] = sigmaX[4][ii] + dt*accglob[1];
    sigmaXplus[5][ii] = sigmaX[5][ii] + dt*accglob[2];
    sigmaXplus[6][ii] = sigmaX[6][ii] + dt*omegaglob[0];
    sigmaXplus[7][ii] = sigmaX[7][ii] + dt*omegaglob[1];
    sigmaXplus[8][ii] = sigmaX[8][ii] + dt*omegaglob[2];
    // Predict
    xpred[0] += Wm[ii]*sigmaXplus[0][ii];
    xpred[1] += Wm[ii]*sigmaXplus[1][ii];
    xpred[2] += Wm[ii]*sigmaXplus[2][ii];
    xpred[3] += Wm[ii]*sigmaXplus[3][ii];
    xpred[4] += Wm[ii]*sigmaXplus[4][ii];
    xpred[5] += Wm[ii]*sigmaXplus[5][ii];
    xpred[6] += Wm[ii]*sigmaXplus[6][ii];
    xpred[7] += Wm[ii]*sigmaXplus[7][ii];
    xpred[8] += Wm[ii]*sigmaXplus[8][ii];
  }

  // COVARIANCE
  tracePxx = 0.0f;
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      Pxx[ii][jj] = q[ii][jj];
      for (int kk=0; kk<NSIGMA; kk++) {
        Pxx[ii][jj] += Wc[kk] * (sigmaXplus[ii][kk] - xpred[ii]) * (sigmaXplus[jj][kk] - xpred[jj]);
      }
      if (ii == jj) {
        tracePxx += Pxx[ii][jj];
      }
    }
  }
}

// UPDATE
// WITH sufficient measurements
static void estimatorBoldermanDynMeas(void) {
  // Reset ypred
  for (int ii=0; ii<NOUT; ii++) {
    ypred[ii] = 0.0f;
  }
  // REDIFINE POSITION IN DISTANCE -> measurement equations
  for (int ii=0; ii<NSIGMA; ii++) {
    sigmaYplus[0][ii] = sqrtf((sigmaXplus[0][ii]-anchor[0][0])*(sigmaXplus[0][ii]-anchor[0][0]) + (sigmaXplus[1][ii]-anchor[1][0])*(sigmaXplus[1][ii]-anchor[1][0]) + (sigmaXplus[2][ii]-anchor[2][0])*(sigmaXplus[2][ii]-anchor[2][0]));
    sigmaYplus[1][ii] = sqrtf((sigmaXplus[0][ii]-anchor[0][1])*(sigmaXplus[0][ii]-anchor[0][1]) + (sigmaXplus[1][ii]-anchor[1][1])*(sigmaXplus[1][ii]-anchor[1][1]) + (sigmaXplus[2][ii]-anchor[2][1])*(sigmaXplus[2][ii]-anchor[2][1]));
    sigmaYplus[2][ii] = sqrtf((sigmaXplus[0][ii]-anchor[0][2])*(sigmaXplus[0][ii]-anchor[0][2]) + (sigmaXplus[1][ii]-anchor[1][2])*(sigmaXplus[1][ii]-anchor[1][2]) + (sigmaXplus[2][ii]-anchor[2][2])*(sigmaXplus[2][ii]-anchor[2][2]));
    ypred[0] += Wm[ii] * sigmaYplus[0][ii];
    ypred[1] += Wm[ii] * sigmaYplus[1][ii];
    ypred[2] += Wm[ii] * sigmaYplus[2][ii];
  }
  // OUTPUT COVARIANCE -> measurement equations covariance
  for (int ii=0; ii<NOUT; ii++) {
    for (int jj=0; jj<NOUT; jj++) {
      Pyy[ii][jj] = r[ii][jj];
      for (int kk=0; kk<NSIGMA; kk++) {
        Pyy[ii][jj] += Wc[kk] * (sigmaYplus[ii][kk]-ypred[ii]) * (sigmaYplus[jj][kk]-ypred[jj]);
      }
    }
  }
  // CROSS COVARIANCE -> state and output covariance
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<NOUT; jj++) {
      Pxy[ii][jj] = 0.0f;
      for (int kk=0; kk<NSIGMA; kk++) {
        Pxy[ii][jj] += Wc[kk] * (sigmaXplus[ii][kk]-xpred[ii]) * (sigmaYplus[jj][kk]-ypred[jj]);
      }
    }
  }
  // INVERSE of OUTPUT COVARIANCE -> should be invertible
  // CRAMER's rule for computing inverse
  /*
  determinantPyy = Pyy[0][0]*(Pyy[1][1]*Pyy[2][2]-Pyy[1][2]*Pyy[2][1]) - Pyy[0][1]*(Pyy[1][0]*Pyy[2][2]-Pyy[1][2]*Pyy[2][0]) + Pyy[0][2]*(Pyy[1][0]*Pyy[2][1]-Pyy[2][0]*Pyy[1][1]);
  if (determinantPyy>LIMIT_INVERT || determinantPyy<LIMIT_INVERT) {
    Pyyinv[0][0] = (Pyy[1][1]*Pyy[2][2]-Pyy[1][2]*Pyy[2][1])/determinantPyy;
    Pyyinv[0][1] = (Pyy[0][2]*Pyy[2][1]-Pyy[0][1]*Pyy[2][2])/determinantPyy;
    Pyyinv[0][2] = (Pyy[0][1]*Pyy[1][2]-Pyy[1][1]*Pyy[0][2])/determinantPyy;
    Pyyinv[1][0] = (Pyy[2][0]*Pyy[1][2]-Pyy[1][0]*Pyy[2][2])/determinantPyy;
    Pyyinv[1][1] = (Pyy[0][0]*Pyy[2][2]-Pyy[2][0]*Pyy[0][2])/determinantPyy;
    Pyyinv[1][2] = (Pyy[1][0]*Pyy[0][2]-Pyy[0][0]*Pyy[1][2])/determinantPyy;
    Pyyinv[2][0] = (Pyy[1][0]*Pyy[2][1]-Pyy[1][1]*Pyy[2][0])/determinantPyy;
    Pyyinv[2][1] = (Pyy[2][0]*Pyy[0][1]-Pyy[0][0]*Pyy[2][1])/determinantPyy;
    Pyyinv[2][2] = (Pyy[0][0]*Pyy[1][1]-Pyy[1][0]*Pyy[0][1])/determinantPyy;
  } else {
    consolePrintf("Determinant Pyy too small!");
    estimatorBoldermanPredEst();
    resetPxx();
  }
  */
  // GAUSSIAN ELIMINATION
  for (int ii=0; ii<NOUT; ii++) {
    for (int jj=0; jj<NOUT; jj++) {
      Pyyhelp[ii][jj] = Pyy[ii][jj];
    }
  }
  static arm_matrix_instance_f32 Pyym = {NOUT, NOUT, (float *)Pyyhelp};
  static arm_matrix_instance_f32 Pyyinvm = {NOUT, NOUT, (float *)Pyyinv};
  mat_inv(&Pyym, &Pyyinvm);

  // KALMAN GAIN -> Kalman gain is given as K=Pxy*inv(Pyy)
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<NOUT; jj++) {
      K[ii][jj] = 0.0f;
      for (int kk=0; kk<NOUT; kk++) {
        K[ii][jj] += Pxy[ii][kk]*Pyyinv[kk][jj];
      }
    }
  }
  for (int ii=0; ii<NOUT; ii++) {
    for (int jj=0; jj<NOUT; jj++) {
      Pyyinv[ii][jj] = 0.0f;
      Pyy[ii][jj] = 0.0f;
    }
  }
  // COVARIANCE UPDATE -> Pxxnew = Pxxold - K*Pyy*trans(K) = Pxxold - Pxy*inv(Pyy)*Pyy*trans(K) = Pxxold - Pxy*trans(K)
  // deltaPxx = Pxy*trans(K)
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      deltaPxx[ii][jj] = 0.0f;
      for (int kk=0; kk<NOUT; kk++) {
        deltaPxx[ii][jj] += Pxy[ii][kk]*K[jj][kk];
      }
    }
  }
  // PERFORM THE REAL UPDATE -> x = xpred + K(y-ypred) && Pxx = Pxxold - deltaPxx
  for (int ii=0; ii<N; ii++) {
    x[ii] = xpred[ii];
    for (int jj=0; jj<NOUT; jj++) {
      x[ii] += K[ii][jj] * (y[jj]-ypred[jj]);
    }
  }
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      Pxx[ii][jj] -= deltaPxx[ii][jj];
    }
  }
  /*
  for (int ii=0; ii<NOUT; ii++) {
    for (int jj=0; jj<N; jj++) {
      Pxy[jj][ii] = 0.0f;
      K[jj][ii] = 0.0f;
    }
  }
  */
  // ENSURE the angles are given;
  safetyAngleBounds();
}
// WITHOUT sufficient distance measurements
static void estimatorBoldermanPredEst(void) {
  for (int ii=0; ii<N; ii++) {
    x[ii] = xpred[ii];
  }
}

// EXTERNALIZE state
static void estimatorBoldermanStateSave(state_t *state, uint32_t osTick) {
  // POSITION -> GLOBAL FRAME
  state->position = (point_t) {
    .timestamp = osTick,
    .x = x[0],
    .y = x[1],
    .z = x[2]
  };
  // VELOCITY -> GLOBAL FRAME
  state->velocity = (velocity_t) {
    .timestamp = osTick,
    .x = x[3],
    .y = x[4],
    .z = x[5]
  };
  // ACCELERATION -> GLOBAL FRAME IN UNIT G
  //calculateRgb(x[6], x[7], x[8]);
  //accext[0] = Rgb[0][0]*acceleration[0] + Rgb[1][0]*acceleration[1] + Rgb[2][0]*acceleration[2];
  //accext[1] = Rgb[0][1]*acceleration[0] + Rgb[1][1]*acceleration[1] + Rgb[2][1]*acceleration[2];
  //accext[2] = Rgb[0][2]*acceleration[0] + Rgb[1][2]*acceleration[1] + Rgb[2][2]*acceleration[2] - GRAVITY_MAGNITUDE;
  accext[0] = cosf(x[6])*cosf(x[8])*acceleration[0] + (cosf(x[8])*sinf(x[7])*sinf(x[6])-cosf(x[6])*sinf(x[8]))*acceleration[1] + (cosf(x[8])*sinf(x[7])*sinf(x[6])+sinf(x[6])*sinf(x[8]))*acceleration[2];
  accext[1] = cosf(x[6])*sinf(x[8])*acceleration[0] + (sinf(x[8])*sinf(x[7])*sinf(x[6])+cosf(x[6])*cosf(x[8]))*acceleration[1] + (sinf(x[8])*sinf(x[7])*cosf(x[6])-sinf(x[6])*cosf(x[8]))*acceleration[2];
  accext[2] = -sinf(x[7])*acceleration[0] + cosf(x[7])*sinf(x[6])*acceleration[1] + cosf(x[7])*cosf(x[6])*acceleration[2] - GRAVITY_MAGNITUDE;

  state->acc = (acc_t) {
    .timestamp = osTick,
    .x = ((accext[0])/GRAVITY_MAGNITUDE),
    .y = ((accext[1])/GRAVITY_MAGNITUDE),
    .z = ((accext[2])/GRAVITY_MAGNITUDE),
  };
  // ATTITUDE -> GLOBAL FRAME, PITCH INVERTED!
  state->attitude = (attitude_t) {
    .timestamp = osTick,
    .roll = (x[6] * RAD_TO_DEG),
    .pitch = -1.0f*(x[7] * RAD_TO_DEG),
    .yaw = (x[8] * RAD_TO_DEG)
  };
}

/*
bool estimatorBoldermanEnqueueDistance(distanceMeasurement_t *measurement) {
  if (isInit) {
    uint32_t yindex = yCount % 3;
    y[yindex] = measurement->distance;
    anchor[0][yindex] = measurement->x;
    anchor[1][yindex] = measurement->y;
    anchor[2][yindex] = measurement->z;
    yCount++;
  }
  return (pdTRUE);
}
*/
bool estimatorBoldermanEnqueueDistance(distanceMeasurement_t *measurement) {
  // Only do this after the initialization
  if (busy) {
    return (pdTRUE);
  } else {
    if (isInit) {
      // Only use measurement when it is within a reasonable distance
      if (measurement->distance < MAX_DISTANCE && measurement->distance > MIN_DISTANCE) {
        // We have the measurement here, use it accordingly:
        // - if the anchor position is not there yet, put it in
        // - if the anchor position is there, overwrite the distance (more recent measurements are used)
        for (int ii=0; ii<3; ii++) {
          if (measurement->x == anchor[0][ii]
            && measurement->y == anchor[1][ii]
            && measurement->z == anchor[2][ii]) {
              //consolePrintf("Second measurement from same anchor before cleaning \n");
              y[ii] = measurement->distance;
              return (pdTRUE);
          }
        }
        uint32_t yindex = yCount % 3;
        y[yindex] = measurement->distance;
        anchor[0][yindex] = measurement->x;
        anchor[1][yindex] = measurement->y;
        anchor[2][yindex] = measurement->z;
        yCount++;
        return (pdTRUE);
      }
    }
  }
  return (pdTRUE);
}

/*
--------------------------------------------------------------------------------
*/
// INITIALIZATION
// TEST for initialization
bool estimatorBoldermanTest(void)
{
  return isInit;
}
// INITIALIZATION performed
void estimatorBoldermanInit(void)
{
  consolePrintf("Initialization start \n");

  // INITIALIZE all variables
  lastPrediction = xTaskGetTickCount();
  thrustAccumulator = 0.0f;
  accAccumulator = (Axis3f){.axis={0}};
  gyroAccumulator = (Axis3f){.axis={0}};
  magAccumulator = (Axis3f){.axis={0}};
  thrustAccumulatorCount = 0;
  accAccumulatorCount = 0;
  gyroAccumulatorCount = 0;
  magAccumulatorCount = 0;

  // Define the weights
  lambda = alpha*alpha*(NCALC+kappa)-NCALC;
  Wm[0] = lambda/(NCALC+lambda);
  Wc[0] = lambda/(NCALC+lambda) + (1.0f-alpha*alpha+beta);
  for (int ii=1; ii<NSIGMA; ii++) {
    Wm[ii] = 1.0f/(2.0f*(NCALC+lambda));
    Wc[ii] = 1.0f/(2.0f*(NCALC+lambda));
  }
  // Define covariances
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      if (ii == jj) {
        Pxx[ii][jj] = Pxxdiag[ii]*Pxxdiag[ii];
        q[ii][jj] = qdiag[ii]*qdiag[ii];
      } else {
        Pxx[ii][jj] = 0.0f;
        q[ii][jj] = 0.0f;
      }
    }
  }
  for (int ii=0; ii<NOUT; ii++) {
    for (int jj=0; jj<NOUT; jj++) {
      if (ii == jj) {
        r[ii][jj] = rdiag[ii]*rdiag[ii];
      } else {
        r[ii][jj] = 0.0f;
      }
    }
  }
  // isINIT lets system know initialization is done
  isInit = true;
  consolePrintf("Initialization finished \n");
}


/*
--------------------------------------------------------------------------------
*/
// CACLULATION functions
static void resetPxx(void) {
  consolePrintf("Reset Pxx \n");
  for (int ii=0; ii<N; ii++) {
    for (int jj=0; jj<N; jj++) {
      if (ii == jj) {
        Pxx[ii][jj] = Pxxdiag[ii]*Pxxdiag[ii];
      } else {
        Pxx[ii][jj] = 0.0f;
      }
    }
  }
}
/*
// Rotation matrix
static void calculateRgb(float phi, float theta, float psi) {
  Rgb[0][0] = cosf(phi)*cosf(psi);
  Rgb[0][1] = cosf(phi)*sinf(psi);
  Rgb[0][2] = -sinf(theta);
  Rgb[1][0] = cosf(psi)*sinf(theta)*sinf(phi)-cosf(phi)*sinf(psi);
  Rgb[1][1] = sinf(psi)*sinf(theta)*sinf(phi)+cosf(phi)*cosf(psi);
  Rgb[1][2] = cosf(theta)*sinf(phi);
  Rgb[2][0] = cosf(psi)*sinf(theta)*sinf(phi)+sinf(phi)*sinf(psi);
  Rgb[2][1] = sinf(psi)*sinf(theta)*cosf(phi)-sinf(phi)*cosf(psi);
  Rgb[2][2] = cosf(theta)*cosf(phi);
}
// Angular rate matrix
static void calculateWinv(float phi, float theta) {
  // ONLY calculate when it is numerically possible to invert
  if (cosf(theta)>LIMIT_INVERT || cosf(theta)<-LIMIT_INVERT) {
    Winv[0][0] = 1.0f;
    Winv[0][1] = sinf(theta)*sinf(phi)/cosf(theta);
    Winv[0][2] = sinf(theta)*cosf(phi)/cosf(theta);
    Winv[1][0] = 0.0f;
    Winv[1][1] = cosf(theta)*cosf(phi)/cosf(theta);
    Winv[1][2] = -cosf(theta)*sinf(phi)/cosf(theta);
    Winv[2][0] = 0.0f;
    Winv[2][1] = sinf(phi)/cosf(theta);
    Winv[2][2] = cosf(phi)/cosf(theta);
  } else {
    Winv[0][0] = 0.0f;
    Winv[0][1] = 0.0f;
    Winv[0][2] = 0.0f;
    Winv[1][0] = 0.0f;
    Winv[1][1] = 0.0f;
    Winv[1][2] = 0.0f;
    Winv[2][0] = 0.0f;
    Winv[2][1] = 0.0f;
    Winv[2][2] = 0.0f;
  }
}
*/


// MIRROR the angles within the allowed region [-PI, +PI]
static void safetyAngleBounds(void) {
  while (x[6]>PI || x[6]<-PI) {
    if (x[6] > PI) {
      x[6] -= 2*PI;
    } else if (x[6] < -PI) {
      x[6] += 2*PI;
    }
  }
  while (x[7]>PI || x[7]<-PI) {
    if (x[7] > PI) {
      x[7] -= 2*PI;
    } else if (x[7] < -PI) {
      x[7] += 2*PI;
    }
  }
  while (x[8]>PI || x[8]<-PI) {
    if (x[8] > PI) {
      x[8] -= 2*PI;
    } else if (x[8] < -PI) {
      x[8] += 2*PI;
    }
  }
}


// CHOLESKY DECOMPOSITION
// WRITTEN by; MARCUS GREIFF
int assert_element(float val){
  /***********************************************************************
   * Here we check that the diagonal element of R is greater than or
   * or equal to some limit, ensuring that the decomposition is done
   * correctly and allowing detection of non-PSD matrices A
   **********************************************************************/
  float limit = LIMIT_CHOLESKY;
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
int cholesky_decomposition(float (*A)[N], float (*R)[N], int n) {
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

/*
--------------------------------------------------------------------------------
END MODEL
--------------------------------------------------------------------------------
*/


// LOG GROUPS for every variable you want to keep track of
LOG_GROUP_START(BOLDERMAN_state)
  LOG_ADD(LOG_FLOAT, tracePxx, &tracePxx)
  LOG_ADD(LOG_FLOAT, x_est, &x[0])
  LOG_ADD(LOG_FLOAT, y_est, &x[1])
  LOG_ADD(LOG_FLOAT, z_est, &x[2])
  LOG_ADD(LOG_FLOAT, xd_est, &x[3])
  LOG_ADD(LOG_FLOAT, yd_est, &x[4])
  LOG_ADD(LOG_FLOAT, zd_est, &x[5])
  LOG_ADD(LOG_FLOAT, roll_est, &x[6])
  LOG_ADD(LOG_FLOAT, pitch_est, &x[7])
  LOG_ADD(LOG_FLOAT, yaw_est, &x[8])
LOG_GROUP_STOP(BOLDERMAN_state)
LOG_GROUP_START(BOLDERMAN_meas)
  LOG_ADD(LOG_FLOAT, covx, &Pxx[0][0])
  LOG_ADD(LOG_FLOAT, covy, &Pxx[1][1])
  LOG_ADD(LOG_FLOAT, covz, &Pxx[2][2])
  LOG_ADD(LOG_FLOAT, covxd, &Pxx[3][3])
  LOG_ADD(LOG_FLOAT, covyd, &Pxx[4][4])
  LOG_ADD(LOG_FLOAT, covzd, &Pxx[5][5])
  LOG_ADD(LOG_FLOAT, covroll, &Pxx[6][6])
  LOG_ADD(LOG_FLOAT, covpitch, &Pxx[7][7])
  LOG_ADD(LOG_FLOAT, covyaw, &Pxx[8][8])
LOG_GROUP_STOP(BOLDERMAN_meas)
