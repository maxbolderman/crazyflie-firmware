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

// VARIABLES standard
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

// PROTOTYPES of functions used in estimator
static void estimatorBoldermanUpdate(state_t *state, float thrust, Axis3f *acc, Axis3f *gyro, Axis3f *mag, float dt, uint32_t osTick);

/*
--------------------------------------------------------------------------------
START MODEL
--------------------------------------------------------------------------------
*/
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

    // CALL update
    float dt = (float)(osTick-lastPrediction)/configTICK_RATE_HZ;
    estimatorBoldermanUpdate(state, thrustAccumulator, &accAccumulator, &gyroAccumulator, &magAccumulator, dt, osTick);

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
  //estimatorBoldermanPredict(dt);

  // UPDATE
  if (yCount >= 3) {
    // Combine dynamical model with measurements
    //estimatorBoldermanDynMeas();
    // Reset variables
    consolePrintf("Full rank distance measurement \n");
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
    consolePrintf("No full rank distance measurement \n");
    // Only use dynamical model, no full space measurements available
    //estimatorBoldermanPredEst();
  }

  // EXTERNALIZE
  //estimatorBoldermanStateSave(state, osTick);
}


// TEST for initialization
bool estimatorBoldermanTest(void)
{
  return isInit;
}
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

  // INIT let system know initialization is done
  isInit = true;
  consolePrintf("Initialization finished \n");
}



bool estimatorBoldermanEnqueueDistance(distanceMeasurement_t *measurement)
{
  // Only do this after the initialization
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
            consolePrintf("Second measurement from same anchor before cleaning \n");
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
  return (pdTRUE);
}

/*
--------------------------------------------------------------------------------
END MODEL
--------------------------------------------------------------------------------
*/


// LOG GROUPS for every variable you want to keep track of
LOG_GROUP_START(BOLDERMAN_meas)
  LOG_ADD(LOG_FLOAT, acc_x, &acceleration[0])
  LOG_ADD(LOG_FLOAT, acc_y, &acceleration[1])
  LOG_ADD(LOG_FLOAT, acc_z, &acceleration[2])
  LOG_ADD(LOG_FLOAT, gyro_x, &omega[0])
  LOG_ADD(LOG_FLOAT, gyro_y, &omega[1])
  LOG_ADD(LOG_FLOAT, gyro_z, &omega[2])
  LOG_ADD(LOG_FLOAT, dist_1, &y[0])
  LOG_ADD(LOG_FLOAT, dist_2, &y[1])
  LOG_ADD(LOG_FLOAT, dist_3, &y[2])
LOG_GROUP_STOP(BOLDERMAN_meas)
