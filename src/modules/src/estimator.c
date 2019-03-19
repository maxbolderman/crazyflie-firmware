#define DEBUG_MODULE "ESTIMATOR"
#include "debug.h"

#include "cfassert.h"
#include "estimator.h"
#include "estimator_complementary.h"
#include "estimator_kalman.h"
#include "estimator_bolderman.h"

//#define DEFAULT_ESTIMATOR boldermanEstimator
#define DEFAULT_ESTIMATOR boldermanEstimator
static StateEstimatorType currentEstimator = anyEstimator;

// prototype of initEstimator, used to initialize estimator, for corresponding estimator
static void initEstimator();

// Defines the different functions for an estimator (3 possibilities)
typedef struct {
  void (*init)(void);
  bool (*test)(void);
  void (*update)(state_t *state, sensorData_t *sensors, control_t *control, const uint32_t tick);
} EstimatorFcns;

// Defines the estimator functions depending on the estimator used
static EstimatorFcns estimatorFunctions[] = {
  {.init = 0, .test = 0, .update = 0}, // Any
  {.init = estimatorComplementaryInit, .test = estimatorComplementaryTest, .update = estimatorComplementary},
  {.init = estimatorKalmanInit, .test = estimatorKalmanTest, .update = estimatorKalman},
  {.init = estimatorBoldermanInit, .test = estimatorBoldermanTest, .update = estimatorBolderman},
};

// Determines the estimator used for this purpose
void stateEstimatorInit(StateEstimatorType estimator) {
  if (estimator < 0 || estimator >= StateEstimatorTypeCount) {
    return;
  }

  currentEstimator = estimator;

  if (anyEstimator == currentEstimator) {
    currentEstimator = DEFAULT_ESTIMATOR;
  }

  StateEstimatorType forcedEstimator = ESTIMATOR_NAME;
  if (forcedEstimator != anyEstimator) {
    DEBUG_PRINT("Estimator type forced\n");
    currentEstimator = forcedEstimator;
  }

  initEstimator();

  DEBUG_PRINT("Using estimator %d\n", currentEstimator);
}

// Returns the estimator used (for this purpose boldermanEstimator)
StateEstimatorType getStateEstimator(void) {
  return currentEstimator;
}

/* Three functions below are different for the estimator used.
 What these functions are and what they do is defined in their own files
 Note that different fucntions are called depending on currentEstimator */
// Initializes the estimator (calls estimatorBoldermanInit() )
static void initEstimator() {
  estimatorFunctions[currentEstimator].init();
}

// Performs a test, returns boolean (calls estimatorBoldermanTest() )
bool stateEstimatorTest(void) {
  return estimatorFunctions[currentEstimator].test();
}

// Updates (calls estimatorBolderman() for update)
void stateEstimator(state_t *state, sensorData_t *sensors, control_t *control, const uint32_t tick) {
  estimatorFunctions[currentEstimator].update(state, sensors, control, tick);
}
