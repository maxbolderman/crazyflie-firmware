#ifndef __ESTIMATOR_BOLDERMAN_H__
#define __ESTIMATOR_BOLDERMAN_H__

#include <stdint.h>
#include "stabilizer_types.h"

// Defines prototypes for the functions, which are further defined in .c file
void estimatorBoldermanInit(void);
bool estimatorBoldermanTest(void);
void estimatorBolderman(state_t *state, sensorData_t *sensors, control_t *control, const uint32_t tick);
bool estimatorBoldermanEnqueueDistance(distanceMeasurement_t *measurement);
//bool can_enqueue_distance_measurement(void);
//void has_enqueued_distance_measurement(void);

#endif // __ESTIMATOR_BOLDERMAN_H__
