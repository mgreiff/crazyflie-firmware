/**
 *    ||          ____  _ __
 * +------+      / __ )(_) /_______________ _____  ___
 * | 0xBC |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * +------+    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *  ||  ||    /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie Firmware
 *
 * Copyright (C) 2011-2016 Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */
#include <math.h>

#include "FreeRTOS.h"
#include "task.h"

#include "system.h"
#include "log.h"
#include "param.h"

#include "stabilizer.h"

#include "sensors.h"
#include "ext_position.h"
#include "sitaw.h"

#define SETPOINT_TYPE_M
#define CONTROL_TYPE_M
#define POWER_DISTRIBUTION_TYPE_M

// Setpoint
#include "M_sequenceCommander.h"  // New commander
#include "M_flatnessGenerator.h"  // Flatness generator
//#include "commander.h"            // Old commander

// Controller
#include "M_control.h"            // New
#include "controller.h"           // Old

// Power distribution
#include "M_powerDistribution.h"  // New
#include "power_distribution.h"   // Old

// State estimator
#ifdef ESTIMATOR_TYPE_kalman
#include "estimator_kalman.h"
#else
#include "estimator.h"
#endif

static bool isInit;

// State variables for the stabilizer
static setpoint_t setpoint;
static sensorData_t sensorData;
static state_t state;
static control_t control;

// State variables for the M-control scheme.
static M_setpoint_t M_setpoint;
static M_mode_t M_mode;
static M_state_t M_state;
static M_control_t M_controlsignal;

static void stabilizerTask(void* param);

void stabilizerInit(void)
{
  if(isInit)
    return;

  sensorsInit();
  stateEstimatorInit();

#if defined(SETPOINT_TYPE_M)
  commanderInit_M();
#endif

#if defined(CONTROL_TYPE_M)
  stateControllerInit_M();
#else
  stateControllerInit();
#endif

#if defined(POWER_DISTRIBUTION_TYPE_M)
  powerDistributionInit_M();
#else
  powerDistributionInit();
#endif

#if defined(SITAW_ENABLED)
  sitAwInit();
#endif

  xTaskCreate(stabilizerTask, STABILIZER_TASK_NAME,
              STABILIZER_TASK_STACKSIZE, NULL, STABILIZER_TASK_PRI, NULL);

  isInit = true;
}

bool stabilizerTest(void)
{
  bool pass = true;

  pass &= sensorsTest();
  pass &= stateEstimatorTest();
#ifndef CONTROL_TYPE_M
  pass &= stateControllerTest();
#endif

#ifndef POWER_DISTRIBUTION_TYPE_M
  pass &= powerDistributionTest();
#endif

  return pass;
}

/* The stabilizer loop runs at 1kHz (stock) or 500Hz (kalman). It is the
 * responsibility of the different functions to run slower by skipping call
 * (ie. returning without modifying the output structure).
 */

static void stabilizerTask(void* param)
{
  uint32_t tick = 0;
  uint32_t lastWakeTime;
  vTaskSetApplicationTaskTag(0, (void*)TASK_STABILIZER_ID_NBR);

  //Wait for the system to be fully started to start stabilization loop
  systemWaitStart();

  // Wait for sensors to be calibrated
  lastWakeTime = xTaskGetTickCount ();
  while(!sensorsAreCalibrated()) {
    vTaskDelayUntil(&lastWakeTime, F2T(RATE_MAIN_LOOP));
  }

  while(1) {
    vTaskDelayUntil(&lastWakeTime, F2T(RATE_MAIN_LOOP));
    
    // State estimation
#if defined(ESTIMATOR_TYPE_kalman)
    stateEstimatorUpdate(&M_state, &sensorData, &M_controlsignal);
#else
    sensorsAcquire(&sensorData, tick);
    stateEstimator(&state, &sensorData, tick);
#endif

    // ~~~ Commander and references ~~~
#if defined(SETPOINT_TYPE_M)
    commanderGetSetpoint_M(&M_setpoint, &M_mode, &M_state, tick);
    flatnessExpansion_M(&M_setpoint, &M_mode);
#else
    getExtPosition(&state);
    commanderGetSetpoint(&setpoint, &state);
#endif

    // ~~~ Controller and situational awareness ~~~
#if defined(CONTROL_TYPE_M)
    if (M_setpoint.position.z > 0.1f) {
      M_mode.controller.current = 4;
      M_mode.flatness.eulerAngles = true;
      M_mode.flatness.quaternion = true;
      M_mode.flatness.any = true;
      M_mode.flatness.bodyacceleration = false;
      M_mode.flatness.bodyrate = true;
      M_mode.flatness.torque = false;
    } else {
      M_mode.controller.current = 0;
      M_mode.flatness.eulerAngles = false;
      M_mode.flatness.quaternion = false;
      M_mode.flatness.any = false;
      M_mode.flatness.bodyacceleration = false;
      M_mode.flatness.bodyrate = false;
      M_mode.flatness.torque = false;
    }
    stateController_M(&M_controlsignal, &M_setpoint, &M_mode, &M_state, tick);
#else
    setpoint.position = M_setpoint.position;
    setpoint.velocity = M_setpoint.position;
    setpoint.attitude.roll = M_setpoint.eulerAngles.pitch;
    setpoint.attitude.pitch = -M_setpoint.eulerAngles.roll;
    setpoint.attitude.yaw = M_setpoint.eulerAngles.yaw;
    
    if (setpoint.position.z > 0.3f){
      setpoint.mode.x = modeAbs;
      setpoint.mode.y = modeAbs;
      setpoint.mode.z = modeAbs;
      setpoint.mode.roll = modeDisable;
      setpoint.mode.pitch = modeDisable;
      setpoint.mode.yaw = modeAbs;
    }else{
      setpoint.mode.x = modeDisable;
      setpoint.mode.y = modeDisable;
      setpoint.mode.z = modeDisable;
      setpoint.mode.roll = modeDisable;
      setpoint.mode.pitch = modeDisable;
      setpoint.mode.yaw = modeDisable;
    }
    sitAwUpdateSetpoint(&setpoint, &sensorData, &state);
    stateController(&control, &sensorData, &state, &setpoint, tick);
#endif


    // ~~~ Power distribution ~~~
#if defined(POWER_DISTRIBUTION_TYPE_M)
    powerDistribution_M(&M_mode, &M_controlsignal);
#else
    powerDistribution(&control);
#endif
    tick++;
  }
}

// Marcus log groups (required for compatibility with the ROS applicaton)
LOG_GROUP_START(reference)
LOG_ADD(LOG_FLOAT, x, &M_setpoint.position.x)
LOG_ADD(LOG_FLOAT, y, &M_setpoint.position.y)
LOG_ADD(LOG_FLOAT, z, &M_setpoint.position.z)
LOG_ADD(LOG_FLOAT, roll, &M_setpoint.eulerAngles.roll)
LOG_ADD(LOG_FLOAT, pitch, &M_setpoint.eulerAngles.pitch)
LOG_ADD(LOG_FLOAT, yaw, &M_setpoint.eulerAngles.yaw)
LOG_GROUP_STOP(reference)

LOG_GROUP_START(measured)
LOG_ADD(LOG_FLOAT, x, &M_state.position.x)
LOG_ADD(LOG_FLOAT, y, &M_state.position.y)
LOG_ADD(LOG_FLOAT, z, &M_state.position.z)
LOG_ADD(LOG_FLOAT, roll, &M_state.eulerAngles.roll)
LOG_ADD(LOG_FLOAT, pitch, &M_state.eulerAngles.pitch)
LOG_ADD(LOG_FLOAT, yaw, &M_state.eulerAngles.yaw)
LOG_GROUP_STOP(measured)

LOG_GROUP_START(acc)
LOG_ADD(LOG_FLOAT, x, &sensorData.acc.x)
LOG_ADD(LOG_FLOAT, y, &sensorData.acc.y)
LOG_ADD(LOG_FLOAT, z, &sensorData.acc.z)
LOG_GROUP_STOP(acc)

LOG_GROUP_START(baro)
LOG_ADD(LOG_FLOAT, asl, &sensorData.baro.asl)
LOG_ADD(LOG_FLOAT, temp, &sensorData.baro.temperature)
LOG_ADD(LOG_FLOAT, pressure, &sensorData.baro.pressure)
LOG_GROUP_STOP(baro)

LOG_GROUP_START(gyro)
LOG_ADD(LOG_FLOAT, x, &sensorData.gyro.x)
LOG_ADD(LOG_FLOAT, y, &sensorData.gyro.y)
LOG_ADD(LOG_FLOAT, z, &sensorData.gyro.z)
LOG_GROUP_STOP(gyro)

LOG_GROUP_START(mag)
LOG_ADD(LOG_FLOAT, x, &sensorData.mag.x)
LOG_ADD(LOG_FLOAT, y, &sensorData.mag.y)
LOG_ADD(LOG_FLOAT, z, &sensorData.mag.z)
LOG_GROUP_STOP(mag)

// Parameters for mode switching
PARAM_GROUP_START(controller)
PARAM_ADD(PARAM_UINT8, mode, &M_mode.controller.current)
PARAM_GROUP_STOP(controller)
