/*
 * Copyright (C) 2014 Freek van Tienen <freek.v.tienen@gmail.com>
 *
 * This file is part of paparazzi.
 *
 * paparazzi is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * paparazzi is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with paparazzi; see the file COPYING.  If not, write to
 * the Free Software Foundation, 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

/** @file modules/loggers/file_logger.c
 *  @brief File logger for Linux based autopilots
 */

#include "file_logger.h"

#include <stdio.h>
#include "std.h"

#include "subsystems/imu.h"
#include "subsystems/actuators/motor_mixing.h"
#include "firmwares/rotorcraft/stabilization.h"
#include "firmwares/rotorcraft/stabilization/stabilization_indi.h"
#include "boards/bebop/actuators.h"


/** Set the default File logger path to the USB drive */
#ifndef FILE_LOGGER_PATH
#define FILE_LOGGER_PATH /data/video/usb
#endif

/** The file pointer */
static FILE *file_logger = NULL;

/** Start the file logger and open a new file */
void file_logger_start(void)
{
  uint32_t counter = 0;
  char filename[512];
//  filename = "WLS_LOGGER"

  // Check for available files
  sprintf(filename, "%s/%05d.csv", STRINGIFY(FILE_LOGGER_PATH), counter);
  while ((file_logger = fopen(filename, "r"))) {
    fclose(file_logger);

    counter++;
    sprintf(filename, "%s/%05d.csv", STRINGIFY(FILE_LOGGER_PATH), counter);
  }

  file_logger = fopen(filename, "w");

/*  if (file_logger != NULL) {
    fprintf(
      file_logger,
 	"counter,  phi, theta, psi, phiref, thetaref, psiref, w_obs_1, w_obs_2, w_obs_3, w_obs_4, w_ref_1, w_ref_2, w_ref_3, w_ref_4, mrate_p, mrate_q, mrate_r, refacc_p, refacc_q, refacc_r, estacc_p, estacc_q, estacc_r, w_wls_1, w_wls_2, w_wls_3, w_wls_4, w_act_1, w_act_2, w_act_3, w_act_4, w_mm_1, w_mm_2, w_mm_3, w_mm_4, vlog_1, vlog_2, vlog_3, vlog_4, vlog_5, cmd_ROLL, cmd_PITCH, cmd_YAW, cmd_THRUST \n"); */

 if (file_logger != NULL) {
    fprintf(
      file_logger,
	"counter \n, G[0][0], G[0][1], G[0][2], G[0][3] \n G[1][0], G[1][1], G[1][2], G[1][3] \n G[2][0], G[2][1], G[2][2], G[2][3] \n G2[0][0], G2[0][1], G2[0][2], G2[0][3] \n");
  }
}

/** Stop the logger an nicely close the file */
void file_logger_stop(void)
{
  if (file_logger != NULL) {
    fclose(file_logger);
    file_logger = NULL;
  }
}

/** Log the values to a csv file */
void file_logger_periodic(void)
{
  if (file_logger == NULL) {
    return;
  }
  static uint32_t counter;
//  struct Int32Quat *quat = stateGetNedToBodyQuat_i();
/*  fprintf(file_logger,"%d,", counter);
  fprintf(file_logger, " %f, %f, %f, %f, %f, %f,", phim, thetam, psim, phiref, thetaref, psiref);
  fprintf(file_logger, " %d, %d, %d, %d, %d, %d, %d, %d,", actuators_bebop.rpm_obs[0], actuators_bebop.rpm_obs[1], actuators_bebop.rpm_obs[2], actuators_bebop.rpm_obs[3], actuators_bebop.rpm_ref[0], actuators_bebop.rpm_ref[1], actuators_bebop.rpm_ref[2], actuators_bebop.rpm_ref[3]);
  fprintf(file_logger, " %f, %f, %f, %f, %f, %f,", pratem, qratem, rratem, prateref, qrateref, rrateref);
  fprintf(file_logger, " %f, %f, %f,", accestlog[0], accestlog[1], accestlog[2]);
  fprintf(file_logger, " %f, %f, %f, %f, %f, %f, %f, %f,", u[0], u[1], u[2], u[3], u_cmd_log[0], u_cmd_log[1], u_cmd_log[2], u_cmd_log[3]);
  fprintf(file_logger, " %d, %d, %d, %d,", motor_mixing.commands[0], motor_mixing.commands[1], motor_mixing.commands[2], motor_mixing.commands[3]);
  fprintf(file_logger, " %f, %f, %f, %f, %f,", vlog[0], vlog[1], vlog[2], vlog[3], vlog[4]);
  fprintf(file_logger, " %f, %f, %f, %f \n", cmd_log[0], cmd_log[1], cmd_log[2], cmd_log[3]); */

  // LOGGING ACTUATOR EFFECTIVENESS
  fprintf(file_logger,"%d \n",
          counter);
  // PRINT DUCKING matrices
  fprintf(file_logger,"%f,%f,%f,%f \n",
          G1wls[0][0],
          G1wls[0][1],
          G1wls[0][2],
          G1wls[0][3]);
 fprintf(file_logger,"%f,%f,%f,%f, \n",
          G1wls[1][0],
          G1wls[1][1],
          G1wls[1][2],
          G1wls[1][3]);
fprintf(file_logger,"%f,%f,%f,%f \n",
          G1wls[2][0],
          G1wls[2][1],
          G1wls[2][2],
          G1wls[2][3]);
fprintf(file_logger,"%f,%f,%f,%f \n",
          G2wls[0],
          G2wls[1],
          G2wls[2],
          G2wls[3]);

  counter++;
}
