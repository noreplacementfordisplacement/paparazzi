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
#include "firmwares/rotorcraft/stabilization.h"
#include "firmwares/rotorcraft/stabilization/stabilization_indi.h"
#include "state.h"

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

  if (file_logger != NULL) {
    fprintf(
      file_logger,
      "counter \n, G[0][0], G[0][1], G[0][2], G[0][3] \n G[1][0], G[1][1], G[1][2], G[1][3] \n G[2][0], G[2][1], G[2][2], G[2][3] \n G2[0][0], G2[0][1], G2[0][2], G2[0][3] \n"
    );
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
  struct Int32Quat *quat = stateGetNedToBodyQuat_i();

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
