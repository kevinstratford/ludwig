/*****************************************************************************
 *
 *  colloid_io_rt.c
 *
 *  Run time colloid I/O settings.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2018 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <string.h>

#include "pe.h"
#include "runtime.h"
#include "colloid_io_rt.h"


/*****************************************************************************
 *
 *  colloid_io_run_time
 *
 *****************************************************************************/

int colloid_io_run_time(pe_t * pe, rt_t * rt, cs_t * cs,
			colloids_info_t * cinfo,
			colloid_io_t ** pcio) {

  int io_grid[3] = {1, 1, 1};
  char tmp[BUFSIZ];

  colloid_io_t * cio = NULL;

  assert(pe);
  assert(rt);
  assert(cs);
  assert(cinfo);

  /* Control user input */
  {
    int io_grid_req[3] = {1, 1, 1};
    int isvalid = 1;
    rt_int_parameter_vector(rt, "default_io_grid", io_grid_req);
    rt_int_parameter_vector(rt, "colloid_io_grid", io_grid_req);

    if (1 > io_grid_req[X] || io_grid_req[X] > 8) isvalid = 0;
    if (1 > io_grid_req[Y] || io_grid_req[Y] > 8) isvalid = 0;
    if (1 > io_grid_req[Z] || io_grid_req[Z] > 8) isvalid = 0;
    if (isvalid == 1) {
      io_grid[X] = io_grid_req[X];
      io_grid[Y] = io_grid_req[Y];
      io_grid[Z] = io_grid_req[Z];
    }
    else {
      pe_info(pe, "Colloid i/o grid out-of-range\n");
      pe_info(pe, "io_grid: %d %d %d\n", io_grid[X], io_grid[Y], io_grid[Z]);
      pe_exit(pe, "Please check and try again\n");
    }
  }

  colloid_io_create(pe, cs, io_grid, cinfo, &cio);
  assert(cio);

  /* Default format to ascii, parallel; then check user input */

  colloid_io_format_input_ascii_set(cio);
  colloid_io_format_output_ascii_set(cio);

  strcpy(tmp, "");
  rt_string_parameter(rt, "colloid_io_format", tmp, BUFSIZ);

  if (strncmp("BINARY", tmp, 5) == 0 || strncmp("binary", tmp, 5) == 0) {
    colloid_io_format_input_binary_set(cio);
    colloid_io_format_output_binary_set(cio);
  }

  rt_string_parameter(rt, "colloid_io_format_input", tmp, BUFSIZ);

  if (strncmp("ASCII",  tmp, 5) == 0 || strncmp("ascii", tmp, 5) == 0) {
    colloid_io_format_input_ascii_set(cio);
  }

  if (strncmp("ASCII_SERIAL", tmp, 12) == 0 ||
      strncmp("ascii_serial", tmp, 12) == 0) {
    colloid_io_format_input_ascii_set(cio);
    colloid_io_format_input_serial_set(cio);
  }

  if (strncmp("BINARY", tmp, 6) == 0 || strncmp("binary", tmp, 6) == 0) {
    colloid_io_format_input_binary_set(cio);
  }

  if (strncmp("BINARY_SERIAL", tmp, 13) == 0 ||
      strncmp("binary_serial", tmp, 13) == 0) {
    colloid_io_format_input_binary_set(cio);
    colloid_io_format_input_serial_set(cio);
  }

  rt_string_parameter(rt, "colloid_io_format_output", tmp, BUFSIZ);

  if (strncmp("ASCII",  tmp, 5) == 0 || strncmp("ascii", tmp, 5) == 0) {
    colloid_io_format_output_ascii_set(cio);
  }

  if (strncmp("BINARY", tmp, 6) == 0 || strncmp("binary", tmp, 6) == 0) {
    colloid_io_format_output_binary_set(cio);
  }

  colloid_io_info(cio);

  *pcio = cio;

  return 0;
}

