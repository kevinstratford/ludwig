/*****************************************************************************
 *
 *  colloid_options.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2025 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_COLLOID_OPTIONS_H
#define LUDWIG_COLLOID_OPTIONS_H

#include "colloid_io_options.h"

typedef struct colloid_options_s colloid_options_t;

struct colloid_options_s {

  int    ncell[3];            /* local cell list extent (request) */
  int    nvel;                /* discrete velocities in model */
  double rho0;                /* Colloid density (uniform) FIXME */

  /* i/o details */
  colloid_io_options_t input;  /* input mode, format, ... */
  colloid_io_options_t output; /* output mode, format, ... */
  int                  nfreq;  /* Output frequency */

  /* Additional switches etc */
  int have_colloids;           /* if one or more colloids globally */
};

colloid_options_t colloid_options_default(void);
colloid_options_t colloid_options_have_colloids(int have_colloids);
colloid_options_t colloid_options_ncell(const int ncell[3]);

/* with some control ... */

int colloid_options_ncell_valid(const int nceel[3]);
int colloid_options_valid(const colloid_options_t * options);

int colloid_options_to_json(const colloid_options_t * opts, cJSON ** json);
int colloid_options_from_json(const cJSON * json, colloid_options_t * opts);

#endif
