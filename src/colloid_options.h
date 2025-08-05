/*****************************************************************************
 *
 *  colloid_options.h
 *
 *  This is conflating colloid_io_options_t, and colloid_options_t
 *  FIXME: split these apart.
 *
 *****************************************************************************/

#ifndef LUDWIG_COLLOID_OPTIONS_H
#define LUDWIG_COLLOID_OPTIONS_H

#include "io_options.h"

typedef enum colloid_io_mode_enum {
  COLLOID_IO_MODE_INVALID,
  COLLOID_IO_MODE_ANSI,
  COLLOID_IO_MODE_MPIIO
} colloid_io_mode_enum_t;

int colloid_io_mode_valid(colloid_io_mode_enum_t mode);
const char * colloid_io_mode_to_string(colloid_io_mode_enum_t mode);
colloid_io_mode_enum_t colloid_io_mode_from_string(const char * str);

typedef struct colloid_io_options_s colloid_io_options_t;

struct colloid_io_options_s {
  /* There are some common features with lattice i/o options,
   * but it's not the same thing... */

  colloid_io_mode_enum_t  mode;        /* ansi or  mpio */
  io_record_format_enum_t iorformat;   /* ascii or binary */
  int                     report;      /* verbose report switch */
  int                     iogrid[3];   /* i/o (file) decomposition */
};


/* io first */
colloid_io_options_t colloid_io_options_default(void);
int colloid_io_options_valid(const colloid_io_options_t * opts);

int colloid_io_options_to_json(const colloid_io_options_t * opts,
			       cJSON ** json);
int colloid_io_options_from_json(const cJSON * json,
				 colloid_io_options_t * opts);


typedef struct colloid_options_s colloid_options_t;

struct colloid_options_s {

  int                  ncell[3];      /* local cell list extent (request) */
  int                  nvel;          /* discrete velocities in model */
  double               rho0;          /* Colloid density (uniform) FIXME */

  /* i/o details */
  colloid_io_options_t input;         /* input mode, format, ... */
  colloid_io_options_t output;        /* output mode, format, ... */
  int                  nfreq;         /* Output frequency */

};


colloid_options_t colloid_options_default(void);

/* with some control ... ... */

int colloid_options_ncell_valid(const int nceel[3]);
int colloid_options_valid(const colloid_options_t * options);

int colloid_options_to_json(const colloid_options_t * opts, cJSON ** json);
int colloid_options_from_json(const cJSON * json, colloid_options_t * opts);

#endif
