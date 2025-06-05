/*****************************************************************************
 *
 *  colloid_io_impl_ansi.h
 *
 *  (c) 2025 The University of Edinburgh
 *
 *  Edinburgh Soft MAtter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *****************************************************************************/

#ifndef LUDWIG_COLLOID_IO_IMPL_ANSI_H
#define LUDWIG_COLLOID_IO_IMPL_ANSI_H

#include "colloid_io_impl.h"

typedef struct colloid_io_ansi_s colloid_io_ansi_t;

struct colloid_io_ansi_s {
  colloid_io_impl_vt_t super;             /* superclass block */
  /* metadata / options block to include cs_t pointer */
  /* io gris; ascii/binary */

  /* State of implementation ... */
};

int colloid_io_ansi_create(object, colloid_io_ansi_t ** io);
int colloid_io_ansi_free(colloid_io_ansi_t ** io);

int colloid_io_ansi_initialise(const objext, colloid_io_ansi_t * io);
int colloid_io_ansi_finalise(colloid_io_ansi_t * io);

int colloid_io_ansi_read(colloid_io_ansi_t * io, const char * filename);
int colloid_io_ansi_write(colloid_io_impl_t * io, const char * filename);

#endif
