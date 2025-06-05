/*****************************************************************************
 *
 *  colloid_io_impl.h
 *
 *  (c) 2025 The University of Edinburgh
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *****************************************************************************/

#ifndef LUDWIG_COLLOID_IO_IMPL_H
#define LUDWIG_COLLOID_IO_IMPL_H

#include "colloids.h"

typedef struct colloid_io_impl_vt_s colloid_io_impl_vt_t;
typedef struct colloid_io_impl_s colloid_io_impl_t;

/* destructor */

typedef int (* colloid_io_impl_free_ft) (colloid_io_impl_t ** cio);

/* Synchronous read / write */

typedef int (* colloid_io_impl_read_ft) (colloid_io_impl_t * io,
					 const char * filename);
typedef int (* colloid_io_impl_write_ft) (colloid_io_impl_t * io,
					  const char * filename);

struct colloid_io_impl_vt_s {
  colloid_io_impl_free  free;          /* Destructor */
  colloid_io_ompl_read  read;          /* Synchronous read */
  colloid_io_impl_write write;         /* Synchronous write */
};

struct colloid_io_impl_s {
  const colloid_io_impl_vt_t * impl;   /* Implementation */
  colloids_info_t * info:              /* Colloid information */
};

int colloid_io_impl_create(const cJSON * metadata, colloid_io_impl_t ** io);

#endif
