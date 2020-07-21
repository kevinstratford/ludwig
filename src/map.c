/*****************************************************************************
 *
 *  map.c
 *
 *  Map of fluid lattice sites containing status (fluid, solid, etc)
 *  and with space for additional data such as wetting parameters.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2012-2020 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "pe.h"
#include "coords.h"
#include "coords_field.h"
#include "map_s.h"

static int map_read(FILE * fp, int index, void * self);
static int map_write(FILE * fp, int index, void * self);
static int map_read_ascii(FILE * fp, int index, void * self);
static int map_write_ascii(FILE * fp, int index, void * self);

/*****************************************************************************
 *
 *  map_create
 *
 *  Allocate the map object, the status field, and any additional
 *  data required.
 *
 *****************************************************************************/

__host__ int map_create(pe_t * pe, cs_t * cs, int ndata, map_t ** pobj) {

  int nsites;
  int nhalo;
  int ndevice;
  map_t * obj = NULL;

  assert(pe);
  assert(cs);
  assert(ndata >= 0);
  assert(pobj);

  cs_nsites(cs, &nsites);
  cs_nhalo(cs, &nhalo);

  obj = (map_t*) calloc(1, sizeof(map_t));
  assert(obj);
  if (obj == NULL) pe_fatal(pe, "calloc(map_t) failed\n");

  obj->status = (char*) calloc(nsites, sizeof(char));
  assert(obj->status);
  if (obj->status == NULL) pe_fatal(pe, "calloc(map->status) failed\n");

  obj->pe = pe;
  obj->cs = cs;
  obj->nsite = nsites;
  obj->ndata = ndata;

  /* Could be zero-sized array */

  obj->data = (double*) calloc((size_t) ndata*nsites, sizeof(double));
  assert(obj->data);
  if (ndata > 0 && obj->data == NULL) pe_fatal(pe, "calloc(map->data) failed\n");

  coords_field_init_mpi_indexed(cs, nhalo, 1, MPI_CHAR, obj->halostatus);
  if (obj->ndata) {
    coords_field_init_mpi_indexed(cs, nhalo, obj->ndata, MPI_DOUBLE,
				  obj->halodata);
  }


  /* Allocate target copy of structure (or alias) */

  tdpGetDeviceCount(&ndevice);

  if (ndevice == 0) {
    obj->target = obj;
  }
  else {
    char * tmp = NULL;
    double * dtmp = NULL;

    tdpMalloc((void **) &obj->target, sizeof(map_t));
    tdpMemset(obj->target, 0, sizeof(map_t));
    tdpMalloc((void **) &tmp, nsites*sizeof(char));
    tdpMemset(tmp, 0, nsites*sizeof(char));

    tdpMemcpy(&obj->target->status, &tmp, sizeof(char *),
	      tdpMemcpyHostToDevice);
    tdpMemcpy(&obj->target->is_porous_media, &obj->is_porous_media,
	      sizeof(int), tdpMemcpyHostToDevice); 
    tdpMemcpy(&obj->target->nsite, &obj->nsite, sizeof(int),
	      tdpMemcpyHostToDevice);
    tdpMemcpy(&obj->target->ndata, &obj->ndata, sizeof(int),
	      tdpMemcpyHostToDevice); 

    /* Data */
    if (obj->ndata > 0) {
      size_t nsz = (size_t) obj->ndata*nsites*sizeof(double);
      tdpAssert(tdpMalloc((void **) &dtmp, nsz));
      tdpAssert(tdpMemset(dtmp, 0, nsz));
      tdpAssert(tdpMemcpy(&obj->target->data, &dtmp, sizeof(double *),
		          tdpMemcpyHostToDevice));
    }
  }

  *pobj = obj;

  return 0;
}

/*****************************************************************************
 *
 *  map_free
 *
 *****************************************************************************/

__host__ int map_free(map_t * obj) {

  int ndevice;
  char * tmp = NULL;
  double * dtmp = NULL;

  assert(obj);

  tdpGetDeviceCount(&ndevice);

  if (ndevice > 0) {
    tdpAssert(tdpMemcpy(&tmp, &obj->target->status, sizeof(char *),
			tdpMemcpyDeviceToHost)); 
    tdpAssert(tdpFree(tmp));

    if (obj->ndata > 0) {
      tdpAssert(tdpMemcpy(&dtmp, &obj->target->data, sizeof(double *),
			  tdpMemcpyDeviceToHost));
      tdpAssert(tdpFree(dtmp));
    }
    tdpFree(obj->target);
  }

  MPI_Type_free(&obj->halostatus[X]);
  MPI_Type_free(&obj->halostatus[Y]);
  MPI_Type_free(&obj->halostatus[Z]);

  if (obj->ndata > 0) {
    MPI_Type_free(&obj->halodata[X]);
    MPI_Type_free(&obj->halodata[Y]);
    MPI_Type_free(&obj->halodata[Z]);
  }

  if (obj->info) io_info_free(obj->info);

  free(obj->data);
  free(obj->status);
  free(obj);

  return 0;
}

/*****************************************************************************
 *
 *  map_memcpy
 *
 *****************************************************************************/

__host__ int map_memcpy(map_t * map, tdpMemcpyKind flag) {

  int ndevice;
  char * tmp;

  assert(map);

  tdpGetDeviceCount(&ndevice);

  if (ndevice == 0) {
    /* Ensure we alias */
    assert(map->target == map);
  }
  else {
    tdpMemcpy(&tmp, &map->target->status, sizeof(char *),
	      tdpMemcpyDeviceToHost);

    switch (flag) {
    case tdpMemcpyHostToDevice:
      tdpMemcpy(tmp, map->status, map->nsite*sizeof(char), flag);
      break;
    case tdpMemcpyDeviceToHost:
      tdpMemcpy(map->status, tmp, map->nsite*sizeof(char), flag);
      break;
    default:
      pe_fatal(map->pe, "Bad flag in map_memcpy()\n");
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  map_init_io_info
 *
 *****************************************************************************/

__host__
int map_init_io_info(map_t * obj, int grid[3], int form_in, int form_out) {

  io_info_arg_t args;
  size_t sz;

  assert(obj);
  assert(obj->info == NULL);

  args.grid[X] = grid[X];
  args.grid[Y] = grid[Y];
  args.grid[Z] = grid[Z];
  
  io_info_create(obj->pe, obj->cs, &args, &obj->info);
  if (obj->info == NULL) pe_fatal(obj->pe, "io_info_create(map) failed\n");

  io_info_set_name(obj->info, "map");
  io_info_write_set(obj->info, IO_FORMAT_BINARY, map_write);
  io_info_write_set(obj->info, IO_FORMAT_ASCII, map_write_ascii);
  io_info_read_set(obj->info, IO_FORMAT_BINARY, map_read);
  io_info_read_set(obj->info, IO_FORMAT_ASCII, map_read_ascii);

  sz = sizeof(char) + obj->ndata*sizeof(double);
  io_info_set_bytesize(obj->info, IO_FORMAT_BINARY, sz);
  io_info_set_bytesize(obj->info, IO_FORMAT_ASCII, 2 + 23*obj->ndata + 1);

  io_info_format_set(obj->info, form_in, form_out);
  io_info_metadata_filestub_set(obj->info, "map");

  return 0;
}

/*****************************************************************************
 *
 *  map_io_info
 *
 *****************************************************************************/

__host__ int map_io_info(map_t * obj, io_info_t ** info) {

  assert(obj);
  assert(info);

  *info = obj->info;

  return 0;
}

/*****************************************************************************
 *
 *  map_halo
 *
 *  This should only be required when porous media read from file.
 *
 *****************************************************************************/

__host__ int map_halo(map_t * obj) {

  int nhalo;

  assert(obj);

  cs_nhalo(obj->cs, &nhalo);

  coords_field_halo_rank1(obj->cs, obj->nsite, nhalo, 1, obj->status,
			  MPI_CHAR);
  if (obj->ndata) {
    coords_field_halo_rank1(obj->cs, obj->nsite, nhalo, obj->ndata, obj->data,
			    MPI_DOUBLE);
  }

  return 0;
}

/*****************************************************************************
 *
 *  map_status
 *
 *  Return map status as an integer.
 *
 *****************************************************************************/

__host__ __device__
int map_status(map_t * obj, int index, int * status) {

  assert(obj);
  assert(status);

  *status = (int) obj->status[addr_rank0(obj->nsite, index)];

  return 0;
}

/*****************************************************************************
 *
 *  map_status_set
 *
 *****************************************************************************/

__host__ __device__
int map_status_set(map_t * obj, int index, int status) {

  assert(obj);
  assert(status >= 0);
  assert(status < MAP_STATUS_MAX);

  obj->status[addr_rank0(obj->nsite, index)] = status;

  return 0;
}

/*****************************************************************************
 *
 *  map_ndata
 *
 *****************************************************************************/

__host__ __device__
int map_ndata(map_t * obj, int * ndata) {

  assert(obj);
  assert(ndata);

  *ndata = obj->ndata;

  return 0;
}

/*****************************************************************************
 *
 *  map_data
 *
 *****************************************************************************/

__host__ __device__
int map_data(map_t * obj, int index, double * data) {

  int n;

  assert(obj);
  assert(data);

  for (n = 0; n < obj->ndata; n++) {
    data[n] = obj->data[addr_rank1(obj->nsite, obj->ndata, index, n)];
  }

  return 0;
}

/*****************************************************************************
 *
 *  map_data_set
 *
 *****************************************************************************/

__host__ __device__
int map_data_set(map_t * obj, int index, double * data) {

  int n;

  assert(obj);
  assert(data);

  for (n = 0; n < obj->ndata; n++) {
    obj->data[addr_rank1(obj->nsite, obj->ndata, index, n)] = data[n];
  }

  return 0;
}

/*****************************************************************************
 *
 *  map_pm
 *
 *****************************************************************************/

__host__ int map_pm(map_t * obj, int * flag) {

  assert(obj);
  assert(flag);

  *flag = obj->is_porous_media;

  return 0;
}

/*****************************************************************************
 *
 *  map_pm_set
 *
 *****************************************************************************/

__host__ int map_pm_set(map_t * obj, int flag) {

  assert(obj);

  obj->is_porous_media = flag;

  return 0;
}

/*****************************************************************************
 *
 *  map_volume_allreduce
 *
 *****************************************************************************/

__host__ int map_volume_allreduce(map_t * obj, int status, int * volume) {

  int vol_local;
  MPI_Comm comm;

  assert(obj);
  assert(volume);

  map_volume_local(obj, status, &vol_local);

  cs_cart_comm(obj->cs, &comm);
  MPI_Allreduce(&vol_local, volume, 1, MPI_INT, MPI_SUM, comm);

  return 0;
}

/*****************************************************************************
 *
 *  map_volume_local
 *
 *****************************************************************************/

int map_volume_local(map_t * obj, int status_wanted, int * volume) {

  int nlocal[3];
  int ic, jc, kc, index;
  int vol_local;
  int status;

  assert(obj);
  assert(volume);

  cs_nlocal(obj->cs, nlocal);
  vol_local = 0;

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = cs_index(obj->cs, ic, jc, kc);
	map_status(obj, index, &status);
	if (status == status_wanted) vol_local += 1;
      }
    }
  }

  *volume = vol_local;

  return 0;
}

/*****************************************************************************
 *
 *  map_write
 *
 *****************************************************************************/

static int map_write(FILE * fp, int index, void * self) {

  int n, nw;
  int indexf;
  map_t * obj = (map_t*) self;

  assert(fp);
  assert(obj);

  indexf = addr_rank0(obj->nsite, index);
  nw = fwrite(&obj->status[indexf], sizeof(char), 1, fp);
  if (nw != 1) pe_fatal(obj->pe, "fwrite(map->status) failed\n");

  for (n = 0; n < obj->ndata; n++) {
    indexf = addr_rank1(obj->nsite, obj->ndata, index, n);
    nw = fwrite(&obj->data[indexf], sizeof(double), 1, fp);
    if (nw != 1) pe_fatal(obj->pe, "fwrite(map->data) failed\n");
  }

  return 0;
}

/*****************************************************************************
 *
 *  map_read
 *
 *****************************************************************************/

static int map_read(FILE * fp, int index, void * self) {

  int n, nr;
  int indexf;
  map_t * obj = (map_t*) self;

  assert(fp);
  assert(obj);

  indexf = addr_rank0(obj->nsite, index);
  nr = fread(&obj->status[indexf], sizeof(char), 1, fp);
  if (nr != 1) pe_fatal(obj->pe, "fread(map->status) failed");

  for (n = 0; n < obj->ndata; n++) {
    indexf = addr_rank1(obj->nsite, obj->ndata, index, n);
    nr = fread(&obj->data[indexf], sizeof(double), 1, fp);
    if (nr != 1) pe_fatal(obj->pe, "fread(map->data) failed\n");
  }

  return 0;
}

/*****************************************************************************
 *
 *  map_write_ascii
 *
 *****************************************************************************/

static int map_write_ascii(FILE * fp, int index, void * self) {

  int n, nw;
  int indexf;
  int status;
  map_t * obj = (map_t *) self;

  assert(fp);
  assert(obj);

  status = obj->status[addr_rank0(obj->nsite, index)];
  assert(status < 99); /* Fixed format check. */

  nw = fprintf(fp, "%2d", status);
  if (nw != 2) pe_fatal(obj->pe, "fprintf(map->status) failed\n");

  for (n = 0; n < obj->ndata; n++) {
    indexf = addr_rank1(obj->nsite, obj->ndata, index, n);
    fprintf(fp, " %22.15e", obj->data[indexf]);
  }

  nw = fprintf(fp, "\n");
  if (nw != 1) pe_fatal(obj->pe, "fprintf(map) failed\n");

  return 0;
}

/*****************************************************************************
 *
 *  map_read_ascii
 *
 *****************************************************************************/

static int map_read_ascii(FILE * fp, int index, void * self) {

  int n, nr;
  int indexf;
  int status;
  map_t * obj = (map_t*) self;

  assert(fp);
  assert(obj);

  nr = fscanf(fp, "%2d", &status);
  if (nr != 1) pe_fatal(obj->pe, "fscanf(map->status) failed\n");
  obj->status[addr_rank0(obj->nsite, index)] = status;

  for (n = 0; n < obj->ndata; n++) {
    indexf = addr_rank1(obj->nsite, obj->ndata, index, n);
    nr = fscanf(fp, " %le", obj->data + indexf);
    if (nr != 1) pe_fatal(obj->pe, "fscanf(map->data[%d]) failed\n", n);
  }

  return 0;
}
