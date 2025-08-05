/*****************************************************************************
 *
 *  colloid_io_impl_mpio.c
 *
 *  MPI/IO implementation.
 *
 *  (c) 2025 The University of Edinvburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <errno.h>
#include <stdlib.h>

#include "colloid_state_io.h"
#include "colloid_io_impl_mpio.h"
#include "util_fopen.h"

/* Function table */
static colloid_io_impl_vt_t vt_ = {
  (colloid_io_impl_free_ft)  colloid_io_mpio_free,
  (colloid_io_impl_read_ft)  colloid_io_mpio_read,
  (colloid_io_impl_write_ft) colloid_io_mpio_write
};

typedef struct sort_by_key_s {
  int ilocal;       /* a local index to locate colloid in memory */
  int index;        /* unique colloid index is the sort key */
} sort_by_key_t;

int colloid_io_mpio_read_ascii(colloid_io_mpio_t * io, FILE * fp);
int colloid_io_mpio_read_binary(colloid_io_mpio_t * io, FILE * fp);

/*****************************************************************************
 *
 *  compare_by_key
 *
 *  Comparison function for array of sort_by_key_t for qsort() etc.
 *  The sort is on the global id.
 *
 *****************************************************************************/

int compare_by_key(const void * aarg, const void * barg) {

  sort_by_key_t * a = (sort_by_key_t *) aarg;
  sort_by_key_t * b = (sort_by_key_t *) barg;

  int icomp = 0;
  if (a->index > b->index) icomp = +1;
  if (a->index < b->index) icomp = -1;

  return icomp;
}

/*****************************************************************************
 *
 * colloid_io_mpio_create
 *
 *****************************************************************************/

int colloid_io_mpio_create(const colloids_info_t * info,
			   colloid_io_mpio_t ** io) {
  int ifail = 0;
  colloid_io_mpio_t * io_inst = NULL;

  io_inst = (colloid_io_mpio_t *) calloc(1, sizeof(colloid_io_mpio_t));
  if (io_inst == NULL) goto err;

  ifail = colloid_io_mpio_initialise(info, io_inst);
  if (ifail != 0) goto err;

  *io = io_inst;

  return 0;

 err:
  free(io_inst);

  return -1;
}

/*****************************************************************************
 *
 *  colloid_io_mpio_free
 *
 *****************************************************************************/

void colloid_io_mpio_free(colloid_io_mpio_t ** io) {

  assert(io);

  if (*io) {
    colloid_io_mpio_finalise(*io);
    free(*io);
    *io = NULL;
  }

  return;
}

/*****************************************************************************
 *
 *  colloid_io_mpio_initialise
 *
 *****************************************************************************/

int colloid_io_mpio_initialise(const colloids_info_t * info,
			       colloid_io_mpio_t * io) {
  assert(info);
  assert(io);

  *io = (colloid_io_mpio_t) {0};

  io->super.impl = &vt_;
  io->info = (colloids_info_t *) info;

  /* Options! */
  /* FIXME currently in lieu of options */
  {
    int iogrid[3] = {1,1,1};
    int ifail = io_subfile_create(io->info->cs, iogrid, &io->subfile);
    if (ifail != 0) goto err;
  }

  {
    int rank = -1;
    MPI_Comm parent = io->info->cs->commcart;
    MPI_Comm_rank(parent, &rank);
    MPI_Comm_split(parent, io->subfile.index, rank, &io->comm);
  }

  return 0;

 err:
  *io = (colloid_io_mpio_t) {0};

  return -1;
}

/*****************************************************************************
 *
 *  colloid_io_mpio_finalise
 *
 *****************************************************************************/

int colloid_io_mpio_finalise(colloid_io_mpio_t * io) {

  assert(io);

  MPI_Comm_free(&io->comm);

  *io = (colloid_io_mpio_t) {0};
  io->comm = MPI_COMM_NULL;

  return 0;
}

/*****************************************************************************
 *
 *  colloid_io_mpio_read
 *
 *  Read is a case of reading all colloids in the file, and adding
 *  locally only those with relevant positions.
 *
 *  A standard C fread per process. One could try MPI_File_read().
 *
 *****************************************************************************/

int colloid_io_mpio_read(colloid_io_mpio_t * io, const char * filename) {

  int ifail = 0;

  assert(io);
  assert(filename);

  FILE * fp = util_fopen(filename, "r");

  if (fp == NULL) ifail = errno;
  if (ifail != 0) goto err;

  {
    int asc = 1; /* FIXME from options please */
    int bin = 0; /* FIXME */
    if (asc) ifail = colloid_io_mpio_read_ascii(io, fp);
    if (bin) ifail = colloid_io_mpio_read_binary(io, fp);

    fclose(fp);
  }

 err:
  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_mpio_write
 *
 *  The write will generate a file with colloids sorted into unique id
 *  order. This is used to provide te file offset for MPI/IO.
 *
 *  1. Form an array of size nlocal to hold colloid ids in sorted
 *     order.
 *  2. Write a contiguous buffer of colloid data in the sorted order.
 *     The buffer is of type MPI_CHAR and size per element depends on
 *     the ascii or binary option.
 *  3. Form an MPI_Type_indexed_block() with the sorted ids to degine
 *     the file 'filetype'
 *  4. set view MPI_File_set_view(fh, disp, etype, filetpye, datarep, info)
 *  5. MPI_File_write_all(fh, ...)
 *
 *****************************************************************************/

int colloid_io_mpio_write(colloid_io_mpio_t * io, const char * filename) {

  int ifail = 0;

  /* 1. Form the sorted list of colloid ids (local). */

  /* For example, locally we might have colloids with unique ids
   * 10, 4, and 12 in that order (from local pointer list) . So
   * the local indices would be [0, 1, 2]. We would need to write
   * these local indices in the order [1, 0, 2], ie., ids [4, 10, 12]. */

  /* The corresponsing global displacements are the (index - 1) in
   * the range [0, ntotal - 1]. */

  int nlocal = 0;
  int ntotal = 0;
  int nbufsz = COLLOID_BUFSZ;
  int * idisp = NULL;
  char * buf = NULL;
  colloid_t * pc = NULL;
  sort_by_key_t * sort = NULL;

  assert(io);
  assert(filename);

  colloids_info_ntotal(io->info, &ntotal);
  colloids_info_nlocal(io->info, &nlocal);
  if (ntotal == 0) return 0;

  assert(nlocal > 0); /* FIXME No zero sized buffers please */

  sort = (sort_by_key_t *) calloc(nlocal, sizeof(sort_by_key_t));

  colloids_info_local_head(io->info, &pc);

  for (int n = 0; pc; pc = pc->nextlocal, n++) {
    sort[nlocal - (1 + n)].ilocal = n;
    sort[nlocal - (1 + n)].index  = pc->s.index;
  }

  qsort(sort, nlocal, sizeof(sort_by_key_t), compare_by_key);

  /* Write the ordered particle data to the contiguous buffer */

  colloids_info_local_head(io->info, &pc);

  for (int ibuf = 0; pc; pc = pc->nextlocal, ibuf++) {
    int icol = sort[ibuf].ilocal;
    colloid_state_io_write_buf_ascii(&pc->s, buf + icol*nbufsz);
  }

  /* Displacements are global indices - 1 */

  idisp = (int *) calloc(nlocal, sizeof(int));

  for (int ibuf = 0; ibuf < nlocal; ibuf++) {
    int icol = sort[ibuf].ilocal;
    idisp[icol] = sort[icol].index - 1;
  }

  {
    MPI_Datatype etype = MPI_DATATYPE_NULL;
    MPI_Datatype filetype = MPI_DATATYPE_NULL;

    MPI_Type_contiguous(nbufsz, MPI_CHAR, &etype);
    MPI_Type_commit(&etype);
    MPI_Type_create_indexed_block(nlocal, 1, idisp, etype, &filetype);
    MPI_Type_commit(&filetype);

    /* File */
    {
      MPI_File fh = MPI_FILE_NULL;
      MPI_Info info = MPI_INFO_NULL;
      MPI_Offset hdisp = 23*sizeof(char); /* ASCII ONLY */

      int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
      int rank = 0;


      ifail = MPI_File_open(io->comm, filename, amode, info, &fh);
      if (ifail != MPI_SUCCESS) goto err_file;

      MPI_Comm_rank(io->comm, &rank);

      if (rank == 0) {
	/* Write header */
	/* ascii header is 23 char for total number (total number in file) */
	char hdr[BUFSIZ] = {0};
	snprintf(hdr, 23 + 1, "%22d\n", ntotal);
	ifail = MPI_File_write_at(fh, 0, hdr, 23, MPI_CHAR, MPI_STATUS_IGNORE);
	if (ifail != MPI_SUCCESS) goto err_file;
      }

      ifail = MPI_File_set_view(fh, hdisp, etype, filetype, "native", info);
      if (ifail != MPI_SUCCESS) goto err_file;
      ifail = MPI_File_write_all(fh, buf, nlocal, etype, MPI_STATUS_IGNORE);
      if (ifail != MPI_SUCCESS) goto err_file;
      ifail = MPI_File_close(&fh);
    }
  err_file:
    
    MPI_Type_free(&filetype);
    MPI_Type_free(&etype);
  }

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_mpio_read_ascii
 *
 *****************************************************************************/

int colloid_io_mpio_read_ascii(colloid_io_mpio_t * io, FILE * fp) {

  int ifail = 0;
  int nread = 0; 

  assert(io);
  assert(fp);

  int nr = fscanf(fp, "%22d\n", &nread);
  if (nr != 1) goto err;

  for (int n = 0; n < nread; n++) {
    char buf[COLLOID_BUFSZ + 1] = {0}; /* Make sure we have a \0 */
    colloid_state_t s = {0};

    nr = fread(buf, COLLOID_BUFSZ*sizeof(char), 1, fp);
    if (nr != 1) goto err;
    ifail = colloid_state_io_read_buf_ascii(&s, buf);
    if (ifail == 0) goto err;
    /* Add if local and assign the state */
    {
      colloid_t * pc = NULL;
      colloids_info_add_local(io->info, s.index, s.r, &pc);
      if (pc) pc->s = s;
    }
  }

 err:

  /* This is collective so cannot by-pass on error */
  {
    int ntotal = 0;
    int nlocal = 0;
    colloids_info_nlocal(io->info, &nlocal);
    MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, io->comm);
    if (ntotal != nread) ifail = -1;
  }
   

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_mpio_read_binary
 *
 *****************************************************************************/

int colloid_io_mpio_read_binary(colloid_io_mpio_t * io, FILE * fp) {

  int ifail = 0;
  int nread = 0; 

  assert(io);
  assert(fp);

  /* Header */
  int nr = fread(&nread, sizeof(int), 1, fp);
  if (nr != 1) goto err;

  /* Data: one at a time at the moment */
  for (int n = 0; n < nread; n++) {
    const size_t nbufsz = sizeof(colloid_state_t);
    char buf[nbufsz + 1] = {}; /* Make sure we have a \0 */
    colloid_state_t s = {0};

    nr = fread(buf, nbufsz, 1, fp);
    if (nr != 1) goto err;
    ifail = colloid_state_io_read_buf(&s, buf);
    if (ifail == 0) goto err;
    /* Add if local and assign the state */
    {
      colloid_t * pc = NULL;
      colloids_info_add_local(io->info, s.index, s.r, &pc);
      if (pc) pc->s = s;
    }
  }

 err:

  /* This is collective so cannot by-pass on error */
  /* The only check we can make here is that the global total is correct
   * compared with the number declared in the file */
  /* FIXME: THIS IS NO GOOD IF NFILES > 1 (need the xcomm) */
  {
    int nlocal = 0;
    int ntotal = 0;
    colloids_info_nlocal(io->info, &nlocal);
    MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, io->comm);
    if (ntotal != nread) ifail = -1;
  }
  /* If the total is corrent, then we can set the info total */

  return ifail;
}
