/*****************************************************************************
 *
 *  colloid_io_impl_ansi.c
 *
 *  Original "ANSI" implementation for colloid input and output.
 *
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2025 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "colloid_io_impl_ansi.h"
#include "colloid_array_util.h"
#include "util_fopen.h"

/* Function table */
static colloid_io_impl_vt_t vt_ = {
    (colloid_io_impl_free_ft)  colloid_io_ansi_free,
    (colloid_io_impl_read_ft)  colloid_io_ansi_read,
    (colloid_io_impl_write_ft) colloid_io_ansi_write
};

static int colloid_io_ansi_check_read(colloid_io_ansi_t * io, int ngroup);

int colloid_io_ansi_write_ascii(FILE * fp, const colloid_array_t * buf);
int colloid_io_ansi_write_binary(FILE * fp, const colloid_array_t * buf);

int colloid_io_ansi_read_ascii(colloid_io_ansi_t * io, FILE * fp);
int colloid_io_ansi_read_binary(colloid_io_ansi_t * io, FILE * fp);

/*****************************************************************************
 *
 *  colloid_io_ansi_create
 *
 *****************************************************************************/

int colloid_io_ansi_create(const colloids_info_t * info,
                           colloid_io_ansi_t ** io) {

  int                 ifail   = 0;
  colloid_io_ansi_t * io_inst = NULL;

  io_inst = (colloid_io_ansi_t *) calloc(1, sizeof(colloid_io_ansi_t));
  if (io_inst == NULL) goto err;

  ifail = colloid_io_ansi_initialise(info, io_inst);
  if (ifail != 0) goto err;

  *io = io_inst;

  return 0;

 err:
  free(io_inst);

  return -1;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_free
 *
 *****************************************************************************/

void colloid_io_ansi_free(colloid_io_ansi_t ** io) {

  assert(io);

  if (*io) {
    colloid_io_ansi_finalise(*io);
    free(*io);
    *io = NULL;
  }

  return;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_initialise
 *
 *****************************************************************************/

int colloid_io_ansi_initialise(const colloids_info_t * info,
                               colloid_io_ansi_t * io) {
  assert(info);
  assert(io);

  *io = (colloid_io_ansi_t) {0};

  io->super.impl = &vt_;
  io->info       = (colloids_info_t *) info;

  /* FIXME currently in lieu of options. Input or output?  */
  {
    int iogrid[3] = {1, 1, 1};
    int ifail     = io_subfile_create(io->info->cs, iogrid, &io->subfile);
    if (ifail != 0) goto err;
  }

  {
    int      rank   = -1;
    MPI_Comm parent = io->info->cs->commcart;
    MPI_Comm_rank(parent, &rank);
    MPI_Comm_split(parent, io->subfile.index, rank, &io->comm);
  }

  return 0;

 err:
  *io = (colloid_io_ansi_t) {0};

  return -1;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_finalise
 *
 *****************************************************************************/

int colloid_io_ansi_finalise(colloid_io_ansi_t * io) {

  assert(io);

  MPI_Comm_free(&io->comm);

  *io      = (colloid_io_ansi_t) {0};
  io->comm = MPI_COMM_NULL;

  return 0;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_write
 *
 *  Write information on the colloids in the local domain proper (not
 *  including the halos) to the specified file.
 *
 *  In parallel, one file per io group is used. Processes within a
 *  group aggregate data to root, and root writes. So this does not
 *  scale very well.
 *
 *****************************************************************************/

int colloid_io_ansi_write(colloid_io_ansi_t * io, const char * filename) {

  int   ifail  = 0;
  int   myrank = -1;
  int   nranks = -1;
  int   nlocal = 0;
  int   ntotal = 0;
  int * nclist = NULL;
  int * displ  = NULL;

  /* consolidated buffers for particle data */
  colloid_array_t cbuf = {0};
  colloid_array_t rbuf = {0};

  assert(io);

  colloids_info_ntotal(io->info, &ntotal);
  if (ntotal == 0) return 0; /* All ranks must agree */

  MPI_Comm_rank(io->comm, &myrank);
  MPI_Comm_size(io->comm, &nranks);

  /* Gather a list of colloid numbers from each rank in the group */

  nclist = (int *) calloc(nranks, sizeof(int));
  assert(nclist);
  if (nclist == NULL) goto err;

  colloids_info_nlocal(io->info, &nlocal);
  assert(nlocal > 0); /* Assume all files have at least 1 colloid */

  MPI_Gather(&nlocal, 1, MPI_INT, nclist, 1, MPI_INT, 0, io->comm);

  /* Allocate local buffer, pack. */

  ifail = colloid_array_alloc(0, nlocal, &cbuf);
  if (ifail != 0) goto err;

  {
    colloid_t * pc = NULL;
    colloids_info_local_head(io->info, &pc);

    for (int n = 0; pc; pc = pc->nextlocal, n += 1) {
      memcpy(&cbuf.data[n], &pc->s, sizeof(colloid_state_t));
    }
  }

  if (myrank == 0) {

    /* Work out displacements and the total */
    /* Allocate total receive buffer */

    displ = (int *) malloc(nranks*sizeof(int));
    assert(displ);
    if (displ == NULL) goto err;

    displ[0] = 0;
    for (int n = 1; n < nranks; n++) {
      displ[n] = displ[n - 1] + nclist[n - 1];
    }

    ntotal = 0;
    for (int n = 0; n < nranks; n++) {
      ntotal    += nclist[n];
      displ[n]  *= sizeof(colloid_state_t); /* to bytes */
      nclist[n] *= sizeof(colloid_state_t); /* ditto */
    }

    ifail = colloid_array_alloc(0, ntotal, &rbuf);
    if (ifail != 0) goto err;
  }

  MPI_Gatherv(cbuf.data, nlocal*sizeof(colloid_state_t), MPI_BYTE, rbuf.data,
              nclist, displ, MPI_BYTE, 0, io->comm);

  if (myrank == 0) {
    int asc = (io->info->options.output.iorformat == IO_RECORD_ASCII);
    int bin = (io->info->options.output.iorformat == IO_RECORD_BINARY);

    FILE * fp = util_fopen(filename, "w");

    if (fp == NULL) ifail = errno;
    if (ifail != 0) goto err;

    if (asc) ifail = colloid_io_ansi_write_ascii(fp, &rbuf);
    if (bin) ifail = colloid_io_ansi_write_binary(fp, &rbuf);

    /* Keep ifail from write function, and assume fclose() is ok... */
    fclose(fp);
  }

 err:

  colloid_array_free(&rbuf);
  colloid_array_free(&cbuf);

  free(displ);
  free(nclist);

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_write_ascii
 *
 *  Write nc states to file. Return number of bytes written.
 *
 *****************************************************************************/

int colloid_io_ansi_write_ascii(FILE * fp, const colloid_array_t * buf) {

  int ifail = 0;

  assert(fp);
  assert(buf);

  fprintf(fp, "%22d\n", buf->ntotal); /* HEADER */

  for (int n = 0; n < buf->ntotal; n++) {
    ifail += colloid_state_write_ascii(&buf->data[n], fp);
  }

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_write_binary
 *
 *  Write nc states to file.
 *
 *****************************************************************************/

int colloid_io_ansi_write_binary(FILE * fp, const colloid_array_t * buf) {

  size_t nw = 0; /* No. of items written */

  assert(fp);
  assert(buf);

  /* Header */
  nw = fwrite(&buf->ntotal, sizeof(int), 1, fp);
  if (nw != 1) goto err;

  /* Data */
  nw = fwrite(buf->data, sizeof(colloid_state_t), buf->ntotal, fp);
  if (nw != (size_t) buf->ntotal) goto err;

  return 0;

 err:
  return errno;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_read
 *
 *  This is the driver routine to read colloid information from file.
 *
 *  The read is a free-for-all in which all processes in the group
 *  read the entire file. Each colloid is then added or discarded
 *  based on its position.
 *
 *****************************************************************************/

int colloid_io_ansi_read(colloid_io_ansi_t * io, const char * filename) {

  assert(io);
  assert(filename);

  int    ifail = 0;
  FILE * fp    = util_fopen(filename, "r");

  if (fp == NULL) ifail = errno;
  if (ifail != 0) goto err;

  {
    int asc = (io->info->options.input.iorformat == IO_RECORD_ASCII);
    int bin = (io->info->options.input.iorformat == IO_RECORD_BINARY);

    if (asc) ifail = colloid_io_ansi_read_ascii(io, fp);
    if (bin) ifail = colloid_io_ansi_read_binary(io, fp);

    fclose(fp);
    if (ifail != 0) goto err; /* keep ifail from read function */
  }

  /* Check data for old-style "type" entries */
  /* This is scheduled for removal */
  {
    /* Check for old 'type' component */
    int    ntype         = 0;
    int    ntype_updates = colloids_type_check(io->info);
    pe_t * pe            = io->info->pe;
    MPI_Allreduce(&ntype_updates, &ntype, 1, MPI_INT, MPI_SUM, io->comm);
    if (ntype > 0) {
      pe_info(pe, "One or more colloids were updated as this looks\n");
      pe_info(pe, "like an old format input file. See note at\n");
      pe_info(pe, "https://ludwig.epcc.ed.ac.uk/outputs/colloid.html\n");
    }
    {
      /* Any ellipsoids must be checked. */
      int nbad = colloids_ellipsoid_abc_check(io->info);
      if (nbad > 0) {
        pe_exit(pe, "One or more ellipses in input file fail checks\n");
      }
    }
  }

 err:
  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_read_ascii
 *
 *****************************************************************************/

int colloid_io_ansi_read_ascii(colloid_io_ansi_t * io, FILE * fp) {

  int ifail = 0;
  int nread = 0; /* Header; number of colloids in this file */

  assert(io);
  assert(fp);

  int nr = fscanf(fp, "%22d\n", &nread);
  if (nr != 1) {
    goto err;
  }

  for (int n = 0; n < nread; n++) {
    colloid_state_t s     = {0};
    colloid_t *     pc    = NULL;
    int             ifail = colloid_state_read_ascii(&s, fp);
    if (ifail != 0) goto err;
    colloids_info_add_local(io->info, s.index, s.r, &pc);
    if (pc) pc->s = s;
  }

 err:

  ifail = colloid_io_ansi_check_read(io, nread);

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_read_binary
 *
 *  (colloid_io_ansi_t * io, const char * filename) better?
 *
 *****************************************************************************/

int colloid_io_ansi_read_binary(colloid_io_ansi_t * io, FILE * fp) {

  int ifail = 0;
  int nread = 0; /* Header: number of colloid entries in this file. */

  assert(io);
  assert(fp);

  size_t nr = fread(&nread, sizeof(int), 1, fp); /* HEADER */
  if (nr != 1) goto err;

  for (int n = 0; n < nread; n++) {
    colloid_state_t s     = {0};
    colloid_t *     pc    = NULL;
    int             ifail = colloid_state_read_binary(&s, fp);
    if (ifail != 0) goto err;

    colloids_info_add_local(io->info, s.index, s.r, &pc);
    if (pc) pc->s = s;
  }

  /* Collective, so cannot be by-passed. Fall through... */

 err:

  /* FIXME: This produces a message to stdout. Avoid the message? */
  ifail = colloid_io_ansi_check_read(io, nread);

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_ansi_check_read
 *
 *  Check the number of colloids in the list is consistent with that
 *  in the file (ngroup), and set the total.
 *
 *  If we haven't lost any particles, we can set the global total and
 *  proceed.
 *
 *****************************************************************************/

static int colloid_io_ansi_check_read(colloid_io_ansi_t * io, int ngroup) {

  int myrank = -1;
  int nlocal = 0;
  int ntotal = 0;

  assert(io);

  MPI_Comm_rank(io->comm, &myrank);

  colloids_info_nlocal(io->info, &nlocal);

  MPI_Reduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, 0, io->comm);

  if (myrank == 0) {
    if (ntotal != ngroup) {
      pe_t * pe = io->info->pe;
      pe_verbose(pe, "Colloid I/O group %d\n", io->subfile.index);
      pe_verbose(pe, "Colloids in file: %d Got %d\n", ngroup, ntotal);
      pe_fatal(pe,   "Total number of colloids not consistent with file\n");
    }
  }

  colloids_info_ntotal_set(io->info);
  colloids_info_ntotal(io->info, &ntotal);
  pe_info(io->info->pe, "Read a total of %d colloids from file\n", ntotal);

  return 0;
}
