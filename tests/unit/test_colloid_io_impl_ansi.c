/*****************************************************************************
 *
 *  test_colloid_io_impl_ansi.c
 *
 *****************************************************************************/

#include <assert.h>

#include "colloid_io_impl_ansi.h"

int test_colloid_array_initialise(void);

int test_colloid_io_ansi_initialise(pe_t * pe);
int test_colloid_io_ansi_create(pe_t * pe);

int test_colloid_io_ansi_read(pe_t * pe);
int test_colloid_io_ansi_write(pe_t * pe);

/*****************************************************************************
 *
 *  test_colloid_io_impl_ansi_suite
 *
 *****************************************************************************/

int test_colloid_io_impl_ansi_suite(void) {

  int ifail = 0;
  pe_t * pe = NULL;

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);

  /* If struct changes, the tests need updating... */
  assert(sizeof(colloid_io_ansi_t) == 80);

  test_colloid_array_initialise();

  test_colloid_io_ansi_initialise(pe);
  test_colloid_io_ansi_create(pe);

  test_colloid_io_ansi_read(pe);
  test_colloid_io_ansi_write(pe);

  pe_info(pe, "%-9s %s\n", "PASS", __FILE__);
  pe_free(pe);

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_array_initialise
 *
 *****************************************************************************/

int test_colloid_array_initialise(void) {

  int ifail = 0;

  {
    int ntotal = 13;
    colloid_array_t array = {0};

    ifail = colloid_array_initialise(ntotal, &array);
    assert(ifail == 0);
    assert(array.ntotal == ntotal);
    assert(array.data);
    assert(array.data[0].index == 0);

    ifail = colloid_array_finalise(&array);
    assert(ifail == 0);
    assert(array.ntotal == 0);
    assert(array.data   == NULL);
  }

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_io_ansi_initialise
 *
 *****************************************************************************/

int test_colloid_io_ansi_initialise(pe_t * pe) {

  int ifail = 0;
  cs_t * cs = NULL;

  cs_create(pe, &cs);
  cs_init(cs);

  {
    int ncell[3] = {8, 8, 8};
    colloids_info_t * info = NULL;
    colloid_io_ansi_t io = {0};

    colloids_info_create(pe, cs, ncell, &info);
    ifail = colloid_io_ansi_initialise(info, &io);
    assert(ifail == 0);

    /* subfile exsits ... */
    assert(io.subfile.nfile == 1);
    assert(io.subfile.index == 0);

    /* new communicator exsits */
    assert(io.comm != MPI_COMM_NULL);
    assert(io.comm != cs->commcart);

    colloid_io_ansi_finalise(&io);
    colloids_info_free(info);
    assert(io.comm == MPI_COMM_NULL);
  }

  cs_free(cs);

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_io_ansi_create
 *
 *****************************************************************************/

int test_colloid_io_ansi_create(pe_t * pe) {

  int ifail = 0;
  cs_t * cs = NULL;

  cs_create(pe, &cs);
  cs_init(cs);

  {
    int ncell[3] = {8, 8, 8};
    colloids_info_t * info = NULL;
    colloid_io_ansi_t * io = NULL;

    colloids_info_create(pe, cs, ncell, &info);
    ifail = colloid_io_ansi_create(info, &io);
    assert(ifail == 0);
    assert(io->comm != MPI_COMM_NULL);

    colloid_io_ansi_free(&io);
    colloids_info_free(info);
    assert(io == NULL);
  }

  cs_free(cs);

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_io_ansi_read
 *
 *****************************************************************************/

int test_colloid_io_ansi_read(pe_t * pe) {

  int ifail = 0;
  cs_t * cs = NULL;

  cs_create(pe, &cs);
  cs_init(cs);

  /* ASCII; from existing file */
  {
    int ncell[3] = {3, 3, 3};
    colloids_info_t * info = NULL;
    colloid_io_ansi_t * io = NULL;

    colloids_info_create(pe, cs, ncell, &info);
    ifail = colloid_io_ansi_create(info, &io);

    ifail = colloid_io_ansi_read(io, "colloid-ansi-ascii.001-001");
    ifail = colloid_io_ansi_write(io, "colloid-tmp.dat");

    colloid_io_ansi_free(&io);
    colloids_info_free(info);
  }

  /* Binary; requires an option to control format */
  {
    printf("Binary read pending\n");
    assert(0);
  }

  cs_free(cs);

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_Io_ansi_write
 *
 *****************************************************************************/

int test_colloid_io_ansi_write(pe_t * pe) {

  int ifail = 0;

  printf("ansi write pending\n");

  return ifail;
}
