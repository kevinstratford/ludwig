/*****************************************************************************
 *
 *  test_colloids_halo.c
 *
 *  Halo swap test.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2017 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <float.h>

#include "pe.h"
#include "coords.h"
#include "colloids.h"
#include "colloids_halo.h"
#include "tests.h"

int test_colloids_halo111(pe_t * pe, cs_t * cs);
int test_colloids_halo211(pe_t * pe, cs_t * cs);
int test_colloids_halo_repeat(pe_t * pe, cs_t * cs);
static void test_position(cs_t * cs, const double r1[3], const double r2[3]);

/*****************************************************************************
 *
 *  test_colloids_halo_suite
 *
 *****************************************************************************/

int test_colloids_halo_suite(void) {

  int ntotal[3] = {1024, 1024, 1024};
  pe_t * pe = NULL;
  cs_t * cs = NULL;

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);
  cs_create(pe, &cs);
  cs_ntotal_set(cs, ntotal);
  cs_init(cs);

  test_colloids_halo111(pe, cs);
  test_colloids_halo211(pe, cs);
  test_colloids_halo_repeat(pe, cs);

  pe_info(pe, "PASS     ./unit/test_colloids_halo\n");
  cs_free(cs);
  pe_free(pe);

  return 0;
}

/*****************************************************************************
 *
 *  test_colloids_halo111
 *
 *  Please one colloid in the bottom left corner locally and swap.
 *  So the local cell is {1, 1, 1}.
 *
 *****************************************************************************/

int test_colloids_halo111(pe_t * pe, cs_t * cs) {

  int ncell[3] = {2, 2, 2};
  int ntotal[3];
  int noffset[3];
  int mpi_cartsz[3];
  int ncount[2];
  int index;
  int ncolloid;
  double r0[3];
  double r1[3];
  double lmin[3];

  colloid_t * pc;
  colloids_info_t * cinfo = NULL;
  colloid_halo_t * halo = NULL;

  assert(pe);
  assert(cs);

  cs_lmin(cs, lmin);
  cs_ntotal(cs, ntotal);
  cs_cartsz(cs, mpi_cartsz);

  colloids_info_create(pe, cs, ncell, &cinfo);
  assert(cinfo);

  colloids_halo_create(cinfo, &halo);
  assert(halo);

  cs_nlocal_offset(cs, noffset);

  r0[X] = lmin[X] + 1.0*(noffset[X] + 1);
  r0[Y] = lmin[Y] + 1.0*(noffset[Y] + 1);
  r0[Z] = lmin[Z] + 1.0*(noffset[Z] + 1);

  index = 1 + pe_mpi_rank(pe);
  colloids_info_add_local(cinfo, index, r0, &pc);
  assert(pc);

  colloids_halo_send_count(halo, X, ncount);
  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 1);

  colloids_halo_send_count(halo, Y, ncount);
  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 1);

  colloids_halo_send_count(halo, Z, ncount);
  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 1);

  colloids_halo_dim(halo, X);

  /* All process should now have one particle in upper x halo region,
   * and the send count in Y (back) should be 2 */

  colloids_info_cell_count(cinfo, ncell[X] + 1, 1, 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nlocal(cinfo, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nallocated(cinfo, &ncolloid);
  test_assert(ncolloid == 2);

  colloids_info_cell_list_head(cinfo, ncell[X] + 1, 1, 1, &pc);
  test_assert(pc != NULL);


  r1[X] = r0[X] + 1.0*ntotal[X]/mpi_cartsz[X];
  r1[Y] = r0[Y];
  r1[Z] = r0[Z];

  test_position(cs, r1, pc->s.r);

  colloids_halo_send_count(halo, Y, ncount);

  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 2);

  colloids_halo_dim(halo, Y);

  /* There should now be two additional images, and the send count
   * goes up to 4 */

  colloids_info_cell_count(cinfo, 1, ncell[Y] + 1, 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_cell_count(cinfo, ncell[X] + 1, ncell[Y] + 1, 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nlocal(cinfo, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nallocated(cinfo, &ncolloid);
  test_assert(ncolloid == 4);

  colloids_info_cell_list_head(cinfo, 1, ncell[Y] + 1, 1, &pc);
  test_assert(pc != NULL);


  r1[X] = r0[X];
  r1[Y] = r0[Y] + 1.0*ntotal[Y]/mpi_cartsz[Y];
  r1[Z] = r0[Z];

  test_position(cs, r1, pc->s.r);

  colloids_halo_send_count(halo, Z, ncount);
  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 4);

  colloids_halo_dim(halo, Z);

  /* We should end up with eight in total, seven of which are
   * periodic images. */

  colloids_info_cell_count(cinfo, 1, 1, ncell[Z] + 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_cell_count(cinfo, ncell[X] + 1, 1, ncell[Z] + 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_cell_count(cinfo, 1, ncell[Y] + 1, ncell[Z] + 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_cell_count(cinfo, ncell[X] + 1, ncell[Y] + 1, ncell[Z] + 1,
			   &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nlocal(cinfo, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nallocated(cinfo, &ncolloid);
  test_assert(ncolloid == 8);

  colloids_halo_free(halo);
  colloids_info_free(cinfo);

  return 0;
}

/*****************************************************************************
 *
 *  test_colloids_halo211
 *
 *  The local cell of the real particle is {2, 1, 1}.
 *
 *****************************************************************************/

int test_colloids_halo211(pe_t * pe, cs_t * cs) {

  int ncell[3] = {2, 2, 2};
  int ntotal[3];
  int noffset[3];
  int mpi_cartsz[3];
  int ncount[2];
  int index;
  int ncolloid;
  double r0[3];
  double r1[3];
  double lmin[3];
  double lcell[3];

  colloid_t * pc = NULL;
  colloid_halo_t * halo = NULL;
  colloids_info_t * cinfo = NULL;

  assert(pe);
  assert(cs);

  cs_lmin(cs, lmin);
  cs_ntotal(cs, ntotal);
  cs_cartsz(cs, mpi_cartsz);

  colloids_info_create(pe, cs, ncell, &cinfo);
  assert(cinfo);

  colloids_halo_create(cinfo, &halo);
  assert(halo);

  cs_nlocal_offset(cs, noffset);
  colloids_info_lcell(cinfo, lcell);

  r0[X] = lmin[X] + 1.0*noffset[X] + lcell[X];
  r0[Y] = lmin[Y] + 1.0*(noffset[Y] + 1);
  r0[Z] = lmin[Z] + 1.0*(noffset[Z] + 1);

  index = 1 + pe_mpi_rank(pe);
  colloids_info_add_local(cinfo, index, r0, &pc);
  assert(pc);

  colloids_halo_send_count(halo, X, ncount);
  test_assert(ncount[FORWARD] == 1);
  test_assert(ncount[BACKWARD] == 0);

  colloids_halo_send_count(halo, Y, ncount);
  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 1);

  colloids_halo_send_count(halo, Z, ncount);
  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 1);

  colloids_halo_dim(halo, X);

  /* All process should now have one particle in lower x halo region,
   * and the send count in Y (back) should be 2 */

  colloids_info_cell_count(cinfo, 0, 1, 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nlocal(cinfo, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nallocated(cinfo, &ncolloid);
  test_assert(ncolloid == 2);

  colloids_info_cell_list_head(cinfo, 0, 1, 1, &pc);
  test_assert(pc != NULL);

  r1[X] = r0[X] - 1.0*ntotal[X]/mpi_cartsz[X];
  r1[Y] = r0[Y];
  r1[Z] = r0[Z];
  test_position(cs, r1, pc->s.r);

  colloids_halo_send_count(halo, Y, ncount);
  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 2);

  colloids_halo_dim(halo, Y);

  /* There should now be two additional images, and the send count
   * goes up to 4 */

  colloids_info_cell_count(cinfo, 2, ncell[Y] + 1, 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_cell_count(cinfo, 0, ncell[Y] + 1, 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nlocal(cinfo, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nallocated(cinfo, &ncolloid);
  test_assert(ncolloid == 4);

  colloids_info_cell_list_head(cinfo, 2, ncell[Y] + 1, 1, &pc);
  test_assert(pc != NULL);

  r1[X] = r0[X];
  r1[Y] = r0[Y] + 1.0*ntotal[Y]/mpi_cartsz[Y];
  r1[Z] = r0[Z];
  test_position(cs, r1, pc->s.r);

  colloids_halo_send_count(halo, Z, ncount);
  test_assert(ncount[FORWARD] == 0);
  test_assert(ncount[BACKWARD] == 4);

  colloids_halo_dim(halo, Z);

  /* We should end up with eight in total (locally) */

  colloids_info_cell_count(cinfo, 2, 1, ncell[Z] + 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_cell_count(cinfo, 0, 1, ncell[Z] + 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_cell_count(cinfo, 2, ncell[Y] + 1, ncell[Z] + 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_cell_count(cinfo, 0, ncell[Y] + 1, ncell[Z] + 1, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nlocal(cinfo, &ncolloid);
  test_assert(ncolloid == 1);

  colloids_info_nallocated(cinfo, &ncolloid);
  test_assert(ncolloid == 8);

  colloids_halo_free(halo);
  colloids_info_free(cinfo);

  return 0;
}

/*****************************************************************************
 *
 *  test_colloids_halo_repeat
 *
 *  Make sure repeat halo swap doesn't multiply particles.
 *
 *****************************************************************************/

int test_colloids_halo_repeat(pe_t * pe, cs_t * cs) {

  int ncell[3] = {2, 2, 2};
  int noffset[3];
  int index;
  int ncolloid;
  double r0[3];
  double lmin[3];

  colloid_t * pc = NULL;
  colloids_info_t * cinfo = NULL;

  assert(pe);
  assert(cs);

  cs_lmin(cs, lmin);

  colloids_info_create(pe, cs, ncell, &cinfo);
  assert(cinfo);

  cs_nlocal_offset(cs, noffset);

  r0[X] = lmin[X] + 1.0*(noffset[X] + 1);
  r0[Y] = lmin[Y] + 1.0*(noffset[Y] + 1);
  r0[Z] = lmin[Z] + 1.0*(noffset[Z] + 1);

  index = 1 + pe_mpi_rank(pe);
  colloids_info_add_local(cinfo, index, r0, &pc);
  assert(pc);

  index = 1 + pe_mpi_size(pe) + pe_mpi_rank(pe);
  colloids_info_add_local(cinfo, index, r0, &pc);
  index = 1 + 2*pe_mpi_size(pe) + pe_mpi_rank(pe);
  colloids_info_add_local(cinfo, index, r0, &pc);

  colloids_halo_state(cinfo);
  colloids_halo_state(cinfo);

  colloids_info_nlocal(cinfo, &ncolloid);
  test_assert(ncolloid == 3);

  colloids_info_nallocated(cinfo, &ncolloid);
  test_assert(ncolloid == 24);

  colloids_info_free(cinfo);

  return 0;
}

/*****************************************************************************
 *
 *  test_position
 *
 *  The periodic halo swap contains a factor (1.0 - DBL_EPSILON)*L
 *  which must be allowed for in the tolerance, hence the slightly
 *  odd-looking formulation below.
 *
 *****************************************************************************/

void test_position(cs_t * cs, const double r1[3], const double r2[3]) {

  double tolerance;
  double ltot[3];

  assert(cs);
  cs_ltot(cs, ltot);

  tolerance = (1.0 + DBL_EPSILON)*DBL_EPSILON*ltot[X];
  test_assert(fabs(r1[X] - r2[X]) < tolerance);
  tolerance = (1.0 + DBL_EPSILON)*DBL_EPSILON*ltot[Y];
  test_assert(fabs(r1[Y] - r2[Y]) < tolerance);
  tolerance = (1.0 + DBL_EPSILON)*DBL_EPSILON*ltot[Z];
  test_assert(fabs(r1[Z] - r2[Z]) < tolerance);

  return;
}
