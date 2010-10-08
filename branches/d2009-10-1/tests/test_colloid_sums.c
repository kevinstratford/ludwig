/*****************************************************************************
 *
 *  test_colloid_sums.c
 *
 *  Test of the various sum routines.
 *
 *  $Id: test_colloid_sums.c,v 1.1.2.3 2010-09-23 17:12:02 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010 The University of Edinburgh
 *
 *****************************************************************************/

#include <math.h>

#include "pe.h"
#include "coords.h"
#include "colloids.h"
#include "colloids_halo.h"
#include "colloid_sums.h"
#include "tests.h"

static int dim_; /* Current direction */

static void test_colloid_sums_1d(void);
static void test_colloid_sums_reference_set(colloid_t * cref, int seed);
static void test_colloid_sums_copy(colloid_t ref, colloid_t * pc);
static void test_colloid_sums_assert(colloid_t c1, colloid_t * c2);
static void test_colloid_sums_edge(const int ncell[3], const double r0[3]);

/*****************************************************************************
 *
 *  main
 *
 *****************************************************************************/

int main(int argc, char ** argv) {

  pe_init(argc, argv);
  info("Checking colloid sum messages\n");

  test_colloid_sums_1d();

  info("Tests complete\n");
  pe_finalise();

  return 0;
}

/*****************************************************************************
 *
 *  test_colloid_sums_1d
 *
 *  Place one colloid locally at the centre of the domain.
 *  In each direction in turn, we put it close to the (periodic) boundary.
 *
 *****************************************************************************/

static void test_colloid_sums_1d(void) {

  int ntotal[3] = {1024, 512, 256};
  int nlocal[3];
  int ncell[3];

  double r0[3];

  coords_ntotal_set(ntotal);
  coords_init();

  /* We use cells > 2 to prevent copies in directions other than dim_
   * when a call to colloids_halo_state() is made. */

  coords_nlocal(nlocal);

  dim_ = X;
  ncell[X] = 4;
  ncell[Y] = 3;
  ncell[Z] = 3;

  r0[X] = Lmin(X) + 0.5;
  r0[Y] = Lmin(Y) + 0.5*nlocal[Y];
  r0[Z] = Lmin(Z) + 0.5*nlocal[Z];

  test_colloid_sums_edge(ncell, r0);

  ncell[X] = 2;
  ncell[Y] = 4;
  ncell[Z] = 3;

  test_colloid_sums_edge(ncell, r0);

  dim_ = Y;
  r0[X] = Lmin(X) + 0.5*nlocal[X];
  r0[Y] = Lmin(Y) + 0.5;
  r0[Z] = Lmin(Z) + 0.5*nlocal[Z];

  test_colloid_sums_edge(ncell, r0);

  dim_ = Z;
  r0[X] = Lmin(X) + 0.5*nlocal[X];
  r0[Y] = Lmin(Y) + 0.5*nlocal[Y];
  r0[Z] = Lmin(Z) + 0.5;

  test_colloid_sums_edge(ncell, r0);

  coords_finish();

  return;
}

/*****************************************************************************
 *
 *  test_colloid_sums_edge
 *
 *  Place a single particle at r0 to test the communication.
 *
 *****************************************************************************/

static void test_colloid_sums_edge(const int ncell[3], const double r0[3]) {

  int index;
  int ic, jc, kc;

  colloid_t * pc;
  colloid_t   cref1;   /* All ranks get the same reference colloids */
  colloid_t   cref2;

  test_colloid_sums_reference_set(&cref1, 1);
  test_colloid_sums_reference_set(&cref2, 2);
  colloids_cell_ncell_set(ncell);

  colloids_init();

  /* This must work in parallel to initialise only a single particle
   * which only gets swapped in the x-direction. */

  index = 1;
  pc = colloid_add_local(index, r0);
  if (pc) {
    test_colloid_sums_copy(cref1, pc);
  }

  index = 2;
  pc = colloid_add_local(index, r0);
  if (pc) {
    test_colloid_sums_copy(cref2, pc);
  }

  colloids_halo_state();
  colloid_sums_dim(X, COLLOID_SUM_STRUCTURE);
  colloid_sums_dim(X, COLLOID_SUM_DYNAMICS);
  colloid_sums_dim(X, COLLOID_SUM_ACTIVE);

  if (dim_ == Y || dim_ == Z) {
    colloid_sums_dim(Y, COLLOID_SUM_STRUCTURE);
    colloid_sums_dim(Y, COLLOID_SUM_DYNAMICS);
    colloid_sums_dim(Y, COLLOID_SUM_ACTIVE);
  }

  if (dim_ == Z) {
    colloid_sums_dim(Z, COLLOID_SUM_STRUCTURE);
    colloid_sums_dim(Z, COLLOID_SUM_DYNAMICS);
    colloid_sums_dim(Z, COLLOID_SUM_ACTIVE);
  }

  /* Everywhere check colloid index = 1 has the correct sum */

  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

	pc = colloids_cell_list(ic, jc, kc);

	if (pc) {
	  /* Check the totals */
	  if (pc->s.index == 1) test_colloid_sums_assert(cref1, pc);
	  if (pc->s.index == 2) test_colloid_sums_assert(cref2, pc);
	}
	/* Next cell */
      }
    }
  }

  /* Finish */

  colloids_finish();

  return;
}

/*****************************************************************************
 *
 *  test_colloid_sums_reference_set
 *
 *****************************************************************************/

static void test_colloid_sums_reference_set(colloid_t * pc, int seed) {

  int ia;
  int ivalue;

  ivalue = seed;

  /* STURCTURE */
  /* Note, we haven't included s.deltaphi as this is part of the
   * state, which is involved in the halo swap, as well as the sum. */

  pc->sumw = 1.0*ivalue++;
  pc->deltam = 1.0*ivalue++;

  for (ia = 0; ia < 3; ia++) {
    pc->cbar[ia] = 1.0*ivalue++;
    pc->rxcbar[ia] = 1.0*ivalue++;
  }
 
  /* DYNAMICS */

  for (ia = 0; ia < 3; ia++) {
    pc->f0[ia] = 1.0*ivalue++;
    pc->t0[ia] = 1.0*ivalue++;
    pc->force[ia] = 1.0*ivalue++;
    pc->torque[ia] = 1.0*ivalue++;
  }

  pc->sump = 1.0*ivalue++;

  for (ia = 0; ia < 21; ia++) {
    pc->zeta[ia] = 1.0*ivalue++;
  }

  /* ACTIVE */

  for (ia = 0; ia < 3; ia++) {
    pc->fc0[ia] = 1.0*ivalue++;
    pc->tc0[ia] = 1.0*ivalue++;
  }

  return;
}

/*****************************************************************************
 *
 *  test_colloid_sums_copy
 *
 *  Copy the sum information from reference to pc.
 *
 *****************************************************************************/

static void test_colloid_sums_copy(colloid_t ref, colloid_t * pc) {

  int ia;

  pc->sumw = ref.sumw;
  pc->sump = ref.sump;

  for (ia = 0; ia < 3; ia++) {
    pc->cbar[ia] = ref.cbar[ia];
    pc->rxcbar[ia] = ref.rxcbar[ia];
    pc->f0[ia] = ref.f0[ia];
    pc->t0[ia] = ref.t0[ia];
    pc->force[ia] = ref.force[ia];
    pc->torque[ia] = ref.torque[ia];
    pc->fc0[ia] = ref.fc0[ia];
    pc->tc0[ia] = ref.tc0[ia];
  }

  for (ia = 0; ia < 21; ia++) {
    pc->zeta[ia] = ref.zeta[ia];
  }

  return;
}

/*****************************************************************************
 *
 *  test_colloid_sums_assert
 *
 *  Assert that two colloids hold the same summed quantities.
 *
 *****************************************************************************/

static void test_colloid_sums_assert(colloid_t c1, colloid_t * c2) {

  int ia;

  /* STRUCTURE */

  test_assert(fabs(c1.sumw - c2->sumw) < TEST_DOUBLE_TOLERANCE);

  for (ia = 0; ia < 3; ia++) {
    test_assert(fabs(c1.cbar[ia] - c2->cbar[ia]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(c1.rxcbar[ia] - c2->rxcbar[ia]) < TEST_DOUBLE_TOLERANCE);
  }

  /* DYNAMICS */

  test_assert(fabs(c1.sump - c2->sump) < TEST_DOUBLE_TOLERANCE);

  for (ia = 0; ia < 3; ia++) {
    test_assert(fabs(c1.f0[ia] - c2->f0[ia]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(c1.t0[ia] - c2->t0[ia]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(c1.force[ia] - c2->force[ia]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(c1.torque[ia] - c2->torque[ia]) < TEST_DOUBLE_TOLERANCE);
  }

  for (ia = 0; ia < 21; ia++) {
    test_assert(fabs(c1.zeta[ia] - c2->zeta[ia]) < TEST_DOUBLE_TOLERANCE);
  }

  /* ACTIVE */

  for (ia = 0; ia < 3; ia++) {
    test_assert(fabs(c1.fc0[ia] - c2->fc0[ia]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(c1.tc0[ia] - c2->tc0[ia]) < TEST_DOUBLE_TOLERANCE);
  }

  return;
}