/*****************************************************************************
 *
 *  test_build.c
 *
 *  Test colloid build process, and integrity of links.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2013-2017 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

#include "pe.h"
#include "coords.h"
#include "colloids_halo.h"
#include "colloid_sums.h"
#include "build.h"
#include "tests.h"

static int test_build_links_model_c1(pe_t * pe, cs_t * cs, int nvel,
				     double a0, double r0[3]);
static int test_build_links_model_c2(pe_t * pe, cs_t * cs, int nvel,
				     double a0, double r0[3]);
static int test_build_rebuild_c1(pe_t * pe, cs_t * cs, int nvel,
				 double a0, double r0[3]);

/*****************************************************************************
 *
 *  test_build_suite
 *
 *****************************************************************************/

int test_build_suite(void) {

  int nvel = NVEL;
  double a0;
  double r0[3];
  double delta = 1.0;        /* A small lattice offset */
  double lmin[3];
  double ltot[3];

  pe_t * pe = NULL;
  cs_t * cs = NULL;

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);
  cs_create(pe, &cs);
  cs_init(cs);

  cs_lmin(cs, lmin);
  cs_ltot(cs, ltot);

  a0 = 2.3;
  r0[X] = 0.5*ltot[X]; r0[Y] = 0.5*ltot[Y]; r0[Z] = 0.5*ltot[Z];
  test_build_links_model_c1(pe, cs, nvel, a0, r0);
  test_build_links_model_c2(pe, cs, nvel, a0, r0);
  test_build_rebuild_c1(pe, cs, nvel, a0, r0);

  a0 = 4.77;
  r0[X] = lmin[X] + delta; r0[Y] = 0.5*ltot[Y]; r0[Z] = 0.5*ltot[Z];
  test_build_links_model_c1(pe, cs, nvel, a0, r0);
  test_build_links_model_c2(pe, cs, nvel, a0, r0);
  test_build_rebuild_c1(pe, cs, nvel, a0, r0);

  a0 = 3.84;
  r0[X] = ltot[X]; r0[Y] = ltot[Y]; r0[Z] = ltot[Z];
  test_build_links_model_c1(pe, cs, nvel, a0, r0);
  test_build_links_model_c2(pe, cs, nvel, a0, r0);
  test_build_rebuild_c1(pe, cs, nvel, a0, r0);

  /* Some known cases: place the colloid in the centre and test only
   * in serial, as there is no quick way to compute in parallel. */

  cs_free(cs);
  pe_info(pe, "PASS     ./unit/test_build\n");
  pe_free(pe);

  return 0;
}

/*****************************************************************************
 *
 *  test_build_links_model_c1
 *
 *  Plave a colloid, a0 = 2.3, in the centre of the system, and
 *  examine the link integrity.
 *
 *  We expect:
 *     \sum_b w_b c_b_alpha = 0         for all alpha;
 *     \sum_b w_b r_b x c_b_alpha = 0   ditto.
 *
 *  independent of the model.
 *
 *  Owner tests the result.
 *
 *****************************************************************************/

static int test_build_links_model_c1(pe_t * pe, cs_t * cs, int nvel,
				     double a0, double r0[3]) {

  int ncolloid;
  int ncell[3] = {2, 2, 2};

  map_t * map = NULL;
  lb_model_t model = {};
  colloid_t * pc = NULL;
  colloids_info_t * cinfo = NULL;

  assert(pe);
  assert(cs);

  colloids_info_create(pe, cs, ncell, &cinfo);
  colloids_info_map_init(cinfo);
  map_create(pe, cs, 0, &map);
  lb_model_create(nvel, &model);

  /* Place the single colloid and construct the links */

  colloids_info_add_local(cinfo, 1, r0, &pc);
  if (pc) pc->s.a0 = a0;
  colloids_info_ntotal_set(cinfo);

  colloids_halo_state(cinfo);
  build_update_map(cs, cinfo, map);
  build_update_links(cs, cinfo, NULL, map, &model);

  colloid_sums_halo(cinfo, COLLOID_SUM_STRUCTURE);
  colloids_info_ntotal(cinfo, &ncolloid);
  assert(ncolloid == 1);

  /* The owner checks the details */

  if (pc) {

    /* DBL_EPSILON is too tight for these tests; FLT_EPSILON will catch
     * missing links which will be O(1). Could accumulate these sums as
     * integers to avoid this. */
    assert(fabs(pc->cbar[X] - 0.0) < FLT_EPSILON);
    assert(fabs(pc->cbar[Y] - 0.0) < FLT_EPSILON);
    assert(fabs(pc->cbar[Z] - 0.0) < FLT_EPSILON);
    assert(fabs(pc->rxcbar[X] - 0.0) < FLT_EPSILON);
    assert(fabs(pc->rxcbar[Y] - 0.0) < FLT_EPSILON);
    assert(fabs(pc->rxcbar[Z] - 0.0) < FLT_EPSILON);
    assert(fabs(pc->deltam - 0.0) < DBL_EPSILON);
  }

  lb_model_free(&model);
  map_free(map);
  colloids_info_free(cinfo);

  return 0;
}

/*****************************************************************************
 *
 * As for _c1() above, but everyone checks the result.
 *
 *****************************************************************************/

static int test_build_links_model_c2(pe_t * pe, cs_t * cs, int nvel,
				     double a0, double r0[3]) {

  int ic, jc, kc;
  int ncolloid;
  int ncell[3] = {4, 4, 4};

  map_t * map = NULL;
  lb_model_t model = {};
  colloid_t * pc = NULL;
  colloids_info_t * cinfo = NULL;

  assert(pe);
  assert(cs);

  colloids_info_create(pe, cs, ncell, &cinfo);
  colloids_info_map_init(cinfo);
  map_create(pe, cs, 0, &map);
  lb_model_create(nvel, &model);

  /* Place the single colloid and construct the links */

  colloids_info_add_local(cinfo, 1, r0, &pc);
  if (pc) pc->s.a0 = a0;
  colloids_info_ntotal_set(cinfo);

  colloids_halo_state(cinfo);
  build_update_map(cs, cinfo, map);
  build_update_links(cs, cinfo, NULL, map, &model);

  colloid_sums_halo(cinfo, COLLOID_SUM_STRUCTURE);
  colloids_info_ntotal(cinfo, &ncolloid);
  assert(ncolloid == 1);

  /* Everyone checks the details */

  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

	colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

	if (pc) {

	  /* DBL_EPSILON is too tight for these tests; FLT_EPSILON will catch
	   * missing links which will be O(1). Could accumulate these sums as
	   * integers to avoid this. */
	  assert(fabs(pc->cbar[X] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->cbar[Y] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->cbar[Z] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->rxcbar[X] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->rxcbar[Y] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->rxcbar[Z] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->deltam - 0.0) < DBL_EPSILON);
	}

      }
    }
  }

  lb_model_free(&model);
  map_free(map);
  colloids_info_free(cinfo);

  return 0;
}

/*****************************************************************************
 *
 *  test_build_rebuild_c1
 *
 *****************************************************************************/

static int test_build_rebuild_c1(pe_t * pe, cs_t * cs, int nvel,
				 double a0, double r0[3]) {

  int ic, jc, kc;
  int ncolloid;
  int ncell[3] = {2, 2, 2};

  map_t * map = NULL;
  lb_model_t model = {};
  colloid_t * pc = NULL;
  colloids_info_t * cinfo = NULL;

  assert(pe);
  assert(cs);

  colloids_info_create(pe, cs, ncell, &cinfo);
  colloids_info_map_init(cinfo);
  map_create(pe, cs, 0, &map);
  lb_model_create(nvel, &model);

  /* Place the single colloid and construct the links */

  colloids_info_add_local(cinfo, 1, r0, &pc);
  if (pc) {
    pc->s.a0 = a0;
    /* All copies must have this, applied later */
    pc->s.dr[X] = 0.0;
    pc->s.dr[Y] = 0.0;
    pc->s.dr[Z] = 0.0;
  }

  colloids_info_ntotal_set(cinfo);

  colloids_halo_state(cinfo);
  build_update_map(cs, cinfo, map);
  build_update_links(cs, cinfo, NULL, map, &model);

  colloids_info_ntotal(cinfo, &ncolloid);
  assert(ncolloid == 1);

  /* Move, rebuild, and check */

  colloids_info_position_update(cinfo);
  colloids_info_update_cell_list(cinfo);
  colloids_halo_state(cinfo);

  build_update_map(cs, cinfo, map);
  build_update_links(cs, cinfo, NULL, map, &model);
  colloid_sums_halo(cinfo, COLLOID_SUM_STRUCTURE);

  for (ic = 0; ic <= ncell[X] + 1; ic++) {
    for (jc = 0; jc <= ncell[Y] + 1; jc++) {
      for (kc = 0; kc <= ncell[Z] + 1; kc++) {

	colloids_info_cell_list_head(cinfo, ic, jc, kc, &pc);

	if (pc) {

	  /* DBL_EPSILON is too tight for these tests; FLT_EPSILON will catch
	   * missing links which will be O(1). Could accumulate these sums as
	   * integers to avoid this. */
	  assert(fabs(pc->cbar[X] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->cbar[Y] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->cbar[Z] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->rxcbar[X] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->rxcbar[Y] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->rxcbar[Z] - 0.0) < FLT_EPSILON);
	  assert(fabs(pc->deltam - 0.0) < DBL_EPSILON);
	}

      }
    }
  }

  lb_model_free(&model);
  map_free(map);
  colloids_info_free(cinfo);

  return 0;
}
