/*****************************************************************************
 *
 *  test_halo.c
 *
 *  This is a more rigourous test of the halo swap code for the
 *  distributions than appears in test model.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2017 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "pe.h"
#include "coords.h"
#include "model.h"
#include "control.h"
#include "tests.h"

int do_test_const_blocks(void);
int do_test_halo_null(pe_t * pe, cs_t * cs, lb_halo_enum_t halo);
int do_test_halo(pe_t * pe, cs_t * cs, int dim, const lb_halo_enum_t halo);

/*****************************************************************************
 *
 *  test_halo_suite
 *
 *****************************************************************************/

int test_halo_suite(void) {

  pe_t * pe = NULL;
  cs_t * cs = NULL;

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);
  cs_create(pe, &cs);
  cs_init(cs);

  do_test_const_blocks();

  do_test_halo_null(pe, cs, LB_HALO_FULL);
  do_test_halo_null(pe, cs, LB_HALO_REDUCED);

  do_test_halo(pe, cs, X, LB_HALO_FULL);
  do_test_halo(pe, cs, Y, LB_HALO_FULL);
  do_test_halo(pe, cs, Z, LB_HALO_FULL);

  if (pe_mpi_size(pe) == 1) {
    do_test_halo(pe, cs, X, LB_HALO_REDUCED);
    do_test_halo(pe, cs, Y, LB_HALO_REDUCED);
    do_test_halo(pe, cs, Z, LB_HALO_REDUCED);
  }


  pe_info(pe, "PASS     ./unit/test_halo\n");
  cs_free(cs);
  pe_free(pe);

  return 0;
}

int do_test_const_blocks(void) {

#ifdef TEST_TO_BE_REMOVED_WITH_GLOBAL_SYMBOLS
  int i, k;

  for (i = 0; i < CVXBLOCK; i++) {
    for (k = 0; k < xblocklen_cv[i]; k++) {
      test_assert(cv[xdisp_fwd_cv[i] + k][X] == +1);
      test_assert(cv[xdisp_bwd_cv[i] + k][X] == -1);
    }
  }

  for (i = 0; i < CVYBLOCK; i++) {
    for (k = 0; k < yblocklen_cv[i]; k++) {
      test_assert(cv[ydisp_fwd_cv[i] + k][Y] == +1);
      test_assert(cv[ydisp_bwd_cv[i] + k][Y] == -1);
    }
  }

  for (i = 0; i < CVZBLOCK; i++) {
    for (k = 0; k < zblocklen_cv[i]; k++) {
      test_assert(cv[zdisp_fwd_cv[i] + k][Z] == +1);
      test_assert(cv[zdisp_bwd_cv[i] + k][Z] == -1);
    }
  }
#endif
  return 0;
}

/*****************************************************************************
 *
 *  test_halo_null
 *
 *  Null halo test. Make sure no halo information appears in the
 *  domain proper. This works for both full and reduced halos.
 *
 *****************************************************************************/

int do_test_halo_null(pe_t * pe, cs_t * cs, lb_halo_enum_t halo) {

  int nlocal[3], n[3];
  int index, nd, p;
  int ndist = 2;
  int rank;
  int nhalo;
  int nextra;
  double f_actual;

  MPI_Comm comm = MPI_COMM_WORLD;
  lb_t * lb = NULL;

  assert(pe);
  assert(cs);

  cs_nhalo(cs, &nhalo);
  nextra = nhalo - 1;

  MPI_Comm_rank(comm, &rank);

  lb_create(pe, cs, &lb);
  lb_ndist_set(lb, ndist);
  lb_init(lb);
  lb_halo_set(lb, halo);

  cs_nlocal(cs, nlocal);

  /* Set entire distribution (all sites including halos) to 1.0 */

  for (n[X] = 1 - nextra; n[X] <= nlocal[X] + nextra; n[X]++) {
    for (n[Y] = 1 - nextra; n[Y] <= nlocal[Y] + nextra; n[Y]++) {
      for (n[Z] = 1 - nextra; n[Z] <= nlocal[Z] + nextra; n[Z]++) {

	index = cs_index(cs, n[X], n[Y], n[Z]);

	for (nd = 0; nd < ndist; nd++) {
	  for (p = 0; p < lb->model.nvel; p++) {
	    lb_f_set(lb, index, p, nd, 1.0);
	  }
	}

      }
    }
  }

  /* Zero interior */

  for (n[X] = 1; n[X] <= nlocal[X]; n[X]++) {
    for (n[Y] = 1; n[Y] <= nlocal[Y]; n[Y]++) {
      for (n[Z] = 1; n[Z] <= nlocal[Z]; n[Z]++) {

	index = cs_index(cs, n[X], n[Y], n[Z]);

	for (nd = 0; nd < ndist; nd++) {
	  for (p = 0; p < lb->model.nvel; p++) {
	    lb_f_set(lb, index, p, nd, 0.0);
	  }
	}

      }
    }
  }

  /* Swap */

  lb_halo(lb);

  /* Check everywhere in the interior still zero */

  for (n[X] = 1; n[X] <= nlocal[X]; n[X]++) {
    for (n[Y] = 1; n[Y] <= nlocal[Y]; n[Y]++) {
      for (n[Z] = 1; n[Z] <= nlocal[Z]; n[Z]++) {

	index = cs_index(cs, n[X], n[Y], n[Z]);

	for (nd = 0; nd < ndist; nd++) {
	  for (p = 0; p < lb->model.nvel; p++) {
	    lb_f(lb, index, p, nd, &f_actual);

	    /* everything should still be zero inside the lattice */
	    test_assert(fabs(f_actual - 0.0) < DBL_EPSILON);
	  }
	}

      }
    }
  }

  lb_free(lb);

  return 0;
}

/*****************************************************************************
 *
 *  do_test_halo
 *
 *  Test the halo swap for the distributions for coordinate direction dim.
 *
 *  Note that the reduced halo swaps are only meaningful in
 *  parallel. They will automatically work in serial.
 *
 *****************************************************************************/

int do_test_halo(pe_t * pe, cs_t * cs, int dim, lb_halo_enum_t halo) {

  int nhalo;
  int nlocal[3], n[3];
  int offset[3];
  int mpi_cartsz[3];
  int mpi_cartcoords[3];
  int nd;
  int ndist = 2;
  int nextra;
  int index, p, d;

  double ltot[3];
  double f_expect, f_actual;
  lb_t * lb = NULL;

  assert(pe);
  assert(cs);
  assert(dim == X || dim == Y || dim == Z);


  lb_create(pe, cs, &lb);
  lb_ndist_set(lb, ndist);
  lb_init(lb);
  lb_halo_set(lb, halo);

  cs_nhalo(cs, &nhalo);
  nextra = nhalo;

  cs_ltot(cs, ltot);
  cs_nlocal(cs, nlocal);
  cs_nlocal_offset(cs, offset);
  cs_cartsz(cs, mpi_cartsz);
  cs_cart_coords(cs, mpi_cartcoords);

  /* Zero entire distribution (all sites including halos) */

  for (n[X] = 1 - nextra; n[X] <= nlocal[X] + nextra; n[X]++) {
    for (n[Y] = 1 - nextra; n[Y] <= nlocal[Y] + nextra; n[Y]++) {
      for (n[Z] = 1 - nextra; n[Z] <= nlocal[Z] + nextra; n[Z]++) {

	index = cs_index(cs, n[X], n[Y], n[Z]);

	for (nd = 0; nd < ndist; nd++) {
	  for (p = 0; p < lb->model.nvel; p++) {
	    lb_f_set(lb, index, p, nd, -1.0);
	  }
	}

      }
    }
  }

  /* Set the interior sites to get swapped with a value related to
   * absolute position */

  for (n[X] = 1; n[X] <= nlocal[X]; n[X]++) {
    for (n[Y] = 1; n[Y] <= nlocal[Y]; n[Y]++) {
      for (n[Z] = 1; n[Z] <= nlocal[Z]; n[Z]++) {

	index = cs_index(cs, n[X], n[Y], n[Z]);

	if (n[X] <= nhalo || n[X] > nlocal[X] - nhalo ||
	    n[Y] <= nhalo || n[Y] > nlocal[Y] - nhalo ||
	    n[Z] <= nhalo || n[Z] > nlocal[Z] - nhalo) {

	  for (nd = 0; nd < ndist; nd++) {
	    for (p = 0; p < lb->model.nvel; p++) {
	      lb_f_set(lb, index, p, nd, 1.0*(offset[dim] + n[dim]));
	    }
	  }
	}

      }
    }
  }

  lb_memcpy(lb, tdpMemcpyHostToDevice);
  lb_halo(lb);
  lb_memcpy(lb, tdpMemcpyDeviceToHost);

  /* Check the results (all sites for distribution halo).
   * The halo regions should contain a copy of the above, while the
   * interior sites are unchanged */

  /* Note the distribution halo swaps are always width 1, irrespective
   * of nhalo */

  for (n[X] = 0; n[X] <= nlocal[X] + 1; n[X]++) {
    for (n[Y] = 0; n[Y] <= nlocal[Y] + 1; n[Y]++) {
      for (n[Z] = 0; n[Z] <= nlocal[Z] + 1; n[Z]++) {

	index = cs_index(cs, n[X], n[Y], n[Z]);

	for (nd = 0; nd < ndist; nd++) {
	  for (d = 0; d < 3; d++) {

	    /* 'Left' side */
	    if (dim == d && n[d] == 0) {

	      f_expect = offset[dim];
	      if (mpi_cartcoords[dim] == 0) f_expect = ltot[dim];

	      for (p = 0; p < lb->model.nvel; p++) {
		lb_f(lb, index, p, nd, &f_actual);
		test_assert(fabs(f_actual-f_expect) < DBL_EPSILON);
	      }
	    }

	    /* 'Right' side */
	    if (dim == d && n[d] == nlocal[d] + 1) {

	      f_expect = offset[dim] + nlocal[dim] + 1.0;
	      if (mpi_cartcoords[dim] == mpi_cartsz[dim] - 1) f_expect = 1.0;

	      for (p = 0; p < lb->model.nvel; p++) {
		lb_f(lb, index, p, nd, &f_actual);
		test_assert(fabs(f_actual-f_expect) < DBL_EPSILON);
	      }
	    }
	  }
	}
	/* Next site */
      }
    }
  }

  lb_free(lb);

  return 0;
}
