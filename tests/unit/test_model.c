/*****************************************************************************
 *
 *  test_model.c
 *
 *  Unit test for the currently compiled model (D3Q15 or D3Q19).
 *
 *  Edinburgh Soft Matter and Statistical Physics Group
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2021 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "pe.h"
#include "coords.h"
#include "util.h"
#include "lb_model_s.h"
#include "tests.h"

static void test_model_constants(void);
static void test_model_velocity_set(void);

int do_test_model_distributions(pe_t * pe, cs_t * cs);
int do_test_model_halo_swap(pe_t * pe, cs_t * cs);
int do_test_model_reduced_halo_swap(pe_t * pe, cs_t * cs);
int do_test_lb_model_io(pe_t * pe, cs_t * cs);
int do_test_d3q19_ghosts(void);
static  int test_model_is_domain(cs_t * cs, int ic, int jc, int kc);

/*****************************************************************************
 *
 *  test_model_suite
 *
 *****************************************************************************/

int test_model_suite(void) {

  pe_t * pe = NULL;
  cs_t * cs = NULL;

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);
  cs_create(pe, &cs);
  cs_init(cs);

  /* Test model structure (coordinate-independent stuff) */

  test_model_constants();
  test_model_velocity_set();

  /* Now test actual distributions */

  do_test_model_distributions(pe, cs);
  do_test_model_halo_swap(pe, cs);
  if (DATA_MODEL == DATA_MODEL_AOS && NSIMDVL == 1) {
    do_test_model_reduced_halo_swap(pe, cs);
  }
  do_test_lb_model_io(pe, cs);
  do_test_d3q19_ghosts();

  pe_info(pe, "PASS     ./unit/test_model\n");
  cs_free(cs);
  pe_free(pe);

  return 0;
}

/*****************************************************************************
 *
 *  test_model_constants
 *
 *  Check the various constants associated with the reduced halo swap.
 *
 *****************************************************************************/

static void test_model_constants(void) {

  int i, k, p;

  for (i = 0; i < CVXBLOCK; i++) {
    for (k = 0; k < xblocklen_cv[i]; k++) {
      p = xdisp_fwd_cv[i] + k;
      test_assert(p >= 0 && p < NVEL);
      test_assert(cv[p][X] == +1);
      p = xdisp_bwd_cv[i] + k;
      test_assert(p >= 0 && p < NVEL);
      test_assert(cv[p][X] == -1);
    }
  }

  for (i = 0; i < CVYBLOCK; i++) {
    for (k = 0; k < yblocklen_cv[i]; k++) {
      p = ydisp_fwd_cv[i] + k;
      test_assert(p >= 0 && p < NVEL);
      test_assert(cv[p][Y] == +1);
      p = ydisp_bwd_cv[i] + k;
      test_assert(p >= 0 && p < NVEL);
      test_assert(cv[p][Y] == -1);
    }
  }

  for (i = 0; i < CVZBLOCK; i++) {
    for (k = 0; k < zblocklen_cv[i]; k++) {
      p = zdisp_fwd_cv[i] + k;
      test_assert(p >= 0 && p < NVEL);
      test_assert(cv[p][Z] == +1);
      p = zdisp_bwd_cv[i] + k;
      test_assert(p >= 0 && p < NVEL);
      test_assert(cv[p][Z] == -1);
    }
  }

  return;
}

/*****************************************************************************
 *
 *  test_model_velocity_set
 *
 *  Check the velocities, kinetic projector, tables of eigenvectors
 *  etc etc are all consistent for the current model.
 *
 *****************************************************************************/

static void test_model_velocity_set(void) {

  int i, j, k, p;
  double dij;
  double sum, sumx, sumy, sumz;

  LB_CS2_DOUBLE(cs2);
  LB_RCS2_DOUBLE(rcs2);
  KRONECKER_DELTA_CHAR(d_);

  /* Number of hydrodynamic modes: */

  test_assert(NHYDRO == (1 + NDIM + NDIM*(NDIM+1)/2));

  /* Speed of sound */

  test_assert(fabs(rcs2 - 3.0) < TEST_DOUBLE_TOLERANCE);

  /* Kronecker delta */

  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      if (i == j) {
	test_assert(fabs(d_[i][j] - 1.0) < TEST_DOUBLE_TOLERANCE);
      }
      else {
	test_assert(fabs(d_[i][j] - 0.0) < TEST_DOUBLE_TOLERANCE);
      }
    }
  }

  /* Checking cv[0] */

  test_assert(cv[0][X] == 0);
  test_assert(cv[0][Y] == 0);
  test_assert(cv[0][Z] == 0);

  /* Checking cv[p][X] = -cv[NVEL-p][X] (p != 0)  */

  for (p = 1; p < NVEL; p++) {
    test_assert(cv[p][X] == -cv[NVEL-p][X]);
    test_assert(cv[p][Y] == -cv[NVEL-p][Y]);
    test_assert(cv[p][Z] == -cv[NVEL-p][Z]);
  }

  /* Sum of quadrature weights, velcoities */

  sum = 0.0; sumx = 0.0; sumy = 0.0; sumz = 0.0;

  for (p = 0; p < NVEL; p++) {
    sum += wv[p];
    sumx += wv[p]*cv[p][X];
    sumy += wv[p]*cv[p][Y];
    sumz += wv[p]*cv[p][Z];
  }

  test_assert(fabs(sum - 1.0) < TEST_DOUBLE_TOLERANCE);
  test_assert(fabs(sumx) < TEST_DOUBLE_TOLERANCE);
  test_assert(fabs(sumy) < TEST_DOUBLE_TOLERANCE);
  test_assert(fabs(sumz) < TEST_DOUBLE_TOLERANCE);

  /* Quadratic terms = cs^2 d_ij */

  /* Checking wv[p]*cv[p][i]*cv[p][j] */

  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      sum = 0.0;
      for (p = 0; p < NVEL; p++) {
	sum += wv[p]*cv[p][i]*cv[p][j];
      }
      test_assert(fabs(sum - d_[i][j]/rcs2) < TEST_DOUBLE_TOLERANCE);
    }
  }

  /* Checking wv[p]*q_[p][i][j]... */

  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      sum = 0.0;
      for (p = 0; p < NVEL; p++) {
	sum += wv[p]*(cv[p][i]*cv[p][j] - cs2*d_[i][j]);
      }
      test_assert(fabs(sum - 0.0) < TEST_DOUBLE_TOLERANCE);
    }
  }

  /* Checking wv[p]*cv[p][i]*q_[p][j][k]... */

  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      for (k = 0; k < NDIM; k++) {
	sum = 0.0;
	for (p = 0; p < NVEL; p++) {
	  sum += wv[p]*cv[p][i]*(cv[p][i]*cv[p][j] - cs2*d_[i][j]);
	}
	test_assert(fabs(sum - 0.0) < TEST_DOUBLE_TOLERANCE);
      }
    }
  }

  /* Check ma_ against rho, cv conserved quantities */

  for (p = 0; p < NVEL; p++) {
    test_assert(fabs(ma_[0][p] - 1.0) < TEST_DOUBLE_TOLERANCE);
    for (i = 0; i < NDIM; i++) {
      test_assert(fabs(ma_[1+i][p] - cv[p][i]) < TEST_DOUBLE_TOLERANCE);
    }
  }

  /* Check ma_ against stress modes */

  for (p = 0; p < NVEL; p++) {
    k = 0;
    for (i = 0; i < NDIM; i++) {
      for (j = i; j < NDIM; j++) {
	double q_ij = cv[p][i]*cv[p][j] - cs2*d_[i][j];
	test_assert(fabs(ma_[1 + NDIM + k++][p] - q_ij)
		    < TEST_DOUBLE_TOLERANCE);
      }
    }
  }

  /* Checking normalisers norm_[i]*wv[p]*ma_[i][p]*ma_[j][p] = dij. */

  for (i = 0; i < NVEL; i++) {
    for (j = 0; j < NVEL; j++) {
      dij = (double) (i == j);
      sum = 0.0;
      for (p = 0; p < NVEL; p++) {
	sum += norm_[i]*wv[p]*ma_[i][p]*ma_[j][p];
      }
      test_assert(fabs(sum - dij) < TEST_DOUBLE_TOLERANCE);
    }
  }

  /* Checking ma_[i][p]*mi_[p][j] = dij ... */

  for (i = 0; i < NVEL; i++) {
    for (j = 0; j < NVEL; j++) {
      dij = (double) (i == j);
      sum = 0.0;
      for (p = 0; p < NVEL; p++) {
	sum += ma_[i][p]*mi_[p][j];
      }
      test_assert(fabs(sum - dij) < TEST_DOUBLE_TOLERANCE);
    }
  }

  return;
}

/*****************************************************************************
 *
 *  do_test_model_distributions
 *
 *  Test the distribution interface.
 *
 *****************************************************************************/

int do_test_model_distributions(pe_t * pe, cs_t * cs) {

  int i, n, p;
  int index = 1;
  int ndist = 2;
  double fvalue, fvalue_expected;
  double u[3];

  lb_t * lb;

  assert(pe);
  assert(cs);

  /* Tests of the basic distribution functions. */

  lb_create(pe, cs, &lb);
  assert(lb);
  lb_ndist(lb, &n);
  assert(n == 1); /* Default */

  lb_ndist_set(lb, ndist);
  lb_init(lb);

  /* Report the number of distributions */

  lb_ndist(lb, &n);
  assert(n == ndist);

  for (n = 0; n < ndist; n++) {
    for (p = 0; p < NVEL; p++) {
      fvalue_expected = 0.01*n + wv[p];
      lb_f_set(lb, index, p, n, fvalue_expected);
      lb_f(lb, index, p, n, &fvalue);
      assert(fabs(fvalue - fvalue_expected) < DBL_EPSILON);
    }

    /* Check zeroth moment... */

    fvalue_expected = 0.01*n*NVEL + 1.0;
    lb_0th_moment(lb, index, (lb_dist_enum_t) n, &fvalue);
    assert(fabs(fvalue - fvalue_expected) <= DBL_EPSILON);

    /* Check first moment... */

    lb_1st_moment(lb, index, (n == 0) ? LB_RHO : LB_PHI, u);

    for (i = 0; i < NDIM; i++) {
      assert(fabs(u[i] - 0.0) < DBL_EPSILON);
    }
  }

  lb_free(lb);

  return 0;
}

/*****************************************************************************
 *
 *  do_test_model_halo_swap
 *
 *  Test full halo swap.
 *
 *****************************************************************************/

int do_test_model_halo_swap(pe_t * pe, cs_t * cs) {

  int i, j, k, p;
  int n, ndist = 2;
  int index, nlocal[3];
  const int nextra = 1;  /* Distribution halo width always 1 */
  double f_expect;
  double f_actual;

  lb_t * lb = NULL;

  assert(pe);
  assert(cs);

  lb_create(pe, cs, &lb);
  assert(lb);
  lb_ndist_set(lb, ndist);
  lb_init(lb);

  cs_nlocal(cs, nlocal);

  /* The test relies on a uniform decomposition in parallel:
   *
   * f[0] or f[X] is set to local x index,
   * f[1] or f[Y] is set to local y index
   * f[2] or f[Z] is set to local z index
   * remainder are set to velocity index. */

  for (i = 1; i <= nlocal[X]; i++) {
    for (j = 1; j <= nlocal[Y]; j++) {
      for (k = 1; k <= nlocal[Z]; k++) {

	index = cs_index(cs, i, j, k);

	for (n = 0; n < ndist; n++) {
	  lb_f_set(lb, index, X, n, (double) (i));
	  lb_f_set(lb, index, Y, n, (double) (j));
	  lb_f_set(lb, index, Z, n, (double) (k));

	  for (p = 3; p < NVEL; p++) {
	    lb_f_set(lb, index, p, n, (double) p);
	  }
	}
      }
    }
  }

  lb_memcpy(lb, tdpMemcpyHostToDevice);
  lb_halo(lb);
  lb_memcpy(lb, tdpMemcpyDeviceToHost);

  /* Test all the sites not in the interior */

  for (i = 1 - nextra; i <= nlocal[X] + nextra; i++) {
    if (i >= 1 && i <= nlocal[X]) continue;
    for (j = 1 - nextra; j <= nlocal[Y] + nextra; j++) {
      if (j >= 1 && j <= nlocal[Y]) continue;
      for (k = 1 - nextra; k <= nlocal[Z] + nextra; k++) {
	if (k >= 1 && k <= nlocal[Z]) continue;

	index = cs_index(cs, i, j, k);

	for (n = 0; n < ndist; n++) {

	  f_expect = 1.0*abs(i - nlocal[X]);
	  lb_f(lb, index, X, n, &f_actual);
	  test_assert(fabs(f_actual - f_expect) < DBL_EPSILON);

	  f_expect = 1.0*abs(j - nlocal[Y]);
	  lb_f(lb, index, Y, n, &f_actual);
	  test_assert(fabs(f_actual - f_expect) < DBL_EPSILON);

	  f_expect = 1.0*abs(k - nlocal[Z]);
	  lb_f(lb, index, Z, n, &f_actual);
	  test_assert(fabs(f_actual - f_expect) < DBL_EPSILON);

	  for (p = 3; p < NVEL; p++) {
	    lb_f(lb, index, p, n, &f_actual);
	    f_expect = (double) p;
	    test_assert(fabs(f_actual - f_expect) < DBL_EPSILON);
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
 *  do_test_model_reduced_halo_swap
 *
 *****************************************************************************/

int do_test_model_reduced_halo_swap(pe_t * pe, cs_t * cs) {  

  int i, j, k, p;
  int icdt, jcdt, kcdt;
  int index, nlocal[3];
  int n, ndist = 2;
  const int nextra = 1;

  double f_expect;
  double f_actual;

  lb_t * lb = NULL;

  assert(pe);
  assert(cs);

  lb_create(pe, cs, &lb);
  assert(lb);
  lb_ndist_set(lb, ndist);
  lb_init(lb);
  lb_halo_set(lb, LB_HALO_REDUCED);

  cs_nlocal(cs, nlocal);

  /* Set everything which is NOT in a halo */

  for (i = 1; i <= nlocal[X]; i++) {
    for (j = 1; j <= nlocal[Y]; j++) {
      for (k = 1; k <= nlocal[Z]; k++) {
	index = cs_index(cs, i, j, k);
	for (n = 0; n < ndist; n++) {
	  for (p = 0; p < NVEL; p++) {
	    f_expect = 1.0*(n*NVEL + p);
	    lb_f_set(lb, index, p, n, f_expect);
	  }
	}
      }
    }
  }

  lb_halo_via_struct(lb);

  /* Now check that the interior sites are unchanged */

  for (i = 1; i <= nlocal[X]; i++) {
    for (j = 1; j <= nlocal[Y]; j++) {
      for (k = 1; k <= nlocal[Z]; k++) {
	index = cs_index(cs, i, j, k);
	for (n = 0; n < ndist; n++) {
	  for (p = 0; p < NVEL; p++) {
	    lb_f(lb, index, p, n, &f_actual);
	    f_expect = 1.0*(n*NVEL +  p);
	    test_assert(fabs(f_expect - f_actual) < DBL_EPSILON);
	  }
	}
      }
    }
  }

  /* Also check the halos sites. The key test of the reduced halo
   * swap is that distributions for which r + c_i dt takes us into
   * the domain proper must be correct. */

  for (i = 1 - nextra; i <= nlocal[X] + nextra; i++) {
    if (i >= 1 && i <= nlocal[X]) continue;
    for (j = 1 - nextra; j <= nlocal[Y] + nextra; j++) {
      if (j >= 1 && j <= nlocal[Y]) continue;
      for (k = 1 - nextra; k <= nlocal[Z] + nextra; k++) {
	if (k >= 1 && k <= nlocal[Z]) continue;

	index = cs_index(cs, i, j, k);

	for (n = 0; n < ndist; n++) {
	  for (p = 0; p < NVEL; p++) {

	    lb_f(lb, index, p, n, &f_actual);
	    f_expect = 1.0*(n*NVEL + p);

	    icdt = i + cv[p][X];
	    jcdt = j + cv[p][Y];
	    kcdt = k + cv[p][Z];

	    if (test_model_is_domain(cs, icdt, jcdt, kcdt)) {
	      test_assert(fabs(f_actual - f_expect) < DBL_EPSILON);
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

/*****************************************************************************
 *
 *  test_model_is_domain
 *
 *  Is (ic, jc, kc) in the domain proper?
 *
 *****************************************************************************/

static int test_model_is_domain(cs_t * cs, int ic, int jc, int kc) {

  int nlocal[3];
  int iam = 1;

  assert(cs);

  cs_nlocal(cs, nlocal);

  if (ic < 1) iam = 0;
  if (jc < 1) iam = 0;
  if (kc < 1) iam = 0;
  if (ic > nlocal[X]) iam = 0;
  if (jc > nlocal[Y]) iam = 0;
  if (kc > nlocal[Z]) iam = 0;

  return iam;
}

/*****************************************************************************
 *
 *  do_test_lb_model_io
 *
 *****************************************************************************/

int do_test_lb_model_io(pe_t * pe, cs_t * cs) {

  int ndist = 2;
  lb_t * lbrd = NULL;
  lb_t * lbwr = NULL;

  assert(pe);
  assert(cs);

  lb_create_ndist(pe, cs, ndist, &lbrd);
  lb_create_ndist(pe, cs, ndist, &lbwr);

  lb_init(lbwr);
  lb_init(lbrd);

  /* Write */

  /* Read */

  /* Compare */

  lb_free(lbwr);
  lb_free(lbrd);

  return 0;
}

/*****************************************************************************
 *
 *  test_d3q19_ghosts
 *
 *  In comparison with Chen and Ladd (2007) we have
 *
 *   chi1  (2cs^2 - 3)(3c_z^2 - cs^2)           mode[10]
 *   chi2  (2cs^2 - 3)(c_y^2 - c_x^2)           mode[14]
 *   chi3  3cs^4 - 6cs^2 + 1                    mode[18]
 *
 *   jchi1 is in fact rho chi3 cv 
 *   jchi1[X] (3*cs^4 - 6cs^2 + 1) cx           mode[11]
 *   jchi1[Y] (3*cs^4 - 6cs^2 + 1) cy           mode[12]
 *   jchi1[Z] (3*cs^4 - 6cs^2 + 1) cz           mode[13]
 *
 *   jchi2 is rho chi2 cv
 *   jchi2[X]  (2cs^2 - 3)(c_y^2 - c_x^2) c_x   mode[15]
 *   jchi2[Y]  (2cs^2 - 3)(c_y^2 - c_x^2) c_y   mode[16]
 *   jchi2[Z]  (2cs^2 - 3)(c_y^2 - c_x^2) c_z   mode[17]
 *
 *   The expressions for the ghost currents appearing in Chun and Ladd
 *   are not quite consistent; the reason for this is unclear.
 *   Note that c_x and c_z are transposed in chi1 and chi2 cf Chun and Ladd.
 *
 *****************************************************************************/

int do_test_d3q19_ghosts(void) {

  int p;
  double rho = 1.0;
  double cs2, chi1, chi2, chi3;
  double jchi1[3], jchi2[3];

  if (NVEL != 19) return 0;

  for (p = 0; p < NVEL; p++) {

    cs2 = cv[p][X]*cv[p][X] + cv[p][Y]*cv[p][Y] + cv[p][Z]*cv[p][Z];
    chi1 = (2.0*cs2 - 3.0)*(3.0*cv[p][Z]*cv[p][Z] - cs2);
    chi2 = (2.0*cs2 - 3.0)*(cv[p][Y]*cv[p][Y] - cv[p][X]*cv[p][X]);
    chi3 = 3.0*cs2*cs2 - 6.0*cs2 + 1;

    jchi1[X] = rho*chi1*cv[p][X];
    jchi1[Y] = rho*chi1*cv[p][Y];
    jchi1[Z] = rho*chi1*cv[p][Z];

    jchi2[X] = rho*chi2*cv[p][X];
    jchi2[Y] = rho*chi2*cv[p][Y];
    jchi2[Z] = rho*chi2*cv[p][Z];

    test_assert(fabs(ma_[10][p] - chi1) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(ma_[14][p] - chi2) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(ma_[18][p] - chi3) < TEST_DOUBLE_TOLERANCE);

    test_assert(fabs(ma_[11][p] - jchi1[X]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(ma_[12][p] - jchi1[Y]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(ma_[13][p] - jchi1[Z]) < TEST_DOUBLE_TOLERANCE);

    test_assert(fabs(ma_[15][p] - jchi2[X]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(ma_[16][p] - jchi2[Y]) < TEST_DOUBLE_TOLERANCE);
    test_assert(fabs(ma_[17][p] - jchi2[Z]) < TEST_DOUBLE_TOLERANCE);
  }

  return 0;
}
