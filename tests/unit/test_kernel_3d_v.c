/*****************************************************************************
 *
 *  test_kernel_3d_v.c
 *
 *  Kernel helper for vectorised kernels.
 *
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2024 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>

#include "kernel_3d_v.h"

int test_kernel_3d_v(pe_t * pe);
int test_kernel_3d_v_coords(pe_t * pe);
int test_kernel_3d_v_mask(pe_t * pe);
int test_kernel_3d_v_cs_index(pe_t * pe);

__global__ void test_kernel_3d_v_coords_kernel(kernel_3d_v_t k3v, cs_t * cs);
__global__ void test_kernel_3d_v_mask_kernel(kernel_3d_v_t k3v);
__global__ void test_kernel_3d_v_cs_index_kernel(kernel_3d_v_t k3v, cs_t * cs);

int util_imin(int a, int b) {return a < b ? a : b;}

/*****************************************************************************
 *
 *  test_kernel_3d_v_suite
 *
 *****************************************************************************/

int test_kernel_3d_v_suite(void) {

  pe_t * pe = NULL;

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);

  test_kernel_3d_v(pe);
  test_kernel_3d_v_coords(pe);
  test_kernel_3d_v_mask(pe);
  test_kernel_3d_v_cs_index(pe);

  /* Hard limit for pass-by-vlaue arguments; must not exceed 64 bytes */
  assert(sizeof(kernel_3d_v_t) <= 64);

  pe_info(pe, "%-9s %s\n", "PASS", __FILE__);
  pe_free(pe);

  return 0;
}

/*****************************************************************************
 *
 *  test_kernel_3d_v
 *
 *****************************************************************************/

int test_kernel_3d_v(pe_t * pe) {

  int ifail = 0;

  int nhalo = 2;
  int ntotal[3] = {64, 32, 16};
  int nlocal[3] = {0};
  cs_t * cs = NULL;

  cs_create(pe, &cs);
  cs_ntotal_set(cs, ntotal);
  cs_nhalo_set(cs, nhalo);
  cs_init(cs);
  cs_nlocal(cs, nlocal);

  /* Vector length 1 */
  {
    cs_limits_t lim = {1, nlocal[X], 1, nlocal[Y], 1, nlocal[Z]};
    kernel_3d_v_t k3v = kernel_3d_v(cs, lim, 1);

    assert(k3v.nhalo       == nhalo);
    assert(k3v.nlocal[X]   == nlocal[X]);
    assert(k3v.nlocal[Y]   == nlocal[Y]);
    assert(k3v.nlocal[Z]   == nlocal[Z]);
    assert(k3v.kindex0     == cs_index(cs, 1, 1-nhalo, 1-nhalo));
    assert(k3v.nklocal[X]  == nlocal[X]);
    assert(k3v.nklocal[Y]  == nlocal[Y] + 2*nhalo);
    assert(k3v.nklocal[Z]  == nlocal[Z] + 2*nhalo);
    assert(k3v.kiterations == k3v.nklocal[X]*k3v.nklocal[Y]*k3v.nklocal[Z]);
    assert(k3v.lim.imin    == lim.imin);
    assert(k3v.lim.imax    == lim.imax);
    assert(k3v.lim.jmin    == lim.jmin);
    assert(k3v.lim.jmax    == lim.jmax);
    assert(k3v.lim.kmin    == lim.kmin);
    assert(k3v.lim.kmax    == lim.kmax);

    if (k3v.nhalo != nhalo) ifail = -1;
  }

  /* Vector length 2; only kindex0 is relevant */
  {
    int nsimdvl = 2;
    cs_limits_t lim = {1, nlocal[X], 1, nlocal[Y], 1, nlocal[Z]};
    kernel_3d_v_t k3v = kernel_3d_v(cs, lim, nsimdvl);

    assert(k3v.kindex0 <= cs_index(cs, 1, 1-nhalo, 1-nhalo));
    assert(k3v.kindex0 % nsimdvl == 0);

    if (k3v.kindex0 % nsimdvl != 0) ifail = -1;
  }

  /* Vector length 4; ditto ... */
  {
    int nsimdvl = 4;
    cs_limits_t lim = {1, nlocal[X], 1, nlocal[Y], 1, nlocal[Z]};
    kernel_3d_v_t k3v = kernel_3d_v(cs, lim, nsimdvl);

    assert(k3v.kindex0 <= cs_index(cs, 1, 1-nhalo, 1-nhalo));
    assert(k3v.kindex0 % nsimdvl == 0);

    if (k3v.kindex0 % nsimdvl != 0) ifail = -1;
  }

  /* Compile time vector length */
  {
    cs_limits_t lim = {1, nlocal[X], 1, nlocal[Y], 1, nlocal[Z]};
    kernel_3d_v_t k3v = kernel_3d_v(cs, lim, NSIMDVL);

    assert(k3v.kindex0 <= cs_index(cs, 1, 1-nhalo, 1-nhalo));
    assert(k3v.kindex0 % NSIMDVL == 0);

    if (k3v.kindex0 % NSIMDVL != 0) ifail = -1;
  }

  cs_free(cs);

  return ifail;
}

/*****************************************************************************
 *
 *  test_kernel_3d_v_coords
 *
 *****************************************************************************/

int test_kernel_3d_v_coords(pe_t * pe) {

  int ifail = 0;

  int nhalo = 2;
  int ntotal[3] = {64, 32, 16};
  int nlocal[3] = {0};
  cs_t * cs = NULL;

  cs_create(pe, &cs);
  cs_ntotal_set(cs, ntotal);
  cs_nhalo_set(cs, nhalo);
  cs_init(cs);
  cs_nlocal(cs, nlocal);

  /* Kernel. The requested vector length must be available as arguments
   * are hardwired as NSIMDVL. */
  {
    int nsimdvl = util_imin(4, NSIMDVL);
    dim3 nblk = {};
    dim3 ntpb = {};
    cs_limits_t lim = {0, nlocal[X], 0, nlocal[Y], 0, nlocal[Z]};
    kernel_3d_v_t k3v = kernel_3d_v(cs, lim, nsimdvl);

    kernel_3d_v_exec_conf(&k3v, &nblk, &ntpb);

    tdpLaunchKernel(test_kernel_3d_v_coords_kernel, nblk, ntpb, 0, 0,
		    k3v, cs->target);

    tdpAssert( tdpPeekAtLastError() );
    tdpAssert( tdpDeviceSynchronize() );
  }

  cs_free(cs);

  return ifail;
}

/*****************************************************************************
 *
 *  test_kernel_3d_v_mask
 *
 *****************************************************************************/

int test_kernel_3d_v_mask(pe_t * pe) {

  int ifail = 0;

  int nhalo = 2;
  int ntotal[3] = {64, 32, 16};
  int nlocal[3] = {0};
  cs_t * cs = NULL;

  cs_create(pe, &cs);
  cs_ntotal_set(cs, ntotal);
  cs_nhalo_set(cs, nhalo);
  cs_init(cs);
  cs_nlocal(cs, nlocal);

  /* Kernel. The requested vector length must be available as arguments
   * are hardwired as NSIMDVL. */
  {
    int nsimdvl = util_imin(4, NSIMDVL);
    dim3 nblk = {};
    dim3 ntpb = {};
    cs_limits_t lim = {0, nlocal[X], 0, nlocal[Y], 0, nlocal[Z]};
    kernel_3d_v_t k3v = kernel_3d_v(cs, lim, nsimdvl);

    kernel_3d_v_exec_conf(&k3v, &nblk, &ntpb);

    tdpLaunchKernel(test_kernel_3d_v_mask_kernel, nblk, ntpb, 0, 0, k3v);

    tdpAssert( tdpPeekAtLastError() );
    tdpAssert( tdpDeviceSynchronize() );
  }

  cs_free(cs);

  return ifail;
}

/*****************************************************************************
 *
 *  test_kernel_3d_v_cs_index
 *
 *****************************************************************************/

int test_kernel_3d_v_cs_index(pe_t * pe) {

  int ifail = 0;

  int nhalo = 2;
  int ntotal[3] = {64, 32, 16};
  int nlocal[3] = {0};
  cs_t * cs = NULL;

  cs_create(pe, &cs);
  cs_ntotal_set(cs, ntotal);
  cs_nhalo_set(cs, nhalo);
  cs_init(cs);
  cs_nlocal(cs, nlocal);

  /* Kernel. The requested vector length must be available as arguments
   * are hardwired as NSIMDVL. */
  {
    int nsimdvl = util_imin(4, NSIMDVL);
    dim3 nblk = {};
    dim3 ntpb = {};
    cs_limits_t lim = {0, nlocal[X], 0, nlocal[Y], 0, nlocal[Z]};
    kernel_3d_v_t k3v = kernel_3d_v(cs, lim, nsimdvl);

    kernel_3d_v_exec_conf(&k3v, &nblk, &ntpb);

    tdpLaunchKernel(test_kernel_3d_v_cs_index_kernel, nblk, ntpb, 0, 0,
		    k3v, cs->target);

    tdpAssert( tdpPeekAtLastError() );
    tdpAssert( tdpDeviceSynchronize() );
  }

  cs_free(cs);

  return ifail;
}

/*****************************************************************************
 *
 *  test_kernel_3d_v_coords_kernel
 *
 *****************************************************************************/

__global__ void test_kernel_3d_v_coords_kernel(kernel_3d_v_t k3v, cs_t * cs) {

  int kindex = 0;

  assert(k3v.nsimdvl == NSIMDVL);

  for_simt_parallel(kindex, k3v.kiterations, k3v.nsimdvl) {
    int icv[NSIMDVL];
    int jcv[NSIMDVL];
    int kcv[NSIMDVL];
    kernel_3d_v_coords(&k3v, kindex, icv, jcv, kcv);
    for (int iv = 0; iv < k3v.nsimdvl; iv++) {
      int ic = icv[iv];
      int jc = jcv[iv];
      int kc = kcv[iv];
      int index = cs_index(cs, ic, jc, kc);
      assert((k3v.kindex0 + kindex + iv) == index);
      if ((k3v.kindex0 + kindex + iv)    != index) printf("Fail\n");
    }
  }

  return;
}

/*****************************************************************************
 *
 *  test_kernel_3d_v_mask_kernel
 *
 *****************************************************************************/

__global__ void test_kernel_3d_v_mask_kernel(kernel_3d_v_t k3v) {

  int kindex = 0;

  assert(k3v.nsimdvl == NSIMDVL);

  for_simt_parallel(kindex, k3v.kiterations, k3v.nsimdvl) {
    int icv[NSIMDVL];
    int jcv[NSIMDVL];
    int kcv[NSIMDVL];
    int maskv[NSIMDVL];

    kernel_3d_v_coords(&k3v, kindex, icv, jcv, kcv);
    kernel_3d_v_mask(&k3v, icv, jcv, kcv, maskv);

    for (int iv = 0; iv < k3v.nsimdvl; iv++) {
      int ic = icv[iv];
      int jc = jcv[iv];
      int kc = kcv[iv];
      int mask = maskv[iv];
      if (mask) {
	/* we should be inside the kernel domain */
	int inx = (k3v.lim.imin <= ic && ic <= k3v.lim.imax);
	int iny = (k3v.lim.jmin <= jc && jc <= k3v.lim.jmax);
	int inz = (k3v.lim.kmin <= kc && kc <= k3v.lim.kmax);
	assert(inx && iny && inz);
	if ((inx && iny && inz) == 0) printf("Fail within mask\n");
      }
      else {
	/* we should be outside ... */
	int outx = (k3v.lim.imin > ic || ic > k3v.lim.imax);
	int outy = (k3v.lim.jmin > jc || jc > k3v.lim.jmax);
	int outz = (k3v.lim.kmin > kc || kc > k3v.lim.kmax);
	assert(outx == 0);
	assert(outy || outz);
	if (outx != 0) printf("Fail outx\n");
	if ((outy || outz) == 0) printf("Fail without mask\n");
      }
    }
  }

  return;
}

/*****************************************************************************
 *
 *  test_kernel_3d_v_cs_index_kernel
 *
 *****************************************************************************/

__global__ void test_kernel_3d_v_cs_index_kernel(kernel_3d_v_t k3v, cs_t *cs) {

  int kindex = 0;

  assert(k3v.nsimdvl == NSIMDVL);

  for_simt_parallel(kindex, k3v.kiterations, k3v.nsimdvl) {
    int icv[NSIMDVL];
    int jcv[NSIMDVL];
    int kcv[NSIMDVL];
    int indexv[NSIMDVL];

    kernel_3d_v_coords(&k3v, kindex, icv, jcv, kcv);
    kernel_3d_v_cs_index(&k3v, icv, jcv, kcv, indexv);

    for (int iv = 0; iv < k3v.nsimdvl; iv++) {
      int ic = icv[iv];
      int jc = jcv[iv];
      int kc = kcv[iv];
      int index = cs_index(cs, ic, jc, kc);
      assert(indexv[iv] == index);
      if (indexv[iv] != index) printf("Fail %2d %2d %2d\n", ic, jc, kc);
    }
  }

  return;
}
