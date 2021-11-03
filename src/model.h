/*****************************************************************************
 *
 *  model.h
 *
 *  This includes the appropriate lattice Boltzmann
 *  model d2q9, d3q15, or d3q19 (see Makefile).
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2020 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_MODEL_H
#define LUDWIG_MODEL_H

/* Here is the choice of model (and it's here only). */

#ifdef _D2Q9_
#include "d2q9.h"

enum {NDIM     = NDIM9,
      NVEL     = NVEL9,
      CVXBLOCK = CVXBLOCK9,
      CVYBLOCK = CVYBLOCK9,
      CVZBLOCK = CVZBLOCK9};

#define LB_CV(cv) LB_CV_D2Q9(cv)
#define LB_NORMALISERS_DOUBLE(n) LB_NORMALISERS_D2Q9(n)

#endif

#ifdef _D3Q15_
#include "d3q15.h"

enum {NDIM     = NDIM15,
      NVEL     = NVEL15,
      CVXBLOCK = CVXBLOCK15,
      CVYBLOCK = CVYBLOCK15,
      CVZBLOCK = CVZBLOCK15};

#define LB_CV(cv) LB_CV_D3Q15(cv)
#define LB_NORMALISERS_DOUBLE(n) LB_NORMALISERS_D3Q15(n)

#endif

#ifdef _D3Q19_

#include "d3q19.h"

enum {NDIM     = NDIM19,
      NVEL     = NVEL19,
      CVXBLOCK = CVXBLOCK19,
      CVYBLOCK = CVYBLOCK19,
      CVZBLOCK = CVZBLOCK19};

#define LB_CV(cv) LB_CV_D3Q19(cv)
#define LB_NORMALISERS_DOUBLE(n) LB_NORMALISERS_D3Q19(n)

#endif

#include "pe.h"
#include "coords.h"
#include "io_harness.h"
#include "memory.h"

/* Number of hydrodynamic modes */
enum {NHYDRO = 1 + NDIM + NDIM*(NDIM+1)/2};

/* Labels to locate relaxation times in array[NVEL] */
/* Bulk viscosity is XX in stress */
/* Shear is XY in stress */

enum {LB_TAU_BULK = 1 + NDIM + XX, LB_TAU_SHEAR = 1 + NDIM + XY};

#define LB_CS2_DOUBLE(cs2)   const double cs2 = (1.0/3.0)
#define LB_RCS2_DOUBLE(rcs2) const double rcs2 = 3.0

typedef enum lb_dist_enum_type{LB_RHO = 0, LB_PHI = 1} lb_dist_enum_t;
typedef enum lb_mode_enum_type{LB_GHOST_ON = 0, LB_GHOST_OFF = 1} lb_mode_enum_t;

typedef struct lb_data_s lb_t;
typedef enum lb_halo_enum_type {LB_HALO_FULL,
				LB_HALO_REDUCED,
				LB_HALO_HOST,
				LB_HALO_TARGET} lb_halo_enum_t;

typedef enum {LB_RELAXATION_M10, LB_RELAXATION_BGK, LB_RELAXATION_TRT}
  lb_relaxation_enum_t;

__host__ int lb_create_ndist(pe_t * pe, cs_t * cs, int ndist, lb_t ** lb);
__host__ int lb_create(pe_t * pe, cs_t * cs, lb_t ** lb);
__host__ int lb_init(lb_t * lb);
__host__ int lb_free(lb_t * lb);
__host__ int lb_memcpy(lb_t * lb, tdpMemcpyKind flag);
__host__ int lb_collide_param_commit(lb_t * lb);
__host__ int lb_halo(lb_t * lb);
__host__ int lb_halo_swap(lb_t * lb, lb_halo_enum_t flag);
__host__ int lb_halo_via_copy(lb_t * lb);
__host__ int lb_halo_via_struct(lb_t * lb);
__host__ int lb_halo_set(lb_t * lb, lb_halo_enum_t halo);
__host__ int lb_io_info(lb_t * lb, io_info_t ** io_info);
__host__ int lb_io_info_set(lb_t * lb, io_info_t * io_info, int fin, int fout);
__host__ int lb_io_rho_set(lb_t *lb, io_info_t * io_rho, int fin, int fout);

__host__ int lb_io_info_commit(lb_t * lb, io_info_args_t args);

__host__ __device__ int lb_nvel(lb_t * lb, int * nvel);
__host__ __device__ int lb_ndim(lb_t * lb, int * ndim);
__host__ __device__ int lb_ndist(lb_t * lb, int * ndist);
__host__ __device__ int lb_f(lb_t * lb, int index, int p, int n, double * f);
__host__ __device__ int lb_f_set(lb_t * lb, int index, int p, int n, double f);
__host__ __device__ int lb_0th_moment(lb_t * lb, int index, lb_dist_enum_t nd,
				      double * rho);
__host__ __device__ int lb_f_index(lb_t * lb, int index, int n, double f[NVEL]);
__host__ __device__ int lb_f_index_set(lb_t * lb, int index, int n,
				       double f[NVEL]);
__host__ __device__ int lb_f_multi_index_part(lb_t * lb, int index, int n,
					      double f[NVEL][NSIMDVL], int nv);
__host__ __device__ int lb_f_multi_index_set_part(lb_t * lb, int index, int n,
						  double f[NVEL][NSIMDVL],
						  int nv);
__host__ __device__ int lb_f_multi_index(lb_t * lb, int index, int n,
					 double f[NVEL][NSIMDVL]);
__host__ __device__ int lb_f_multi_index_set(lb_t * lb, int index, int n,
					     double f[NVEL][NSIMDVL]);

/* These  could be __host__ __device__ pending removal of
 * static constants */

__host__ int lb_nblock(lb_t * lb, int dim, int * nblock);
__host__ int lb_init_rest_f(lb_t * lb, double rho0);
__host__ int lb_order(lb_t * lb);
__host__ int lb_ndist_set(lb_t * lb, int ndist);
__host__ int lb_1st_moment(lb_t * lb, int index, lb_dist_enum_t nd, double g[3]);
__host__ int lb_2nd_moment(lb_t * lb, int index, lb_dist_enum_t nd, double s[3][3]);
__host__ int lb_0th_moment_equilib_set(lb_t * lb, int index, int n, double rho);
__host__ int lb_1st_moment_equilib_set(lb_t * lb, int index, double rho, double u[3]);

#endif
