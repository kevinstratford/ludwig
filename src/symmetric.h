/****************************************************************************
 *
 *  symmetric.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group
 *  and Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2021 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 ****************************************************************************/

#ifndef LUDWIG_SYMMETRIC_H
#define LUDWIG_SYMMETRIC_H

#include "memory.h"
#include "free_energy.h"
#include "field.h"
#include "field_grad.h"

typedef struct fe_symm_param_s fe_symm_param_t;
typedef struct fe_symm_s fe_symm_t;

/* Free energy structure */

struct fe_symm_s {
  fe_t super;
  pe_t * pe;                   /* Parallel environment */
  cs_t * cs;                   /* Coordinate system */
  fe_symm_param_t * param;     /* Parameters */
  field_t * phi;               /* Scalar order parameter or composition */
  field_grad_t * dphi;         /* Gradients thereof */
  fe_symm_t * target;          /* Target copy */
};

/* Parameters */

struct fe_symm_param_s {
  double a;                    /* Bulk quadratic term */
  double b;                    /* Bulk quartic term */
  double kappa;                /* Interfacial penalty term */
  double c;                    /* Surface wetting parameter C (uniform) */
  double h;                    /* Surface wetting parameter H (ditto) */
};

__host__ int fe_symm_create(pe_t * pe, cs_t * cs, field_t * f,
			    field_grad_t * grd, fe_symm_t ** p);
__host__ int fe_symm_free(fe_symm_t * fe);
__host__ int fe_symm_param_set(fe_symm_t * fe, fe_symm_param_t values);
__host__ int fe_symm_target(fe_symm_t * fe, fe_t ** target);

__host__ __device__ int fe_symm_param(fe_symm_t * fe, fe_symm_param_t * values);
__host__ __device__ int fe_symm_interfacial_tension(fe_symm_t * fe, double * s);
__host__ __device__ int fe_symm_interfacial_width(fe_symm_t * fe, double * xi);
__host__ __device__ int fe_symm_fed(fe_symm_t * fe, int index, double * fed);
__host__ __device__ int fe_symm_mu(fe_symm_t * fe, int index, double * mu);

__host__ __device__ int fe_symm_str(fe_symm_t * fe, int index, double s[3][3]);
__host__ __device__ void fe_symm_str_v(fe_symm_t * fe, int index,
				       double s[3][3][NSIMDVL]);

/* Some additional host-only utilities */
__host__ int fe_symm_theta_to_h(double theta, double * h);
__host__ int fe_symm_h_to_costheta(double h, double * costheta);


#endif

