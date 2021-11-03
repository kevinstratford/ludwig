/****************************************************************************
 *
 *  fe_surfactant.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group
 *  and Edinburgh Parallel Computing Centre
 *
 *  (c) 2009-2019 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 ****************************************************************************/

#ifndef LUDWIG_FE_SURFACTANT_H
#define LUDWIG_FE_SURFACTANT_H

#include "memory.h"
#include "free_energy.h"
#include "field.h"
#include "field_grad.h"

typedef struct fe_surfactant_s fe_surf_t;
typedef struct fe_surfactant_param_s fe_surf_param_t;

struct fe_surfactant_param_s {
  double a;              /* Symmetric a */
  double b;              /* Symmetric b */
  double kappa;          /* Symmetric kappa */

  double kt;             /* Surfactant kT */
  double epsilon;        /* Surfactant epsilon */
  double beta;           /* Frumpkin isotherm */
  double w;              /* Surfactant w */
};

struct fe_surfactant_s {
  fe_t super;                      /* "Superclass" block */
  pe_t * pe;                       /* Parallel environment */
  cs_t * cs;                       /* Coordinate system */
  fe_surf_param_t * param;         /* Parameters */
  field_t * phi;                   /* Single field with {phi,psi} */
  field_grad_t * dphi;             /* gradients thereof */
  fe_surf_t * target;              /* Device copy */
};

__host__ int fe_surf_create(pe_t * pe, cs_t * cs, field_t * phi,
			    field_grad_t * dphi, fe_surf_param_t param,
			    fe_surf_t ** fe);
__host__ int fe_surf_free(fe_surf_t * fe);
__host__ int fe_surf_info(fe_surf_t * fe);
__host__ int fe_surf_param_set(fe_surf_t * fe, fe_surf_param_t vals);
__host__ int fe_surf_sigma(fe_surf_t * fe, double * sigma);
__host__ int fe_surf_xi0(fe_surf_t * fe,  double * xi0);
__host__ int fe_surf_langmuir_isotherm(fe_surf_t * fe, double * psi_c);
__host__ int fe_surf_target(fe_surf_t * fe, fe_t ** target);

__host__ int fe_surf_param(fe_surf_t * fe, fe_surf_param_t * param);
__host__ int fe_surf_fed(fe_surf_t * fe, int index, double * fed);
__host__ int fe_surf_mu(fe_surf_t * fe, int index, double * mu);
__host__ int fe_surf_str(fe_surf_t * fe, int index, double s[3][3]);
__host__ int fe_surf_str_v(fe_surf_t * fe, int index, double s[3][3][NSIMDVL]);

#endif
