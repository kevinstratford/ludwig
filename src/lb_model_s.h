/*****************************************************************************
 *
 *  lb_model_s.h
 *
 *  LB model data structure implementation.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2014-2020 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LB_MODEL_S_H
#define LB_MODEL_S_H

#include "model.h"
#include "halo_swap.h"
#include "io_harness.h"
#include "stdint.h"

typedef struct lb_collide_param_s lb_collide_param_t;

struct lb_collide_param_s {
  int8_t isghost;                      /* switch for ghost modes */
  int8_t cv[NVEL][3];
  int nsite;
  int ndist;
  double rho0;
  double eta_shear;
  double var_shear;
  double eta_bulk;
  double var_bulk;
  double var_noise[NVEL];
  double rtau[NVEL];
  double wv[NVEL];
  double q[NVEL][3][3];
  double ma[NVEL][NVEL];
  double mi[NVEL][NVEL];
};

struct lb_data_s {

  int ndist;             /* Number of distributions (default one) */
  int nsite;             /* Number of lattice sites (local) */
  int model;             /* MODEL or MODEL_R */

  pe_t * pe;             /* parallel environment */
  cs_t * cs;             /* coordinate system */
  halo_swap_t * halo;    /* halo swap driver */
  io_info_t * io_info;   /* Distributions */ 
  io_info_t * io_rho;    /* Fluid density (here; could be hydrodynamics...) */

  double * f;            /* Distributions */
  double * fprime;       /* used in propagation only */

  lb_collide_param_t * param;   /* Collision parameters REFACTOR THIS */
  lb_relaxation_enum_t nrelax;  /* Relaxation scheme */

  /* MPI data types for halo swaps; these are comupted at runtime
   * to conform to the model selected at compile time */

  MPI_Datatype plane_xy_full;
  MPI_Datatype plane_xz_full;
  MPI_Datatype plane_yz_full;
  MPI_Datatype plane_xy_reduced[2];
  MPI_Datatype plane_xz_reduced[2];
  MPI_Datatype plane_yz_reduced[2];
  MPI_Datatype plane_xy[2];
  MPI_Datatype plane_xz[2];
  MPI_Datatype plane_yz[2];
  MPI_Datatype site_x[2];
  MPI_Datatype site_y[2];
  MPI_Datatype site_z[2];

  lb_t * target;              /* copy of this structure on target */ 
};

/* Data storage: A rank two object */

#define LB_ADDR(nsites, ndist, nvel, index, n, p) \
  addr_rank2(nsites, ndist, nvel, index, n, p)

#endif
