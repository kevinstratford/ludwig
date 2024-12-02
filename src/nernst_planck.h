/*****************************************************************************
 *
 *  nernst_planck.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2012-2022 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_NERNST_PLANCK_H
#define LUDWIG_NERNST_PLANCK_H

#include "psi.h"
#include "free_energy.h"
#include "hydro.h"
#include "map.h"
#include "colloids.h"

int nernst_planck_driver(psi_t * psi, fe_t * fe, map_t * map);
int nernst_planck_driver_d3qx(psi_t * psi, fe_t * fe, hydro_t * hydro,
			      map_t * map, colloids_info_t * cinfo);
int nernst_planck_adjust_multistep(psi_t * psi);

int nernst_planck_maxacc(double * acc);
int nernst_planck_maxacc_set(double acc);

#endif
