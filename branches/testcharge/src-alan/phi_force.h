/*****************************************************************************
 *
 *  phi_force.h
 *
 *  $Id$
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef PHI_FORCE_H
#define PHI_FORCE_H

#include "field.h"
#include "hydro.h"

int phi_force_calculation(field_t * phi, hydro_t * hydro, double dt);
int phi_force_required(int * flag);
int phi_force_required_set(const int flag);
int phi_force_divergence_set(const int flag);

#endif