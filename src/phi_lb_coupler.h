/*****************************************************************************
 *
 *  phi_lb_coupler.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2009-2022 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  Alan Gray (alang@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef PHI_LB_COUPLER_H
#define PHI_LB_COUPLER_H

#include "field.h"
#include "lb_data.h"

__host__ int phi_lb_to_field(field_t * phi, lb_t * lb);
__host__ int phi_lb_from_field(field_t * phi, lb_t * lb);

#endif
