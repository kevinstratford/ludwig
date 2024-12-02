/*****************************************************************************
 *
 *  gradient_d3q27.h
 *
 *  Only the _d2 routine is actually available at the moment.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2022 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_GRADIENT_D3Q27_H
#define LUDWIG_GRADIENT_D3Q27_H

#include "field_grad.h"

__host__ int gradient_d3q27_d2(field_grad_t * fg);
__host__ int gradient_d3q27_d4(field_grad_t * fg);
__host__ int gradient_d3q27_ab(field_grad_t * fg);

#endif
