/*****************************************************************************
 *
 *  d3q15.h
 *
 *  D3Q15 definitions.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group
 *  and Edinburgh Parallel Computing Centre
 *
 *  (c) 2008-2021 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_D3Q15_MODEL_H
#define LUDWIG_D3Q15_MODEL_H

enum {NDIM15 = 3};
enum {NVEL15 = 15};
enum {CVXBLOCK15 = 1};
enum {CVYBLOCK15 = 3};
enum {CVZBLOCK15 = 5};

#include "lb_d3q15.h"

#ifdef _D3Q15_

extern const int xblocklen_cv[CVXBLOCK15];
extern const int xdisp_fwd_cv[CVXBLOCK15];
extern const int xdisp_bwd_cv[CVXBLOCK15];

extern const int yblocklen_cv[CVYBLOCK15];
extern const int ydisp_fwd_cv[CVYBLOCK15];
extern const int ydisp_bwd_cv[CVYBLOCK15];

extern const int zblocklen_cv[CVZBLOCK15];
extern const int zdisp_fwd_cv[CVZBLOCK15];
extern const int zdisp_bwd_cv[CVZBLOCK15];

#endif

#endif

