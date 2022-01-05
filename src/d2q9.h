/*****************************************************************************
 *
 *  d2q9.h
 *
 *  D2Q9 definitions.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2021 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_D2Q9_MODEL_H
#define LUDWIG_D2Q9_MODEL_H

enum {NDIM9 = 2};
enum {NVEL9 = 9};
enum {CVXBLOCK9 = 1};
enum {CVYBLOCK9 = 3};
enum {CVZBLOCK9 = 1};

#include "lb_d2q9.h"

#ifdef _D2Q9_

extern const int xblocklen_cv[CVXBLOCK9];
extern const int xdisp_fwd_cv[CVXBLOCK9];
extern const int xdisp_bwd_cv[CVXBLOCK9];

extern const int yblocklen_cv[CVYBLOCK9];
extern const int ydisp_fwd_cv[CVYBLOCK9];
extern const int ydisp_bwd_cv[CVYBLOCK9];

extern const int zblocklen_cv[CVZBLOCK9];
extern const int zdisp_fwd_cv[CVZBLOCK9];
extern const int zdisp_bwd_cv[CVZBLOCK9];

#endif

#endif
