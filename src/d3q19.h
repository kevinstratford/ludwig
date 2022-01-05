/*****************************************************************************
 *
 *  d3q19.h
 *
 *  D3Q19 definitions.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel computing Centre
 *
 *  (c) 2008-2021 The University of Edinburgh
 *
 *  Contributing authors:
 *    Kevin Stratford (kevin@epcc.ed.ac.uk)
 *    Ronojoy Adhikari computed this D3Q19 basis.
 *
 *****************************************************************************/

#ifndef LUDWIG_D3Q19_MODEL_H
#define LUDWIG_D3Q19_MODEL_H

enum {NDIM19 = 3};
enum {NVEL19 = 19};
enum {CVXBLOCK19 = 1};
enum {CVYBLOCK19 = 3};
enum {CVZBLOCK19 = 5};

#include "lb_d3q19.h"

#ifdef _D3Q19_

extern const int xblocklen_cv[CVXBLOCK19];
extern const int xdisp_fwd_cv[CVXBLOCK19];
extern const int xdisp_bwd_cv[CVXBLOCK19];

extern const int yblocklen_cv[CVYBLOCK19];
extern const int ydisp_fwd_cv[CVYBLOCK19];
extern const int ydisp_bwd_cv[CVYBLOCK19];

extern const int zblocklen_cv[CVZBLOCK19];
extern const int zdisp_fwd_cv[CVZBLOCK19];
extern const int zdisp_bwd_cv[CVZBLOCK19];

#endif

#endif
