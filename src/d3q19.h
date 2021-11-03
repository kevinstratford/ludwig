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

#ifdef _D3Q19_

extern const    int cv[NVEL19][3];
extern const double wv[NVEL19];
extern const double norm_[NVEL19];
extern const double ma_[NVEL19][NVEL19];
extern const double mi_[NVEL19][NVEL19];

extern const int xblocklen_cv[CVXBLOCK19];
extern const int xdisp_fwd_cv[CVXBLOCK19];
extern const int xdisp_bwd_cv[CVXBLOCK19];

extern const int yblocklen_cv[CVYBLOCK19];
extern const int ydisp_fwd_cv[CVYBLOCK19];
extern const int ydisp_bwd_cv[CVYBLOCK19];

extern const int zblocklen_cv[CVZBLOCK19];
extern const int zdisp_fwd_cv[CVZBLOCK19];
extern const int zdisp_bwd_cv[CVZBLOCK19];

#define LB_CV_D3Q19(cv) const int8_t cv[NVEL19][3] = {                \
                           { 0,  0,  0},                              \
			   { 1,  1,  0}, { 1,  0,  1}, { 1,  0,  0},  \
			   { 1,  0, -1}, { 1, -1,  0}, { 0,  1,  1},  \
			   { 0,  1,  0}, { 0,  1, -1}, { 0,  0,  1},  \
			   { 0,  0, -1}, { 0, -1,  1}, { 0, -1,  0},  \
			   { 0, -1, -1}, {-1,  1,  0}, {-1,  0,  1},  \
			   {-1,  0,  0}, {-1,  0, -1}, {-1, -1,  0}};

#define LB_NORMALISERS_D3Q19(norm) const double norm[NVEL19] = { \
    1.0, \
    3.0, 3.0, 3.0, \
    9.0/2.0, 9.0, 9.0, 9.0/2.0, 9.0, 9.0/2.0, \
    3.0/4.0, 3.0/2.0, 3.0/2.0, 3.0/2.0, \
    9.0/4.0, 9.0/2.0, 9.0/2.0, 9.0/2.0, \
    1.0/2.0};

#endif

#endif
