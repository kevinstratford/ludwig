/*****************************************************************************
 *
 *  coords.h
 *
 *  $Id: coords.h,v 1.4 2010-10-15 12:40:02 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef COORDS_H
#define COORDS_H

#include <mpi.h>

#include "targetDP.h"

#define NSYMM 6      /* Elements for general symmetric tensor */

enum cartesian_directions {X, Y, Z};
enum cartesian_neighbours {FORWARD, BACKWARD};
enum upper_triangle {XX, XY, XZ, YY, YZ, ZZ};

HOST void   coords_init(void);
HOST void   coords_finish(void);
HOST void   coords_info(void);
HOST int    N_total(const int);
HOST int    is_periodic(const int);
HOST double L(const int);
HOST double Lmin(const int);
HOST int    cart_rank(void);
HOST int    cart_size(const int);
HOST int    cart_coords(const int);
HOST int    cart_neighb(const int direction, const int dimension);

HOST MPI_Comm cart_comm(void);

HOST void   coords_nlocal(int n[3]);
HOST void   coords_nlocal_offset(int n[3]);
HOST void   coords_nhalo_set(const int nhalo);
HOST int    coords_nhalo(void);
HOST int    coords_ntotal(int ntotal[3]);
HOST void   coords_ntotal_set(const int n[3]);
HOST void   coords_decomposition_set(const int p[3]);
HOST void   coords_reorder_set(const int);
HOST void   coords_periodicity_set(const int p[3]);
HOST int    coords_nsites(void);
HOST int    coords_index(const int ic, const int jc, const int kc);
HOST void   coords_minimum_distance(const double r1[3], const double r2[3],
			       double r12[3]);
HOST void   coords_index_to_ijk(const int index, int coords[3]);
HOST int    coords_strides(int * xs, int * ys, int * zs);

void coords_active_region_radius_set(const double r);
double coords_active_region(const int index);
int    coords_periodic_comm(MPI_Comm * comm);
int    coords_cart_shift(MPI_Comm comm, int dim, int direction, int * rank);
#endif