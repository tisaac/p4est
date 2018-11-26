/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#if !defined(P4_TO_P8)
#include <p4est_neigh.h>
#else
#include <p8est_neigh.h>
#endif

typedef enum
{
  P4EST_NEIGH_DEFAULT = -1,
  P4EST_NEIGH_BASIC = 0,
  P4EST_NEIGH_MPI,
  P4EST_NEIGH_TREE,
  P4EST_NEIGH_NUM_METHODS,
}
p4est_neigh_method_t;

struct p4est_neigh_s
{
  sc_MPI_Comm comm;
  int size;
  int rank;
  int n_neigh;
  int *neigh_procs;
  p4est_neigh_method_t method;
};

p4est_neigh_t *
p4est_neigh_new (p4est_t *p4est, int n_neigh, const int *neigh_procs)
{
  int mpiret;
  p4est_neigh_t *neigh;

  neigh = P4EST_ALLOC_ZERO (p4est_neigh_t, 1);

  P4EST_ASSERT (n_neigh >= 0);

#if !defined(P4EST_ENABLE_MPI)
  P4EST_ASSERT (n_neigh == 0 || n_neigh == 1);
  P4EST_ASSERT (n_neigh == 0 || neigh_procs[0] == 0);
#endif

  mpiret = sc_MPI_Comm_dup (p4est->mpicomm, &(neigh->comm));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (neigh->comm, &(neigh->size));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (neigh->comm, &(neigh->rank));
  SC_CHECK_MPI (mpiret);

  neigh->n_neigh = n_neigh;
  neigh->neigh_procs = P4EST_ALLOC (int, n_neigh);
  memcpy (neigh->neigh_procs, neigh_procs, n_neigh);
  neigh->method = P4EST_NEIGH_BASIC;

  return neigh;
}

void
p4est_neigh_destroy (p4est_neigh_t *neigh)
{
  int mpiret;

  mpiret = sc_MPI_Comm_free(&(neigh->comm));
  SC_CHECK_MPI (mpiret);
  P4EST_FREE (neigh->neigh_procs);

  P4EST_FREE (neigh);
}

void
p4est_neigh_get_procs (p4est_neigh_t *neigh,
                       int *n_neigh,
                       const int **neigh_procs)
{
  if (n_neigh) {
    *n_neigh = neigh->n_neigh;
  }
  if (neigh_procs) {
    *neigh_procs = neigh->neigh_procs;
  }
}

void
p4est_neigh_allgather (p4est_neigh_t *neigh,
                       sc_array_t *send_array,
                       sc_array_t *recv_array,
                       int ordered)
{
}

void
p4est_neigh_allgatherv (p4est_neigh_t *neigh,
                        sc_array_t *send_array,
                        sc_array_t *send_offsets,
                        sc_array_t *recv_buf,
                        sc_array_t *recv_offsets)
{
}

void
p4est_neigh_alltoall (p4est_neigh_t *neigh,
                      sc_array_t *send_array,
                      sc_array_t *recv_array,
                      int ordered)
{
}

void
p4est_neigh_alltoallv (p4est_neigh_t *neigh,
                       sc_array_t *send_array,
                       sc_array_t *send_offsets,
                       sc_array_t *recv_buf,
                       sc_array_t *recv_offsets)
{
}
