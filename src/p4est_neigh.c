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
#include <p4est_bits.h>
#else
#include <p8est_neigh.h>
#include <p8est_bits.h>
#endif

typedef enum
{
  P4EST_NEIGH_DEFAULT = -1,
  P4EST_NEIGH_BASIC = 0,
  P4EST_NEIGH_MPI,
  P4EST_NEIGH_INSUL,
  P4EST_NEIGH_NUM_METHODS,
}
p4est_neigh_method_t;

enum
{
  P4EST_NEIGH_BASIC_ALL,
  P4EST_NEIGH_BASIC_ALLV,
};

struct p4est_neigh_s
{
  sc_MPI_Comm comm;
  int size;
  int rank;
  int n_neigh;
  int *neigh_procs;
  int *neigh_procs_sorted;
  int *neigh_perm;
  p4est_neigh_method_t method;
};

struct p4est_neigh_req_s
{
  sc_array_t *reqs;
  sc_array_t *msgs;
};

typedef struct p4est_neigh_msg_s
{
  int         source;
  size_t      count;
#if defined(P4EST_ENABLE_MPIMPROBE)
  MPI_Message msg;
#endif
}
p4est_neigh_msg_t;

#if defined(P4EST_ENABLE_MPIMPROBE)
static int
p4est_neigh_msg_cmp (const void *A, const void *B)
{
  const p4est_neigh_msg_t *msgA = (const p4est_neigh_msg_t *) A;
  const p4est_neigh_msg_t *msgB = (const p4est_neigh_msg_t *) B;

  return (msgA->source - msgB->source);
}
#endif

p4est_neigh_t *
p4est_neigh_new (p4est_t *p4est, int n_neigh, const int *neigh_procs)
{
  int mpiret;
  p4est_neigh_t *neigh;

  neigh = P4EST_ALLOC_ZERO (p4est_neigh_t, 1);
  neigh->method = P4EST_NEIGH_BASIC;

  P4EST_ASSERT (n_neigh >= 0);

#if !defined(P4EST_ENABLE_MPI)
  P4EST_ASSERT (n_neigh == 0 || n_neigh == 1);
  P4EST_ASSERT (n_neigh == 0 || neigh_procs[0] == 0);
#endif

  if (neigh->method == P4EST_NEIGH_MPI) {
    /* TODO: create distributed graph MPI_Dist_graph_create_adjacent */
  }
  else {
    mpiret = sc_MPI_Comm_dup (p4est->mpicomm, &(neigh->comm));
    SC_CHECK_MPI (mpiret);
  }
  mpiret = sc_MPI_Comm_size (neigh->comm, &(neigh->size));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (neigh->comm, &(neigh->rank));
  SC_CHECK_MPI (mpiret);

  neigh->n_neigh = n_neigh;
  neigh->neigh_procs = P4EST_ALLOC (int, n_neigh);
  memcpy (neigh->neigh_procs, neigh_procs, n_neigh * sizeof(int));

  {
    sc_array_t neigh_array;
    int        is_sorted;

    sc_array_init_data (&neigh_array, neigh->neigh_procs, sizeof(int),
                        n_neigh);

    is_sorted = sc_array_is_sorted (&neigh_array, sc_int_compare);

    if (!is_sorted) {
      int *neigh_sorter = P4EST_ALLOC (int, 2 * n_neigh), nq;
      neigh->neigh_procs_sorted = P4EST_ALLOC (int, n_neigh);
      neigh->neigh_perm = P4EST_ALLOC (int, n_neigh);

      for (nq = 0; nq < n_neigh; nq++) {
        neigh_sorter[2 * nq] = neigh->neigh_procs[nq];
        neigh_sorter[2 * nq + 1] = nq;
      }
      memcpy (neigh->neigh_procs_sorted, neigh->neigh_procs, n_neigh * sizeof(int));
      sc_array_init_data (&neigh_array, neigh_sorter, 2 * sizeof(int),
                          n_neigh);
      sc_array_sort (&neigh_array, sc_int_compare);
      for (nq = 0; nq < n_neigh; nq++) {
        int orig = neigh_sorter[2 * nq + 1];

        neigh->neigh_procs_sorted[nq] = neigh_sorter[2 * nq];
        neigh->neigh_perm[orig] = nq;
      }
      P4EST_FREE (neigh_sorter);
    }
    else {
      neigh->neigh_procs_sorted = neigh->neigh_procs;
    }
  }

#if defined(P4EST_ENABLE_DEBUG)
  {
    int *neigh_in, *neigh_out, nq, diff;

    neigh_in = P4EST_ALLOC_ZERO(int, neigh->size);
    neigh_out = P4EST_ALLOC(int, neigh->size);

    for (nq = 0; nq < neigh->n_neigh; nq++) {
      neigh_in[neigh->neigh_procs[nq]] = 1;
    }

    mpiret = sc_MPI_Alltoall (neigh_in, 1, MPI_INT, neigh_out, 1, MPI_INT,
                              neigh->comm);
    SC_CHECK_MPI (mpiret);

    diff = memcmp (neigh_in, neigh_out, neigh->n_neigh * sizeof(int));

    P4EST_ASSERT (!diff);
    P4EST_FREE (neigh_in);
    P4EST_FREE (neigh_out);
  }
#endif

  return neigh;
}

void
p4est_neigh_destroy (p4est_neigh_t *neigh)
{
  int mpiret;

  mpiret = sc_MPI_Comm_free(&(neigh->comm));
  SC_CHECK_MPI (mpiret);
  P4EST_FREE (neigh->neigh_perm);
  if (neigh->neigh_procs_sorted != neigh->neigh_procs) {
    P4EST_FREE (neigh->neigh_procs_sorted);
  }
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

/* == BASIC == */

static void
p4est_neigh_all_basic_begin (p4est_neigh_t *neigh,
                             sc_array_t *send_array,
                             sc_array_t *recv_array,
                             int ordered, int alltoall,
                             p4est_neigh_req_t **neigh_req_p)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  int *neigh_procs = neigh->neigh_procs;
  size_t nreqs = n_neigh * 2, z;
  sc_MPI_Request * reqs = NULL;
  p4est_neigh_req_t *neigh_req;
  int nq;
  int mpiret;

  neigh_req = P4EST_ALLOC (p4est_neigh_req_t, 1);
  *neigh_req_p = neigh_req;
  neigh_req->reqs = sc_array_new_count (sizeof (sc_MPI_Request *), (size_t) nreqs);
  if (!n_neigh) {
    return;
  }
  reqs = (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0);
  for (z = 0; z < nreqs; z++) {
    reqs[z] = sc_MPI_REQUEST_NULL;
  }
  for (nq = 0; nq < n_neigh; nq++) {
    size_t idx = alltoall ? nq: 0;
    mpiret = sc_MPI_Isend (sc_array_index (send_array, idx),
                           (int) send_array->elem_size,
                           sc_MPI_BYTE,
                           neigh_procs[nq],
                           P4EST_NEIGH_BASIC_ALL,
                           comm,
                           &reqs[2 * nq]);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Irecv (sc_array_index (recv_array, nq),
                           (int) recv_array->elem_size,
                           sc_MPI_BYTE,
                           neigh_procs[nq],
                           P4EST_NEIGH_BASIC_ALL,
                           comm,
                           &reqs[2 * nq + 1]);
    SC_CHECK_MPI (mpiret);
  }
}

static void
p4est_neigh_all_basic_end (p4est_neigh_t *neigh,
                           sc_array_t *send_array,
                           sc_array_t *recv_array,
                           int ordered, int alltoall,
                           p4est_neigh_req_t **neigh_req_p)
{
  p4est_neigh_req_t *neigh_req = *neigh_req_p;
  int mpiret, nreqs = neigh_req->reqs->elem_count;

  if (nreqs) {
    mpiret = sc_MPI_Waitall (nreqs, (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0), sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_destroy (neigh_req->reqs);
  P4EST_FREE (neigh_req);
  *neigh_req_p = NULL;
}

static void
p4est_neigh_all_basic (p4est_neigh_t *neigh,
                       sc_array_t *send_array,
                       sc_array_t *recv_array,
                       int ordered, int alltoall)
{
  p4est_neigh_req_t *req;

  p4est_neigh_all_basic_begin (neigh, send_array, recv_array, ordered, alltoall, &req);
  p4est_neigh_all_basic_end (neigh, send_array, recv_array, ordered, alltoall, &req);
}

static void
p4est_neigh_allv_basic_middle_known (p4est_neigh_t *neigh,
                                     sc_array_t *recv_array,
                                     sc_array_t *recv_offsets,
                                     p4est_neigh_req_t **neigh_req_p);

#if defined(P4EST_ENABLE_MPIMPROBE)

static void
p4est_neigh_allv_basic_mprobe_middle_unknown (p4est_neigh_t *neigh,
                                              sc_array_t *recv_array,
                                              sc_array_t *recv_offsets,
                                              p4est_neigh_req_t **neigh_req_p);

static void
p4est_neigh_allv_basic_mprobe_begin (p4est_neigh_t *neigh,
                                     sc_array_t *send_array,
                                     sc_array_t *send_offsets,
                                     sc_array_t *recv_array,
                                     sc_array_t *recv_offsets,
                                     int alltoall,
                                     p4est_neigh_req_t **neigh_req_p)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  int *neigh_procs = neigh->neigh_procs;
  size_t nmsgs = n_neigh;
  size_t nreqs = n_neigh * 2, z;
  sc_MPI_Request * reqs = NULL;
  size_t *offsets = NULL;
  int nq;
  int mpiret;
  size_t send_count = 0;
  size_t elem_size = send_array->elem_size;
  char *send_buf = NULL;
  p4est_neigh_req_t *neigh_req;

  neigh_req = P4EST_ALLOC (p4est_neigh_req_t, 1);
  *neigh_req_p = neigh_req;
  neigh_req->reqs = sc_array_new_count (sizeof (sc_MPI_Request), (size_t) nreqs);
  neigh_req->msgs = sc_array_new_count (sizeof (p4est_neigh_msg_t), (size_t) nmsgs);

  if (!n_neigh) {
    return;
  }
  reqs = (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0);
  for (z = 0; z < nreqs; z++) {
    reqs[z] = sc_MPI_REQUEST_NULL;
  }
  if (alltoall) {
    offsets = sc_array_index (send_offsets, 0);
    send_buf = send_array->elem_count ? sc_array_index (send_array, 0) : NULL;
  }
  if (n_neigh && send_array->elem_count && !alltoall) {
    send_count = send_array->elem_count;
    send_buf = sc_array_index (send_array, 0);
  }
  /* send messages */
  for (nq = 0; nq < neigh->n_neigh; nq++) {
    char *buf = alltoall ? &send_buf[offsets[nq] * elem_size] : send_buf;
    size_t count = alltoall ? (offsets[nq + 1] - offsets[nq]) * elem_size : send_count * elem_size;

    mpiret = sc_MPI_Isend (buf,
                           count,
                           sc_MPI_BYTE,
                           neigh_procs[nq],
                           P4EST_NEIGH_BASIC_ALLV,
                           comm,
                           &reqs[2 * nq]);
    SC_CHECK_MPI (mpiret);
  }

  if (recv_offsets == NULL || !recv_offsets->elem_count) {
    p4est_neigh_allv_basic_mprobe_middle_unknown (neigh, recv_array, recv_offsets, neigh_req_p);
  }
  else {
    p4est_neigh_allv_basic_middle_known (neigh, recv_array, recv_offsets, neigh_req_p);
  }
}

static void
p4est_neigh_allx_basic_mprobe_begin (p4est_neigh_t *neigh,
                                     sc_array_t **send_arrays,
                                     sc_array_t *recv_array,
                                     sc_array_t *recv_offsets,
                                     p4est_neigh_req_t **neigh_req_p)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  int *neigh_procs = neigh->neigh_procs;
  size_t nmsgs = n_neigh;
  size_t nreqs = n_neigh * 2, z;
  sc_MPI_Request * reqs = NULL;
  int nq;
  int mpiret;
  size_t elem_size = recv_array->elem_size;
  p4est_neigh_req_t *neigh_req;

  neigh_req = P4EST_ALLOC (p4est_neigh_req_t, 1);
  *neigh_req_p = neigh_req;
  neigh_req->reqs = sc_array_new_count (sizeof (sc_MPI_Request), (size_t) nreqs);
  neigh_req->msgs = sc_array_new_count (sizeof (p4est_neigh_msg_t), (size_t) nmsgs);

  if (!n_neigh) {
    return;
  }
  reqs = (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0);
  for (z = 0; z < nreqs; z++) {
    reqs[z] = sc_MPI_REQUEST_NULL;
  }
  /* send messages */
  for (nq = 0; nq < neigh->n_neigh; nq++) {
    char *buf = send_arrays[nq]->array;
    size_t count = send_arrays[nq]->elem_count * elem_size;

    mpiret = sc_MPI_Isend (buf,
                           count,
                           sc_MPI_BYTE,
                           neigh_procs[nq],
                           P4EST_NEIGH_BASIC_ALLV,
                           comm,
                           &reqs[2 * nq]);
    SC_CHECK_MPI (mpiret);
  }

  if (recv_offsets == NULL || !recv_offsets->elem_count) {
    p4est_neigh_allv_basic_mprobe_middle_unknown (neigh, recv_array, recv_offsets, neigh_req_p);
  }
  else {
    p4est_neigh_allv_basic_middle_known (neigh, recv_array, recv_offsets, neigh_req_p);
  }
}

static void
p4est_neigh_allv_basic_mprobe_middle_unknown (p4est_neigh_t *neigh,
                                              sc_array_t *recv_array,
                                              sc_array_t *recv_offsets,
                                              p4est_neigh_req_t **neigh_req_p)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  size_t nmsgs = n_neigh;
  sc_MPI_Request * reqs = NULL;
  size_t *offsets = NULL;
  p4est_neigh_msg_t *msgs = NULL;
  int nq;
  int mpiret;
  int sorted = (recv_offsets != NULL);
  size_t off;
  size_t recv_count;
  size_t elem_size = recv_array->elem_size;
  char *recv_buf = NULL;
  int *neigh_perm = neigh->neigh_perm;
  p4est_neigh_req_t *neigh_req = *neigh_req_p;

  reqs = (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0);
  msgs = (p4est_neigh_msg_t *) sc_array_index (neigh_req->msgs, 0);

  offsets = recv_offsets ? (size_t *) sc_array_push_count (recv_offsets, n_neigh + 1) : NULL;
  /* get sizes with mprobe */
  for (nq = 0, recv_count = 0; nq < n_neigh; nq++) {
    MPI_Status status;
    int        count;

    mpiret = MPI_Mprobe (sc_MPI_ANY_SOURCE,
                         P4EST_NEIGH_BASIC_ALLV,
                         comm,
                         &(msgs[nq].msg),
                         &status);
    SC_CHECK_MPI (mpiret);
    msgs[nq].source = status.MPI_SOURCE;
    mpiret = MPI_Get_count (&status, MPI_BYTE, &count);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (count % elem_size == 0);
    count = count / elem_size;
    msgs[nq].count = (size_t) count;
    recv_count += (size_t) count;
  }
  off = recv_array->elem_count;
  (void) sc_array_push_count (recv_array, recv_count);
  if (recv_array->elem_count) {
    recv_buf = recv_array->array;
  }
  if (sorted) {
    sc_array_t msg_sorter;

    sc_array_init_data (&msg_sorter, msgs, sizeof (p4est_neigh_msg_t),
                        nmsgs);
    sc_array_sort (neigh_req->msgs, p4est_neigh_msg_cmp);
  }
  for (nq = 0; nq < n_neigh; nq++) {
    int q = (sorted && neigh_perm) ? neigh_perm[nq] : nq;
    p4est_neigh_msg_t *msg = &msgs[q];

    P4EST_ASSERT (!sorted || msg->source == neigh->neigh_procs[nq]);
    mpiret = MPI_Imrecv (&recv_buf[off * elem_size],
                         msg->count * elem_size, sc_MPI_BYTE, &(msg->msg),
                         &reqs[2 * nq + 1]);
    SC_CHECK_MPI (mpiret);
    if (offsets) {
      offsets[nq] = off;
    }
    off += msg->count;
  }
  if (offsets) {
    offsets[n_neigh] = off;
  }
  P4EST_ASSERT (off == recv_array->elem_count);
}

static void
p4est_neigh_allv_basic_mprobe_end (p4est_neigh_t *neigh,
                                   sc_array_t *send_array,
                                   sc_array_t *send_offsets,
                                   sc_array_t *recv_array,
                                   sc_array_t *recv_offsets,
                                   int alltoall,
                                   p4est_neigh_req_t **neigh_req_p)
{
  p4est_neigh_req_t *neigh_req = *neigh_req_p;
  int mpiret, nreqs = (int) neigh_req->reqs->elem_count;

  if (nreqs) {
    mpiret = sc_MPI_Waitall (nreqs, (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0), sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_destroy (neigh_req->reqs);
  sc_array_destroy (neigh_req->msgs);
  P4EST_FREE (neigh_req);
  *neigh_req_p = NULL;
}

#else

static void
p4est_neigh_allv_basic_counts_begin (p4est_neigh_t *neigh,
                                     sc_array_t *send_array,
                                     sc_array_t *send_offsets,
                                     sc_array_t *recv_array,
                                     sc_array_t *recv_offsets,
                                     int alltoall,
                                     p4est_neigh_req_t **neigh_req_p)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  int *neigh_procs = neigh->neigh_procs;
  size_t nreqs = n_neigh * 2, z;
  sc_MPI_Request * reqs;
  size_t *soffsets = NULL;
  size_t *roffsets = NULL;
  int nq;
  int mpiret;
  int sorted = (recv_offsets != NULL);
  size_t off;
  size_t send_count = 0, recv_count;
  size_t elem_size = send_array->elem_size;
  char *send_buf = NULL;
  char *recv_buf = NULL;
  sc_array_t *send_counts, *recv_counts;
  p4est_neigh_req_t *neigh_req, *req_inner;

  send_counts = sc_array_new_size (sizeof (size_t), alltoall ? (size_t) n_neigh : 1);
  recv_counts = sc_array_new_size (sizeof (size_t), (size_t) n_neigh);
  if (alltoall) {
    for (nq = 0; nq < n_neigh; nq++) {
      size_t offstart = *((size_t *) sc_array_index (send_offsets, nq));
      size_t offend = *((size_t *) sc_array_index (send_offsets, nq+1));

      *((size_t *) sc_array_index (send_counts, nq)) = offend - offstart;
    }
  }
  else {
    *((size_t *) sc_array_index (send_counts, 0)) = send_array->elem_count;
  }

  p4est_neigh_all_basic_begin (neigh, send_counts, recv_counts, 1, alltoall, &req_inner);
  p4est_neigh_all_basic_end (neigh, send_counts, recv_counts, 1, alltoall, &req_inner);

  for (nq = 0, recv_count = 0; nq < n_neigh; nq++) {
    recv_count += *((size_t *) sc_array_index (recv_counts, nq));
  }

  neigh_req = P4EST_ALLOC (p4est_neigh_req_t, 1);
  *neigh_req_p = neigh_req;

  neigh_req->reqs = sc_array_new_count (sizeof (sc_MPI_Request *), (size_t) nreqs);
  if (!n_neigh) {
    sc_array_destroy (recv_counts);
    sc_array_destroy (send_counts);
    return;
  }

  reqs = (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0);
  for (z = 0; z < nreqs; z++) {
    reqs[z] = sc_MPI_REQUEST_NULL;
  }
  if (alltoall) {
    soffsets = sc_array_index (send_offsets, 0);
  }
  if (sorted) {
    roffsets = sc_array_index (recv_offsets, 0);
  }
  if (send_array->elem_count) {
    send_buf = send_array->array;
  }
  off = recv_array->elem_count;
  (void) sc_array_push_count (recv_array, recv_count);
  if (recv_array->elem_count) {
    recv_buf = recv_array->array;
  }
  /* send / recv messages */
  for (nq = 0; nq < neigh->n_neigh; nq++) {
    char *sbuf = alltoall ? &send_buf[soffsets[nq] * elem_size] : send_buf;
    size_t scount = alltoall ? (soffsets[nq + 1] - soffsets[nq]) * elem_size : send_count * elem_size;
    char *rbuf = &recv_buf[off * elem_size];
    size_t rcount = *((size_t *) sc_array_index (recv_counts, nq));

    if (scount) {
      mpiret = sc_MPI_Isend (sbuf,
                             scount,
                             sc_MPI_BYTE,
                             neigh_procs[nq],
                             P4EST_NEIGH_BASIC_ALLV,
                             comm,
                             &reqs[2 * nq]);
      SC_CHECK_MPI (mpiret);
    }
    if (rcount) {
      mpiret = sc_MPI_Irecv (rbuf,
                             rcount,
                             sc_MPI_BYTE,
                             neigh_procs[nq],
                             P4EST_NEIGH_BASIC_ALLV,
                             comm,
                             &reqs[2 * nq + 1]);
      SC_CHECK_MPI (mpiret);

      if (recv_offsets) {
        roffsets[nq] = off;
      }
    }
    off += rcount;
  }
  if (recv_offsets) {
    roffsets[nq] = off;
  }
  P4EST_ASSERT (off == recv_array->elem_count);

  sc_array_destroy (recv_counts);
  sc_array_destroy (send_counts);
}

static void
p4est_neigh_allx_basic_counts_begin (p4est_neigh_t *neigh,
                                     sc_array_t **send_arrays,
                                     sc_array_t *recv_array,
                                     sc_array_t *recv_offsets,
                                     p4est_neigh_req_t **neigh_req_p)
{
  sc_MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh;
  int *neigh_procs = neigh->neigh_procs;
  size_t nreqs = n_neigh * 2, z;
  sc_MPI_Request * reqs;
  size_t *roffsets = NULL;
  int nq;
  int mpiret;
  int sorted = (recv_offsets != NULL);
  size_t off;
  size_t recv_count;
  size_t elem_size = recv_array->elem_size;
  char *recv_buf = NULL;
  sc_array_t *send_counts, *recv_counts;
  p4est_neigh_req_t *neigh_req, *req_inner;

  send_counts = sc_array_new_size (sizeof (size_t), (size_t) n_neigh);
  recv_counts = sc_array_new_size (sizeof (size_t), (size_t) n_neigh);
  for (nq = 0; nq < n_neigh; nq++) {
    *((size_t *) sc_array_index (send_counts, nq)) = send_arrays[nq]->elem_count;
  }

  p4est_neigh_all_basic_begin (neigh, send_counts, recv_counts, 1, 1, &req_inner);
  p4est_neigh_all_basic_end (neigh, send_counts, recv_counts, 1, 1, &req_inner);

  for (nq = 0, recv_count = 0; nq < n_neigh; nq++) {
    recv_count += *((size_t *) sc_array_index (recv_counts, nq));
  }

  neigh_req = P4EST_ALLOC (p4est_neigh_req_t, 1);
  *neigh_req_p = neigh_req;

  neigh_req->reqs = sc_array_new_count (sizeof (sc_MPI_Request *), (size_t) nreqs);
  if (!n_neigh) {
    sc_array_destroy (recv_counts);
    sc_array_destroy (send_counts);
    return;
  }

  reqs = (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0);
  for (z = 0; z < nreqs; z++) {
    reqs[z] = sc_MPI_REQUEST_NULL;
  }
  if (sorted) {
    roffsets = sc_array_index (recv_offsets, 0);
  }
  off = recv_array->elem_count;
  (void) sc_array_push_count (recv_array, recv_count);
  if (recv_array->elem_count) {
    recv_buf = recv_array->array;
  }
  /* send / recv messages */
  for (nq = 0; nq < neigh->n_neigh; nq++) {
    char *sbuf = send_arrays[nq]->array;
    size_t scount = send_arrays[nq]->elem_count * elem_size;
    char *rbuf = &recv_buf[off * elem_size];
    size_t rcount = *((size_t *) sc_array_index (recv_counts, nq));

    if (scount) {
      mpiret = sc_MPI_Isend (sbuf,
                             scount,
                             sc_MPI_BYTE,
                             neigh_procs[nq],
                             P4EST_NEIGH_BASIC_ALLV,
                             comm,
                             &reqs[2 * nq]);
      SC_CHECK_MPI (mpiret);
    }
    if (rcount) {
      mpiret = sc_MPI_Irecv (rbuf,
                             rcount,
                             sc_MPI_BYTE,
                             neigh_procs[nq],
                             P4EST_NEIGH_BASIC_ALLV,
                             comm,
                             &reqs[2 * nq + 1]);
      SC_CHECK_MPI (mpiret);

      if (recv_offsets) {
        roffsets[nq] = off;
      }
    }
    off += rcount;
  }
  if (recv_offsets) {
    roffsets[nq] = off;
  }
  P4EST_ASSERT (off == recv_array->elem_count);

  sc_array_destroy (recv_counts);
  sc_array_destroy (send_counts);
}

static void
p4est_neigh_allv_basic_counts_end (p4est_neigh_t *neigh,
                                   sc_array_t *send_array,
                                   sc_array_t *send_offsets,
                                   sc_array_t *recv_array,
                                   sc_array_t *recv_offsets,
                                   int alltoall,
                                   p4est_neigh_req_t **neigh_req_p)
{
  p4est_neigh_req_t *neigh_req = *neigh_req_p;
  int mpiret, nreqs = (int) neigh_req->reqs->elem_count;

  if (nreqs) {
    mpiret = sc_MPI_Waitall (nreqs, (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0), sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_destroy (neigh_req->reqs);
  P4EST_FREE (neigh_req);
  *neigh_req_p = NULL;
}
#endif

static void
p4est_neigh_allv_basic_middle_known (p4est_neigh_t *neigh,
                                     sc_array_t *recv_array,
                                     sc_array_t *recv_offsets,
                                     p4est_neigh_req_t **neigh_req_p)
{
  MPI_Comm comm = neigh->comm;
  int n_neigh = neigh->n_neigh, nq;
  const int *neigh_procs = neigh->neigh_procs;
  char *recv_buf = recv_array->array;
  size_t *roffsets = (size_t *) sc_array_index (recv_offsets, 0);
  size_t elem_size = recv_array->elem_size;
  sc_MPI_Request *reqs = NULL;
  p4est_neigh_req_t *neigh_req = *neigh_req_p;

  reqs = (sc_MPI_Request *) sc_array_index (neigh_req->reqs, 0);

  for (nq = 0; nq < n_neigh; nq++) {
    size_t offstart = roffsets[nq];
    size_t offend = roffsets[nq+1];
    int count = (offend - offstart) * elem_size;
    char *buf = &recv_buf[offstart * elem_size];
    int mpiret;

    mpiret = sc_MPI_Irecv (buf,
                           count,
                           sc_MPI_BYTE,
                           neigh_procs[nq],
                           P4EST_NEIGH_BASIC_ALLV,
                           comm,
                           &reqs[2 * nq + 1]);
    SC_CHECK_MPI (mpiret);
  }
}

static void
p4est_neigh_allv_basic_begin (p4est_neigh_t *neigh,
                              sc_array_t *send_array,
                              sc_array_t *send_offsets,
                              sc_array_t *recv_array,
                              sc_array_t *recv_offsets,
                              int alltoall,
                              p4est_neigh_req_t **req)
{
#if defined(P4EST_ENABLE_MPIMPROBE)
  p4est_neigh_allv_basic_mprobe_begin (neigh, send_array, send_offsets, recv_array, recv_offsets, alltoall, req);
#else
  p4est_neigh_allv_basic_counts_begin (neigh, send_array, send_offsets, recv_array, recv_offsets, alltoall, req);
#endif
}

static void
p4est_neigh_allx_basic_begin (p4est_neigh_t *neigh,
                              sc_array_t **send_arrays,
                              sc_array_t *recv_array,
                              sc_array_t *recv_offsets,
                              p4est_neigh_req_t **req)
{
#if defined(P4EST_ENABLE_MPIMPROBE)
  p4est_neigh_allx_basic_mprobe_begin (neigh, send_arrays, recv_array, recv_offsets, req);
#else
  p4est_neigh_allx_basic_counts_begin (neigh, send_arrays, recv_array, recv_offsets, req);
#endif
}

static void
p4est_neigh_allv_basic_end (p4est_neigh_t *neigh,
                            sc_array_t *send_array,
                            sc_array_t *send_offsets,
                            sc_array_t *recv_array,
                            sc_array_t *recv_offsets,
                            int alltoall,
                            p4est_neigh_req_t **req)
{
#if defined(P4EST_ENABLE_MPIMPROBE)
  p4est_neigh_allv_basic_mprobe_end (neigh, send_array, send_offsets, recv_array, recv_offsets, alltoall, req);
#else
  p4est_neigh_allv_basic_counts_end (neigh, send_array, send_offsets, recv_array, recv_offsets, alltoall, req);
#endif
}

static void
p4est_neigh_allx_basic_end (p4est_neigh_t *neigh,
                            sc_array_t **send_array,
                            sc_array_t *recv_array,
                            sc_array_t *recv_offsets,
                            p4est_neigh_req_t **req)
{
  p4est_neigh_allv_basic_end (neigh, NULL, NULL, NULL, NULL, 1, req);
}

static void
p4est_neigh_allv_basic (p4est_neigh_t *neigh,
                        sc_array_t *send_array,
                        sc_array_t *send_offsets,
                        sc_array_t *recv_array,
                        sc_array_t *recv_offsets,
                        int alltoall)
{
  p4est_neigh_req_t *req;

  p4est_neigh_allv_basic_begin (neigh, send_array, send_offsets, recv_array, recv_offsets, alltoall, &req);
  p4est_neigh_allv_basic_end (neigh, send_array, send_offsets, recv_array, recv_offsets, alltoall, &req);
}

static void
p4est_neigh_allx_basic (p4est_neigh_t *neigh,
                        sc_array_t **send_arrays,
                        sc_array_t *recv_array,
                        sc_array_t *recv_offsets)
{
  p4est_neigh_req_t *req;

  p4est_neigh_allx_basic_begin (neigh, send_arrays, recv_array, recv_offsets, &req);
  p4est_neigh_allx_basic_end (neigh, send_arrays, recv_array, recv_offsets, &req);
}

/* == INSUL == */

typedef struct
{
  p4est_topidx_t      nt;
  int                 nl;
  int                 u[P4EST_UTRANSFORM];
}
p4est_tree_neigh_t;

typedef struct
{
  int                 offset[P4EST_INSUL + 1];
  sc_array_t          tnarray;
}
p4est_tree_neigh_info_t;

const static int    p4est_insul_lookup[P4EST_INSUL] = {
#ifndef P4_TO_P8
  0, 2, 1, 0, -1, 1, 2, 3, 3,
#else
  0, 0, 1, 4, 4, 5, 2, 1, 3,
  8, 2, 9, 0, -1, 1, 10, 3, 11,
  4, 2, 5, 6, 5, 7, 6, 3, 7,
#endif
};

const static int    p4est_face_to_insul[P4EST_FACES] = {
#ifndef P4_TO_P8
  3, 5, 1, 7,
#else
  12, 14, 10, 16, 4, 22,
#endif
};

#ifdef P4_TO_P8
const static int    p8est_edge_to_insul[P8EST_EDGES] =
  { 1, 7, 19, 25, 3, 5, 21, 23, 9, 11, 15, 17 };
#endif
const static int    p4est_corner_to_insul[P4EST_CHILDREN] = {
#ifndef P4_TO_P8
  0, 2, 6, 8,
#else
  0, 2, 6, 8, 18, 20, 24, 26,
#endif
};

/* Get and store all of the inter-tree transforms for a tree */
static void
p4est_tree_get_conn_info (p4est_t * p4est, p4est_topidx_t t,
                          p4est_tree_neigh_info_t * info)
{
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_topidx_t      nt;
  p4est_quadrant_t    root, neigh;
  int                 i, j, k, l, nl;
  int                 kstart, kend;
  size_t              s, nneigh;
  p4est_corner_info_t ci;
  sc_array_t         *cta = &ci.corner_transforms;
#ifdef P4_TO_P8
  p8est_edge_info_t   ei;
  sc_array_t         *eta = &ei.edge_transforms;
#endif

#ifdef P4_TO_P8
  kstart = 0;
  kend = 3;
#else
  kstart = 1;
  kend = 2;
#endif

  sc_array_init (cta, sizeof (p4est_corner_transform_t));
#ifdef P4_TO_P8
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif

  sc_array_init_size (&(info->tnarray), sizeof (p4est_tree_neigh_t),
                      P4EST_INSUL);
  sc_array_truncate (&(info->tnarray));

  memset (&root, 0, sizeof (p4est_quadrant_t));

  /* loop over the tree insulation layer,
   * cache neighboring trees and transforms */
  for (k = kstart, l = 0; k < kend; k++) {
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++, l++) {
        int                 type = (k != 1) + (j != 1) + (i != 1);
        int                 id = p4est_insul_lookup[l];

        P4EST_ASSERT (l < P4EST_INSUL);

        info->offset[l] = (int) info->tnarray.elem_count;

        if (!type) {
          p4est_tree_neigh_t *tn = sc_array_push (&info->tnarray);

          tn->nt = t;
          tn->nl = P4EST_INSUL / 2;
        }
        else if (type == 1) {
          /* outside face */
          nt =
            p4est_quadrant_face_neighbor_extra (&root, t, id, &neigh, &nl,
                                                conn);
          nl = nl % P4EST_FACES;
          if (nt >= 0) {
            int                 f[P4EST_FTRANSFORM];
            p4est_tree_neigh_t *tn =
              (p4est_tree_neigh_t *) sc_array_push (&(info->tnarray));

            tn->nt = nt;
            tn->nl = p4est_face_to_insul[nl];
            p4est_find_face_transform (conn, t, id, f);
            p4est_face_transform_to_utransform (f, &(tn->u[0]));
          }
        }
#ifdef P4_TO_P8
        else if (type == 2) {
          /* outside edge */
          p4est_tree_neigh_t *tn;

          p8est_find_edge_transform (conn, t, id, &ei);
          nneigh = eta->elem_count;
          tn =
            (p4est_tree_neigh_t *) sc_array_push_count (&(info->tnarray),
                                                        nneigh);
          for (s = 0; s < nneigh; s++) {
            p8est_edge_transform_t *et = p8est_edge_array_index (eta, s);

            tn[s].nt = et->ntree;
            tn[s].nl = p8est_edge_to_insul[et->nedge];
            p8est_edge_transform_to_utransform (et, id, &(tn[s].u[0]));
          }
          sc_array_truncate (eta);
        }
#endif
        else {
          /* outside corner */
          p4est_tree_neigh_t *tn;

          p4est_find_corner_transform (conn, t, id, &ci);
          nneigh = cta->elem_count;
          tn =
            (p4est_tree_neigh_t *) sc_array_push_count (&(info->tnarray),
                                                        nneigh);
          for (s = 0; s < nneigh; s++) {
            p4est_corner_transform_t *ct = p4est_corner_array_index (cta, s);

            tn[s].nt = ct->ntree;
            tn[s].nl = p4est_corner_to_insul[ct->ncorner];
            p4est_corner_transform_to_utransform (ct, id, &(tn[s].u[0]));
          }
        }
      }
    }
  }
  info->offset[P4EST_INSUL] = info->tnarray.elem_count;

  sc_array_reset (cta);
#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
  return;
}

/* Which insulation neighbor of the root does \a q belong to? */
static int
p4est_root_insul (p4est_quadrant_t * q)
{
  p4est_qcoord_t      x[P4EST_DIM];
  int                 i;
  int                 total;

  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  x[0] = (q->x >> P4EST_MAXLEVEL);
  x[1] = (q->y >> P4EST_MAXLEVEL);
#ifdef P4_TO_P8
  x[2] = (q->z >> P4EST_MAXLEVEL);
#endif
  for (i = P4EST_DIM - 1, total = 0; i >= 0; i--) {
    total *= 3;
    P4EST_ASSERT (-1 <= x[i] && x[i] <= 1);
    total += x[i] + 1;
  }
  P4EST_ASSERT (0 <= total && total < P4EST_INSUL);
  return total;
}

/* Get the coarsest forest within the range */
static void
p4est_neigh_insul_get_coarsest (p4est_t *p4est, sc_array_t * quads)
{
  int rank = p4est->mpirank;
  p4est_quadrant_t s = p4est->global_first_position[rank];
  p4est_quadrant_t e = p4est->global_first_position[rank+1];
  p4est_quadrant_t p;
  int8_t lolimit;
  int sid, eid;

  P4EST_ASSERT (quads->elem_size == sizeof (p4est_quadrant_t));

  if (p4est_quadrant_is_equal_piggy (&s, &e)) {
    return;
  }
  /* transform e from first quad of next rank to last quad of this rank */
  if (e.x == 0 && e.y == 0 &&
#if defined(P4_TO_P8)
      e.z == 0 &&
#endif
      1) {
    e.p.which_tree--;
    e.x = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
    e.y = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
#if defined(P4_TO_P8)
    e.z = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
#endif
  }
  else {
    uint64_t morton = p4est_quadrant_linear_id (&e, P4EST_QMAXLEVEL);
    p4est_quadrant_set_morton (&e, P4EST_QMAXLEVEL, morton-1);
  }
  P4EST_ASSERT (p4est_quadrant_compare_piggy (&s, &e) <= 0);

  /* find the coarsest level for quadrants that can be in quads */
  if (s.p.which_tree == e.p.which_tree) {
    p4est_quadrant_t a, as, ae;

    P4EST_QUADRANT_INIT (&a);
    p4est_nearest_common_ancestor (&s, &e, &a);
    p4est_quadrant_first_descendant(&a, &as, P4EST_QMAXLEVEL);
    p4est_quadrant_last_descendant(&a, &ae, P4EST_QMAXLEVEL);
    if (p4est_quadrant_is_equal (&as, &s) &&
        p4est_quadrant_is_equal (&ae, &e)) {
      /* a is equal to the range */
      p4est_quadrant_t *q = p4est_quadrant_array_push (quads);

      *q = a;
      q->p.which_tree = s.p.which_tree;
      return;
    }
    lolimit = a.level + 1;
  }
  else {
    lolimit = 0;
  }

  /* add quads on the upswing, starting from s */
  {
    int8_t hilimit = P4EST_QMAXLEVEL;
    p4est_quadrant_t *q, *qsib;
    int i;

    /* get coarsest quad that starts as s */
    while ((sid = p4est_quadrant_ancestor_id (&s, hilimit)) == 0 && hilimit > lolimit) {
      hilimit--;
    }
    q = p4est_quadrant_array_push (quads);
    P4EST_QUADRANT_INIT (q);
    p4est_quadrant_ancestor (&s, hilimit, q);
    q->p.which_tree = s.p.which_tree;
    p = *q;
    /* add trailing sibling for all levels past lolimit */
    while (hilimit > lolimit) {
      p4est_quadrant_t np;

      qsib = (p4est_quadrant_t *) sc_array_push_count (quads, P4EST_CHILDREN - (sid + 1));
      for (i = sid + 1; i < P4EST_CHILDREN; i++) {
        P4EST_QUADRANT_INIT (&qsib[i - (sid + 1)]);
        p4est_quadrant_sibling (&p, &qsib[i - (sid + 1)], i);
        qsib[i - (sid + 1)].p.which_tree = p.p.which_tree;
      }
      np = p;
      p4est_quadrant_parent (&np, &p);
      sid = p4est_quadrant_child_id (&p);
      hilimit--;
    }
    P4EST_ASSERT (p.level == lolimit);
    P4EST_ASSERT (sid == p4est_quadrant_child_id (&p));
  }
  eid = p4est_quadrant_ancestor_id (&e, lolimit);

  /* add lolimit quads */
  if (s.p.which_tree == e.p.which_tree) {
    int i;
    p4est_quadrant_t np;
    p4est_quadrant_t *qsib = (p4est_quadrant_t *) sc_array_push_count (quads, eid - (sid + 1));

    P4EST_ASSERT (sid < eid);
    for (i = sid + 1; i < eid; i++) {
      P4EST_QUADRANT_INIT (&qsib[i - (sid + 1)]);
      p4est_quadrant_sibling (&p, &qsib[i - (sid + 1)], i);
      qsib[i - (sid + 1)].p.which_tree = s.p.which_tree;
    }
    np = p;
    p4est_quadrant_sibling (&np, &p, eid);
  }
  else {
    p4est_topidx_t i;

    P4EST_ASSERT (lolimit == 0);
    p4est_quadrant_t *qsib = (p4est_quadrant_t *) sc_array_push_count (quads, e.p.which_tree - (s.p.which_tree + 1));
    for (i = 0; i < e.p.which_tree - (s.p.which_tree + 1); i++) {
      P4EST_QUADRANT_INIT (&qsib[i]);
      qsib[i].x = 0;
      qsib[i].y = 0;
#if defined(P4_TO_P8)
      qsib[i].z = 0;
#endif
      qsib[i].level = 0;
      qsib[i].p.which_tree = i + s.p.which_tree + 1;
    }
    p.x = 0;
    p.y = 0;
#if defined(P4_TO_P8)
    p.z = 0;
#endif
    p.level = 0;
    p.p.which_tree = e.p.which_tree;
  }

  P4EST_ASSERT (p.level == lolimit);
  P4EST_ASSERT (p.p.which_tree == s.p.which_tree);
  P4EST_ASSERT (eid == p4est_quadrant_child_id (&p));

  /* add quads on the upswing, starting with the children of p */
  {
    int8_t l, hilimit = P4EST_QMAXLEVEL;
    p4est_quadrant_t *q;

    /* level of coarsest quad that ends at the same location as e */
    while (hilimit > lolimit && p4est_quadrant_ancestor_id (&e, hilimit) == P4EST_CHILDREN - 1) {
      hilimit--;
    }
    /* add leading children of p that don't overlap e */
    for (l = lolimit + 1; l <= hilimit; l++) {
      int id = p4est_quadrant_ancestor_id (&e, l);
      int i;
      p4est_quadrant_t np;

      np = p;
      p4est_quadrant_first_descendant (&np, &p, l);
      p4est_quadrant_t *qsib = (p4est_quadrant_t *) sc_array_push_count (quads, id);

      for (i = 0; i < id; i++) {
        P4EST_QUADRANT_INIT (&qsib[i]);
        p4est_quadrant_sibling (&p, &qsib[i], i);
        qsib[i].p.which_tree = e.p.which_tree;
      }
    }
    /* add the coarsest quad that ends at the same location as e */
    q = p4est_quadrant_array_push (quads);
    P4EST_QUADRANT_INIT (q);
    p4est_quadrant_ancestor (&e, hilimit, q);
    q->p.which_tree = e.p.which_tree;
  }
  P4EST_ASSERT (sc_array_is_sorted (quads, p4est_quadrant_compare_piggy));
}

/* add quad/insul quad pairs to the list for each quad */
static void
p4est_neigh_insul_get_pairs (p4est_t *p4est, sc_array_t *quads, sc_array_t
                             *insul_quads, p4est_tree_neigh_info_t *info)
{
  size_t nquads = quads->elem_count, z;

  for (z = 0; z < nquads; z++) {
    p4est_quadrant_t *q = p4est_quadrant_array_index (quads, z);
    p4est_topidx_t this_tree = q->p.which_tree;
    int l, m, n;
    int nstart, nend;
    p4est_tree_neigh_info_t *tinfo = &info[this_tree - p4est->first_local_tree];
    p4est_qcoord_t h = P4EST_QUADRANT_LEN (q->level);

#if !defined(P4_TO_P8)
    nstart = 1;
    nend = 2;
#else
    nstart = 0;
    nend = 3;
#endif
    for (n = nstart; n < nend; n++) {
      for (m = 0; m < 3; m++) {
        for (l = 0; l < 3; l++) {
          p4est_quadrant_t nh = *q;
          int root_i;

          nh.x += (l - 1) * h;
          nh.y += (m - 1) * h;
#if defined(P4_TO_P8)
          nh.z += (n - 1) * h;
#endif
          root_i = p4est_root_insul (&nh);
          if (root_i == P4EST_INSUL / 2) {
            p4est_quadrant_t * nq = p4est_quadrant_array_push (insul_quads);

            nq[0] = nh;
            nq[0].p.which_tree = this_tree;
            nq[1] = nh;
            nq[1].p.which_tree = this_tree;
          }
          else {
            size_t s, s_start, s_end;

            s_start = tinfo->offset[root_i];
            s_end = tinfo->offset[root_i + 1];
            for (s = s_start; s < s_end; s++) {
              p4est_quadrant_t   *nq = p4est_quadrant_array_push (insul_quads);
              p4est_tree_neigh_t *tn =
                (p4est_tree_neigh_t *) sc_array_index (&(tinfo->tnarray), s);

              p4est_quadrant_utransform (&nh, &nq[0], &(tn->u[0]), 0);
              nq[0].p.which_tree = tn->nt;
              nq[1] = nh;
              nq[1].p.which_tree = this_tree;
            }
          }
        }
      }
    }
  }
}

static void
p4est_neigh_insul_get_overlaps (sc_array_t *quads, sc_array_t *insul_quads, p4est_tree_neigh_info_t *info)
{
  size_t nquads = quads->elem_count, z;
  p4est_quadrant_t *s, *e;
  p4est_quadrant_t a;

  if (!nquads) {
    return;
  }
  s = p4est_quadrant_array_index (quads, 0);
  e = p4est_quadrant_array_index (quads, nquads - 1);
}

/* == INTERFACE == */

/* != 0: early exit */
static int
p4est_neigh_allgather_setup (p4est_neigh_t *neigh,
                             sc_array_t **send_array_p,
                             sc_array_t **recv_array_p,
                             sc_array_t *send_array_view,
                             sc_array_t *recv_array_view)
{
  size_t total_size;

  sc_array_t *send_array = *send_array_p;
  sc_array_t *recv_array = *recv_array_p;
  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (recv_array->elem_count == neigh->n_neigh * send_array->elem_count);

  /* element count should be the same on each process, so
   * exiting with zero count will be done by each process. */
  total_size = send_array->elem_size * send_array->elem_count;
  if (!total_size) {
    return 1;
  }

  /* reshape the arrays so that they have one element per send */
  if (send_array->elem_count > 1) {
    sc_array_init_data (send_array_view, send_array->array, send_array->elem_size * send_array->elem_count, 1);
    sc_array_init_data (recv_array_view, recv_array->array, send_array->elem_size * send_array->elem_count, neigh->n_neigh);
    *send_array_p = send_array = send_array_view;
    *recv_array_p = recv_array = recv_array_view;
  }

  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (send_array->elem_count == 1);
  P4EST_ASSERT (recv_array->elem_count == neigh->n_neigh);

  return 0;
}

void
p4est_neigh_allgather (p4est_neigh_t *neigh,
                       sc_array_t *send_array,
                       sc_array_t *recv_array,
                       int ordered)
{
  int early_exit;
  sc_array_t send_array_view, recv_array_view;

  early_exit = p4est_neigh_allgather_setup (neigh, &send_array, &recv_array,
                                            &send_array_view,
                                            &recv_array_view);
  if (early_exit) {
    return;
  }

  p4est_neigh_all_basic (neigh, send_array, recv_array, ordered, 0);
}

void
p4est_neigh_iallgather_begin (p4est_neigh_t *neigh,
                              sc_array_t *send_array,
                              sc_array_t *recv_array,
                              int ordered,
                              p4est_neigh_req_t **req)
{
  int early_exit;
  sc_array_t send_array_view, recv_array_view;

  early_exit = p4est_neigh_allgather_setup (neigh, &send_array, &recv_array,
                                            &send_array_view,
                                            &recv_array_view);
  if (early_exit) {
    return;
  }

  p4est_neigh_all_basic_begin (neigh, send_array, recv_array, ordered, 0, req);
}

void
p4est_neigh_iallgather_end (p4est_neigh_t *neigh,
                            sc_array_t *send_array,
                            sc_array_t *recv_array,
                            int ordered,
                            p4est_neigh_req_t **req)
{
  int early_exit;
  sc_array_t send_array_view, recv_array_view;

  early_exit = p4est_neigh_allgather_setup (neigh, &send_array, &recv_array,
                                            &send_array_view,
                                            &recv_array_view);
  if (early_exit) {
    return;
  }

  p4est_neigh_all_basic_end (neigh, send_array, recv_array, ordered, 0, req);
}

static void
p4est_neigh_recv_array_offsets_check (p4est_neigh_t *neigh,
                                      sc_array_t *recv_array,
                                      sc_array_t *recv_offsets)
{
  if (recv_offsets == NULL) {
    /* the user does not know what the incoming sizes are (they will be
     * appended to recv_array), and does not care about the order in which
     * they will be appended */
    P4EST_ASSERT (SC_ARRAY_IS_OWNER (recv_array));
  }
  else {
    size_t elem_count;

    /* the array should contain size_t's */
    P4EST_ASSERT (recv_offsets->elem_size == sizeof (size_t));

    elem_count = recv_offsets->elem_count;
    P4EST_ASSERT (elem_count == 0 || elem_count == neigh->n_neigh + 1);
    if (!elem_count) {
      /* the user does not know what the incoming sizes are (they will be
       * appended to recv_array), but does care about the order in which
       * they will be appended */
      P4EST_ASSERT (SC_ARRAY_IS_OWNER (recv_offsets));
      P4EST_ASSERT (SC_ARRAY_IS_OWNER (recv_array));
    }
    else {
#if defined(P4EST_ENABLE_DEBUG)
      int nq;
      /* the user does know what the incoming sizes are,
       * and has provided offsets for them into the array */

      for (nq = 0; nq < neigh->n_neigh; nq++) {
        size_t offstart = *((size_t *) sc_array_index (recv_offsets, nq));
        size_t offend = *((size_t *) sc_array_index (recv_offsets, nq+1));

        P4EST_ASSERT (offend >= offstart);
        P4EST_ASSERT (offend <= recv_array->elem_count);
      }
#endif
    }
  }
}

/* != 0: early exit */
static int
p4est_neigh_allgatherv_setup (p4est_neigh_t *neigh,
                              sc_array_t *send_array,
                              sc_array_t *recv_array,
                              sc_array_t *recv_offsets)
{
  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  p4est_neigh_recv_array_offsets_check (neigh, recv_array, recv_offsets);

  if (!send_array->elem_size) {
    return 1;
  }
  return 0;
}

void
p4est_neigh_allgatherv (p4est_neigh_t *neigh,
                        sc_array_t *send_array,
                        sc_array_t *recv_array,
                        sc_array_t *recv_offsets)
{
  int early_exit;

  early_exit = p4est_neigh_allgatherv_setup (neigh, send_array, recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allv_basic (neigh, send_array, NULL, recv_array, recv_offsets, 0);
}

void
p4est_neigh_iallgatherv_begin (p4est_neigh_t *neigh,
                               sc_array_t *send_array,
                               sc_array_t *recv_array,
                               sc_array_t *recv_offsets,
                               p4est_neigh_req_t **req)
{
  int early_exit;

  early_exit = p4est_neigh_allgatherv_setup (neigh, send_array, recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allv_basic_begin (neigh, send_array, NULL, recv_array, recv_offsets, 0, req);
}

void
p4est_neigh_iallgatherv_end (p4est_neigh_t *neigh,
                             sc_array_t *send_array,
                             sc_array_t *recv_array,
                             sc_array_t *recv_offsets,
                             p4est_neigh_req_t **req)
{
  int early_exit;

  early_exit = p4est_neigh_allgatherv_setup (neigh, send_array, recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allv_basic_end (neigh, send_array, NULL, recv_array, recv_offsets, 0, req);
}

/* != 0: early exit */
static int
p4est_neigh_alltoall_setup (p4est_neigh_t *neigh,
                            sc_array_t **send_array_p,
                            sc_array_t **recv_array_p,
                            sc_array_t *send_array_view,
                            sc_array_t *recv_array_view)
{
  sc_array_t *send_array = *send_array_p;
  sc_array_t *recv_array = *recv_array_p;

  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (recv_array->elem_count == send_array->elem_count);
  P4EST_ASSERT (neigh->n_neigh == 0 || send_array->elem_count % neigh->n_neigh == 0);

  /* The element size should be the same on each process: if it is zero here,
   * it is zero on all processes, there is no communication necessary */
  if (!send_array->elem_size) {
    return 1;
  }

  /* Resize the arrays so that there is one outgoing message for each neighbor */
  if (send_array->elem_count > neigh->n_neigh) {
    size_t mult = send_array->elem_count / neigh->n_neigh;

    sc_array_init_data (send_array_view, send_array->array, send_array->elem_size * mult, neigh->n_neigh);
    sc_array_init_data (recv_array_view, recv_array->array, send_array->elem_size * mult, neigh->n_neigh);
    *send_array_p = send_array = send_array_view;
    *recv_array_p = recv_array = recv_array_view;
  }

  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (recv_array->elem_count == send_array->elem_count);
  P4EST_ASSERT (send_array->elem_count == neigh->n_neigh);

  return 0;
}

void
p4est_neigh_alltoall (p4est_neigh_t *neigh,
                      sc_array_t *send_array,
                      sc_array_t *recv_array,
                      int ordered)
{
  sc_array_t send_array_view, recv_array_view;
  int early_exit;

  early_exit = p4est_neigh_alltoall_setup (neigh, &send_array, &recv_array, &send_array_view, &recv_array_view);

  if (early_exit) {
    return;
  }

  p4est_neigh_all_basic (neigh, send_array, recv_array, ordered, 1);
}

void
p4est_neigh_ialltoall_begin (p4est_neigh_t *neigh,
                             sc_array_t *send_array,
                             sc_array_t *recv_array,
                             int ordered,
                             p4est_neigh_req_t **req)
{
  sc_array_t send_array_view, recv_array_view;
  int early_exit;

  early_exit = p4est_neigh_alltoall_setup (neigh, &send_array, &recv_array, &send_array_view, &recv_array_view);

  if (early_exit) {
    return;
  }

  p4est_neigh_all_basic_begin (neigh, send_array, recv_array, ordered, 1, req);
}

void
p4est_neigh_ialltoall_end (p4est_neigh_t *neigh,
                           sc_array_t *send_array,
                           sc_array_t *recv_array,
                           int ordered,
                           p4est_neigh_req_t **req)
{
  sc_array_t send_array_view, recv_array_view;
  int early_exit;

  early_exit = p4est_neigh_alltoall_setup (neigh, &send_array, &recv_array, &send_array_view, &recv_array_view);

  if (early_exit) {
    return;
  }

  p4est_neigh_all_basic_end (neigh, send_array, recv_array, ordered, 1, req);
}

/* != 0: early exit */
static int
p4est_neigh_alltoallv_setup (p4est_neigh_t *neigh,
                             sc_array_t *send_array,
                             sc_array_t *send_offsets,
                             sc_array_t *recv_array,
                             sc_array_t *recv_offsets)
{
  P4EST_ASSERT (send_array->elem_size == recv_array->elem_size);
  P4EST_ASSERT (send_offsets->elem_size == sizeof (size_t));
  P4EST_ASSERT (send_offsets->elem_count == neigh->n_neigh + 1);
  p4est_neigh_recv_array_offsets_check (neigh, recv_array, recv_offsets);

  if (!send_array->elem_size) {
    return 1;
  }
  return 0;
}

void
p4est_neigh_alltoallv (p4est_neigh_t *neigh,
                       sc_array_t *send_array,
                       sc_array_t *send_offsets,
                       sc_array_t *recv_array,
                       sc_array_t *recv_offsets)
{
  int early_exit;

  early_exit = p4est_neigh_alltoallv_setup (neigh,
                                            send_array, send_offsets,
                                            recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allv_basic (neigh, send_array, send_offsets, recv_array, recv_offsets, 1);
}

void
p4est_neigh_ialltoallv_begin (p4est_neigh_t *neigh,
                              sc_array_t *send_array,
                              sc_array_t *send_offsets,
                              sc_array_t *recv_array,
                              sc_array_t *recv_offsets,
                              p4est_neigh_req_t **req)
{
  int early_exit;

  early_exit = p4est_neigh_alltoallv_setup (neigh,
                                            send_array, send_offsets,
                                            recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allv_basic_begin (neigh, send_array, send_offsets, recv_array, recv_offsets, 1, req);
}

void
p4est_neigh_ialltoallv_end (p4est_neigh_t *neigh,
                            sc_array_t *send_array,
                            sc_array_t *send_offsets,
                            sc_array_t *recv_array,
                            sc_array_t *recv_offsets,
                            p4est_neigh_req_t **req)
{
  int early_exit;

  early_exit = p4est_neigh_alltoallv_setup (neigh,
                                            send_array, send_offsets,
                                            recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allv_basic_end (neigh, send_array, send_offsets, recv_array, recv_offsets, 1, req);
}

static int
p4est_neigh_alltoallx_setup (p4est_neigh_t *neigh,
                             sc_array_t **send_arrays,
                             sc_array_t *recv_array,
                             sc_array_t *recv_offsets)
{
#if defined(P4EST_ENABLE_DEBUG)
  int i;

  for (i = 0; i < neigh->n_neigh; i++) {
    P4EST_ASSERT (send_arrays[i]->elem_size == recv_array->elem_size);
  }
#endif
  p4est_neigh_recv_array_offsets_check (neigh, recv_array, recv_offsets);

  if (!recv_array->elem_size) {
    return 1;
  }
  return 0;
}

void
p4est_neigh_alltoallx (p4est_neigh_t *neigh,
                       sc_array_t **send_arrays,
                       sc_array_t *recv_array,
                       sc_array_t *recv_offsets)
{
  int early_exit;

  early_exit = p4est_neigh_alltoallx_setup (neigh, send_arrays, recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allx_basic (neigh, send_arrays, recv_array, recv_offsets);
}

void
p4est_neigh_ialltoallx_begin (p4est_neigh_t *neigh,
                              sc_array_t **send_arrays,
                              sc_array_t *recv_array,
                              sc_array_t *recv_offsets,
                              p4est_neigh_req_t **req)
{
  int early_exit;

  early_exit = p4est_neigh_alltoallx_setup (neigh, send_arrays, recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allx_basic_begin (neigh, send_arrays, recv_array, recv_offsets, req);
}

void
p4est_neigh_ialltoallx_end (p4est_neigh_t *neigh,
                            sc_array_t **send_arrays,
                            sc_array_t *recv_array,
                            sc_array_t *recv_offsets,
                            p4est_neigh_req_t **req)
{
  int early_exit;

  early_exit = p4est_neigh_alltoallx_setup (neigh, send_arrays, recv_array, recv_offsets);

  if (early_exit) {
    return;
  }

  p4est_neigh_allx_basic_end (neigh, send_arrays, recv_array, recv_offsets, req);
}
