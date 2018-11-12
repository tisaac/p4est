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

#ifdef P4_TO_P8
#include <p8est_fused.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_algorithms.h>
#include <p8est_search.h>
#else
#include <p4est_fused.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_algorithms.h>
#include <p4est_search.h>
#endif /* !P4_TO_P8 */

typedef struct
{
  int                 counter;
  int8_t             *refine_flags;
}
refine_loop_t;

typedef struct
{
  p4est_quadrant_t    last_processed;
  int                 counter_in;
  int                 counter_out;
  const int8_t       *refine_flags;
  int8_t             *rflags_copy;
}
coarsen_loop_t;

typedef struct
{
  refine_loop_t      *refine_loop;
  coarsen_loop_t     *coarsen_loop;
  const char         *viz_name;
}
fusion_ctx_t;

static inline fusion_ctx_t *
p4est_get_fusion_ctx (p4est_t * p4est)
{
  return (fusion_ctx_t *) p4est->user_pointer;
}

static int
refine_in_loop (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quad)
{
  int                 flag;
  fusion_ctx_t       *ctx = p4est_get_fusion_ctx (p4est);
  refine_loop_t      *loop_ctx = ctx->refine_loop;

  flag = (loop_ctx->refine_flags[loop_ctx->counter++]);
  if (flag == P4EST_FUSED_REFINE) {
    return 1;
  }
  return 0;
}

static int
coarsen_in_loop (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quads[])
{
  int                 i;
  fusion_ctx_t       *ctx = p4est_get_fusion_ctx (p4est);
  coarsen_loop_t     *loop_ctx = ctx->coarsen_loop;
  p4est_quadrant_t    testq;
  int                 flag;

  for (i = 0; i < P4EST_CHILDREN; i++) {

    if (quads[i] == NULL) {
      /* this coarsen callback is not being called on a family: coarsening
       * isn't possible, but we need to advance our counters */
      break;
    }
    testq = *quads[i];
    testq.p.which_tree = which_tree;
    if (!loop_ctx->counter_in
        || (p4est_quadrant_compare_piggy (&(loop_ctx->last_processed), &testq)
            < 0)) {
      flag = loop_ctx->refine_flags[loop_ctx->counter_in + i];
      if (flag != P4EST_FUSED_COARSEN) {
        /* One of the siblings does not want to coarsen, coarsening is not
         * possible */
        break;
      }
    }
  }
  if (i == P4EST_CHILDREN) {
    /* all agree to coarsen */
    /* advance the input counter by P4EST_CHILDREN */
    loop_ctx->counter_in += P4EST_CHILDREN;
    testq = *quads[P4EST_CHILDREN - 1];
    testq.p.which_tree = which_tree;
    loop_ctx->last_processed = testq;
    /* the new forest will have one parent in its place, that we do not want
     * to refine */
    loop_ctx->rflags_copy[loop_ctx->counter_out++] = P4EST_FUSED_KEEP;
    return 1;
  }

  /* do not coarsen, copy flags */
  for (i = 0; i < P4EST_CHILDREN; i++) {
    if (quads[i] == NULL) {
      /* no more quadrants */
      break;
    }
    testq = *quads[i];
    testq.p.which_tree = which_tree;
    if (!loop_ctx->counter_in
        || (p4est_quadrant_compare_piggy (&(loop_ctx->last_processed), &testq)
            < 0)) {
      flag = loop_ctx->refine_flags[loop_ctx->counter_in++];
      loop_ctx->last_processed = *quads[i];
      loop_ctx->rflags_copy[loop_ctx->counter_out++] =
        (flag == P4EST_FUSED_REFINE) ? P4EST_FUSED_REFINE : P4EST_FUSED_KEEP;
    }
  }

  return 0;
}

void
p4est_adapt_fused_reference (p4est_t * p4est,
                             const int8_t * adapt_flag,
                             int copy_data,
                             p4est_balance_obj_t * bobj,
                             int repartition,
                             int partition_for_coarsening,
                             int ghost_layer_width,
                             p4est_connect_type_t
                             ghost_type,
                             p4est_weight_t weight_fn,
                             p4est_init_t init_fn,
                             p4est_replace_t replace_fn,
                             p4est_t ** p4est_out, p4est_ghost_t ** ghost_out)
{
  refine_loop_t       refine_ctx;
  coarsen_loop_t      coarsen_ctx;
  fusion_ctx_t        fusion_ctx;
  void               *orig_ctx;
  int8_t             *aflag_copy;
  int                 has_bobj = (bobj != NULL);

  if (*p4est_out != p4est) {
    *p4est_out = p4est_copy (p4est, copy_data);
  }
  orig_ctx = (*p4est_out)->user_pointer;
  (*p4est_out)->user_pointer = (void *) &fusion_ctx;

  aflag_copy = P4EST_ALLOC (int8_t, (size_t) p4est->local_num_quadrants);

  coarsen_ctx.counter_in = 0;
  coarsen_ctx.counter_out = 0;
  coarsen_ctx.refine_flags = adapt_flag;
  coarsen_ctx.rflags_copy = aflag_copy;

  fusion_ctx.coarsen_loop = &coarsen_ctx;
  /* coarsen families where every child's adapt_flag is P4EST_FUSED_COARSEN */
  p4est_coarsen_ext (*p4est_out, 0, 1, coarsen_in_loop, init_fn, replace_fn);

  refine_ctx.counter = 0;
  refine_ctx.refine_flags = aflag_copy;

  fusion_ctx.refine_loop = &refine_ctx;
  /* refine all quadrants marked P4EST_FUSED_REFINE */
  p4est_refine_ext (*p4est_out, 0, P4EST_QMAXLEVEL, refine_in_loop, init_fn,
                    *replace_fn);

  P4EST_FREE (aflag_copy);
  if (!has_bobj) {
    bobj = p4est_balance_obj_new (p4est->mpicomm);
    p4est_balance_obj_set_init (bobj, init_fn);
    p4est_balance_obj_set_replace (bobj, replace_fn);
  }
  p4est_balance_obj (bobj, *p4est_out);
  if (!has_bobj) {
    p4est_balance_obj_destroy (bobj);
  }

  if (repartition) {
    p4est_partition (*p4est_out, partition_for_coarsening, weight_fn);
  }
  if (ghost_layer_width > 0) {
    int                 i;

    *ghost_out = p4est_ghost_new (*p4est_out, ghost_type);
    for (i = 1; i < ghost_layer_width; i++) {
      p4est_ghost_expand (*p4est_out, *ghost_out);
    }
  }

  (*p4est_out)->user_pointer = orig_ctx;
}

static void
p4est_adapt_fused_partition_ghost (p4est_t * p4est, int repartition,
                                   int partition_for_coarsening,
                                   int ghost_layer_width,
                                   p4est_connect_type_t
                                   ghost_type,
                                   p4est_weight_t weight_fn,
                                   p4est_ghost_t ** ghost_out)
{
  if (repartition) {
    p4est_partition (p4est, partition_for_coarsening, weight_fn);
  }
  if (ghost_layer_width > 0) {
    int                 i;

    *ghost_out = p4est_ghost_new (p4est, ghost_type);
    for (i = 1; i < ghost_layer_width; i++) {
      p4est_ghost_expand (p4est, *ghost_out);
    }
  }
}

#if 0
static int
p4est_connectivity_get_max_neighbors (p4est_connectivity_t * conn)
{
  int                 max_corners[P4EST_CHILDREN];
  p4est_topidx_t      num_trees = conn->num_trees, t;
  int                 i;
  p4est_corner_info_t ci;
  sc_array_t         *cta;
  int                 max;
#ifdef P4_TO_P8
  int                 max_edges[P8EST_EDGES];
  p8est_edge_info_t   ei;
  sc_array_t         *eta;

  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));

  for (i = 0; i < P4EST_CHILDREN; i++) {
    max_corners[i] = 1;
  }
#ifdef P4_TO_P8
  for (i = 0; i < P8EST_EDGES; i++) {
    max_edges[i] = 1;
  }
#endif
  for (t = 0; t < num_trees; t++) {
#ifdef P4_TO_P8
    for (i = 0; i < P8EST_EDGES; i++) {
      p8est_find_edge_transform (conn, t, i, &ei);
      max_edges[i] = SC_MAX (max_edges[i], (int) eta->elem_count);
    }
#endif
    for (i = 0; i < P4EST_CHILDREN; i++) {
      p4est_find_corner_transform (conn, t, i, &ci);
      max_corners[i] = SC_MAX (max_corners[i], (int) cta->elem_count);
    }
  }
#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
  sc_array_reset (cta);
  max = P4EST_FACES;
  for (i = 0; i < P4EST_CHILDREN; i++) {
    max += max_corners[i];
  }
#ifdef P4_TO_P8
  for (i = 0; i < P8EST_EDGES; i++) {
    max += max_edges[i];
  }
#endif

  return max;
}
#endif

static void
compute_reduction_neighbors_add (p4est_t * p4est, p4est_quadrant_t * neigh,
                                 p4est_topidx_t nt,
                                 int32_t touch, int single_proc,
                                 int *in_all_procs,
                                 int *out_all_procs,
                                 sc_array_t * in_procs,
                                 sc_array_t * out_procs)
{
  int                 p, np, nq;
  p4est_quadrant_t    first_desc, last_desc;
#ifdef P4_TO_P8
  int32_t             corner_mask =
    (int32_t) (P4EST_CHILDREN - 1) << (P8EST_EDGES + P4EST_FACES);
#else
  int32_t             corner_mask =
    (int32_t) (P4EST_CHILDREN - 1) << P4EST_FACES;
#endif
  int32_t             rb;

  p4est_quadrant_first_descendant (neigh, &first_desc, P4EST_QMAXLEVEL);
  p4est_quadrant_last_descendant (neigh, &last_desc, P4EST_QMAXLEVEL);
  np = p4est_comm_find_owner (p4est, nt, &first_desc, p4est->mpirank);
  if (single_proc) {
    p4est_quadrant_t   *next_first;

    next_first = &(p4est->global_first_position[np + 1]);
    if (next_first->p.which_tree > nt
        || p4est_quadrant_compare (&last_desc, next_first) < 0) {
      if (!in_all_procs[np]) {
        in_all_procs[np] = 1;
        *((int *) sc_array_push (in_procs)) = np;
      }
    }
    return;
  }
  nq = p4est_comm_find_owner (p4est, nt, &last_desc, np);
  if (np == nq) {
    if (!in_all_procs[np]) {
      in_all_procs[np] = 1;
      *((int *) sc_array_push (in_procs)) = np;
    }
    if (!out_all_procs[np]) {
      out_all_procs[np] = 1;
      *((int *) sc_array_push (out_procs)) = np;
    }
    return;
  }
  for (p = np; p <= nq; p++) {
    p4est_quadrant_t   *lq, *uq;
    p4est_quadrant_t    temp;

    lq = &(p4est->global_first_position[p]);
    uq = &(p4est->global_first_position[p + 1]);
    if (p4est_quadrant_is_equal_piggy (lq, uq)) {
      continue;
    }
    if (!out_all_procs[p]) {
      out_all_procs[p] = 1;
      *((int *) sc_array_push (out_procs)) = p;
    }
    if (p == np) {
      lq = NULL;
    }
    if (p == nq) {
      uq = NULL;
    }
    else {
      uint64_t            next_lid, uid;

      next_lid = p4est_quadrant_linear_id (uq, P4EST_QMAXLEVEL);
      uid = next_lid - 1;
      uq = &temp;
      p4est_quadrant_set_morton (uq, P4EST_QMAXLEVEL, uid);
    }
    if (touch & corner_mask) {
      p4est_quadrant_t    cdesc;
#ifdef P4_TO_P8
      int                 corner = (int) touch >> (P8EST_EDGES + P4EST_FACES);
#else
      int                 corner = (int) touch >> P4EST_FACES;
#endif
      P4EST_ASSERT (1 <= corner && corner < (1 << P4EST_CHILDREN));
      corner = SC_LOG2_32 (corner);
      P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);

      p4est_quadrant_corner_descendant (neigh, &cdesc, corner,
                                        P4EST_QMAXLEVEL);
      if ((lq == NULL || p4est_quadrant_compare (lq, &cdesc) <= 0)
          && (uq == NULL || p4est_quadrant_compare (&cdesc, uq) <= 0)) {
        if (!in_all_procs[p]) {
          in_all_procs[p] = 1;
          *((int *) sc_array_push (in_procs)) = p;
        }
      }
    }
    else {
      rb = p4est_find_range_boundaries (lq, uq, (int) neigh->level,
#ifdef P4_TO_P8
                                        NULL,
#endif
                                        NULL, NULL);
      if (rb & touch) {
        if (!in_all_procs[p]) {
          in_all_procs[p] = 1;
          *((int *) sc_array_push (in_procs)) = p;
        }
      }
    }
  }
}

static void
compute_reduction_neighbors (p4est_t * p4est, p4est_quadrant_t * q,
                             int single_proc,
                             int *in_all_procs,
                             int *out_all_procs,
                             sc_array_t * in_procs, sc_array_t * out_procs)
{
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_topidx_t      which_tree = q->p.which_tree, nt;
  int                 i;
  int32_t             touch;
  sc_array_t          quads;
  sc_array_t          trees;
  sc_array_t          points;

  sc_array_init (&quads, sizeof (p4est_quadrant_t));
  sc_array_init (&trees, sizeof (p4est_topidx_t));
  sc_array_init (&points, sizeof (int));

  for (i = 0; i < P4EST_FACES; i++) {
    p4est_quadrant_t    neigh;
    int                 nf = -1;

    nt =
      p4est_quadrant_face_neighbor_extra (q, which_tree, i, &neigh, &nf,
                                          conn);
    if (nt < 0) {
      continue;
    }
    touch = ((int32_t) 1 << nf);

    compute_reduction_neighbors_add (p4est, &neigh, nt, touch, single_proc,
                                     in_all_procs, out_all_procs, in_procs,
                                     out_procs);

  }
#ifdef P4_TO_P8
  for (i = 0; i < P8EST_EDGES; i++) {
    int                 neneigh, j;

    p8est_quadrant_edge_neighbor_extra (q, which_tree, i, &quads, &trees,
                                        &points, conn);
    neneigh = (int) quads.elem_count;
    for (j = 0; j < neneigh; j++) {
      p4est_quadrant_t   *neigh =
        p4est_quadrant_array_index (&quads, (size_t) j);
      int                 ne;

      nt = *((p4est_topidx_t *) sc_array_index (&trees, (size_t) j));
      if (nt < 0) {
        continue;
      }
      ne = *((int *) sc_array_index_int (&points, j));
      touch = ((int32_t) 1 << (ne + P4EST_FACES));

      compute_reduction_neighbors_add (p4est, neigh, nt, touch, single_proc,
                                       in_all_procs, out_all_procs, in_procs,
                                       out_procs);
    }
    sc_array_truncate (&trees);
    sc_array_truncate (&quads);
    sc_array_truncate (&points);
  }
#endif
  for (i = 0; i < P4EST_CHILDREN; i++) {
    int                 ncneigh, j;

    p4est_quadrant_corner_neighbor_extra (q, which_tree, i, &quads, &trees,
                                          &points, conn);
    ncneigh = (int) quads.elem_count;
    for (j = 0; j < ncneigh; j++) {
      p4est_quadrant_t   *neigh =
        p4est_quadrant_array_index (&quads, (size_t) j);
      int                 nc;

      nt = *((p4est_topidx_t *) sc_array_index (&trees, (size_t) j));
      if (nt < 0) {
        continue;
      }
      nc = *((int *) sc_array_index_int (&points, j));
#ifdef P4_TO_P8
      touch = (1 << (nc + P8EST_EDGES + P4EST_FACES));
#else
      touch = (1 << (nc + P4EST_FACES));
#endif

      compute_reduction_neighbors_add (p4est, neigh, nt, touch, single_proc,
                                       in_all_procs, out_all_procs, in_procs,
                                       out_procs);
    }
    sc_array_truncate (&trees);
    sc_array_truncate (&quads);
    sc_array_truncate (&points);
  }
  sc_array_reset (&trees);
  sc_array_reset (&quads);
  sc_array_reset (&points);
}

static void
p4est_adapt_fused_compute_insulation_comm (p4est_t * p4est,
                                           sc_array_t * out_procs,
                                           sc_array_t * in_procs)
{
  int                 mpisize = p4est->mpisize;
  int                 mpirank = p4est->mpirank;
  p4est_quadrant_t    first_desc = p4est->global_first_position[mpirank];
  p4est_quadrant_t    next_desc = p4est->global_first_position[mpirank + 1];
  p4est_quadrant_t    last_desc;
  int                 maxlevel[2];
  p4est_topidx_t      which_tree[2], t;
  p4est_quadrant_t    ancestors[2][P4EST_MAXLEVEL];
  int                 s;
  int                 alevel;
  int                 l;
  int                *out_all = NULL, *in_all = NULL;

  P4EST_GLOBAL_VERBOSE ("Into compute_reduction_network\n");

  out_all = P4EST_ALLOC_ZERO (int, mpisize);
  in_all = P4EST_ALLOC_ZERO (int, mpisize);
  if (!p4est_quadrant_is_equal_piggy (&first_desc, &next_desc)) {
#if 0
    maxsize = p4est_connectivity_get_max_neighbors (p4est->connectivity);
    P4EST_GLOBAL_VERBOSEF ("Maximum number of quadrant neighbors: %d\n",
                           maxsize);
#endif

    /* mark remote in */
    if (next_desc.x == 0 && next_desc.y == 0 &&
#ifdef P4_TO_P8
        next_desc.z == 0 &&
#endif
        1) {
      p4est_quadrant_t    root;

      memset (&root, 0, sizeof (p4est_quadrant_t));
      p4est_quadrant_last_descendant (&root, &last_desc, P4EST_QMAXLEVEL);
      last_desc.p.which_tree = next_desc.p.which_tree - 1;
    }
    else {
      uint64_t            id;

      id = p4est_quadrant_linear_id (&next_desc, P4EST_QMAXLEVEL);
      p4est_quadrant_set_morton (&last_desc, P4EST_QMAXLEVEL, id - 1);
      last_desc.p.which_tree = next_desc.p.which_tree;
    }
    which_tree[0] = first_desc.p.which_tree;
    which_tree[1] = last_desc.p.which_tree;
    ancestors[0][P4EST_QMAXLEVEL] = first_desc;
    ancestors[1][P4EST_QMAXLEVEL] = last_desc;
    if (which_tree[0] < which_tree[1]) {
      alevel = -1;
    }
    else {
      p4est_quadrant_t    a;

      p4est_nearest_common_ancestor (&first_desc, &last_desc, &a);
      alevel = a.level;
    }

    for (s = 0; s < 2; s++) {
      int                 mask = (s == 0) ? 0 : P4EST_CHILDREN - 1;
      int                 limit = s ? alevel + 1 : 0;

      /* compute all ancestors of first and last quadrants */
      for (l = P4EST_QMAXLEVEL - 1; l >= limit; l--) {
        p4est_quadrant_parent (&ancestors[s][l + 1], &ancestors[s][l]);
        ancestors[s][l].p.which_tree = which_tree[s];
      }
      if (s) {
        for (l = 0; l <= alevel; l++) {
          ancestors[1][l] = ancestors[0][l];
        }
      }
      for (l = P4EST_QMAXLEVEL; l > SC_MAX (0, alevel); l--) {
        int                 cid = p4est_quadrant_child_id (&ancestors[s][l]);

        if (cid ^ mask) {
          break;
        }
      }
      /* maxlevel[s] is the coarsest refinement of that endpoint
       * that is entirely in the process */
      maxlevel[s] = l;
    }
    if (which_tree[0] == which_tree[1]
        && (maxlevel[0] == alevel || maxlevel[1] == alevel)) {
      if (maxlevel[0] != alevel) {
        maxlevel[1] = alevel + 1;
      }
      else if (maxlevel[1] != alevel) {
        maxlevel[0] = alevel + 1;
      }
    }
    for (s = 0; s < 2; s++) {
      int                 limit = s ? alevel + 1 : 0;
      int                 stride = s ? -1 : 1;
      p4est_quadrant_t    quad = ancestors[s][maxlevel[s]];
      p4est_quadrant_t   *stop =
        (which_tree[0] ==
         which_tree[1]) ? &ancestors[s ^ 1][P4EST_QMAXLEVEL] : NULL;

      for (l = limit; l < maxlevel[s]; l++) {
        compute_reduction_neighbors (p4est, &ancestors[s][l], 1, in_all,
                                     out_all, in_procs, out_procs);
      }
      for (;;) {
        p4est_quadrant_t    temp;
        int                 cid;

        quad.p.which_tree = which_tree[s];
        compute_reduction_neighbors (p4est, &quad, 0, in_all, out_all,
                                     in_procs, out_procs);
        cid = p4est_quadrant_child_id (&quad) + stride;
        while (quad.level > SC_MAX (0, alevel)
               && (cid < 0 || cid >= P4EST_CHILDREN)) {

          p4est_quadrant_parent (&quad, &temp);
          quad = temp;
          cid = p4est_quadrant_child_id (&quad) + stride;
        }
        if (quad.level == SC_MAX (0, alevel)) {
          break;
        }
        p4est_quadrant_sibling (&quad, &temp, cid);
        quad = temp;
        if (stop && (p4est_quadrant_overlaps (&quad, stop))) {
          break;
        }
      }
    }
    for (t = which_tree[0] + 1; t < which_tree[1]; t++) {
      p4est_quadrant_t    root;

      root.x = 0;
      root.y = 0;
#ifdef P4_TO_P8
      root.z = 0;
#endif
      root.level = 0;
      root.p.which_tree = t;
      compute_reduction_neighbors (p4est, &root, 0, in_all, out_all, in_procs,
                                   out_procs);
    }
    sc_array_sort (in_procs, sc_int_compare);
    sc_array_sort (out_procs, sc_int_compare);
    sc_array_uniq (in_procs, sc_int_compare);
    sc_array_uniq (out_procs, sc_int_compare);
  }
#ifdef P4EST_ENABLE_DEBUG
  {
    int                *out_all2 = P4EST_ALLOC (int, mpisize);
    int                 mpiret, i;

    mpiret =
      sc_MPI_Alltoall (in_all, 1, sc_MPI_INT, out_all2, 1, sc_MPI_INT,
                       p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    for (i = 0; i < mpisize; i++) {
      P4EST_ASSERT (out_all2[i] == out_all[i]);
    }
    mpiret = sc_MPI_Barrier (p4est->mpicomm);

    P4EST_FREE (out_all2);
  }
#endif
  P4EST_FREE (in_all);
  P4EST_FREE (out_all);
  P4EST_GLOBAL_VERBOSE ("Done compute_reduction_network\n");
}

static void
p4est_superset_callback (sc_array_t * receivers, sc_array_t * extra_receivers,
                         sc_array_t * super_senders, sc_notify_t * notify,
                         void *ctx)
{
  int                 i, j, uniq;
  int                *extra;
  int                *recv;
  int                 num_recv, num_all;
  sc_array_t         *recv_sort;
#ifdef P4EST_ENABLE_DEBUG
  sc_array_t         *extra_copy;
#endif
  sc_flopinfo_t       snap;

  p4est_t            *p4est = (p4est_t *) ctx;

  P4EST_FUNC_SNAP (p4est, &snap);

  if (!sc_array_is_sorted (receivers, sc_int_compare)) {
    recv_sort = sc_array_new_count (sizeof (int), receivers->elem_count);
    sc_array_copy (recv_sort, receivers);
    sc_array_sort (recv_sort, sc_int_compare);
    sc_array_uniq (recv_sort, sc_int_compare);
    P4EST_ASSERT (recv_sort->elem_count == receivers->elem_count);
  }
  else {
    recv_sort = receivers;
  }
  p4est_adapt_fused_compute_insulation_comm (p4est, extra_receivers,
                                             super_senders);
#ifdef P4EST_ENABLE_DEBUG
  extra_copy = sc_array_new_count (sizeof (int), extra_receivers->elem_count);
  sc_array_copy (extra_copy, extra_receivers);
#endif
  num_all = (int) extra_receivers->elem_count;
  num_recv = (int) recv_sort->elem_count;
  extra = (int *) extra_receivers->array;
  recv = (int *) recv_sort->array;
  for (i = 0, j = 0, uniq = 0; i < num_all && j < num_recv;) {
    if (extra[i] < recv[j]) {
      extra[uniq++] = extra[i++];
    }
    else if (extra[i] == recv[j]) {
      i++;
      j++;
    }
    else {
      j++;
    }
  }
  for (; i < num_all; i++) {
    extra[uniq++] = extra[i];
  }
  sc_array_resize (extra_receivers, uniq);
#ifdef P4EST_ENABLE_DEBUG
  for (i = 0; i < num_all; i++) {
    ssize_t             r, e;
    int                 j = *((int *) sc_array_index_int (extra_copy, i));

    r = sc_array_bsearch (recv_sort, &j, sc_int_compare);
    e = sc_array_bsearch (extra_receivers, &j, sc_int_compare);
    P4EST_ASSERT (r < 0 || e < 0);
    P4EST_ASSERT (r >= 0 || e >= 0);
  }
  sc_array_destroy (extra_copy);
#endif
  if (recv_sort != receivers) {
    sc_array_destroy (recv_sort);
  }
  P4EST_FUNC_SHOT (p4est, &snap);
}

void
p4est_adapt_fused (p4est_t * p4est,
                   const int8_t * adapt_flag,
                   int copy_data,
                   p4est_balance_obj_t * bobj,
                   int repartition,
                   int partition_for_coarsening,
                   int ghost_layer_width,
                   p4est_connect_type_t
                   ghost_type,
                   p4est_weight_t weight_fn,
                   p4est_init_t init_fn,
                   p4est_replace_t replace_fn,
                   p4est_t ** p4est_out, p4est_ghost_t ** ghost_out)
{
  p4est_inspect_t     inspect;
  p4est_inspect_t    *inspect_orig;
  sc_notify_t        *notify;
  sc_notify_type_t    type;
  sc_flopinfo_t       snap;
  int                 has_bobj = (bobj != NULL);

  P4EST_FUNC_SNAP (p4est, &snap);

  if (*p4est_out != p4est) {
    *p4est_out = p4est_copy (p4est, copy_data);
  }
  inspect_orig = p4est->inspect;
  if (!inspect_orig) {
    memset (&inspect, 0, sizeof (p4est_inspect_t));
    (*p4est_out)->inspect = &inspect;
  }
  else {
    (*p4est_out)->inspect = inspect_orig;
  }
  if (!has_bobj) {
    bobj = p4est_balance_obj_new (p4est->mpicomm);
    p4est_balance_obj_set_init (bobj, init_fn);
    p4est_balance_obj_set_replace (bobj, replace_fn);
  }
  notify = p4est_balance_obj_get_notify (bobj);
  if (notify) {
    type = sc_notify_get_type (notify);
    if (type == SC_NOTIFY_SUPERSET) {
      sc_notify_superset_set_callback (notify,
                                       p4est_superset_callback,
                                       (void *) p4est);
    }
  }
  p4est_balance_obj_set_adapt_flags (bobj, adapt_flag);
  p4est_balance_obj (bobj, *p4est_out);
  p4est_balance_obj_set_adapt_flags (bobj, NULL);
  if (!has_bobj) {
    p4est_balance_obj_destroy (bobj);
  }
  p4est_adapt_fused_partition_ghost (*p4est_out, repartition,
                                     partition_for_coarsening,
                                     ghost_layer_width, ghost_type, weight_fn,
                                     ghost_out);
  P4EST_FUNC_SHOT (p4est, &snap);
}
