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
#else
#include <p4est_fused.h>
#include <p4est_bits.h>
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
                             p4est_connect_type_t
                             balance_type,
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
  p4est_balance_ext (*p4est_out, balance_type, init_fn, replace_fn);
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

static int
p4est_gloidx_interval_compar (const void *key, const void *base)
{
  p4est_gloidx_t      my_first = *((const p4est_gloidx_t *) key);
  const p4est_gloidx_t *array = (const p4est_gloidx_t *) base;
  p4est_gloidx_t      this_offset = array[0];
  p4est_gloidx_t      next_offset = array[1];

  if (my_first < this_offset)
    return -1;
  if (my_first >= next_offset)
    return 1;
  return 0;
}

static p4est_gloidx_t *
p4est_fused_prefix_sum (p4est_locidx_t * counts, int mpisize)
{
  int i;
  p4est_gloidx_t offset;
  p4est_gloidx_t *offsets = P4EST_ALLOC(p4est_gloidx_t, mpisize + 1);

  for (i = 0, offset = 0; i < mpisize; i++) {
    p4est_locidx_t count = counts[i];

    offsets[i] = offset;
    offset += count;
  }
  offsets[mpisize] = offset;

  return offsets;
}

static int
p4est_fused_key_locate (const p4est_gloidx_t *array, p4est_gloidx_t key, int mpisize)
{
  ptrdiff_t           diff;
  p4est_gloidx_t     *array_loc;

  array_loc = (p4est_gloidx_t *) bsearch ((void *) &key,
                                          (void *) array,
                                          (size_t) mpisize,
                                          sizeof (p4est_gloidx_t),
                                          p4est_gloidx_interval_compar);
  P4EST_ASSERT (array_loc != NULL);
  diff = array_loc - array;
  P4EST_ASSERT (diff >= 0 && diff < (ptrdiff_t) mpisize);
  P4EST_ASSERT (array[diff] <= key && key < array[diff + 1]);
  return (int) diff;
}

static void
p4est_fused_overlap_compute (const p4est_gloidx_t * offsets,
                             const p4est_gloidx_t key[2],
                             int size,
                             int *p_first, /* the lowest offset range whose first position is >= key[0] */
                             int *p_last)  /* the offsets range that contains key[1] - 1 */
{
  int first, last;

  P4EST_ASSERT (key[1] > key[0]);

  first = last = p4est_fused_key_locate (offsets, key[0], size);
  /* We could search for the last key, but we're optimizing for the case when
   * $(last - first) \in O(1)$ */
  while (offsets[last + 1] < key[1]) {
    last++;
  }
  if (offsets[first] < key[0]) {
    first++;
  }
  else {
    while (first > 0 && offsets[first - 1] >= key[0]) {
      first--;
    }
  }
  P4EST_ASSERT (offsets[first] >= key[0]);
  P4EST_ASSERT (first == 0 || offsets[first - 1] < key[0]);
  P4EST_ASSERT (offsets[last] < key[1]);
  P4EST_ASSERT (offsets[last + 1] >= key[1]);
  *p_first = first;
  *p_last = last;
}

#if 0
{
  int                 i, j;
  int                 mpiret;
  int                 size = *p_last + 1 - *p_first;
  int                 QMAXLEVEL;
  p4est_topidx_t      which_tree;
  p4est_tree_t        tree;
  const p4est_topidx_t first_local_tree = p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = p4est->last_local_tree;
  p4est_gloidx_t     *post_offsets;
  p4est_gloidx_t     *pre_offsets;
  p4est_gloidx_t     *local_tree_last_quad_index;
  p4est_quadrant_t   *recv_buf;
  p4est_quadrant_t    my_post_first_pos;
  p4est_gloidx_t      my_pre_first, my_pre_last, my_post_first, my_post_last;
  p4est_topidx_t      flt, llt;
  MPI_Request         recv_req;
  MPI_Request        *send_req;

  pre_offsets = p4est->global_first_quadrant;

  my_pre_first = pre_offsets[mpirank];
  my_pre_last = pre_offsets[mpirank + 1] - 1;
  my_post_first = post_offsets[mpirank];
  my_post_last = post_offsets[mpirank + 1] - 1;

  /* Determine which proc's current range contains my new start */
  pre_first = p4est_fused_key_locate (pre_offsets, my_post_first, size);

  mpiret = sc_MPI_Irecv(&my_post_first_pos, sizeof (p4est_quadrant_t), MPI_BYTE, pre_first,
                        P4EST_COMM_FUSED_PART_1, p4est->mpicomm, &recv_req);
  SC_CHECK_MPI (mpiret);

  /* some of the post first positions are in my range: find out which ones
   * those are and send them to their new owners */

  /* TODO: make sure we don't update the global_first_quadrant before this point */
  flt = p4est->first_local_tree;
  llt = p4est->last_local_tree;
  if (my_pre_last >= my_pre_first) { /* if I have currently have a non-empty partition */
    int                 post_first, post_last, nsends, pre_first, p;
    p4est_topidx_t      t;

    /* find the first process whose post-partition starts at or after the
     * start of my current partition */
    post_first = p4est_fused_key_locate (post_offsets, my_pre_first, mpisize);
    if (post_offsets[post_first] < my_pre_first) {
      post_first++;
      P4EST_ASSERT (post_offsets[post_first] >= my_pre_first);
    } else {
      P4EST_ASSERT (post_offsets[post_first] == my_pre_first);
      while (post_first > 0 && post_offsets[post_first - 1] == my_pre_first) {
        post_first--;
      }
    }
    P4EST_ASSERT (post_first == 0 || post_offsets[post_first - 1] < my_pre_first);
    *p_first = post_first;
    /* find the last process whose post-partition starts at or before the
     * end of my current partition */
    *p_last = post_last = p4est_fused_key_locate (post_offsets, my_pre_last, mpisize);
    P4EST_ASSERT (post_offsets[post_last + 1] > my_pre_last);

    nsends = SC_MIN (post_last + 1 - post_first, 0);

    send_req = P4EST_ALLOC (MPI_Request, nsends);

    for (p = post_first, t = flt; p <= post_last; p++) {
      P4EST_ASSERT (post_offsets[p] >= my_pre_first && post_offsets[p] <= my_pre_last);

      if (p && post_offsets[p] == post_offsets[p - 1]) {
        post_first_positions[p] = post_first_positions[p - 1];
      } else {
        p4est_tree_t *tree;
        p4est_locidx_t local_offset = post_offsets[p] - my_pre_first;
        p4est_locidx_t tree_offset;
        p4est_quadran_t *q;

        tree = p4est_tree_array_index (p4est->trees, t);
        tree_offset = local_offset - tree->quadrants_offset;
        while (t < llt && tree_offset < 0) {
          t++;
          tree = p4est_tree_array_index (p4est->trees, t);
          tree_offset = local_offset - tree->quadrants_offset;
        }
        P4EST_ASSERT (t <= llt && tree_offset >= 0 && (size_t) tree_offset < tree->quadrants->elem_count);
        q = p4est_quadrant_array_index (&tree->quadrants, tree_offset);
        p4est_quadrant_first_descendant (q, &post_first_positions[p], P4EST_QMAXLEVEL);
      }
      mpiret = sc_MPI_Isend(&post_first_position[p], sizeof (p4est_quadrant_t), MPI_BYTE,
                            p, P4EST_COMM_FUSED_PART_1, p4est->mpicomm, &send_req[p - post_first]);
      SC_CHECK_MPI (mpiret);
    }
  }
  *post_first_locations = post_first_positions;
  *send_req_p = send_req;
  *post_offsets_p = post_offsets;
}
#endif

/* TODO: put this in an internal header file */
p4est_locidx_t     *p4est_partition_compute (p4est_t * p4est,
                                             int partition_for_coarsening,
                                             p4est_gloidx_t
                                             global_num_quadrants,
                                             p4est_weight_t weight_fn);

static void
p4est_adapt_fused_partition_ghost (p4est_t * p4est,
                                   p4est_ghost_t ghost,
                                   int repartition,
                                   int partition_for_coarsening,
                                   int ghost_layer_width,
                                   p4est_connect_type_t
                                   ghost_type,
                                   p4est_weight_t weight_fn,
                                   p4est_ghost_t ** ghost_out)
{
  if (ghost == NULL) { /* if we don't have a ghost layer, use the existing methods */
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
  else {
    p4est_locidx_t     *post_num_quads_in_proc;
    p4est_gloidx_t     *post_offsets;
    p4est_quadrant_t   *post_first_positions;
    p4est_ghost_t      *pre_ghost = *ghost_out;
    int                 mpisize = p4est->mpisize;
    int                 mpirank = p4est->mpirank;
    p4est_gloidx_t      my_pre_first, next_pre_first;
    const p4est_gloidx_t *pre_offsets;

    post_num_quads_in_proc =
                            p4est_partition_compute (p4est, partition_for_coarsening, -1, weight_fn);
    post_offsets = p4est_fused_prefix_sum (post_num_quads_in_proc, mpisize);
    P4EST_FREE (post_num_quads_in_proc);
    pre_offsets = p4est->global_first_quadrant;

    post_first_positions = P4EST_ALLOC_ZERO (p4est_quadrant_t, mpisize + 1);

    my_pre_first = pre_offsets[mpirank];
    next_pre_first = pre_offsets[mpirank+1];

    if (my_pre_first < next_pre_first) {
      p4est_topidx-t      llt = p4est->last_local_tree;
      p4est_topidx_t      t = p4est->first_local_tree;
      p4est_tree_t       *tree;
      int                 p, p_first, p_last;

      /* From the number of quadrants in each proc in the new partition, compute
       * the first quadrant for procs that I know */
      p4est_fused_overlap_compute (post_offsets, &pre_offsets[mpirank], post_offsets, &p_first, &p_last);
      tree = p4est_tree_array_index (p4est->trees, t);
      for (p = p_first; p <= p_last; p++) {
        p4est_gloidx_t p_offset = post_offsets[p];

        P4EST_ASSERT (p_offset >= pre_offsets[mpirank] && p_offset < pre_offsets[mpirank + 1]);

        if (p && p_offset == post_offsets[p - 1]) {
          post_first_positions[p] = post_first_positions[p - 1];
        } else {
          p4est_locidx_t local_offset = (p4est_locidx_t) p_offset - my_pre_first;
          p4est_locidx_t tree_offset;
          p4est_quadrant_t *q;

          tree_offset = local_offset - tree->quadrants_offset;
          while (t < llt && tree_offset < 0) {
            t++;
            tree++;
            tree_offset = local_offset - tree->quadrants_offset;
          }
          P4EST_ASSERT (t <= llt && tree_offset >= 0 && (size_t) tree_offset < tree->quadrants->elem_count);
          q = p4est_quadrant_array_index (&tree->quadrants, tree_offset);
          p4est_quadrant_first_descendant (q, &post_first_positions[p], P4EST_QMAXLEVEL);
          post_first_positions[p].p.which_tree = t;
        }
      }
    }
    ierr = sc_MPI_Allreduce
  }
}

void
p4est_adapt_fused (p4est_t * p4est,
                   p4est_ghost_t * ghost,
                   const int8_t * adapt_flag,
                   int copy_data,
                   p4est_connect_type_t
                   balance_type,
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
  p4est_balance_ext (*p4est_out, balance_type, init_fn, replace_fn);
  p4est_adapt_fused_partition_ghost (*p4est_out, ghost, repartition,
                                     partition_for_coarsening,
                                     ghost_layer_width, ghost_type, weight_fn,
                                     ghost_out);

  (*p4est_out)->user_pointer = orig_ctx;
}
