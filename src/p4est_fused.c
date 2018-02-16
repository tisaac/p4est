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

void
p4est_adapt_fused (p4est_t * p4est,
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
