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

static void
p4est_fused_overlap_compute (p4est_t *p4est, p4est_locidx_t *post_num_quads_in_proc,
                            int *p_first, /* the id of the first mpi rank that overlaps my current domain */
                            int *p_last,  /* the id of the last mpi rank that overlaps my current domain */
                            p4est_quadrant_t **post_first_locations)
                            /* size: | p_last + 1 - p_first |: for each p in
                             * [p_last, ..., p_first], a quadrant of size
                             * QMAXLEVEL (smallest possible, see
                             * p4est_quadrant_first_descendant()) that is the first
                             * location for that process in my domain in the
                             * post repartitioning distribution (mimicking the global_first_position array) */
{
	int i,j;
	int size = *p_last + 1 - *p_first;
	int QMAXLEVEL;
	p4est_topidx_t which_tree;
	p4est_tree_t tree;
	const p4est_topidx_t first_local_tree = p4est->first_local_tree;
	const p4est_topidx_t last_local_tree = p4est->last_local_tree;
	p4est_gloidx_t *local_tree_last_quad_index;

	sc_array_t		*trees = p4est->trees;

	local_tree_last_quad_index = P4EST_ALLOC (p4est_gloidx_t, trees->elem_count);



	tree = p4est_tree_array_index(trees,which_tree);
	
	for (which_tree = first_local_tree + 1;
		 which_tree <= last_local_tree; ++which_tree){
	    tree = p4est_tree_array_index (trees, which_tree);
    	local_tree_last_quad_index[which_tree] = tree->quadrants.elem_count
	      + local_tree_last_quad_index[which_tree - 1];
	}

	for (i = 0; i < *post_num_quads_in_proc; i++){
		for ( j = 0 ; j < size; j++){
			if (quadran[i] is in proc_j){
				first = first_descendant(QMAXLEVEL);
			}
		p4est_quadrant_first_descendant(*p4est_quadrant_t input, *p4est_quadrant_t output, QMAXLEVEL);
		}
	}
}

/* TODO: put this in an internal header file */
p4est_locidx_t * p4est_partition_compute(p4est_t *p4est, int partition_for_coarsening,
                                         p4est_gloidx_t global_num_quadrants, p4est_weight_t weight_fn);

static void
p4est_adapt_fused_partition_ghost (p4est_t * p4est, int repartition,
                                   int partition_for_coarsening,
                                   int ghost_layer_width,
                                   p4est_connect_type_t
                                   ghost_type,
                                   p4est_weight_t weight_fn,
                                   p4est_ghost_t ** ghost_out)
{
#if 0
  if (repartition) {
    p4est_partition (p4est, partition_for_coarsening, weight_fn);
  }
#else
  p4est_locidx_t *post_num_quads_in_proc;
  int             p_first, p_last;
  p4est_quadrant_t *post_first_locations;

  post_num_quads_in_proc = p4est_partition_compute (p4est, partition_for_coarsening,
                                               -1, weight_fn);

  /* From the number of quadrants in each proc in the new partition, compute
   * the first quadrant */

  p4est_fused_overlap_compute (p4est, post_num_quads_in_proc, &p_first, &p_last,
                               &post_first_locations);


  P4EST_FREE (post_num_quads_in_proc);
#endif
  if (ghost_layer_width > 0) {
    int                 i;

    *ghost_out = p4est_ghost_new (p4est, ghost_type);
    for (i = 1; i < ghost_layer_width; i++) {
      p4est_ghost_expand (p4est, *ghost_out);
    }
  }
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
  p4est_adapt_fused_partition_ghost (*p4est_out, repartition,
                                     partition_for_coarsening,
                                     ghost_layer_width, ghost_type, weight_fn,
                                     ghost_out);

  (*p4est_out)->user_pointer = orig_ctx;
}
