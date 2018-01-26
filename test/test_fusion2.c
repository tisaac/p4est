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
#ifndef P4_TO_P8
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#else
#include <p8est_bits.h>
#include <p8est.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#endif
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

#ifndef P4_TO_P8
static int          refine_level = 5;
#else
static int          refine_level = 4;
#endif

typedef struct
{
  double              time;
  double              radius;
  double              x0[P4EST_DIM];
  double              velocity[P4EST_DIM];
  int                 max_level;
}
fusion_sphere_t;

typedef struct
{
  fusion_sphere_t    *sphere;
  int				  context; // Whether refine or coarsen
}
test_fusion_t;

enum
{
  TIMINGS_INIT,
  TIMINGS_NUM_STATS
};

/* Given in (x,y,[z],t), an initial position, a velocity, and a radius,
 * answer the question: is this point inside the ball? */
static int
in_sphere (const double coordinate[P4EST_DIM], double time, double radius,
           const double x0[P4EST_DIM], const double velocity[P4EST_DIM])
{
  double              xTime[P4EST_DIM]; /* The center at the current time */
  int                 i;
  double              dist = 0.;

  for (i = 0; i < P4EST_DIM; i++) {
    xTime[i] = x0[i] + time * velocity[i];
    dist += (coordinate[i] - xTime[i]) * (coordinate[i] - xTime[i]);
  }

  if (dist < radius * radius) {
    return 1;
  }
  return 0;
}

static int
quadrant_on_sphere_boundary (p4est_t * p4est, p4est_locidx_t which_tree,
                             p4est_quadrant_t * q, double time, double radius,
                             const double x0[P4EST_DIM],
                             const double velocity[P4EST_DIM])
{
  int                 i;
  int                 inside = -1;
  p4est_connectivity_t *conn = p4est->connectivity;

  p4est_qcoord_t      h = P4EST_QUADRANT_LEN (q->level);

  for (i = 0; i < P4EST_CHILDREN; i++) {
    p4est_qcoord_t      tree_coord[P4EST_DIM];
    double              coordinate[3];
    int                 d, this_inside;

    tree_coord[0] = q->x;
    tree_coord[1] = q->y;
#ifdef P4_TO_P8
    tree_coord[2] = q->z;
#endif

    for (d = 0; d < P4EST_DIM; d++) {
      if (i & (1 << d)) {
        tree_coord[d] += h;
      }
    }

#ifndef P4_TO_P8
    p4est_qcoord_to_vertex (conn, which_tree,
                            tree_coord[0], tree_coord[1], coordinate);
#else
    p8est_qcoord_to_vertex (conn, which_tree,
                            tree_coord[0], tree_coord[1], tree_coord[2],
                            coordinate);
#endif

    this_inside = in_sphere (coordinate, time, radius, x0, velocity);
    if (i) {
      inside = this_inside;
    }
    else {
      if (this_inside != inside) {
        return 1;
      }
    }
  }
  return 0;
}


static int
refine_sphere_boundary (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant)
{
  fusion_sphere_t    *sphere = (fusion_sphere_t *) p4est->user_pointer;

  if (quadrant->level >= sphere->max_level) {
    return 0;
  }
  return quadrant_on_sphere_boundary (p4est, which_tree, quadrant,
                                      sphere->time, sphere->radius,
                                      sphere->x0, sphere->velocity);
}


static int
coarsen_sphere_int_ext (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant)
{
  fusion_sphere_t    *sphere = (fusion_sphere_t *) p4est->user_pointer;

  if (quadrant->level =<  1) {
    return 0;
  }
  return quadrant_on_sphere_boundary (p4est, which_tree, quadrant,
                                      sphere->time, sphere->radius,
                                      sphere->x0, sphere->velocity);
}


static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  /* for now, just call refine_sphere_boundary:
   * TODO: allow for other types of refinement */

  return refine_on_sphere_boundary (p4est, which_tree, quadrant);
}

typedef struct
{
  int                 counter;
  int                *refine_flags;
}
refine_loop_t;

/* create an array of ints with values meaning:
 *   0 : keep
 *   1 : refine
 *   2 : coarsen
 */
enum
{
  FUSION_KEEP = 0,
  FUSION_REFINE,
  FUSION_COARSEN
};

static int
refine_in_loop (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quad)
{
  int                 flag;
  refine_loop_t      *loop_ctx = (refine_loop_t *) (p4est->user_pointer);

  flag = (loop_ctx->refine_flags[loop_ctx->counter++]);
  if (flag == FUSION_REFINE) {
    return 1;
  }
  return 0;
}

static void
mark_leaves (p4est_t * p4est, int *refine_flags, test_fusion_t * ctx)
{
  p4est_locidx_t      i;
  p4est_locidx_t      num_local = p4est->local_num_quadrants;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_topidx_t      num_trees = conn->num_trees, t;

  /* TODO: replace with a meaningful (or more than one meaningful, controlled
   * by ctx) refinment pattern */

  p4est->user_pointer = (void *) &(ctx->sphere);
  for (t = 0, i = 0; t < num_trees; t++) {
    p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
    sc_array_t         *quadrants = &tree.quadrants;
    p4est_locidx_t      num_quadrants =
      (p4est_locidx_t) quadrants->elem_count;
    p4est_locidx_t      j;

    for (j = 0; j < num_quadrants; j++, i++) {
      p4est_quadrant_t   *q =
        p4est_quadrant_array_index (quadrants, (size_t) j);
	  int				  lv = q->level;
      int                 on_boundary;
	  int 				  int_ext;
	  int				  c_count;
		if (ctx->context == 0){
	      on_boundary = refine_sphere_boundary (p4est, t, quadrants);
          refine_flags[i] = (on_boundary) ? FUSION_REFINE : FUSION_KEEP;
		}
	    else{
			/* need to make it s.t. it does coarsening check on 
			 * octants on same level. Currently, it is done on
			 * for all quadrants, which needs to be fixed */
	      int_ext = coarsen_sphere_int_ext (p4est, t, quadrants);
		  c_count += int_ext;
		}	
	}
		if (ctx->context != 0){
			if (c_count != 0){
	            refine_flags[i] = FUSION_KEEP;
			}
			else{
	            refine_flags[i] = (on_boundary) ? FUSION_COARSEN : FUSION_KEEP;
			}
		}
  }
}

enum
{
  FUSION_FULL_LOOP,
  FUSION_TIME_COARSEN,
  FUSION_TIME_REFINE,
  FUSION_TIME_BALANCE,
  FUSION_TIME_PARTITION,
  FUSION_TIME_GHOST,
  FUSION_NUM_STATS
};

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t      *ghost;
  int                 i;
  p4est_lnodes_t     *lnodes;   /* runtime option: time lnodes construction */
  int                 first_argc;
  int                 type;
  int                 num_tests = 3;
  sc_statinfo_t       stats[FUSION_NUM_STATS];
  sc_options_t       *opt;
  int                 log_priority = SC_LP_ESSENTIAL;
  fusion_sphere_t     sphere;
  test_fusion_t       ctx;

  /* initialize default values for sphere:
   * TODO: make configurable */
  sphere.time = 0.;
  sphere.radius = 0.25;
  for (i = 0; i < P4EST_DIM; i++) {
    sphere.x0[i] = 0.5;
  }
  sphere.velocity[0] = 0.1;
  sphere.velocity[1] = 0.02;
#ifdef P4_TO_P8
  sphere.velocity[2] = -0.05;
#endif
  sphere.max_level = max_level;

  /* initialize MPI and libsc, p4est packages */
  mpiret = sc_MPI_Init (&argc, &argv);
  mpicomm = sc_MPI_COMM_WORLD;
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'k', "num-tests", &num_tests, num_tests,
                      "The number of instances of timiing the fusion loop");
  sc_options_add_int (opt, 'q', "quiet", &log_priority, log_priority,
                      "Degree of quietude of output");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return 1;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

  sc_set_log_defaults (NULL, NULL, log_priority);
  p4est_init (NULL, log_priority);

  /* set values in ctx here */
  /* random 1 for now */
  ctx.sphere = &sphere;
  ctx.context = 0; // 0 means refine, rest coarsen

#ifndef P4_TO_P8
  conn = p4est_connectivity_new_moebius ();
#else
  conn = p8est_connectivity_new_rotcubes ();
#endif

  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  p4est->user_pointer = (void *) &sphere;
  p4est_refine (p4est, 1, refine_fn, NULL);

  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  p4est_partition (p4est, 0, NULL);

  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  sc_stats_init (&stats[FUSION_FULL_LOOP], "Full loop");
  sc_stats_init (&stats[FUSION_TIME_REFINE], "Local refinement");
  sc_stats_init (&stats[FUSION_TIME_BALANCE], "Balance");
  sc_stats_init (&stats[FUSION_TIME_PARTITION], "Partition");
  sc_stats_init (&stats[FUSION_TIME_GHOST], "Ghost");

  for (i = 0; i <= num_tests; i++) {
    p4est_t            *forest_copy;
    int                *refine_flags, *rflags_copy;
    p4est_ghost_t      *gl_copy;
    refine_loop_t       loop_ctx;
    sc_flopinfo_t       fi_full, snapshot_full;
    sc_flopinfo_t       fi_refine, snapshot_refine;
    sc_flopinfo_t       fi_balance, snapshot_balance;
    sc_flopinfo_t       fi_partition, snapshot_partition;
    sc_flopinfo_t       fi_ghost, snapshot_ghost;
    sc_flopinfo_t       fi_coarsen, snapshot_coarsen;

    if (!i) {
      P4EST_GLOBAL_PRODUCTION ("Timing loop 0 (discarded)\n");
    }
    else {
      P4EST_GLOBAL_PRODUCTIONF ("Timing loop %d\n", i);
    }
    sc_log_indent_push_count (p4est_package_id, 2);

    forest_copy = p4est_copy (p4est, 0 /* do not copy data */ );

    /* predefine which leaves we want to refine and coarsen */

    refine_flags = P4EST_ALLOC (int, p4est->local_num_quadrants);

    sphere.time = 1.;
	ctx.context = 0;
    mark_leaves (forest_copy, refine_flags, &ctx);

    /* Once the leaves have been marked, we have sufficient information to
     * complete a refinement cycle: start the timing at this point */
    sc_flops_snap (&fi_full, &snapshot_full);

    /* start the timing of one instance of the timing cycle */
    /* see sc_flops_snap() / sc_flops_shot() in timings2.c */

    /* non-recursive refinement loop: the callback simply checks the flags
     * that we have defined for which leaves we want to refine */

    /* TODO: conduct coarsening before refinement.  This requires making a
     * copy of refine_flags: creating one that is valid after coarsening */

    sc_flops_snap (&fi_coarsen, &snapshot_coarsen);

    /* coarsen 
	   copy the flags. Then coarsen the interior and exterior of the circle */
    rflags_copy = P4EST_STRDUP (refine_flags);

	ctx.context = 1; // coarsen
	mark_leaves (forest_copy, rflags_copy, &ctx);

    sc_flops_shot (&fi_coarsen, &snapshot_coarsen);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_COARSEN],
                           snapshot_coarsen.iwtime);
    }

    sc_flops_snap (&fi_refine, &snapshot_refine);
    loop_ctx.counter = 0;
    loop_ctx.refine_flags = refine_flags;

    forest_copy->user_pointer = (void *) &loop_ctx;
    p4est_refine (forest_copy, 0 /* non-recursive */ , refine_in_loop, NULL);
    sc_flops_shot (&fi_refine, &snapshot_refine);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_REFINE],
                           snapshot_refine.iwtime);
    }

    sc_flops_shot (&fi_balance, &snapshot_balance);
    p4est_balance (forest_copy, P4EST_CONNECT_FULL, NULL);
    sc_flops_shot (&fi_balance, &snapshot_balance);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_BALANCE],
                           snapshot_balance.iwtime);
    }

    sc_flops_shot (&fi_partition, &snapshot_partition);
    p4est_partition (forest_copy, 0, NULL);
    sc_flops_shot (&fi_partition, &snapshot_partition);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_PARTITION],
                           snapshot_partition.iwtime);
    }

    sc_flops_shot (&fi_ghost, &snapshot_ghost);
    gl_copy = p4est_ghost_new (forest_copy, P4EST_CONNECT_FULL);
    sc_flops_shot (&fi_ghost, &snapshot_ghost);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_TIME_GHOST], snapshot_ghost.iwtime);
    }

    /* end  the timing of one instance of the timing cycle */
    sc_flops_shot (&fi_full, &snapshot_full);
    if (i) {
      sc_stats_accumulate (&stats[FUSION_FULL_LOOP], snapshot_full.iwtime);
    }

    /* clean up */
    P4EST_FREE (refine_flags);
    p4est_ghost_destroy (gl_copy);
    p4est_destroy (forest_copy);

    sc_log_indent_pop_count (p4est_package_id, 2);
  }

  /* accumulate and print statistics */
  sc_stats_compute (mpicomm, FUSION_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL,
                  FUSION_NUM_STATS, stats, 1, 1);

  /* clean up */
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  sc_options_destroy (opt);

  /* exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
