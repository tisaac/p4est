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
#include <p4est_extended.h>
#include <p4est_balance_obj.h>
#include <p4est_vtk.h>
#include <p4est_fused.h>
#include <p4est_algorithms.h>
#else
#include <p8est_bits.h>
#include <p8est.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_extended.h>
#include <p8est_balance_obj.h>
#include <p8est_vtk.h>
#include <p8est_fused.h>
#include <p8est_algorithms.h>
#endif
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <sc_notify.h>

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
  int                 min_level;
}
fusion_sphere_t;

typedef struct
{
  int                 counter;
  int                *refine_flags;
}
refine_loop_t;

typedef struct
{
  p4est_quadrant_t    last_processed;
  int                 counter_in;
  int                 counter_out;
  int                *refine_flags;
  int                *rflags_copy;
}
coarsen_loop_t;

typedef struct
{
  fusion_sphere_t    *sphere;
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
    if (!i) {
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
  fusion_ctx_t       *ctx = p4est_get_fusion_ctx (p4est);
  fusion_sphere_t    *sphere = ctx->sphere;

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
  fusion_ctx_t       *ctx = p4est_get_fusion_ctx (p4est);
  fusion_sphere_t    *sphere = ctx->sphere;

  if (quadrant->level <= SC_MAX (1, sphere->min_level)) {
    return 0;
  }
  return !quadrant_on_sphere_boundary (p4est, which_tree, quadrant,
                                       sphere->time, sphere->radius,
                                       sphere->x0, sphere->velocity);
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  /* for now, just call refine_sphere_boundary:
   * TODO: allow for other types of refinement */

  return refine_sphere_boundary (p4est, which_tree, quadrant);
}

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
  fusion_ctx_t       *ctx = p4est_get_fusion_ctx (p4est);
  refine_loop_t      *loop_ctx = ctx->refine_loop;

  flag = (loop_ctx->refine_flags[loop_ctx->counter++]);
  if (flag == FUSION_REFINE) {
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
      if (flag != FUSION_COARSEN) {
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
    loop_ctx->rflags_copy[loop_ctx->counter_out++] = FUSION_KEEP;
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
        (flag == FUSION_REFINE) ? FUSION_REFINE : FUSION_KEEP;
    }
  }

  return 0;
}

static void
mark_leaves (p4est_t * p4est, int *refine_flags, fusion_ctx_t * ctx)
{
  p4est_locidx_t      i;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_topidx_t      num_trees = conn->num_trees, t;

  /* TODO: replace with a meaningful (or more than one meaningful, controlled
   * by ctx) refinment pattern */

  for (t = 0, i = 0; t < num_trees; t++) {
    p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
    sc_array_t         *quadrants = &(tree->quadrants);
    p4est_locidx_t      num_quadrants =
      (p4est_locidx_t) quadrants->elem_count;
    p4est_locidx_t      j;

    for (j = 0; j < num_quadrants; j++, i++) {
      p4est_quadrant_t   *q =
        p4est_quadrant_array_index (quadrants, (size_t) j);

      refine_flags[i] = FUSION_KEEP;
      if (refine_sphere_boundary (p4est, t, q)) {
        refine_flags[i] = FUSION_REFINE;
      }
      else if (coarsen_sphere_int_ext (p4est, t, q)) {
        refine_flags[i] = FUSION_COARSEN;
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
  FUSION_TIME_OPTIMIZED,
  FUSION_NUM_STATS
};

static double
fusion_compute_h (p4est_t * p4est)
{
  p4est_topidx_t      flt, llt, t;
  p4est_connectivity_t *conn = p4est->connectivity;
  double              h_min = -1., h_min_global;
  int                 numverts = conn->num_vertices, v, d;
  int                 mpierr;
  double              bbox[3][2] = { {0.} };

  flt = p4est->first_local_tree;
  llt = p4est->last_local_tree;

  for (v = 0; v < numverts; v++) {
    double             *vert = &conn->vertices[3 * v];

    if (!v) {
      for (d = 0; d < 3; d++) {
        bbox[d][0] = bbox[d][1] = vert[d];
      }
    }
    else {
      for (d = 0; d < 3; d++) {
        bbox[d][0] = SC_MIN (bbox[d][0], vert[d]);
        bbox[d][1] = SC_MAX (bbox[d][1], vert[d]);
      }
    }
  }
  h_min = sqrt (SC_SQR ((bbox[0][1] - bbox[0][1]))
                + SC_SQR ((bbox[1][1] - bbox[1][0]))
                + SC_SQR ((bbox[2][1] - bbox[2][0])));

  for (t = flt; t <= llt; t++) {
    p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
    sc_array_t         *quadrants = &(tree->quadrants);
    p4est_locidx_t      num_quadrants =
      (p4est_locidx_t) quadrants->elem_count;
    p4est_locidx_t      j;

    for (j = 0; j < num_quadrants; j++) {
      p4est_locidx_t      i;
      p4est_quadrant_t   *q =
        p4est_quadrant_array_index (quadrants, (size_t) j);
      double              first_coord[3] = { 0., 0., 0. };
      p4est_qcoord_t      h = P4EST_QUADRANT_LEN (q->level);

      for (i = 0; i < P4EST_CHILDREN; i++) {
        p4est_qcoord_t      tree_coord[P4EST_DIM];
        double              coordinate[3];
        int                 d;

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
        p4est_qcoord_to_vertex (p4est->connectivity, t,
                                tree_coord[0], tree_coord[1], coordinate);
#else
        p8est_qcoord_to_vertex (p4est->connectivity, t,
                                tree_coord[0], tree_coord[1], tree_coord[2],
                                coordinate);
#endif
        if (!i) {
          for (d = 0; d < P4EST_DIM; d++) {
            first_coord[d] = coordinate[d];
          }
        }
        else {
          double              dist = 0.;

          for (d = 0; d < P4EST_DIM; d++) {
            double              disp = coordinate[d] - first_coord[d];
            dist += disp * disp;
          }
          dist = sqrt (dist);
          if (h_min < 0. || dist < h_min) {
            h_min = dist;
          }
        }
      }
    }
  }
  P4EST_ASSERT (h_min > 0.);
  /* TODO: handle empty processes */
  mpierr =
    sc_MPI_Allreduce (&h_min, &h_min_global, 1, sc_MPI_DOUBLE, sc_MPI_MIN,
                      p4est->mpicomm);
  SC_CHECK_MPI (mpierr);
  P4EST_ASSERT (h_min_global > 0.);

  return h_min_global;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t      *ghost = NULL;
  int                 i;
#if 0
  /* TODO: fusion for finite element nodes as well */
  p4est_lnodes_t     *lnodes;   /* runtime option: time lnodes construction */
#endif
  int                 first_argc;
  int                 num_tests = 3;
  int                 max_level = refine_level;
  int                 min_level = 1;
  int                 skip_reference = 0;
  sc_statinfo_t       stats[FUSION_NUM_STATS];
  sc_options_t       *opt;
  int                 log_priority = SC_LP_ESSENTIAL;
  fusion_ctx_t        ctx;
  fusion_sphere_t     sphere;
  double              velnorm = 0.;
  double              mindist = -1.;
  const char         *out_base_name = NULL;
  const char         *conn_string = NULL;
  const char         *notify_name;
  const char         *balance_name = NULL;
  int                 ntop, nint, nbot;
  int                 balance_sort_root = 0;
  int                 balance_sort_tdeg = 0;
  p4est_inspect_t     inspect;
  p4est_balance_obj_t *bobj;
  sc_notify_t         *notify;
  sc_notify_type_t    notify_type;
  p4est_balance_method_t balance_method;
  sc_statistics_t    *istats;
  p4est_vtk_context_t *vtk;

  /* initialize default values for sphere:
   * TODO: make configurable */
  sphere.time = 0.;
  sphere.radius = 0.25;
  for (i = 0; i < P4EST_DIM; i++) {
    sphere.x0[i] = 0.5;
  }
  sphere.velocity[0] = 0.4;
  sphere.velocity[1] = 0.2;
#ifdef P4_TO_P8
  sphere.velocity[2] = -0.5;
#endif
  sphere.max_level = max_level;
  sphere.min_level = min_level;

  for (i = 0; i < P4EST_DIM; i++) {
    velnorm += sphere.velocity[i] * sphere.velocity[i];
  }
  velnorm = sqrt (velnorm);

  /* initialize MPI and libsc, p4est packages */
  mpiret = sc_MPI_Init (&argc, &argv);
  mpicomm = sc_MPI_COMM_WORLD;
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* grab parameters for notify_nary from command line */
  ntop = sc_notify_nary_ntop_default;
  nint = sc_notify_nary_nint_default;
  nbot = sc_notify_nary_nbot_default;

  /* process command line arguments */
  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'k', "num-tests", &num_tests, num_tests,
                      "The number of instances of timiing the fusion loop");
  sc_options_add_int (opt, 'q', "quiet", &log_priority, log_priority,
                      "Degree of quietude of output");
  sc_options_add_string (opt, 'o', "output", &out_base_name, NULL,
                         "Base name of visualization output");
  sc_options_add_string (opt, 'c', "conn", &conn_string, NULL,
                         "Name of connectivity");
  sc_options_add_int (opt, 'x', "max-level", &sphere.max_level,
                      sphere.max_level, "Maximum refinement level");
  sc_options_add_int (opt, 'i', "min-level", &sphere.min_level,
                      sphere.min_level, "Minimum refinement level");
  sc_options_add_int (opt, '\0', "skip-reference", &skip_reference,
                      skip_reference,
                      "skip timing reference fusion algorithm");
  sc_options_add_int (opt, '\0', "ntop", &ntop, ntop,
                      "wayness of top of sc_notify nary tree");
  sc_options_add_int (opt, '\0', "nbot", &nbot, nbot,
                      "wayness of bot of sc_notify nary tree");
  sc_options_add_int (opt, '\0', "nint", &nint, nint,
                      "wayness of int of sc_notify nary tree");
  sc_options_add_string (opt, 'n', "notify-type", &notify_name, NULL,
                         "Notify algorithm (see sc_notify.h) for type strings");
  sc_options_add_string (opt, '\0', "balance-method", &balance_name,
                         balance_name,
                         "2:1 balance algorithm (see p4est_balance_obj.h for type string)");
  sc_options_add_int (opt, '\0', "balance-sort-root", &balance_sort_root,
                      balance_sort_root,
                      "use fully synchronous sorting in 2:1 balance");
  sc_options_add_int (opt, '\0', "balance-sort-tdeg", &balance_sort_tdeg,
                      balance_sort_tdeg,
                      "tree-branching degree in balance-sort neigbor comm");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return 1;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

  /* set parameters for notify_nary from command line */
  sc_notify_nary_ntop_default = ntop;
  sc_notify_nary_nint_default = nint;
  sc_notify_nary_nbot_default = nbot;

  bobj = p4est_balance_obj_new (mpicomm);

  memset (&inspect, 0, sizeof (p4est_inspect_t));
  notify = NULL;
  if (notify_name) {
    for (i = 0; i < SC_NOTIFY_NUM_TYPES; i++) {
      if (!strcmp (notify_name, sc_notify_type_strings[i])) {
        notify_type = i;
        break;
      }
    }
    notify = sc_notify_new (mpicomm);
    sc_notify_set_type (notify, notify_type);
    p4est_balance_obj_set_notify (bobj, notify);
  }
  if (balance_name) {
    for (i = 0; i < P4EST_BALANCE_NUM_METHODS; i++) {
      if (!strcmp (balance_name, p4est_balance_method_strings[i])) {
        balance_method = i;
        break;
      }
    }
    p4est_balance_obj_set_method (bobj, balance_method);
  }
  if (balance_sort_root) {
    p4est_balance_obj_sort_set_use_root (bobj, 1);
  }
  p4est_balance_obj_sort_set_neigh_tree_degree (bobj, balance_sort_tdeg);

  sc_set_log_defaults (NULL, NULL, log_priority);
  p4est_init (NULL, log_priority);

  /* set values in ctx here */
  /* random 1 for now */
  ctx.sphere = &sphere;
  ctx.viz_name = out_base_name;

  if (conn_string) {
    conn = p4est_connectivity_new_byname (conn_string);
  }
  else {
#ifndef P4_TO_P8
    conn = p4est_connectivity_new_unitsquare ();
#else
    conn = p8est_connectivity_new_unitcube ();
#endif
  }

  p4est = p4est_new_ext (mpicomm, conn, 0 /* min quadrants per proc */ ,
                         sphere.min_level /* min quadrant level */ ,
                         1 /* fill uniform */ ,
                         0 /* data size */ ,
                         NULL, NULL);
  p4est->user_pointer = (void *) &ctx;

  p4est_refine (p4est, 1, refine_fn, NULL);

  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  p4est_partition (p4est, 0, NULL);

  mindist = fusion_compute_h (p4est);

  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  p4est_balance_obj_sort_set_ghost (bobj, ghost);

  sc_stats_init (&stats[FUSION_FULL_LOOP], "Full loop");
  sc_stats_init (&stats[FUSION_TIME_REFINE], "Local refinement");
  sc_stats_init (&stats[FUSION_TIME_COARSEN], "Local coarsening");
  sc_stats_init (&stats[FUSION_TIME_BALANCE], "Balance");
  sc_stats_init (&stats[FUSION_TIME_PARTITION], "Partition");
  sc_stats_init (&stats[FUSION_TIME_GHOST], "Ghost");
  sc_stats_init (&stats[FUSION_TIME_OPTIMIZED], "Optimized");

  istats = sc_statistics_new (mpicomm);

  for (i = 0; i <= num_tests; i++) {
    p4est_t            *forest_copy;
    int                *refine_flags, *rflags_copy;
    p4est_ghost_t      *gl_copy;
    refine_loop_t       ref_loop_ctx;
    coarsen_loop_t      crs_loop_ctx;
    sc_flopinfo_t       fi_full, snapshot_full;
    sc_flopinfo_t       fi_refine, snapshot_refine;
    sc_flopinfo_t       fi_balance, snapshot_balance;
    sc_flopinfo_t       fi_partition, snapshot_partition;
    sc_flopinfo_t       fi_ghost, snapshot_ghost;
    sc_flopinfo_t       fi_coarsen, snapshot_coarsen;
    sc_flopinfo_t       fi_opt, snapshot_opt;

    if (!i) {
      P4EST_GLOBAL_PRODUCTION ("Timing loop 0 (discarded)\n");
    }
    else {
      P4EST_GLOBAL_PRODUCTIONF ("Timing loop %d\n", i);
      inspect.stats = istats;
      p4est_balance_obj_set_stats (bobj, istats);
    }
    sc_log_indent_push_count (p4est_package_id, 2);

    sc_flops_start (&fi_full);
    sc_flops_start (&fi_coarsen);
    sc_flops_start (&fi_refine);
    sc_flops_start (&fi_balance);
    sc_flops_start (&fi_partition);
    sc_flops_start (&fi_ghost);
    sc_flops_start (&fi_opt);

    forest_copy = p4est_copy (p4est, 0 /* do not copy data */ );
    forest_copy->user_pointer = (void *) &ctx;

    /* predefine which leaves we want to refine and coarsen */

    refine_flags = P4EST_ALLOC (int, p4est->local_num_quadrants);

    ctx.sphere->time = 0.5 * mindist / velnorm;
    P4EST_GLOBAL_STATISTICSF ("Velocity: %g, Simulation time: %g\n", velnorm,
                              ctx.sphere->time);
    mark_leaves (forest_copy, refine_flags, &ctx);

    if (!i && ctx.viz_name) {
      char                buffer[BUFSIZ] = { '\0' };

      snprintf (buffer, BUFSIZ, "%s_pre", ctx.viz_name);
      vtk = p4est_vtk_context_new (forest_copy, buffer);
      p4est_vtk_context_set_continuous (vtk, 1);
      p4est_vtk_context_set_scale (vtk, 1.);
      vtk = p4est_vtk_write_header (vtk);
      vtk = p4est_vtk_write_cell_dataf (vtk,
                                  0, /* write tree */
                                  1, /* write level */
                                  2, /* write rank */
                                  0, /* wrap rank */
                                  0, /* no scalars */
                                  0, /* no vectors */
                                  vtk);
      p4est_vtk_write_footer (vtk);
    }

    if (!skip_reference) {
      /* Once the leaves have been marked, we have sufficient information to
       * complete a refinement cycle: start the timing at this point */
      sc_flops_snap (&fi_full, &snapshot_full);

      /* start the timing of one instance of the timing cycle */
      /* see sc_flops_snap() / sc_flops_shot() in timings2.c */

      /* non-recursive refinement loop: the callback simply checks the flags
       * that we have defined for which leaves we want to refine */

      /* TODO: conduct coarsening before refinement.  This requires making a
       * copy of refine_flags: creating one that is valid after coarsening */

      /* make space for copy of flags */

      rflags_copy = P4EST_ALLOC (int, p4est->local_num_quadrants);

      sc_flops_snap (&fi_coarsen, &snapshot_coarsen);

      /* coarsen copy the flags. Then coarsen the interior and exterior of the
       * circle */

      crs_loop_ctx.counter_in = 0;
      crs_loop_ctx.counter_out = 0;
      crs_loop_ctx.refine_flags = refine_flags;
      crs_loop_ctx.rflags_copy = rflags_copy;
      ctx.coarsen_loop = &crs_loop_ctx;
      {
        p4est_locidx_t      n_in = forest_copy->local_num_quadrants;

        p4est_coarsen_ext (forest_copy, 0 /* non-recursive */ ,
                           1 /* callback on ophans */ , coarsen_in_loop, NULL,
                           NULL);
        P4EST_ASSERT (crs_loop_ctx.counter_in == n_in);
        P4EST_ASSERT (crs_loop_ctx.counter_out ==
                      forest_copy->local_num_quadrants);
      }

      if (!i && ctx.viz_name) {
        char                buffer[BUFSIZ] = { '\0' };

        snprintf (buffer, BUFSIZ, "%s_first_coarsen", ctx.viz_name);
        vtk = p4est_vtk_context_new (forest_copy, buffer);
        p4est_vtk_context_set_continuous (vtk, 0);
        vtk = p4est_vtk_write_header (vtk);
        vtk = p4est_vtk_write_cell_dataf (vtk,
                                    0, /* write tree */
                                    1, /* write level */
                                    2, /* write rank */
                                    0, /* wrap rank */
                                    0, /* no scalars */
                                    0, /* no vectors */
                                    vtk);
        p4est_vtk_write_footer (vtk);
      }

      sc_flops_shot (&fi_coarsen, &snapshot_coarsen);
      if (i) {
        sc_stats_accumulate (&stats[FUSION_TIME_COARSEN],
                             snapshot_coarsen.iwtime);
      }

      sc_flops_snap (&fi_refine, &snapshot_refine);
      ref_loop_ctx.counter = 0;
      ref_loop_ctx.refine_flags = rflags_copy;
      ctx.refine_loop = &ref_loop_ctx;
      p4est_refine (forest_copy, 0 /* non-recursive */ , refine_in_loop,
                    NULL);
      sc_flops_shot (&fi_refine, &snapshot_refine);
      if (i) {
        sc_stats_accumulate (&stats[FUSION_TIME_REFINE],
                             snapshot_refine.iwtime);
      }

      if (!i && ctx.viz_name) {
        char                buffer[BUFSIZ] = { '\0' };

        snprintf (buffer, BUFSIZ, "%s_second_refine", ctx.viz_name);
        vtk = p4est_vtk_context_new (forest_copy, buffer);
        p4est_vtk_context_set_continuous (vtk, 0);
        vtk = p4est_vtk_write_header (vtk);
        vtk = p4est_vtk_write_cell_dataf (vtk,
                                    0, /* write tree */
                                    1, /* write level */
                                    2, /* write rank */
                                    0, /* wrap rank */
                                    0, /* no scalars */
                                    0, /* no vectors */
                                    vtk);
        p4est_vtk_write_footer (vtk);
      }

      P4EST_FREE (rflags_copy);

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
        sc_stats_accumulate (&stats[FUSION_TIME_GHOST],
                             snapshot_ghost.iwtime);
      }

      /* end  the timing of one instance of the timing cycle */
      sc_flops_shot (&fi_full, &snapshot_full);
      if (i) {
        sc_stats_accumulate (&stats[FUSION_FULL_LOOP], snapshot_full.iwtime);
      }

      if (!i && ctx.viz_name) {
        char                buffer[BUFSIZ] = { '\0' };

        snprintf (buffer, BUFSIZ, "%s_post", ctx.viz_name);
        vtk = p4est_vtk_context_new (forest_copy, buffer);
        p4est_vtk_context_set_continuous (vtk, 1);
        p4est_vtk_context_set_scale (vtk, 1.);
        vtk = p4est_vtk_write_header (vtk);
        vtk = p4est_vtk_write_cell_dataf (vtk,
                                          0, /* write tree */
                                          1, /* write level */
                                          2, /* write rank */
                                          0, /* wrap rank */
                                          0, /* no scalars */
                                          0, /* no vectors */
                                          vtk);
        p4est_vtk_write_footer (vtk);
      }
    }

    /* Check instrumented version against the reference version */
    {
      p4est_locidx_t      li;
      int8_t             *flags8;
      p4est_t            *forest_copy_ref = NULL;
      p4est_ghost_t      *ghost_copy_ref = NULL;
      p4est_t            *forest_copy_opt = NULL;
      p4est_ghost_t      *ghost_copy_opt = NULL;

      flags8 = P4EST_ALLOC (int8_t, p4est->local_num_quadrants);
      /* (Copy flags to int8_t) */
      for (li = 0; li < p4est->local_num_quadrants; li++) {
        flags8[li] =
          (refine_flags[li] ==
           FUSION_KEEP) ? P4EST_FUSED_KEEP : (refine_flags[li] ==
                                              FUSION_REFINE) ?
          P4EST_FUSED_REFINE : P4EST_FUSED_COARSEN;
      }

      if (!skip_reference) {
        p4est_adapt_fused_reference (p4est, flags8, 0 /* no data */ ,
                                     bobj,
                                     1 /* yes repartition */ ,
                                     0 /* no repartition for coarsening */ ,
                                     1 /* ghost layer width */ ,
                                     P4EST_CONNECT_FULL,
                                     NULL, NULL, NULL, &forest_copy_ref,
                                     &ghost_copy_ref);

        SC_CHECK_ABORT (p4est_is_equal (forest_copy, forest_copy_ref, 0),
                        "Instrumented adaptivity cycle forest different from reference\n");
        SC_CHECK_ABORT (p4est_ghost_is_equal (gl_copy, ghost_copy_ref),
                        "Instrumented adaptivity cycle ghost layer different from reference\n");
      }

      p4est->inspect = &inspect;

      mpiret = sc_MPI_Barrier (mpicomm);
      SC_CHECK_MPI (mpiret);
      sc_flops_snap (&fi_opt, &snapshot_opt);
      p4est_adapt_fused (p4est, flags8, 0 /* no data */ ,
                         bobj, 1 /* yes repartition */ ,
                         0 /* no repartition for coarsening */ ,
                         1 /* ghost layer width */ ,
                         P4EST_CONNECT_FULL,
                         NULL, NULL, NULL, &forest_copy_opt, &ghost_copy_opt);
      sc_flops_shot (&fi_opt, &snapshot_opt);
      if (i) {
        sc_stats_accumulate (&stats[FUSION_TIME_OPTIMIZED],
                             snapshot_opt.iwtime);
      }

      p4est->inspect = NULL;

      if (!skip_reference) {
        SC_CHECK_ABORT (p4est_is_equal (forest_copy_ref, forest_copy_opt, 0),
                        "Optimized adaptivity cycle forest different from reference\n");
        SC_CHECK_ABORT (p4est_ghost_is_equal (ghost_copy_ref, ghost_copy_opt),
                        "Optimized adaptivity cycle ghost layer different from reference\n");
        p4est_ghost_destroy (ghost_copy_ref);
        p4est_destroy (forest_copy_ref);
      }

      p4est_ghost_destroy (ghost_copy_opt);
      p4est_destroy (forest_copy_opt);
      P4EST_FREE (flags8);
    }

    /* clean up */
    P4EST_FREE (refine_flags);
    if (!skip_reference) {
      p4est_ghost_destroy (gl_copy);
    }
    p4est_destroy (forest_copy);

    sc_log_indent_pop_count (p4est_package_id, 2);
  }

  /* accumulate and print statistics */
  sc_stats_compute (mpicomm, FUSION_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL,
                  FUSION_NUM_STATS, stats, 1, 1);

  sc_statistics_compute (istats);
  sc_statistics_print (istats, p4est_package_id, SC_LP_ESSENTIAL, 1, 1);
  sc_statistics_destroy (istats);

  if (notify) {
    sc_notify_destroy (notify);
  }

  /* clean up */
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  p4est_balance_obj_destroy (bobj);

  sc_options_destroy (opt);

  /* exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
