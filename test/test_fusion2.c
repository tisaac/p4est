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
  int dummy;
}
test_fusion_t;

enum
{
  TIMINGS_INIT,
  TIMINGS_NUM_STATS
};

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int refine_in_loop (p4est_t *p4est, p4est_topidx_t which_tree,
                           p4est_quadrant_t *quad);

static void mark_leaves (p4est_t *, int *, test_fusion_t *ctx);

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t      *ghost;
  int                 num_cycles = 2;
  int                 i;
  p4est_lnodes_t     *lnodes;
  int                 first_argc;
  int                 type;
  int                 num_tests = 3;
  sc_statinfo_t		  stats[num_tests+1];
  sc_flopinfo_t		  fi, snapshot;
  sc_options_t       *opt;
  test_fusion_t       ctx;

  /* initialize MPI and libsc, p4est packages */
  mpiret = sc_MPI_Init (&argc, &argv);
  mpicomm = sc_MPI_COMM_WORLD;
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm,1, 1, NULL, SC_LP_DEFAULT);

#ifndef P4EST_ENABLE_DEBUG 
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS); 
#endif 

  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  opt = sc_options_new (argv[0]);

  sc_options_add_int (opt, 'k', "num-tests", &num_tests, num_tests,
                      "The number of instances of timiing the fusion loop");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return 1;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

  /* set values in ctx here */

#ifndef P4_TO_P8
  conn = p4est_connectivity_new_moebius ();
#else
  conn = p8est_connectivity_new_rotcubes ();
#endif

  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  p4est_refine (p4est, 1, refine_fn, NULL);

  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  p4est_partition (p4est, 0, NULL);

  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  for (i = 0; i <= num_tests; i++) {

    sc_flops_snap (&fi, &snapshot);
    p4est_t *forest_copy;
    int     *refine_flags;
    p4est_ghost_t *gl_copy;

    forest_copy = p4est_copy (p4est, 0 /* do not copy data */);

    /* predefine which leaves we want to refine and coarsen */

    /* create an array of ints with values meaning:
     *   0 : keep
     *   1 : refine
     *   2 : coarsen
     */

    refine_flags = P4EST_ALLOC (int, p4est->local_num_quadrants);

    mark_leaves (forest_copy, refine_flags, &ctx);

    /* start the timing of one instance of the timing cycle */
    /* see sc_flops_snap() / sc_flops_shot() in timings2.c */

    /* non-recursive refinement loop: the callback simply checks the flags
     * that we have defined for which leaves we want to refine */


    p4est_refine (forest_copy, 0 /* non-recursive */, refine_in_loop, NULL);

    p4est_balance (forest_copy, P4EST_CONNECT_FULL, NULL);

    p4est_partition (forest_copy, 0, NULL);

    gl_copy = p4est_ghost_new (forest_copy, P4EST_CONNECT_FULL);

    /* end  the timing of one instance of the timing cycle */

    /* clean up */
    P4EST_FREE (refine_flags);
    p4est_ghost_destroy (gl_copy);
    p4est_destroy (forest_copy);

    sc_flops_shot (&fi, &snapshot);
    sc_stats_set1 (&stats[i], snapshot.iwtime, "Refine Loops");
  }

  /* accumulate and print statistics */

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
