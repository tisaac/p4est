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
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_io.h>
#include <p8est_search.h>
#include <p8est_balance_obj.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_io.h>
#include <p4est_search.h>
#include <p4est_balance_obj.h>
#endif /* !P4_TO_P8 */
#include <sc_io.h>
#include <sc_notify.h>
#include <sc_ranges.h>
#include <sc_search.h>
#ifdef P4EST_HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef P4EST_ENABLE_MPIIO
#define P4EST_MPIIO_WRITE
#endif

#ifdef P4EST_HAVE_UNISTD_H
#include <unistd.h>
#endif

#define p4est_num_ranges (25)

#ifndef P4_TO_P8

static int          p4est_uninitialized_key;
void               *P4EST_DATA_UNINITIALIZED = &p4est_uninitialized_key;

#endif /* P4_TO_P8 */

static const size_t number_toread_quadrants = 32;

void
p4est_qcoord_to_vertex (p4est_connectivity_t * connectivity,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y,
#ifdef P4_TO_P8
                        p4est_qcoord_t z,
#endif
                        double vxyz[3])
{
  const double       *vertices = connectivity->vertices;
#ifdef P4EST_ENABLE_DEBUG
  const p4est_topidx_t num_vertices = connectivity->num_vertices;
#endif
  const p4est_topidx_t *vindices;
  int                 xi, yi;
  double              wx[2], wy[2];
#ifdef P4_TO_P8
  int                 zi;
  double              wz[2];
#endif
  double              xfactor, yfactor;
  p4est_topidx_t      vindex;

  P4EST_ASSERT (num_vertices > 0);
  P4EST_ASSERT (vertices != NULL);
  P4EST_ASSERT (treeid >= 0 && treeid < connectivity->num_trees);

  P4EST_ASSERT (connectivity->tree_to_vertex != NULL);
  vindices = connectivity->tree_to_vertex + P4EST_CHILDREN * treeid;

  vxyz[0] = vxyz[1] = vxyz[2] = 0.;

  P4EST_ASSERT (x >= 0 && x <= P4EST_ROOT_LEN);
  wx[1] = (double) x / (double) P4EST_ROOT_LEN;
  wx[0] = 1. - wx[1];

  P4EST_ASSERT (y >= 0 && y <= P4EST_ROOT_LEN);
  wy[1] = (double) y / (double) P4EST_ROOT_LEN;
  wy[0] = 1. - wy[1];

#ifdef P4_TO_P8
  P4EST_ASSERT (z >= 0 && z <= P4EST_ROOT_LEN);
  wz[1] = (double) z / (double) P4EST_ROOT_LEN;
  wz[0] = 1. - wz[1];

  for (zi = 0; zi < 2; ++zi) {
#endif
    for (yi = 0; yi < 2; ++yi) {
#ifdef P4_TO_P8
      yfactor = wz[zi] * wy[yi];
#else
      yfactor = wy[yi];
#endif
      for (xi = 0; xi < 2; ++xi) {
        xfactor = yfactor * wx[xi];

        vindex = *vindices++;
        P4EST_ASSERT (vindex >= 0 && vindex < num_vertices);

        vxyz[0] += xfactor * vertices[3 * vindex + 0];
        vxyz[1] += xfactor * vertices[3 * vindex + 1];
        vxyz[2] += xfactor * vertices[3 * vindex + 2];
      }
    }
#ifdef P4_TO_P8
  }
#endif
}

size_t
p4est_memory_used (p4est_t * p4est)
{
  int                 mpisize;
  size_t              size;
  p4est_topidx_t      nt;
  p4est_tree_t       *tree;

  /* do not assert p4est_is_valid since it is collective */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->connectivity != NULL);
  P4EST_ASSERT (p4est->trees != NULL);

  mpisize = p4est->mpisize;
  size = sizeof (p4est_t) +
    (mpisize + 1) * (sizeof (p4est_gloidx_t) + sizeof (p4est_quadrant_t));

  size += sc_array_memory_used (p4est->trees, 1);
  for (nt = 0; nt < p4est->connectivity->num_trees; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    size += sc_array_memory_used (&tree->quadrants, 0);
  }

  if (p4est->data_size > 0) {
    P4EST_ASSERT (p4est->user_data_pool != NULL);
    size += sc_mempool_memory_used (p4est->user_data_pool);
  }
  P4EST_ASSERT (p4est->quadrant_pool != NULL);
  size += sc_mempool_memory_used (p4est->quadrant_pool);

  return size;
}

long
p4est_revision (p4est_t * p4est)
{
  /* do not assert p4est_is_valid since it is collective */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->revision >= 0);

  return p4est->revision;
}

p4est_t            *
p4est_new (sc_MPI_Comm mpicomm, p4est_connectivity_t * connectivity,
           size_t data_size, p4est_init_t init_fn, void *user_pointer)
{
  return p4est_new_ext (mpicomm, connectivity, 0, 0, 1,
                        data_size, init_fn, user_pointer);
}

p4est_t            *
p4est_new_ext (sc_MPI_Comm mpicomm, p4est_connectivity_t * connectivity,
               p4est_locidx_t min_quadrants, int min_level, int fill_uniform,
               size_t data_size, p4est_init_t init_fn, void *user_pointer)
{
  int                 num_procs, rank;
  int                 i, must_remove_last_quadrant;
  int                 level;
  uint64_t            first_morton, last_morton, miu, count;
  p4est_topidx_t      jt, num_trees;
  p4est_gloidx_t      tree_num_quadrants, global_num_quadrants;
  p4est_gloidx_t      first_tree, first_quadrant, first_tree_quadrant;
  p4est_gloidx_t      last_tree, last_quadrant, last_tree_quadrant;
  p4est_gloidx_t      quadrant_index;
  p4est_t            *p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_quadrant_t    a, b, c;
  p4est_quadrant_t   *global_first_position;
  sc_array_t         *tquadrants;

  P4EST_GLOBAL_PRODUCTIONF
    ("Into " P4EST_STRING
     "_new with min quadrants %lld level %d uniform %d\n",
     (long long) min_quadrants, SC_MAX (min_level, 0), fill_uniform);
  p4est_log_indent_push ();

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));
  P4EST_ASSERT (min_level <= P4EST_QMAXLEVEL);

  /* create p4est object and assign some data members */
  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  p4est->data_size = data_size;
  p4est->user_pointer = user_pointer;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  /* set parallel environment */
  p4est_comm_parallel_env_assign (p4est, mpicomm);
  num_procs = p4est->mpisize;
  rank = p4est->mpirank;

  /* allocate memory pools */
  if (p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  else {
    p4est->user_data_pool = NULL;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* determine uniform level of initial tree */
  tree_num_quadrants = 1;
  for (level = 0; level < P4EST_QMAXLEVEL; ++level) {
    if (tree_num_quadrants >=
        (num_procs * (p4est_gloidx_t) min_quadrants + (num_trees - 1))
        / num_trees) {
      break;
    }
    tree_num_quadrants *= P4EST_CHILDREN;
    P4EST_ASSERT (tree_num_quadrants > 0);
  }
  for (; level < min_level; ++level) {
    tree_num_quadrants *= P4EST_CHILDREN;
    P4EST_ASSERT (tree_num_quadrants > 0);
  }
  P4EST_ASSERT (level <= P4EST_QMAXLEVEL
                && tree_num_quadrants <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);

  /* compute global number of quadrants */
  global_num_quadrants = tree_num_quadrants * num_trees;
  P4EST_GLOBAL_PRODUCTIONF ("New " P4EST_STRING
                            " with %lld trees on %d processors\n",
                            (long long) num_trees, num_procs);
  P4EST_GLOBAL_INFOF ("Initial level %d potential global quadrants"
                      " %lld per tree %lld\n",
                      level, (long long) global_num_quadrants,
                      (long long) tree_num_quadrants);

  /* compute index of first tree for this processor */
  first_quadrant =
    p4est_partition_cut_gloidx (global_num_quadrants, rank, num_procs);
  first_tree = first_quadrant / tree_num_quadrants;
  first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
  last_quadrant = p4est_partition_cut_gloidx (global_num_quadrants, rank + 1,
                                              num_procs) - 1;
  P4EST_VERBOSEF
    ("first tree %lld first quadrant %lld global quadrant %lld\n",
     (long long) first_tree, (long long) first_tree_quadrant,
     (long long) first_quadrant);
  P4EST_ASSERT (first_tree_quadrant < tree_num_quadrants);

  /* compute index of last tree for this processor */
  if (first_quadrant <= last_quadrant) {
    last_tree = last_quadrant / tree_num_quadrants;
    last_tree_quadrant = last_quadrant - last_tree * tree_num_quadrants;
    P4EST_VERBOSEF
      ("last tree %lld last quadrant %lld global quadrant %lld\n",
       (long long) last_tree, (long long) last_tree_quadrant,
       (long long) last_quadrant);

    /* check ranges of various integers to be 32bit compatible */
    P4EST_ASSERT (first_tree <= last_tree && last_tree < num_trees);
    P4EST_ASSERT (0 <= first_tree_quadrant && 0 <= last_tree_quadrant);
    P4EST_ASSERT (last_tree_quadrant < tree_num_quadrants);
    if (first_tree == last_tree) {
      P4EST_ASSERT (first_tree_quadrant <= last_tree_quadrant);
    }
  }
  else {
    P4EST_VERBOSE ("Empty processor");
    P4EST_ASSERT (0 <= first_tree && 0 <= first_tree_quadrant);
    first_tree = -1;
    last_tree = -2;
    last_tree_quadrant = -1;
  }

  /* allocate trees and quadrants */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    sc_array_init (&tree->quadrants, sizeof (p4est_quadrant_t));
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
    tree->quadrants_offset = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    for (; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = -1;
    }
    tree->maxlevel = 0;
  }
  p4est->local_num_quadrants = 0;
  p4est->global_num_quadrants = 0;

  /* for every locally non-empty tree fill first and last quadrant */
  P4EST_QUADRANT_INIT (&a);
  P4EST_QUADRANT_INIT (&b);
  P4EST_QUADRANT_INIT (&c);
  for (jt = first_tree; jt <= last_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;

    quad = NULL;
    if (!fill_uniform) {        /* fill with coarsest possible quadrants */
      must_remove_last_quadrant = 0;

      /* set morton id of first quadrant and initialize user data */
      if (jt == first_tree) {
        p4est_quadrant_set_morton (&a, level, first_tree_quadrant);
      }
      else {
        p4est_quadrant_set_morton (&a, level, 0);
      }
#ifdef P4_TO_P8
      P4EST_LDEBUGF ("tree %lld first morton 0x%llx 0x%llx 0x%llx\n",
                     (long long) jt, (long long) a.x,
                     (long long) a.y, (long long) a.z);
#else
      P4EST_LDEBUGF ("tree %lld first morton 0x%llx 0x%llx\n",
                     (long long) jt, (long long) a.x, (long long) a.y);
#endif
      p4est_quadrant_first_descendant (&a, &tree->first_desc,
                                       P4EST_QMAXLEVEL);

      /* set morton id of last quadrant */
      if (tree_num_quadrants == 1 ||
          (jt == first_tree
           && first_tree_quadrant == tree_num_quadrants - 1)) {
        /* There is only a in the tree */
        quad = p4est_quadrant_array_push (tquadrants);
        *quad = a;
        p4est_quadrant_init_data (p4est, jt, quad, init_fn);
        tree->maxlevel = a.level;
        tree->quadrants_per_level[a.level] = 1;
      }
      else {
        if (jt == last_tree) {
          if (last_tree_quadrant == tree_num_quadrants - 1) {
            quadrant_index = last_tree_quadrant;
          }
          else {
            quadrant_index = last_tree_quadrant + 1;
            must_remove_last_quadrant = 1;
          }
          p4est_quadrant_set_morton (&b, level, quadrant_index);
        }
        else {
          p4est_quadrant_set_morton (&b, level, tree_num_quadrants - 1);
        }
#ifdef P4_TO_P8
        P4EST_LDEBUGF ("tree %lld last morton 0x%llx 0x%llx 0x%llx\n",
                       (long long) jt, (long long) b.x,
                       (long long) b.y, (long long) b.z);
#else
        P4EST_LDEBUGF ("tree %lld last morton 0x%llx 0x%llx\n",
                       (long long) jt, (long long) b.x, (long long) b.y);
#endif
        /* fill up tree between a and b with coarse quadrants */
        p4est_complete_region (p4est, &a, 1, &b, !must_remove_last_quadrant,
                               tree, jt, init_fn);
        quad = p4est_quadrant_array_index (tquadrants,
                                           tquadrants->elem_count - 1);
      }
    }
    else {                      /* fill tree with quadrants of given level */
      /* determine range of quadrants in this tree */
      first_morton = (uint64_t)
        (jt == first_tree ? first_tree_quadrant : 0);
      last_morton = (uint64_t)
        (jt == last_tree ? last_tree_quadrant : tree_num_quadrants - 1);
      count = last_morton - first_morton + 1;
      P4EST_ASSERT (count > 0);

      /* populate quadrant array in Morton order */
      sc_array_resize (tquadrants, (size_t) count);
      for (miu = 0; miu < count; ++miu) {
        quad = p4est_quadrant_array_index (tquadrants, (size_t) miu);
        p4est_quadrant_set_morton (quad, level, first_morton + miu);
        p4est_quadrant_init_data (p4est, jt, quad, init_fn);
      }

      /* remember first tree position */
      p4est_quadrant_first_descendant (p4est_quadrant_array_index
                                       (tquadrants, 0), &tree->first_desc,
                                       P4EST_QMAXLEVEL);

      /* set tree counters */
      tree->maxlevel = (int8_t) level;
      tree->quadrants_per_level[level] = (p4est_locidx_t) count;
    }

#if 0
    P4EST_VERBOSEF ("tree %lld quadrants %llu\n", (long long) jt,
                    (unsigned long long) tquadrants->elem_count);
#endif

    tree->quadrants_offset = p4est->local_num_quadrants;
    p4est->local_num_quadrants += tquadrants->elem_count;
    p4est_quadrant_last_descendant (quad, &tree->last_desc, P4EST_QMAXLEVEL);
  }
  if (last_tree >= 0) {
    for (; jt < num_trees; ++jt) {
      tree = p4est_tree_array_index (p4est->trees, jt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  /* compute some member variables */
  p4est->first_local_tree = first_tree;
  p4est->last_local_tree = last_tree;
  p4est->global_first_quadrant = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  if (!fill_uniform && level > 0) {
    /* this performs an allgather to count all quadrants */
    p4est_comm_count_quadrants (p4est);
  }
  else {
    /* for a uniform forest we know all global information a priori */
    for (i = 0; i <= num_procs; ++i) {
      p4est->global_first_quadrant[i] =
        p4est_partition_cut_gloidx (global_num_quadrants, i, num_procs);
    }
    p4est->global_num_quadrants = global_num_quadrants;
  }

  /* fill in global partition information */
  global_first_position = P4EST_ALLOC_ZERO (p4est_quadrant_t, num_procs + 1);
  for (i = 0; i <= num_procs; ++i) {
    first_quadrant =
      p4est_partition_cut_gloidx (global_num_quadrants, i, num_procs);
    first_tree = first_quadrant / tree_num_quadrants;
    first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
    p4est_quadrant_set_morton (&c, level, first_tree_quadrant);
    global_first_position[i].x = c.x;
    global_first_position[i].y = c.y;
#ifdef P4_TO_P8
    global_first_position[i].z = c.z;
#endif
    global_first_position[i].level = P4EST_QMAXLEVEL;
    global_first_position[i].p.which_tree = first_tree;
  }
  p4est->global_first_position = global_first_position;

  /* print more statistics */
  P4EST_VERBOSEF ("total local quadrants %lld\n",
                  (long long) p4est->local_num_quadrants);

  P4EST_ASSERT (p4est->revision == 0);
  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_new with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  return p4est;
}

void
p4est_destroy (p4est_t * p4est)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              qz;
#endif
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;

  for (jt = 0; jt < p4est->connectivity->num_trees; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);

#ifdef P4EST_ENABLE_DEBUG
    for (qz = 0; qz < tree->quadrants.elem_count; ++qz) {
      p4est_quadrant_t   *quad =
        p4est_quadrant_array_index (&tree->quadrants, qz);
      p4est_quadrant_free_data (p4est, quad);
    }
#endif

    sc_array_reset (&tree->quadrants);
  }
  sc_array_destroy (p4est->trees);

  if (p4est->user_data_pool != NULL) {
    sc_mempool_destroy (p4est->user_data_pool);
  }
  sc_mempool_destroy (p4est->quadrant_pool);

  p4est_comm_parallel_env_release (p4est);
  P4EST_FREE (p4est->global_first_quadrant);
  P4EST_FREE (p4est->global_first_position);
  P4EST_FREE (p4est);
}

p4est_t            *
p4est_copy (p4est_t * input, int copy_data)
{
  return p4est_copy_ext (input, copy_data, 0 /* don't duplicate MPI comm */ );
}

p4est_t            *
p4est_copy_ext (p4est_t * input, int copy_data, int duplicate_mpicomm)
{
  const p4est_topidx_t num_trees = input->connectivity->num_trees;
  const p4est_topidx_t first_tree = input->first_local_tree;
  const p4est_topidx_t last_tree = input->last_local_tree;
  size_t              icount;
  size_t              zz;
  p4est_topidx_t      jt;
  p4est_t            *p4est;
  p4est_tree_t       *itree, *ptree;
  p4est_quadrant_t   *iq, *pq;
  sc_array_t         *iquadrants, *pquadrants;

  /* create a shallow copy and zero out dependent fields */
  p4est = P4EST_ALLOC (p4est_t, 1);
  memcpy (p4est, input, sizeof (p4est_t));
  p4est->global_first_quadrant = NULL;
  p4est->global_first_position = NULL;
  p4est->trees = NULL;
  p4est->user_data_pool = NULL;
  p4est->quadrant_pool = NULL;

  /* set parallel environment */
  p4est_comm_parallel_env_assign (p4est, input->mpicomm);
  if (duplicate_mpicomm) {
    p4est_comm_parallel_env_duplicate (p4est);
  }

  /* allocate a user data pool if necessary and a quadrant pool */
  if (copy_data && p4est->data_size > 0) {
    p4est->user_data_pool = sc_mempool_new (p4est->data_size);
  }
  else {
    p4est->data_size = 0;
  }
  p4est->quadrant_pool = sc_mempool_new (sizeof (p4est_quadrant_t));

  /* copy quadrants for each tree */
  p4est->trees = sc_array_new (sizeof (p4est_tree_t));
  sc_array_resize (p4est->trees, num_trees);
  for (jt = 0; jt < num_trees; ++jt) {
    itree = p4est_tree_array_index (input->trees, jt);
    ptree = p4est_tree_array_index (p4est->trees, jt);
    memcpy (ptree, itree, sizeof (p4est_tree_t));
    sc_array_init (&ptree->quadrants, sizeof (p4est_quadrant_t));
  }
  for (jt = first_tree; jt <= last_tree; ++jt) {
    itree = p4est_tree_array_index (input->trees, jt);
    iquadrants = &itree->quadrants;
    icount = iquadrants->elem_count;
    ptree = p4est_tree_array_index (p4est->trees, jt);
    pquadrants = &ptree->quadrants;
    sc_array_resize (pquadrants, icount);
    memcpy (pquadrants->array, iquadrants->array,
            icount * sizeof (p4est_quadrant_t));
    if (p4est->data_size > 0) {
      P4EST_ASSERT (copy_data);
      for (zz = 0; zz < icount; ++zz) {
        iq = p4est_quadrant_array_index (iquadrants, zz);
        pq = p4est_quadrant_array_index (pquadrants, zz);
        pq->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
        memcpy (pq->p.user_data, iq->p.user_data, p4est->data_size);
      }
    }
  }

  /* allocate and copy global quadrant count */
  p4est->global_first_quadrant =
    P4EST_ALLOC (p4est_gloidx_t, p4est->mpisize + 1);
  memcpy (p4est->global_first_quadrant, input->global_first_quadrant,
          (p4est->mpisize + 1) * sizeof (p4est_gloidx_t));

  /* allocate and copy global partition information */
  p4est->global_first_position = P4EST_ALLOC (p4est_quadrant_t,
                                              p4est->mpisize + 1);
  memcpy (p4est->global_first_position, input->global_first_position,
          (p4est->mpisize + 1) * sizeof (p4est_quadrant_t));

  /* the copy starts with a revision count of zero */
  p4est->revision = 0;

  /* check for valid p4est and return */
  P4EST_ASSERT (p4est_is_valid (p4est));

  return p4est;
}

void
p4est_reset_data (p4est_t * p4est, size_t data_size,
                  p4est_init_t init_fn, void *user_pointer)
{
  int                 doresize;
  size_t              zz;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;

  doresize = (p4est->data_size != data_size);

  p4est->data_size = data_size;
  p4est->user_pointer = user_pointer;

  if (doresize) {
    if (p4est->user_data_pool != NULL) {
      sc_mempool_destroy (p4est->user_data_pool);
    }
    if (p4est->data_size > 0) {
      p4est->user_data_pool = sc_mempool_new (p4est->data_size);
    }
    else {
      p4est->user_data_pool = NULL;
    }
  }

  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      q = p4est_quadrant_array_index (tquadrants, zz);
      if (doresize) {
        if (p4est->data_size > 0) {
          q->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
        }
        else {
          q->p.user_data = NULL;
        }
      }
      if (init_fn != NULL) {
        init_fn (p4est, jt, q);
      }
    }
  }
}

void
p4est_refine (p4est_t * p4est, int refine_recursive,
              p4est_refine_t refine_fn, p4est_init_t init_fn)
{
  p4est_refine_ext (p4est, refine_recursive, -1, refine_fn, init_fn, NULL);
}

void
p4est_refine_ext (p4est_t * p4est, int refine_recursive, int allowed_level,
                  p4est_refine_t refine_fn, p4est_init_t init_fn,
                  p4est_replace_t replace_fn)
{
  p4est_gloidx_t      old_gnq;

  if (allowed_level < 0) {
    allowed_level = P4EST_QMAXLEVEL;
  }
  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_refine with %lld total quadrants,"
                            " allowed level %d\n",
                            (long long) p4est->global_num_quadrants,
                            allowed_level);
  p4est_log_indent_push ();
  P4EST_ASSERT (p4est_is_valid (p4est));
  /* remember input quadrant count; it will not decrease */
  old_gnq = p4est->global_num_quadrants;
  p4est_refine_ext_dirty (p4est, refine_recursive, allowed_level, refine_fn,
                          init_fn, replace_fn);
  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);
  P4EST_ASSERT (p4est->global_num_quadrants >= old_gnq);
  if (old_gnq != p4est->global_num_quadrants) {
    ++p4est->revision;
  }

  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_refine with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

void
p4est_refine_ext_dirty (p4est_t * p4est, int refine_recursive,
                        int allowed_level, p4est_refine_t refine_fn,
                        p4est_init_t init_fn, p4est_replace_t replace_fn)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              quadrant_pool_size, data_pool_size;
#endif
  int                 firsttime;
  int                 i, maxlevel;
  p4est_topidx_t      nt;
  size_t              incount, current, restpos, movecount;
  sc_list_t          *list;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q, *qalloc, *qpop;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
#ifdef P4_TO_P8
  p4est_quadrant_t   *c4, *c5, *c6, *c7;
#endif
  sc_array_t         *tquadrants;
  p4est_quadrant_t   *family[8];
  sc_flopinfo_t       snap;
  p4est_quadrant_t    parent, *pp = &parent;

  P4EST_FUNC_SNAP (p4est, &snap);

  if (allowed_level < 0) {
    allowed_level = P4EST_QMAXLEVEL;
  }
  P4EST_ASSERT (0 <= allowed_level && allowed_level <= P4EST_QMAXLEVEL);
  P4EST_ASSERT (refine_fn != NULL);

  /*
     q points to a quadrant that is an array member
     qalloc is a quadrant that has been allocated through quadrant_pool
     qpop is a quadrant that has been allocated through quadrant_pool
     never mix these two types of quadrant pointers

     The quadrant->pad8 field of list quadrants is interpreted as boolean
     and set to true for quadrants that have already been refined.
   */
  list = sc_list_new (NULL);
  p4est->local_num_quadrants = 0;

  /* loop over all local trees */
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    tree->quadrants_offset = p4est->local_num_quadrants;
    tquadrants = &tree->quadrants;
#ifdef P4EST_ENABLE_DEBUG
    quadrant_pool_size = p4est->quadrant_pool->elem_count;
    data_pool_size = 0;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }
#endif

    /* initial log message for this tree */
    P4EST_VERBOSEF ("Into refine tree %lld with %llu\n", (long long) nt,
                    (unsigned long long) tquadrants->elem_count);

    /* reset the quadrant counters */
    maxlevel = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }

    /* run through the array to find first quadrant to be refined */
    q = NULL;
    incount = tquadrants->elem_count;
    for (current = 0; current < incount; ++current) {
      q = p4est_quadrant_array_index (tquadrants, current);
      if (refine_fn (p4est, nt, q) && (int) q->level < allowed_level) {
        break;
      }
      maxlevel = SC_MAX (maxlevel, (int) q->level);
      ++tree->quadrants_per_level[q->level];
    }
    if (current == incount) {
      /* no refinement occurs in this tree */
      p4est->local_num_quadrants += incount;
      continue;
    }
    P4EST_ASSERT (q != NULL);

    /* now we have a quadrant to refine, prepend it to the list */
    qalloc = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
    *qalloc = *q;               /* never prepend array members directly */
    qalloc->pad8 = 0;           /* this quadrant has not been refined yet */
    (void) sc_list_prepend (list, qalloc);      /* only new quadrants */

    P4EST_QUADRANT_INIT (&parent);

    /*
       current points to the next array member to write
       restpos points to the next array member to read
     */
    restpos = current + 1;

    /* run through the list and refine recursively */
    firsttime = 1;
    while (list->elem_count > 0) {
      qpop = p4est_quadrant_list_pop (list);
      if (firsttime ||
          ((refine_recursive || !qpop->pad8) &&
           refine_fn (p4est, nt, qpop) &&
           (int) qpop->level < allowed_level)) {
        firsttime = 0;
        sc_array_resize (tquadrants,
                         tquadrants->elem_count + P4EST_CHILDREN - 1);

        if (replace_fn != NULL) {
          /* do not free qpop's data yet: we will do this when the parent
           * is replaced */
          parent = *qpop;
        }
        else {
          p4est_quadrant_free_data (p4est, qpop);
        }
        c0 = qpop;
        c1 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c2 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c3 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);

#ifdef P4_TO_P8
        c4 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c5 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c6 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
        c7 = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);

        p8est_quadrant_children (qpop, c0, c1, c2, c3, c4, c5, c6, c7);
#else
        p4est_quadrant_children (qpop, c0, c1, c2, c3);
#endif
        p4est_quadrant_init_data (p4est, nt, c0, init_fn);
        p4est_quadrant_init_data (p4est, nt, c1, init_fn);
        p4est_quadrant_init_data (p4est, nt, c2, init_fn);
        p4est_quadrant_init_data (p4est, nt, c3, init_fn);
        c0->pad8 = c1->pad8 = c2->pad8 = c3->pad8 = 1;

#ifdef P4_TO_P8
        p4est_quadrant_init_data (p4est, nt, c4, init_fn);
        p4est_quadrant_init_data (p4est, nt, c5, init_fn);
        p4est_quadrant_init_data (p4est, nt, c6, init_fn);
        p4est_quadrant_init_data (p4est, nt, c7, init_fn);
        c4->pad8 = c5->pad8 = c6->pad8 = c7->pad8 = 1;

        (void) sc_list_prepend (list, c7);
        (void) sc_list_prepend (list, c6);
        (void) sc_list_prepend (list, c5);
        (void) sc_list_prepend (list, c4);
#endif
        (void) sc_list_prepend (list, c3);
        (void) sc_list_prepend (list, c2);
        (void) sc_list_prepend (list, c1);
        (void) sc_list_prepend (list, c0);

        if (replace_fn != NULL) {
          /* in family mode we always call the replace callback right
           * away */
          family[0] = c0;
          family[1] = c1;
          family[2] = c2;
          family[3] = c3;
#ifdef P4_TO_P8
          family[4] = c4;
          family[5] = c5;
          family[6] = c6;
          family[7] = c7;
#endif
          replace_fn (p4est, nt, 1, &pp, P4EST_CHILDREN, family);
          p4est_quadrant_free_data (p4est, &parent);
        }
      }
      else {
        /* need to make room in the array to store this new quadrant */
        if (restpos < incount && current == restpos) {
          movecount = SC_MIN (incount - restpos, number_toread_quadrants);
          while (movecount > 0) {
            q = p4est_quadrant_array_index (tquadrants, restpos);
            qalloc = p4est_quadrant_mempool_alloc (p4est->quadrant_pool);
            *qalloc = *q;       /* never append array members directly */
            qalloc->pad8 = 0;   /* has not been refined yet */
            (void) sc_list_append (list, qalloc);       /* only new quadrants */
            --movecount;
            ++restpos;
          }
        }

        /* store new quadrant and update counters */
        q = p4est_quadrant_array_index (tquadrants, current);
        *q = *qpop;
        maxlevel = SC_MAX (maxlevel, (int) qpop->level);
        ++tree->quadrants_per_level[qpop->level];
        ++current;
        sc_mempool_free (p4est->quadrant_pool, qpop);
      }
    }
    tree->maxlevel = (int8_t) maxlevel;
    p4est->local_num_quadrants += tquadrants->elem_count;

    P4EST_ASSERT (restpos == incount);
    P4EST_ASSERT (current == tquadrants->elem_count);
    P4EST_ASSERT (list->first == NULL && list->last == NULL);
    P4EST_ASSERT (quadrant_pool_size == p4est->quadrant_pool->elem_count);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size + tquadrants->elem_count ==
                    p4est->user_data_pool->elem_count + incount);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    P4EST_VERBOSEF ("Done refine tree %lld now %llu\n", (long long) nt,
                    (unsigned long long) tquadrants->elem_count);
  }
  if (p4est->last_local_tree >= 0) {
    for (; nt < p4est->connectivity->num_trees; ++nt) {
      tree = p4est_tree_array_index (p4est->trees, nt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  sc_list_destroy (list);

  P4EST_FUNC_SHOT (p4est, &snap);
}

void
p4est_coarsen (p4est_t * p4est, int coarsen_recursive,
               p4est_coarsen_t coarsen_fn, p4est_init_t init_fn)
{
  p4est_coarsen_ext (p4est, coarsen_recursive, 0, coarsen_fn, init_fn, NULL);
}

void
p4est_coarsen_ext (p4est_t * p4est, int coarsen_recursive,
                   int callback_orphans, p4est_coarsen_t coarsen_fn,
                   p4est_init_t init_fn, p4est_replace_t replace_fn)
{
  p4est_gloidx_t      old_gnq;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_coarsen with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
  p4est_log_indent_push ();
  P4EST_ASSERT (p4est_is_valid (p4est));
  /* remember input quadrant count; it will not increase */
  old_gnq = p4est->global_num_quadrants;
  p4est_coarsen_ext_dirty (p4est, coarsen_recursive, callback_orphans,
                           coarsen_fn, init_fn, replace_fn);
  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);
  P4EST_ASSERT (p4est->global_num_quadrants <= old_gnq);
  if (old_gnq != p4est->global_num_quadrants) {
    ++p4est->revision;
  }

  P4EST_ASSERT (p4est_is_valid (p4est));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_coarsen with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

void
p4est_coarsen_ext_dirty (p4est_t * p4est,
                         int coarsen_recursive, int callback_orphans,
                         p4est_coarsen_t coarsen_fn, p4est_init_t init_fn,
                         p4est_replace_t replace_fn)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              data_pool_size;
#endif
  int                 i, maxlevel;
  int                 isfamily;
  size_t              zz;
  size_t              incount, removed;
  size_t              window, start, length, cidz;
  p4est_locidx_t      num_quadrants, prev_offset;
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *c[P4EST_CHILDREN];
  p4est_quadrant_t   *cfirst, *clast;
  sc_array_t         *tquadrants;
  sc_flopinfo_t       snap;
  p4est_quadrant_t    qtemp;

  P4EST_FUNC_SNAP (p4est, &snap);

  P4EST_ASSERT (coarsen_fn != NULL);

  P4EST_QUADRANT_INIT (&qtemp);

  /* loop over all local trees */
  prev_offset = 0;
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;
#ifdef P4EST_ENABLE_DEBUG
    data_pool_size = 0;
    if (p4est->user_data_pool != NULL) {
      data_pool_size = p4est->user_data_pool->elem_count;
    }
#endif
    removed = 0;

    /* initial log message for this tree */
    P4EST_VERBOSEF ("Into coarsen tree %lld with %llu\n", (long long) jt,
                    (unsigned long long) tquadrants->elem_count);

    /* state information */
    window = 0;                 /* start position of sliding window in array */
    start = 1;                  /* start position of hole in window/array */
    length = 0;                 /* length of hole in window/array */

    /* run through the array and coarsen recursively */
    incount = tquadrants->elem_count;
    while (window + P4EST_CHILDREN + length <= incount) {
      P4EST_ASSERT (window < start);

      cidz = incount;
      isfamily = 1;
      for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
        c[zz] = (window + zz < start) ?
          p4est_quadrant_array_index (tquadrants, window + zz) :
          p4est_quadrant_array_index (tquadrants, window + length + zz);

        if (zz != (size_t) p4est_quadrant_child_id (c[zz])) {
          isfamily = 0;
          if (callback_orphans) {
            c[1] = NULL;
            (void) coarsen_fn (p4est, jt, c);
          }
          break;
        }
      }
      /* in a complete tree, the only way P4EST_CHILDREN consecutive quadrants
       * have the correct consecutive child_id's is if they are, in fact, a
       * family.
       */
      P4EST_ASSERT (!isfamily || p4est_quadrant_is_familypv (c));
      if (isfamily && coarsen_fn (p4est, jt, c)) {
        /* coarsen this family of quadrants */
        if (replace_fn == NULL) {
          for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
            p4est_quadrant_free_data (p4est, c[zz]);
          }
        }
        tree->quadrants_per_level[c[0]->level] -= P4EST_CHILDREN;
        cfirst = c[0];
        if (replace_fn != NULL) {
          qtemp = *(c[0]);
          c[0] = &qtemp;
        }
        p4est_quadrant_parent (c[0], cfirst);
        p4est_quadrant_init_data (p4est, jt, cfirst, init_fn);
        tree->quadrants_per_level[cfirst->level] += 1;
        p4est->local_num_quadrants -= P4EST_CHILDREN - 1;
        removed += P4EST_CHILDREN - 1;

        cidz = (size_t) p4est_quadrant_child_id (cfirst);
        start = window + 1;
        length += P4EST_CHILDREN - 1;

        if (replace_fn != NULL) {
          replace_fn (p4est, jt, P4EST_CHILDREN, c, 1, &cfirst);
          for (zz = 0; zz < P4EST_CHILDREN; zz++) {
            p4est_quadrant_free_data (p4est, c[zz]);
          }
        }
      }

      if (cidz <= window && coarsen_recursive) {
        window -= cidz;
      }
      else {
        ++window;
        if (window == start && start + length < incount) {
          if (length > 0) {
            cfirst = p4est_quadrant_array_index (tquadrants, start);
            clast = p4est_quadrant_array_index (tquadrants, start + length);
            *cfirst = *clast;
          }
          start = window + 1;
        }
      }
    }

    /* adjust final array size */
    if (length > 0) {
      for (zz = start + length; zz < incount; ++zz) {
        cfirst = p4est_quadrant_array_index (tquadrants, zz - length);
        clast = p4est_quadrant_array_index (tquadrants, zz);
        *cfirst = *clast;
      }
      sc_array_resize (tquadrants, incount - length);
    }

    /* call remaining orphans */
    if (callback_orphans) {
      c[1] = NULL;
      for (zz = window; zz < incount - length; ++zz) {
        c[0] = p4est_quadrant_array_index (tquadrants, zz);
        (void) coarsen_fn (p4est, jt, c);
      }
    }

    /* compute maximum level */
    maxlevel = 0;
    num_quadrants = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      P4EST_ASSERT (tree->quadrants_per_level[i] >= 0);
      num_quadrants += tree->quadrants_per_level[i];    /* same type */
      if (tree->quadrants_per_level[i] > 0) {
        maxlevel = i;
      }
    }
    tree->maxlevel = (int8_t) maxlevel;
    tree->quadrants_offset = prev_offset;
    prev_offset += num_quadrants;

    /* do some sanity checks */
    P4EST_ASSERT (num_quadrants == (p4est_locidx_t) tquadrants->elem_count);
    P4EST_ASSERT (tquadrants->elem_count == incount - removed);
    if (p4est->user_data_pool != NULL) {
      P4EST_ASSERT (data_pool_size - removed ==
                    p4est->user_data_pool->elem_count);
    }
    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (p4est_tree_is_complete (tree));

    /* final log message for this tree */
    P4EST_VERBOSEF ("Done coarsen tree %lld now %llu\n", (long long) jt,
                    (unsigned long long) tquadrants->elem_count);
  }
  if (p4est->last_local_tree >= 0) {
    for (; jt < p4est->connectivity->num_trees; ++jt) {
      tree = p4est_tree_array_index (p4est->trees, jt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  P4EST_FUNC_SHOT (p4est, &snap);
}

void
p4est_balance (p4est_t * p4est, p4est_connect_type_t btype,
               p4est_init_t init_fn)
{
  p4est_balance_ext (p4est, btype, init_fn, NULL);
}

extern void p4est_balance_tworound (p4est_t *p4est, p4est_connect_type_t btype,
                                    p4est_init_t init_fn, p4est_replace_t replace_fn);
extern void p4est_balance_sort (p4est_t *p4est, p4est_connect_type_t btype,
                                p4est_init_t init_fn, p4est_replace_t replace_fn);

void
p4est_balance_ext (p4est_t * p4est, p4est_connect_type_t btype,
                   p4est_init_t init_fn, p4est_replace_t replace_fn)
{
  p4est_balance_obj_t *bobj;

  bobj = p4est_balance_obj_new (p4est->mpicomm);

  p4est_balance_obj_set_connect (bobj, btype);
  p4est_balance_obj_set_init (bobj, init_fn);
  p4est_balance_obj_set_replace (bobj, replace_fn);
  p4est_balance_obj (bobj, p4est);
  p4est_balance_obj_destroy (bobj);
}

void
p4est_partition (p4est_t * p4est, int allow_for_coarsening,
                 p4est_weight_t weight_fn)
{
  (void) p4est_partition_ext (p4est, allow_for_coarsening, weight_fn);
}

/* p4est may not be valid yet: we expect
   each process to know its own local number of quadrants, but the global
   total may not be known yet */
p4est_locidx_t     *
p4est_partition_compute (p4est_t * p4est, int partition_for_coarsening,
                         p4est_gloidx_t global_num_quadrants,
                         p4est_weight_t weight_fn)
{
  p4est_gloidx_t      global_shipped = 0;
  sc_flopinfo_t       snap;
#ifdef P4EST_ENABLE_MPI
  int                 mpiret;
  int                 low_source, high_source;
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
  const p4est_locidx_t local_num_quadrants = p4est->local_num_quadrants;
  int                 i, p;
  int                 send_lowest, send_highest;
  int                 num_sends, rcount, base_index;
  size_t              lz;
  ssize_t             lowers;
  p4est_topidx_t      nt;
  p4est_locidx_t      kl, qlocal;
  p4est_locidx_t     *num_quadrants_in_proc;
  p4est_gloidx_t      prev_quadrant, next_quadrant;
  p4est_gloidx_t      send_index, recv_low, recv_high, qcount;
  p4est_gloidx_t     *send_array;
  int64_t             weight, weight_sum;
  int64_t             cut, my_lowcut, my_highcut;
  int64_t            *local_weights;    /* cumulative weights by quadrant */
  int64_t            *global_weight_sums;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  MPI_Request        *send_requests, recv_requests[2];
  MPI_Status          recv_statuses[2];
  p4est_gloidx_t      num_corrected;
#endif /* P4EST_ENABLE_MPI */

  P4EST_FUNC_SNAP (p4est, &snap);

  /* allocate new quadrant distribution counts */
  num_quadrants_in_proc = P4EST_ALLOC (p4est_locidx_t, num_procs);

  if (weight_fn == NULL) {
    /* Divide up the quadrants equally */
    if (global_num_quadrants < 0) {
      /* global count is not known */
      p4est_gloidx_t      glocal_num_quadrants =
        (p4est_gloidx_t) local_num_quadrants;

      mpiret =
        sc_MPI_Allreduce (&glocal_num_quadrants, &global_num_quadrants, 1,
                          P4EST_MPI_GLOIDX, MPI_SUM, p4est->mpicomm);
      SC_CHECK_MPI (mpiret);
      P4EST_ASSERT (global_num_quadrants >= 0);
    }
    for (p = 0, next_quadrant = 0; p < num_procs; ++p) {
      prev_quadrant = next_quadrant;
      next_quadrant =
        p4est_partition_cut_gloidx (global_num_quadrants, p + 1, num_procs);
      qcount = next_quadrant - prev_quadrant;
      P4EST_ASSERT (0 <= qcount
                    && qcount <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);
      num_quadrants_in_proc[p] = (p4est_locidx_t) (qcount);
    }
  }
  else {
    /* do a weighted partition */
    local_weights = P4EST_ALLOC (int64_t, local_num_quadrants + 1);
    global_weight_sums = P4EST_ALLOC (int64_t, num_procs + 1);
    P4EST_VERBOSEF ("local quadrant count %lld\n",
                    (long long) local_num_quadrants);

    /* linearly sum weights across all trees */
    kl = 0;
    local_weights[0] = 0;
    for (nt = first_tree; nt <= last_tree; ++nt) {
      tree = p4est_tree_array_index (p4est->trees, nt);
      for (lz = 0; lz < tree->quadrants.elem_count; ++lz, ++kl) {
        q = p4est_quadrant_array_index (&tree->quadrants, lz);
        weight = (int64_t) weight_fn (p4est, nt, q);
        P4EST_ASSERT (weight >= 0);
        local_weights[kl + 1] = local_weights[kl] + weight;
      }
    }
    P4EST_ASSERT (kl == local_num_quadrants);
    weight_sum = local_weights[local_num_quadrants];
    P4EST_VERBOSEF ("local weight sum %lld\n", (long long) weight_sum);

    /* distribute local weight sums */
    global_weight_sums[0] = 0;
    mpiret = MPI_Allgather (&weight_sum, 1, MPI_LONG_LONG_INT,
                            &global_weight_sums[1], 1, MPI_LONG_LONG_INT,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    /* adjust all arrays to reflect the global weight */
    for (i = 0; i < num_procs; ++i) {
      global_weight_sums[i + 1] += global_weight_sums[i];
    }
    if (rank > 0) {
      weight_sum = global_weight_sums[rank];
      for (kl = 0; kl <= local_num_quadrants; ++kl) {
        local_weights[kl] += weight_sum;
      }
    }
    P4EST_ASSERT (local_weights[0] == global_weight_sums[rank]);
    P4EST_ASSERT (local_weights[local_num_quadrants] ==
                  global_weight_sums[rank + 1]);
    weight_sum = global_weight_sums[num_procs];

    if (rank == 0) {
      for (i = 0; i <= num_procs; ++i) {
        P4EST_GLOBAL_VERBOSEF ("Global weight sum [%d] %lld\n",
                               i, (long long) global_weight_sums[i]);
      }
    }

    /* if all quadrants have zero weight we do nothing */
    if (weight_sum == 0) {
      P4EST_FREE (local_weights);
      P4EST_FREE (global_weight_sums);
      P4EST_FREE (num_quadrants_in_proc);

      /* in particular, there is no need to bumb the revision counter */
      P4EST_ASSERT (global_shipped == 0);

      P4EST_FUNC_SHOT (p4est, &snap);
      return NULL;              /* NULL inidicates no shipping to be done */
    }

    /* determine processor ids to send to */
    send_lowest = num_procs;
    send_highest = 0;
    for (i = 1; i <= num_procs; ++i) {
      cut = p4est_partition_cut_uint64 (weight_sum, i, num_procs);
      if (global_weight_sums[rank] < cut &&
          cut <= global_weight_sums[rank + 1]) {
        send_lowest = SC_MIN (send_lowest, i);
        send_highest = SC_MAX (send_highest, i);
      }
    }
    /*
     * send low cut to send_lowest..send_highest
     * and high cut to send_lowest-1..send_highest-1
     */
    P4EST_LDEBUGF ("my send peers %d %d\n", send_lowest, send_highest);

    num_sends = 2 * (send_highest - send_lowest + 1);
    if (num_sends <= 0) {
      num_sends = 0;
      send_requests = NULL;
      send_array = NULL;
    }
    else {
      send_requests = P4EST_ALLOC (MPI_Request, num_sends);
      send_array = P4EST_ALLOC (p4est_gloidx_t, num_sends);
      lowers = 0;
      for (i = send_lowest; i <= send_highest; ++i) {
        base_index = 2 * (i - send_lowest);
        if (i < num_procs) {
          /* do binary search in the weight array */
          lowers = sc_search_lower_bound64 (p4est_partition_cut_uint64
                                            (weight_sum, i, num_procs),
                                            local_weights,
                                            (size_t) local_num_quadrants + 1,
                                            (size_t) lowers);
          P4EST_ASSERT (lowers > 0
                        && (p4est_locidx_t) lowers <= local_num_quadrants);

          /* send low bound */
          send_index = send_array[base_index + 1] =
            (p4est_gloidx_t) lowers + p4est->global_first_quadrant[rank];
          P4EST_LDEBUGF ("send A %d %d %d index %lld base %d to %d\n",
                         send_lowest, i, send_highest,
                         (long long) send_index, base_index + 1, i);
          mpiret =
            MPI_Isend (&send_array[base_index + 1], 1, P4EST_MPI_GLOIDX, i,
                       P4EST_COMM_PARTITION_WEIGHTED_LOW, p4est->mpicomm,
                       &send_requests[base_index + 1]);
          SC_CHECK_MPI (mpiret);
        }
        else {
          lowers = 0;
          send_index = global_num_quadrants;
          send_requests[base_index + 1] = MPI_REQUEST_NULL;
          send_array[base_index + 1] = -1;
        }

        /* send high bound */
        send_array[base_index] = send_index;
        P4EST_LDEBUGF ("send B %d %d %d index %lld base %d to %d\n",
                       send_lowest, i, send_highest,
                       (long long) send_index, base_index, i - 1);
        mpiret = MPI_Isend (&send_array[base_index], 1, P4EST_MPI_GLOIDX,
                            i - 1, P4EST_COMM_PARTITION_WEIGHTED_HIGH,
                            p4est->mpicomm, &send_requests[base_index]);
        SC_CHECK_MPI (mpiret);
      }
    }

    /* determine processor ids to receive from and post irecv */
    i = 0;
    my_lowcut = p4est_partition_cut_uint64 (weight_sum, rank, num_procs);
    if (my_lowcut == 0) {
      recv_low = 0;
      recv_requests[0] = MPI_REQUEST_NULL;
      low_source = -1;
    }
    else {
      for (; i < num_procs; ++i) {
        if (global_weight_sums[i] < my_lowcut &&
            my_lowcut <= global_weight_sums[i + 1]) {
          P4EST_LDEBUGF ("recv A from %d\n", i);
          mpiret = MPI_Irecv (&recv_low, 1, P4EST_MPI_GLOIDX, i,
                              P4EST_COMM_PARTITION_WEIGHTED_LOW,
                              p4est->mpicomm, &recv_requests[0]);
          SC_CHECK_MPI (mpiret);
          break;
        }
      }
      P4EST_ASSERT (i < num_procs);
      low_source = i;
    }
    my_highcut = p4est_partition_cut_uint64 (weight_sum, rank + 1, num_procs);
    if (my_highcut == 0) {
      recv_high = 0;
      recv_requests[1] = MPI_REQUEST_NULL;
      high_source = -1;
    }
    else {
      for (; i < num_procs; ++i) {
        if (global_weight_sums[i] < my_highcut &&
            my_highcut <= global_weight_sums[i + 1]) {
          P4EST_LDEBUGF ("recv B from %d\n", i);
          mpiret = MPI_Irecv (&recv_high, 1, P4EST_MPI_GLOIDX, i,
                              P4EST_COMM_PARTITION_WEIGHTED_HIGH,
                              p4est->mpicomm, &recv_requests[1]);
          SC_CHECK_MPI (mpiret);
          break;
        }
      }
      P4EST_ASSERT (i < num_procs);
      high_source = i;
    }
    P4EST_LDEBUGF ("my recv peers %d %d cuts %lld %lld\n",
                   low_source, high_source,
                   (long long) my_lowcut, (long long) my_highcut);

    /* free temporary memory */
    P4EST_FREE (local_weights);
    P4EST_FREE (global_weight_sums);

    /* wait for sends and receives to complete */
    if (num_sends > 0) {
      mpiret = MPI_Waitall (num_sends, send_requests, MPI_STATUSES_IGNORE);
      SC_CHECK_MPI (mpiret);
      P4EST_FREE (send_requests);
      P4EST_FREE (send_array);
    }
    mpiret = MPI_Waitall (2, recv_requests, recv_statuses);
    SC_CHECK_MPI (mpiret);
    if (my_lowcut != 0) {
      SC_CHECK_ABORT (recv_statuses[0].MPI_SOURCE == low_source,
                      "Wait low source");
      SC_CHECK_ABORT (recv_statuses[0].MPI_TAG ==
                      P4EST_COMM_PARTITION_WEIGHTED_LOW, "Wait low tag");
      mpiret = MPI_Get_count (&recv_statuses[0], P4EST_MPI_GLOIDX, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount == 1, "Wait low count %d", rcount);
    }
    if (my_highcut != 0) {
      SC_CHECK_ABORT (recv_statuses[1].MPI_SOURCE == high_source,
                      "Wait high source");
      SC_CHECK_ABORT (recv_statuses[1].MPI_TAG ==
                      P4EST_COMM_PARTITION_WEIGHTED_HIGH, "Wait high tag");
      mpiret = MPI_Get_count (&recv_statuses[1], P4EST_MPI_GLOIDX, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount == 1, "Wait high count %d", rcount);
    }

    /* communicate the quadrant ranges */
    qcount = recv_high - recv_low;
    P4EST_LDEBUGF ("weighted partition count %lld\n", (long long) qcount);
    P4EST_ASSERT (qcount >= 0 && qcount <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);
    qlocal = (p4est_locidx_t) qcount;
    mpiret = MPI_Allgather (&qlocal, 1, P4EST_MPI_LOCIDX,
                            num_quadrants_in_proc, 1, P4EST_MPI_LOCIDX,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

#if(0)
    /* run through the count array and repair zero ranges */
    for (i = 0; i < num_procs; ++i) {
      if (num_quadrants_in_proc[i] == 0) {
        for (p = i - 1; p >= 0; --p) {
          P4EST_ASSERT (num_quadrants_in_proc[p] > 0);
          if (num_quadrants_in_proc[p] > 1) {
            --num_quadrants_in_proc[p];
            ++num_quadrants_in_proc[i];
            break;
          }
        }
        if (p < 0) {
          for (p = i + 1; p < num_procs; ++p) {
            P4EST_ASSERT (num_quadrants_in_proc[p] >= 0);
            if (num_quadrants_in_proc[p] > 1) {
              --num_quadrants_in_proc[p];
              ++num_quadrants_in_proc[i];
              break;
            }
          }
          P4EST_ASSERT (p < num_procs);
        }
      }
    }
#endif
  }

  /* correct partition */
  if (partition_for_coarsening) {
    num_corrected =
      p4est_partition_for_coarsening (p4est, num_quadrants_in_proc);
    P4EST_GLOBAL_INFOF
      ("Designated partition for coarsening %lld quadrants moved\n",
       (long long) num_corrected);
  }

  P4EST_FUNC_SHOT (p4est, &snap);
  return num_quadrants_in_proc;
}

p4est_gloidx_t
p4est_partition_ext (p4est_t * p4est, int partition_for_coarsening,
                     p4est_weight_t weight_fn)
{
  p4est_gloidx_t      global_shipped = 0;
  const p4est_gloidx_t global_num_quadrants = p4est->global_num_quadrants;
#ifdef P4EST_ENABLE_MPI
  p4est_locidx_t     *num_quadrants_in_proc;
#endif /* P4EST_ENABLE_MPI */
  sc_flopinfo_t       snap;

  P4EST_FUNC_SNAP (p4est, &snap);

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_GLOBAL_PRODUCTIONF
    ("Into " P4EST_STRING
     "_partition with %lld total quadrants\n",
     (long long) p4est->global_num_quadrants);

  /* this function does nothing in a serial setup */
  if (p4est->mpisize == 1) {
    P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_partition no shipping\n");

    /* in particular, there is no need to bumb the revision counter */
    P4EST_ASSERT (global_shipped == 0);
    return global_shipped;
  }

  p4est_log_indent_push ();

#ifdef P4EST_ENABLE_MPI
  /* allocate new quadrant distribution counts */
  num_quadrants_in_proc =
    p4est_partition_compute (p4est, partition_for_coarsening,
                             p4est->global_num_quadrants, weight_fn);

  if (!num_quadrants_in_proc) {
    P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_partition no shipping\n");
    P4EST_FUNC_SHOT (p4est, &snap);
    return 0;
  }

  /* run the partition algorithm with proper quadrant counts */
  global_shipped = p4est_partition_given (p4est, num_quadrants_in_proc);
  if (global_shipped) {
    /* the partition of the forest has changed somewhere */
    ++p4est->revision;
  }
  P4EST_FREE (num_quadrants_in_proc);

  /* check validity of the p4est */
  P4EST_ASSERT (p4est_is_valid (p4est));
#endif /* P4EST_ENABLE_MPI */

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done " P4EST_STRING "_partition shipped %lld quadrants %.3g%%\n",
     (long long) global_shipped,
     global_shipped * 100. / global_num_quadrants);

  P4EST_FUNC_SHOT (p4est, &snap);
  return global_shipped;
}

p4est_gloidx_t
p4est_partition_for_coarsening (p4est_t * p4est,
                                p4est_locidx_t * num_quadrants_in_proc)
{
#ifdef P4EST_ENABLE_MPI
  int                 num_procs = p4est->mpisize;
  int                 rank = p4est->mpirank;
  int                 mpiret;
  p4est_gloidx_t      global_num_quadrants = p4est->global_num_quadrants;
  int                 i, send_lowest, send_highest, num_sends;
  int                 parent_index;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;
  p4est_locidx_t      num_quadrants_in_tree;
  p4est_topidx_t      it, tree_index;
  p4est_gloidx_t      iq, quad_id_near_cut;
  p4est_gloidx_t      min_quad_id, max_quad_id;
  int8_t              quad_near_cut_level;
  p4est_gloidx_t     *partition_now = p4est->global_first_quadrant;
  p4est_gloidx_t     *partition_new;
  p4est_quadrant_t   *parent_send;
  MPI_Request        *send_requests;
  MPI_Request        *receive_requests;
  int                 receive_lowest, receive_highest, num_receives;
  int                 process_with_cut = -1, process_with_cut_recv_id = -1;
  p4est_quadrant_t   *parent_receive;
  int                *receive_process;
  int                *correction, correction_local = 0;
  int                 current_proc, next_proc;
  p4est_gloidx_t      num_moved_quadrants;
  sc_flopinfo_t       snap;

  P4EST_FUNC_SNAP (p4est, &snap);

  /* create array with first quadrants of new partition */
  partition_new = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  partition_new[0] = 0;
  for (i = 1; i < num_procs; i++) {     /* loop over all processes */
    partition_new[i] = partition_new[i - 1] + num_quadrants_in_proc[i - 1];
  }
  partition_new[num_procs] = global_num_quadrants;

  /* BEGIN: send */
  if (partition_now[rank] < partition_now[rank + 1]) {
    /* if this process has quadrants */
    /* determine number and min/max process ids to send to */
    num_sends = 0;              /* number of sends */
    send_lowest = num_procs;    /* lowest process id */
    send_highest = 0;           /* highest process id */
    for (i = 1; i < num_procs; i++) {
      /* loop over all processes (without first) */
      if (partition_new[i] < partition_new[i + 1] &&
          partition_now[rank] <= partition_new[i] + P4EST_CHILDREN - 2 &&
          partition_new[i] - P4EST_CHILDREN + 1 < partition_now[rank + 1]) {
        /* if this process has relevant quadrants for process `i` */
        num_sends++;
        send_lowest = SC_MIN (send_lowest, i);
        send_highest = SC_MAX (send_highest, i);
      }
    }
  }
  else {
    /* set number of messages to send */
    num_sends = 0;
  }

  if (num_sends > 0) {          /* if this process sends messages */
    /* allocate send messages */
    send_requests = P4EST_ALLOC (MPI_Request, num_sends);
    parent_send = P4EST_ALLOC (p4est_quadrant_t, num_sends);

    /* array index of send messages */
    parent_index = 0;

    for (i = send_lowest; i <= send_highest; i++) {
      /* loop over all process candidates to send to */
      if (!(partition_new[i] < partition_new[i + 1] &&
            partition_now[rank] <= partition_new[i] + P4EST_CHILDREN - 2 &&
            partition_new[i] - P4EST_CHILDREN + 1 <
            partition_now[rank + 1])) {
        /* if this process has no relevant quadrants for process `i` */
        continue;
      }

      /* get nearest quadrant `quad_id_near_cut` to cut `partition_new[i]` */
      if (partition_now[rank] <= partition_new[i] &&
          partition_new[i] < partition_now[rank + 1]) {
        /* if cut is owned by this process */
        quad_id_near_cut = partition_new[i];
      }
      else {
        if (P4EST_GLOIDX_ABS (partition_new[i] - partition_now[rank]) <
            P4EST_GLOIDX_ABS (partition_new[i] - partition_now[rank + 1] +
                              1)) {
          quad_id_near_cut = partition_now[rank];
        }
        else {
          quad_id_near_cut = partition_now[rank + 1] - 1;
        }
      }

      /* get tree `tree` of quadrant `quad_id_near_cut` */
      num_quadrants_in_tree = partition_now[rank + 1] - partition_now[rank];
      tree_index = -1;
      for (it = p4est->first_local_tree; it <= p4est->last_local_tree; it++) {
        /* loop over all local trees */
        tree = p4est_tree_array_index (p4est->trees, it);
        if (tree->quadrants_offset <= quad_id_near_cut - partition_now[rank]) {
          tree_index = it;
        }
        else {
          num_quadrants_in_tree = tree->quadrants_offset;
          break;
        }
      }
      tree = p4est_tree_array_index (p4est->trees, tree_index);
      num_quadrants_in_tree -= tree->quadrants_offset;

      /* get quadrant with index `quad_id_near_cut` */
      q = p4est_quadrant_array_index (&tree->quadrants,
                                      quad_id_near_cut - partition_now[rank] -
                                      tree->quadrants_offset);

      /* get level of quadrant near cut */
      quad_near_cut_level = q->level;

      if (quad_near_cut_level > 0) {
        /* if quadrant near cut is not root of tree, i.e. level is not zero */
        /* get parent of quadrant near cut */
        p4est_quadrant_parent (q, &parent_send[parent_index]);

        /* get min quadrant with same parent */
        min_quad_id = quad_id_near_cut;
        for (iq = quad_id_near_cut;
             iq >= SC_MAX (partition_now[rank] + tree->quadrants_offset,
                           partition_new[i] - P4EST_CHILDREN + 1); iq--) {
          /* loop over eligible quadrants */
          /* get quadrant with index `iq` */
          q = p4est_quadrant_array_index (&tree->quadrants,
                                          iq - partition_now[rank] -
                                          tree->quadrants_offset);

          /* check quadrant `iq` */
          if (q->level == quad_near_cut_level) {        /* if same level */
            if (p4est_quadrant_is_parent (&parent_send[parent_index], q)) {
              /* if same parent */
              min_quad_id = iq;
            }
            else {
              break;
            }
          }
          else {
            break;
          }
        }

        /* get max quadrant with same parent */
        max_quad_id = quad_id_near_cut;
        for (iq = quad_id_near_cut;
             iq <=
             SC_MIN (partition_now[rank] + tree->quadrants_offset +
                     num_quadrants_in_tree - 1,
                     partition_new[i] + P4EST_CHILDREN - 2); iq++) {
          /* loop over eligible quadrants */

          /* get quadrant `iq` */
          q = p4est_quadrant_array_index (&tree->quadrants,
                                          iq - partition_now[rank] -
                                          tree->quadrants_offset);

          /* check quadrant `iq` */
          if (q->level == quad_near_cut_level) {        /* if same level */
            if (p4est_quadrant_is_parent (&parent_send[parent_index], q)) {
              /* if same parent */
              max_quad_id = iq;
            }
            else {
              break;
            }
          }
          else {
            break;
          }
        }

        /* write tree */
        parent_send[parent_index].p.piggy3.which_tree = tree_index;

        if (quad_id_near_cut == partition_new[i]) {
          /* if this process has cut */
          /* encode number of quadrants with same parent before and after
           * `partition_new[i]` into one integer */
          parent_send[parent_index].p.piggy3.local_num =
            (partition_new[i] - min_quad_id) * P4EST_CHILDREN +
            (max_quad_id - partition_new[i]);
        }
        else {
          /* write number of quadrants with same parent */
          parent_send[parent_index].p.piggy3.local_num =
            max_quad_id - min_quad_id + 1;
        }

        /* MPI send: parent */
        mpiret = MPI_Isend (&parent_send[parent_index],
                            sizeof (p4est_quadrant_t), MPI_BYTE,
                            i,
                            P4EST_COMM_PARTITION_CORRECTION, p4est->mpicomm,
                            &send_requests[parent_index]);
        SC_CHECK_MPI (mpiret);
      }
      else {
        /* if quadrant near cut is root of tree, i.e., level is zero,
           set parent as tree root `q` */
        parent_send[parent_index].level = q->level;
        parent_send[parent_index].x = q->x;
        parent_send[parent_index].y = q->y;
#ifdef P4_TO_P8
        parent_send[parent_index].z = q->z;
#endif

        /* write tree */
        parent_send[parent_index].p.piggy3.which_tree = tree_index;

        /* write number of quadrants with same "parent" */
        parent_send[parent_index].p.piggy3.local_num = 1;

        /* MPI send: root of tree */
        mpiret = MPI_Isend (&parent_send[parent_index],
                            sizeof (p4est_quadrant_t), MPI_BYTE,
                            i,
                            P4EST_COMM_PARTITION_CORRECTION, p4est->mpicomm,
                            &send_requests[parent_index]);
        SC_CHECK_MPI (mpiret);
      }

      /* increment parent index */
      parent_index++;
    }
  }
  /* END: send */

  /* BEGIN: receive */
  if (rank != 0 && partition_new[rank] < partition_new[rank + 1]) {
    /* if this process should get quadrants */
    /* determine process ids to receive from */
    num_receives = 0;           /* number of receives */
    receive_lowest = num_procs; /* lowest process id */
    receive_highest = 0;        /* highest process id */
    for (i = 0; i < num_procs; i++) {   /* loop over all processes */
      if (partition_now[i] < partition_now[i + 1] &&
          partition_now[i] <= partition_new[rank] + P4EST_CHILDREN - 2 &&
          partition_new[rank] - P4EST_CHILDREN + 1 < partition_now[i + 1]) {
        /* if process `i` has relevant quadrants for this process */
        num_receives++;
        receive_lowest = SC_MIN (receive_lowest, i);
        receive_highest = SC_MAX (receive_highest, i);

        if (partition_now[i] <= partition_new[rank] &&
            partition_new[rank] < partition_now[i + 1]) {
          /* if cut is owned by process `i` */
          /* process id that sends parent of cut quadrant */
          process_with_cut = i;

          /* array index of receive messages of process with cut  */
          process_with_cut_recv_id = num_receives - 1;
        }
      }
    }
  }
  else {
    /* set number of messages to receive */
    num_receives = 0;

    /* set correction */
    correction_local = 0;
  }

  if (num_receives > 0) {       /* if this process receives messages */
    /* allocate receive messages */
    receive_requests = P4EST_ALLOC (MPI_Request, num_receives);
    parent_receive = P4EST_ALLOC (p4est_quadrant_t, num_receives);
    receive_process = P4EST_ALLOC (int, num_receives);

    /* array index of receive messages */
    parent_index = 0;

    for (i = receive_lowest; i <= receive_highest; i++) {
      /* loop over all process candidates to receive from */
      if (!(partition_now[i] < partition_now[i + 1] &&
            partition_now[i] <= partition_new[rank] + P4EST_CHILDREN - 2 &&
            partition_new[rank] - P4EST_CHILDREN + 1 <
            partition_now[i + 1])) {
        /* if process `i` has no relevant quadrants for this process */
        continue;
      }
      /* store process index */
      receive_process[parent_index] = i;

      /* MPI receive */
      mpiret = MPI_Irecv (&parent_receive[parent_index],
                          sizeof (p4est_quadrant_t), MPI_BYTE,
                          i,
                          P4EST_COMM_PARTITION_CORRECTION, p4est->mpicomm,
                          &receive_requests[parent_index]);
      SC_CHECK_MPI (mpiret);

      /* increment parent index */
      parent_index++;
    }
  }
  /* END: receive */

  /* BEGIN: wait for MPI receive to complete */
  if (num_receives > 0) {
    /* wait for receives to complete */
    mpiret =
      MPI_Waitall (num_receives, receive_requests, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    /* free receive memory */
    P4EST_FREE (receive_requests);
  }
  /* END: wait for MPI recieve to complete */

  /* BEGIN: compute correction with received quadrants */
  if (num_receives > 0) {
    /* if this process received quadrants */
    min_quad_id = partition_new[rank];  /* min quadrant id with same parent */
    max_quad_id = partition_new[rank];  /* max quadrant id with same parent */

    for (i = 0; i < num_receives; i++) {
      /* loop over all received (parent or tree root) quadrants */
      if (parent_receive[i].p.piggy3.which_tree ==
          parent_receive[process_with_cut_recv_id].p.piggy3.which_tree
          &&
          p4est_quadrant_is_equal (&parent_receive[i],
                                   &parent_receive[process_with_cut_recv_id]
          )) {
        /* if trees and parents are equal */
        /* decrease/increase min/max quadrant with same parent */
        if (receive_process[i] < process_with_cut) {
          /* if before process with cut */
          /* decrease min quadrant */
          min_quad_id -= parent_receive[i].p.piggy3.local_num;
        }
        else if (receive_process[i] > process_with_cut) {
          /* if after process with cut */
          /* increase max quadrant */
          max_quad_id += parent_receive[i].p.piggy3.local_num;
        }
        else {
          /* decrease min quadrant */
          min_quad_id -=
            parent_receive[i].p.piggy3.local_num / P4EST_CHILDREN;

          /* increase max quadrant */
          max_quad_id +=
            parent_receive[i].p.piggy3.local_num % P4EST_CHILDREN;
        }
      }
    }

    /* compute correction */
    correction_local =
      (int) p4est_partition_correction (partition_new, num_procs, rank,
                                        min_quad_id, max_quad_id);

    /* free receive memory */
    P4EST_FREE (parent_receive);
    P4EST_FREE (receive_process);
  }
  /* END: compute correction with received parent quadrants */

  /* free memory */
  P4EST_FREE (partition_new);

  /* communicate corrections */
  correction = P4EST_ALLOC (int, num_procs);
  mpiret = MPI_Allgather (&correction_local, 1, MPI_INT,
                          correction, 1, MPI_INT, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* BEGIN: wait for MPI send to complete */
  if (num_sends > 0) {
    /* wait for sends to complete */
    mpiret = MPI_Waitall (num_sends, send_requests, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    /* free send memory */
    P4EST_FREE (parent_send);
    P4EST_FREE (send_requests);
  }
  /* END: wait for MPI send to complete */

  /* correct partition */
  current_proc =
    p4est_next_nonempty_process (0, num_procs, num_quadrants_in_proc);
  next_proc =
    p4est_next_nonempty_process (current_proc + 1, num_procs,
                                 num_quadrants_in_proc);
  num_moved_quadrants = 0;
  while (current_proc < num_procs) {
    /* loop over all non empty processes */
    /* compute correct partition for process `current_proc` */
    if (0 < current_proc && current_proc < num_procs) {
      /* if any process but first */
      num_quadrants_in_proc[current_proc] += correction[current_proc];
      /* input is just a locidx, but the result is gloidx so we cast cleanly */
      num_moved_quadrants += P4EST_GLOIDX_ABS (correction[current_proc]);
    }
    if (next_proc < num_procs) {
      /* if first process or next process is feasible */
      num_quadrants_in_proc[current_proc] -= correction[next_proc];
    }

    /* increase process ids */
    current_proc = next_proc;
    next_proc =
      p4est_next_nonempty_process (next_proc + 1, num_procs,
                                   num_quadrants_in_proc);
  }

  /* free memory */
  P4EST_FREE (correction);

  /* return absolute number of moved quadrants */
  P4EST_FUNC_SHOT (p4est, &snap);
  return num_moved_quadrants;
#else
  P4EST_ASSERT (num_quadrants_in_proc[0] == p4est->local_num_quadrants);

  return 0;
#endif
}

unsigned
p4est_checksum (p4est_t * p4est)
{
#ifdef P4EST_HAVE_ZLIB
  uLong               treecrc, crc;
  size_t              scount, ssum;
  p4est_topidx_t      nt;
  p4est_tree_t       *tree;
  sc_array_t          checkarray;

  P4EST_ASSERT (p4est_is_valid (p4est));

  sc_array_init (&checkarray, 4);
  crc = adler32 (0, Z_NULL, 0);
  ssum = 0;
  for (nt = p4est->first_local_tree; nt <= p4est->last_local_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    treecrc =
      (uLong) p4est_quadrant_checksum (&tree->quadrants, &checkarray, 0);
    scount = 4 * checkarray.elem_count;
    ssum += scount;
    crc = adler32_combine (crc, treecrc, (z_off_t) scount);
  }
  sc_array_reset (&checkarray);
  P4EST_ASSERT ((p4est_locidx_t) ssum ==
                p4est->local_num_quadrants * 4 * (P4EST_DIM + 1));

  return p4est_comm_checksum (p4est, (unsigned) crc, ssum);
#else
  sc_abort_collective
    ("Configure did not find a recent enough zlib.  Abort.\n");

  return 0;
#endif /* !P4EST_HAVE_ZLIB */
}

void
p4est_save (const char *filename, p4est_t * p4est, int save_data)
{
  p4est_save_ext (filename, p4est, save_data, 1);
}

void
p4est_save_ext (const char *filename, p4est_t * p4est,
                int save_data, int save_partition)
{
  const int           headc = 6;
  const int           align = 32;
#ifdef P4EST_ENABLE_MPI
  int                 mpiret;
#ifndef P4EST_MPIIO_WRITE
  MPI_Status          mpistatus;
#endif
#endif
  int                 retval;
  int                 num_procs, save_num_procs, rank;
  int                 i;
  long                fpos = -1, foffset;
  size_t              data_size, qbuf_size, comb_size, head_count;
  size_t              zz, zcount;
  uint64_t           *u64a;
  FILE               *file;
#ifdef P4EST_MPIIO_WRITE
  MPI_File            mpifile;
  MPI_Offset          mpipos;
  MPI_Offset          mpithis;
#else
  long                fthis;
#endif
  p4est_topidx_t      jt, num_trees;
  p4est_gloidx_t     *pertree;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  char               *lbuf, *bp;
  p4est_qcoord_t     *qpos;
  sc_array_t         *tquadrants;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_save %s\n", filename);
  p4est_log_indent_push ();

  P4EST_ASSERT (p4est_connectivity_is_valid (p4est->connectivity));
  P4EST_ASSERT (p4est_is_valid (p4est));

  /* when data is not saved the size is set to zero */
  data_size = save_data ? p4est->data_size : 0;

  /* zero data size is effectively not saved */
  if (data_size == 0) {
    save_data = 0;
  }

  /* other parameters */
  num_trees = p4est->connectivity->num_trees;
  num_procs = p4est->mpisize;
  save_num_procs = save_partition ? num_procs : 1;
  head_count = (size_t) (headc + save_num_procs) + (size_t) num_trees;
  rank = p4est->mpirank;
  qbuf_size = (P4EST_DIM + 1) * sizeof (p4est_qcoord_t);
  comb_size = qbuf_size + data_size;
  pertree = P4EST_ALLOC (p4est_gloidx_t, num_trees + 1);
  p4est_comm_count_pertree (p4est, pertree);

  if (rank == 0) {
    p4est_connectivity_save (filename, p4est->connectivity);

    /* open file after writing connectivity to it */
    file = fopen (filename, "ab");
    SC_CHECK_ABORT (file != NULL, "file open");

    /* align the start of the header */
    fpos = ftell (file);
    SC_CHECK_ABORT (fpos > 0, "first file tell");
    while (fpos % align != 0) {
      retval = fputc ('\0', file);
      SC_CHECK_ABORT (retval == 0, "first file align");
      ++fpos;
    }

    /* write format and partition information */
    u64a = P4EST_ALLOC (uint64_t, head_count);
    u64a[0] = P4EST_ONDISK_FORMAT;
    u64a[1] = (uint64_t) sizeof (p4est_qcoord_t);
    u64a[2] = (uint64_t) sizeof (p4est_quadrant_t);
    u64a[3] = (uint64_t) data_size;
    u64a[4] = (uint64_t) save_data;
    u64a[5] = (uint64_t) save_num_procs;
    if (save_partition) {
      P4EST_ASSERT (save_num_procs == num_procs);
      for (i = 0; i < num_procs; ++i) {
        u64a[headc + i] = (uint64_t) p4est->global_first_quadrant[i + 1];
      }
    }
    else {
      P4EST_ASSERT (save_num_procs == 1);
      u64a[headc] = (uint64_t) p4est->global_first_quadrant[num_procs];
    }
    for (jt = 0; jt < num_trees; ++jt) {
      u64a[headc + save_num_procs + jt] = (uint64_t) pertree[jt + 1];
    }
    sc_fwrite (u64a, sizeof (uint64_t), head_count,
               file, "write header information");
    P4EST_FREE (u64a);
    fpos += head_count * sizeof (uint64_t);

    /* align the start of the quadrants */
    fpos = ftell (file);
    SC_CHECK_ABORT (fpos > 0, "second file tell");
    while (fpos % align != 0) {
      retval = fputc ('\0', file);
      SC_CHECK_ABORT (retval == 0, "second file align");
      ++fpos;
    }

#ifdef P4EST_MPIIO_WRITE
    /* We will close the sequential access to the file */
    /* best attempt to flush file to disk */
    retval = fflush (file);
    SC_CHECK_ABORT (retval == 0, "file flush");
#ifdef P4EST_HAVE_FSYNC
    retval = fsync (fileno (file));
    SC_CHECK_ABORT (retval == 0, "file fsync");
#endif
    retval = fclose (file);
    SC_CHECK_ABORT (retval == 0, "file close");
    file = NULL;
#else
    /* file is still open for sequential write mode */
#endif
  }
  else {
    file = NULL;
  }
  P4EST_FREE (pertree);

#ifndef P4EST_MPIIO_WRITE
  if (rank > 0) {
    /* wait for sequential synchronization */
#ifdef P4EST_ENABLE_MPI
    mpiret = MPI_Recv (&fpos, 1, MPI_LONG, rank - 1, P4EST_COMM_SAVE,
                       p4est->mpicomm, &mpistatus);
    SC_CHECK_MPI (mpiret);
#endif

    /* open file after all previous processors have written to it */
    file = fopen (filename, "rb+");
    SC_CHECK_ABORT (file != NULL, "file open");
  }
#else
  /* Every core opens the file in append mode */
  mpiret = MPI_File_open (p4est->mpicomm, (char *) filename,
                          MPI_MODE_WRONLY | MPI_MODE_APPEND |
                          MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &mpifile);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_File_get_position (mpifile, &mpipos);
  SC_CHECK_MPI (mpiret);
#endif

  if (rank > 0) {
    /* seek to the beginning of this processor's storage */
    foffset = (long) (p4est->global_first_quadrant[rank] * comb_size);

#ifndef P4EST_MPIIO_WRITE
    fthis = fpos + foffset;
    retval = fseek (file, fthis, SEEK_SET);
    SC_CHECK_ABORT (retval == 0, "seek data");
#else
    mpithis = mpipos + (MPI_Offset) foffset;
    mpiret = MPI_File_seek (mpifile, mpithis, MPI_SEEK_SET);
    SC_CHECK_MPI (mpiret);
#endif
  }

  /* write quadrant coordinates and data interleaved */
  for (jt = p4est->first_local_tree; jt <= p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    tquadrants = &tree->quadrants;
    zcount = tquadrants->elem_count;

    /* storage that will be written for this tree */
    bp = lbuf = P4EST_ALLOC (char, comb_size * zcount);
    for (zz = 0; zz < zcount; ++zz) {
      qpos = (p4est_locidx_t *) bp;
      q = p4est_quadrant_array_index (tquadrants, zz);
      *qpos++ = q->x;
      *qpos++ = q->y;
#ifdef P4_TO_P8
      *qpos++ = q->z;
#endif
      *qpos++ = (p4est_qcoord_t) q->level;
      if (save_data) {
        memcpy (qpos, q->p.user_data, data_size);
      }
      bp += comb_size;
    }
#ifndef P4EST_MPIIO_WRITE
    sc_fwrite (lbuf, comb_size, zcount, file, "write quadrants");
#else
    sc_mpi_write (mpifile, lbuf, comb_size * zcount, MPI_BYTE,
                  "write quadrants");
#endif
    P4EST_FREE (lbuf);
  }

#ifndef P4EST_MPIIO_WRITE
  /* best attempt to flush file to disk */
  retval = fflush (file);
  SC_CHECK_ABORT (retval == 0, "file flush");
#ifdef P4EST_HAVE_FSYNC
  retval = fsync (fileno (file));
  SC_CHECK_ABORT (retval == 0, "file fsync");
#endif
  retval = fclose (file);
  SC_CHECK_ABORT (retval == 0, "file close");
  file = NULL;

  /* initiate sequential synchronization */
#ifdef P4EST_ENABLE_MPI
  if (rank < num_procs - 1) {
    mpiret = MPI_Send (&fpos, 1, MPI_LONG, rank + 1, P4EST_COMM_SAVE,
                       p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
#endif
#else
  mpiret = MPI_File_close (&mpifile);
  SC_CHECK_MPI (mpiret);
#endif

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTION ("Done " P4EST_STRING "_save\n");
}

p4est_t            *
p4est_load (const char *filename, sc_MPI_Comm mpicomm, size_t data_size,
            int load_data, void *user_pointer,
            p4est_connectivity_t ** connectivity)
{
  return p4est_load_ext (filename, mpicomm, data_size, load_data,
                         0, 0, user_pointer, connectivity);
}

p4est_t            *
p4est_load_ext (const char *filename, sc_MPI_Comm mpicomm, size_t data_size,
                int load_data, int autopartition, int broadcasthead,
                void *user_pointer, p4est_connectivity_t ** connectivity)
{
  int                 retval;
  p4est_t            *p4est;
  sc_io_source_t     *src;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_load %s\n", filename);
  p4est_log_indent_push ();

  /* open file on all processors */

  src = sc_io_source_new (SC_IO_TYPE_FILENAME, SC_IO_ENCODE_NONE, filename);
  SC_CHECK_ABORT (src != NULL, "file source: possibly file not found");

  p4est = p4est_source_ext (src, mpicomm, data_size, load_data, autopartition,
                            broadcasthead, user_pointer, connectivity);

  retval = sc_io_source_destroy (src);
  SC_CHECK_ABORT (!retval, "source destroy");

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done " P4EST_STRING "_load with %lld total quadrants\n",
     (long long) p4est->global_num_quadrants);

  return p4est;
}

p4est_t            *
p4est_source_ext (sc_io_source_t * src, sc_MPI_Comm mpicomm, size_t data_size,
                  int load_data, int autopartition, int broadcasthead,
                  void *user_pointer, p4est_connectivity_t ** connectivity)
{
  const int           headc = 6;
  const int           align = 32;
  int                 root = 0;
  int                 retval;
  int                 mpiret;
  int                 num_procs, rank;
  int                 save_num_procs;
  int                 save_data;
  int                 i;
  uint64_t           *u64a, u64int;
  size_t              conn_bytes, file_offset;
  size_t              save_data_size;
  size_t              qbuf_size, comb_size, head_count;
  size_t              zz, zcount, zpadding;
  p4est_topidx_t      jt, num_trees;
  p4est_gloidx_t     *gfq;
  p4est_gloidx_t     *pertree;
  p4est_qcoord_t     *qap;
  p4est_connectivity_t *conn;
  p4est_t            *p4est;
  sc_array_t         *qarr, *darr;
  char               *dap, *lbuf;

  /* set some parameters */
  P4EST_ASSERT (src->bytes_out == 0);
  P4EST_ASSERT (connectivity != NULL);
  if (data_size == 0) {
    load_data = 0;
  }
  qbuf_size = (P4EST_DIM + 1) * sizeof (p4est_qcoord_t);

  /* retrieve MPI information */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* the first part of the header determines further offsets */
  save_data_size = (size_t) ULONG_MAX;
  save_num_procs = -1;
  conn = NULL;
  conn_bytes = 0;
  u64a = P4EST_ALLOC (uint64_t, headc + 1);
  if (!broadcasthead || rank == root) {

    /* read the forest connectivity */
    conn = p4est_connectivity_source (src);
    SC_CHECK_ABORT (conn != NULL, "connectivity source");
    zcount = src->bytes_out;
    zpadding = (align - zcount % align) % align;
    retval = sc_io_source_read (src, NULL, zpadding, NULL);
    SC_CHECK_ABORT (!retval, "source padding");
    conn_bytes = src->bytes_out;

    /* read format and some basic partition parameters */
    retval = sc_io_source_read (src, u64a, sizeof (uint64_t) * (size_t) headc,
                                NULL);
    SC_CHECK_ABORT (!retval, "read format");
    SC_CHECK_ABORT (u64a[0] == P4EST_ONDISK_FORMAT, "invalid format");
    SC_CHECK_ABORT (u64a[1] == (uint64_t) sizeof (p4est_qcoord_t),
                    "invalid coordinate size");
    SC_CHECK_ABORT (u64a[2] == (uint64_t) sizeof (p4est_quadrant_t),
                    "invalid quadrant size");
    save_data_size = (size_t) u64a[3];
    save_data = (int) u64a[4];
    if (load_data) {
      SC_CHECK_ABORT (save_data_size == data_size, "invalid data size");
      SC_CHECK_ABORT (save_data, "quadrant data not saved");
    }
    save_num_procs = (int) u64a[5];
    SC_CHECK_ABORT (autopartition || num_procs == save_num_procs,
                    "num procs mismatch");

    /* piggy-back the bytes for the connectivity onto first message */
    u64a[headc + 0] = (uint64_t) conn_bytes;
  }
  if (broadcasthead) {

    /* broadcast connectivity and first part of header */
    conn = p4est_connectivity_bcast (conn, root, mpicomm);
    mpiret = sc_MPI_Bcast (u64a, headc + 1, sc_MPI_LONG_LONG_INT,
                           root, mpicomm);
    SC_CHECK_MPI (mpiret);
    if (rank != root) {

      /* make sure the rest of the processes has the information */
      SC_CHECK_ABORT (u64a[0] == P4EST_ONDISK_FORMAT, "invalid format");
      save_data_size = (size_t) u64a[3];
      save_data = (int) u64a[4];
      save_num_procs = (int) u64a[5];
      conn_bytes = (size_t) u64a[headc + 0];
    }
  }
  P4EST_ASSERT (save_num_procs >= 0);
  P4EST_ASSERT (save_data_size != (size_t) ULONG_MAX);
  *connectivity = conn;
  comb_size = qbuf_size + save_data_size;
  file_offset = conn_bytes + headc * sizeof (uint64_t);

  /* create partition data */
  gfq = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  gfq[0] = 0;
  if (!broadcasthead || rank == root) {
    if (!autopartition) {
      P4EST_ASSERT (num_procs == save_num_procs);
      u64a = P4EST_REALLOC (u64a, uint64_t, num_procs);
      sc_io_source_read (src, u64a, sizeof (uint64_t) * (size_t) num_procs,
                         NULL);
      SC_CHECK_ABORT (!retval, "read quadrant partition");
      for (i = 0; i < num_procs; ++i) {
        gfq[i + 1] = (p4est_gloidx_t) u64a[i];
      }
    }
    else {
      /* ignore saved partition and compute a new uniform one */
      retval = sc_io_source_read
        (src, NULL, (long) ((save_num_procs - 1) * sizeof (uint64_t)), NULL);
      SC_CHECK_ABORT (!retval, "seek over ignored partition");
      retval = sc_io_source_read (src, &u64int, sizeof (uint64_t), NULL);
      SC_CHECK_ABORT (!retval, "read quadrant count");
      for (i = 1; i <= num_procs; ++i) {
        gfq[i] = p4est_partition_cut_uint64 (u64int, i, num_procs);
      }
    }
  }
  if (broadcasthead) {
    mpiret = sc_MPI_Bcast (gfq + 1, num_procs, P4EST_MPI_GLOIDX,
                           root, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  zcount = (size_t) (gfq[rank + 1] - gfq[rank]);
  file_offset += save_num_procs * sizeof (uint64_t);

  /* read pertree data */
  num_trees = conn->num_trees;
  pertree = P4EST_ALLOC (p4est_gloidx_t, num_trees + 1);
  pertree[0] = 0;
  if (!broadcasthead || rank == root) {
    u64a = P4EST_REALLOC (u64a, uint64_t, num_trees);
    retval = sc_io_source_read (src, u64a, sizeof (uint64_t) * (size_t)
                                num_trees, NULL);
    SC_CHECK_ABORT (!retval, "read pertree information");
    for (jt = 0; jt < num_trees; ++jt) {
      pertree[jt + 1] = (p4est_gloidx_t) u64a[jt];
    }
    SC_CHECK_ABORT (gfq[num_procs] == pertree[num_trees], "pertree mismatch");
  }
  if (broadcasthead) {
    mpiret = sc_MPI_Bcast (pertree + 1, num_trees, P4EST_MPI_GLOIDX,
                           root, mpicomm);
    SC_CHECK_MPI (mpiret);
  }
  P4EST_FREE (u64a);
  file_offset += num_trees * sizeof (uint64_t);

  /* seek to the beginning of this processor's storage */
  if (!broadcasthead || rank == root) {
    P4EST_ASSERT (file_offset == src->bytes_out);
    file_offset = 0;
  }
  head_count = (size_t) (headc + save_num_procs) + (size_t) num_trees;
  zpadding = (align - (head_count * sizeof (uint64_t)) % align) % align;
  if (zpadding > 0 || rank > 0) {
    retval = sc_io_source_read (src, NULL, (long)
                                (file_offset + zpadding +
                                 gfq[rank] * comb_size), NULL);
    SC_CHECK_ABORT (!retval, "seek data");
  }

  /* read quadrant coordinates and data interleaved */
  qarr =
    sc_array_new_size (sizeof (p4est_qcoord_t), (P4EST_DIM + 1) * zcount);
  qap = (p4est_qcoord_t *) qarr->array;
  darr = NULL;
  dap = NULL;
  lbuf = NULL;
  if (load_data) {
    P4EST_ASSERT (data_size == save_data_size && data_size > 0);
    darr = sc_array_new_size (data_size, zcount);
    dap = darr->array;
    lbuf = P4EST_ALLOC (char, comb_size);
  }
  for (zz = 0; zz < zcount; ++zz) {
    if (load_data) {
      retval = sc_io_source_read (src, lbuf, comb_size, NULL);
      SC_CHECK_ABORT (!retval, "read quadrant with data");
      memcpy (qap, lbuf, qbuf_size);
      memcpy (dap, lbuf + qbuf_size, data_size);
    }
    else {
      retval = sc_io_source_read (src, qap, qbuf_size, NULL);
      SC_CHECK_ABORT (!retval, "read quadrant with data");
      if (save_data_size > 0) {
        retval = sc_io_source_read (src, NULL, save_data_size, NULL);
        SC_CHECK_ABORT (!retval, "seek over data");
      }
    }
    qap += P4EST_DIM + 1;
    dap += data_size;
  }
  P4EST_FREE (lbuf);

  /* seek every process to the end of the source (in case there is data
   * appended to the end of this source) */
  if (gfq[num_procs] > gfq[rank + 1]) {
    retval = sc_io_source_read
      (src, NULL, (long) (gfq[num_procs] - gfq[rank + 1]) * comb_size, NULL);
    SC_CHECK_ABORT (!retval, "seek to end of data");
  }

  /* create p4est from accumulated information */
  p4est = p4est_inflate (mpicomm, conn, gfq, pertree,
                         qarr, darr, user_pointer);
  sc_array_destroy (qarr);
  if (darr != NULL) {
    sc_array_destroy (darr);
  }
  P4EST_FREE (pertree);
  P4EST_FREE (gfq);

  /* assert that we loaded a valid forest and return */
  SC_CHECK_ABORT (p4est_is_valid (p4est), "invalid forest");

  return p4est;
}
