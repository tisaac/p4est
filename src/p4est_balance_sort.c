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
#include <p8est_balance_obj.h>
#include <p8est_extended.h>
#include <p8est_bits.h>
#include <p8est_search.h>
#include <p8est_communication.h>
#include <p8est_algorithms.h>
#include <p8est_balance_seeds.h>
#else
#include <p4est_balance_obj.h>
#include <p4est_extended.h>
#include <p4est_bits.h>
#include <p4est_search.h>
#include <p4est_communication.h>
#include <p4est_algorithms.h>
#include <p4est_balance_seeds.h>
#endif
#include "p4est_balance_obj_impl.h"

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
static int
p4est_tree_get_conn_info (p4est_balance_obj_t * bobj,
                          p4est_t * p4est, p4est_topidx_t t,
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
  sc_flopinfo_t       snap;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);

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
  P4EST_BAL_FUNC_SHOT (bobj, &snap);
  return 0;
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

static void
p4est_tree_ghost_neigh_add_r (p4est_t * p4est, p4est_topidx_t from_tree,
                              int insul, int insul_max,
                              p4est_topidx_t which_tree, int nlevel,
                              int32_t touch, int p_fi, int p_la,
                              p4est_quadrant_t * fi, p4est_quadrant_t * la,
                              int *proc_hash, sc_array_t * procs)
{
  int32_t             rb;
  p4est_quadrant_t    fi_stack[64];
  p4est_quadrant_t    la_stack[64];
  int                 p_fi_stack[64];
  int                 p_la_stack[64];
  int                 stack_pointer = 0;
  int                 rank = p4est->mpirank;

  fi_stack[0] = *fi;
  la_stack[0] = *la;
  p_fi_stack[0] = p_fi;
  p_la_stack[0] = p_la;

  for (;;) {
    P4EST_ASSERT (stack_pointer < 64);
    rb =
      p4est_find_range_boundaries (&fi_stack[stack_pointer],
                                   &la_stack[stack_pointer], nlevel, NULL,
#ifdef P4_TO_P8
                                   NULL,
#endif
                                   NULL);
    if (!(rb & touch)) {
      /* move down the stack */
      while (stack_pointer > 0
             && p_la_stack[stack_pointer] == p_la_stack[stack_pointer - 1]) {
        stack_pointer--;
      }
      if (!stack_pointer) {
        return;
      }
      p_fi_stack[stack_pointer] = p_fi_stack[stack_pointer - 1];
      p_la_stack[stack_pointer] = p_la_stack[stack_pointer - 1];
      fi_stack[stack_pointer] = fi_stack[stack_pointer - 1];
      la_stack[stack_pointer] = la_stack[stack_pointer - 1];
      continue;
    }
    if (p_fi_stack[stack_pointer] == p_la_stack[stack_pointer]) {       /* (uniquely) add this tree's portion of the process range to the neighborhood */
      int                 p = p_fi_stack[stack_pointer];
      int                 idx = proc_hash[p];

      if (p != rank && !idx) {
        /* process has not been seen before, add it */
        *((int *) sc_array_push (procs)) = p;
        proc_hash[p] = (int) procs->elem_count;
      }
      while (stack_pointer > 0
             && p_la_stack[stack_pointer] == p_la_stack[stack_pointer - 1]) {
        stack_pointer--;
      }
      if (!stack_pointer) {
        return;
      }
      p_fi_stack[stack_pointer] = p_fi_stack[stack_pointer - 1];
      p_la_stack[stack_pointer] = p_la_stack[stack_pointer - 1];
      fi_stack[stack_pointer] = fi_stack[stack_pointer - 1];
      la_stack[stack_pointer] = la_stack[stack_pointer - 1];
      continue;
    }
    /* recurse only if this range of processes is adjacent */
    {
      int                 p = p_fi_stack[stack_pointer];
      int                 q = p_la_stack[stack_pointer];
      int                 num_procs = q + 1 - p;
      int                 mid = p + num_procs / 2;
      uint64_t            id;
      p4est_quadrant_t   *lalo;

      P4EST_ASSERT (stack_pointer < 63);
      p_fi_stack[stack_pointer + 1] = p;
      p_la_stack[stack_pointer + 1] = mid - 1;
      fi_stack[stack_pointer + 1] = fi_stack[stack_pointer];

      lalo = &p4est->global_first_position[mid];
      id = p4est_quadrant_linear_id (lalo, P4EST_QMAXLEVEL);
      p4est_quadrant_set_morton (&la_stack[stack_pointer + 1],
                                 P4EST_QMAXLEVEL, id - 1);
      p_fi_stack[stack_pointer] = mid;
      fi_stack[stack_pointer] = *lalo;

      stack_pointer++;
    }
  }
}

static void
p4est_tree_ghost_neighborhood_add (p4est_t * p4est, int from_tree, int insul,
                                   int insul_max, p4est_topidx_t which_tree,
                                   p4est_quadrant_t * q, p4est_quadrant_t * n,
                                   int *proc_hash, sc_array_t * procs)
{
  p4est_qcoord_t      diff[3];
  int                 ninsul;
  int                 hshift;
  int                 i;
  int                 l = q->level;
  int32_t             touch_table[P4EST_INSUL] = {
#ifndef P4_TO_P8
    0x10, 0x04, 0x20,
    0x01, -1, 0x02,
    0x40, 0x08, 0x80,
#else
    0x0040000, 0x0000040, 0x0080000,
    0x0000400, 0x0000010, 0x0000800,
    0x0100000, 0x0000080, 0x0200000,

    0x0004000, 0x0000004, 0x0008000,
    0x0000001, -1, 0x0000002,
    0x0010000, 0x0000008, 0x0020000,

    0x0400000, 0x0000100, 0x0800000,
    0x0001000, 0x0000020, 0x0002000,
    0x1000000, 0x0000200, 0x2000000,
#endif
  };
  int32_t             touch;
  p4est_quadrant_t    fi, la;
  int                 rank = p4est->mpirank;
  int                 p_fi, p_la;

  P4EST_ASSERT (n->level == l);
  P4EST_ASSERT (p4est_quadrant_is_valid (n));
  P4EST_ASSERT (p4est_quadrant_is_extended (q));

  hshift = P4EST_MAXLEVEL - l;
  diff[0] = (q->x - n->x) >> hshift;
  diff[1] = (q->y - n->y) >> hshift;
#ifdef P4_TO_P8
  diff[2] = (q->z - n->z) >> hshift;
#else
  diff[2] = 0;
#endif
  ninsul = 0;
  for (i = P4EST_DIM - 1, ninsul = 0; i >= 0; i--) {
    ninsul *= 3;
    P4EST_ASSERT (-1 <= diff[i] && diff[i] <= 1);
    ninsul += diff[i] + 1;
  }
  P4EST_ASSERT (ninsul != P4EST_INSUL / 2);
  touch = touch_table[ninsul];

  p4est_quadrant_first_descendant (n, &fi, P4EST_QMAXLEVEL);
  p4est_quadrant_last_descendant (n, &la, P4EST_QMAXLEVEL);
  p_fi = p4est_comm_find_owner (p4est, which_tree, &fi, rank);
  p_la = p4est_comm_find_owner (p4est, which_tree, &la, p_fi);
  if (p_fi == rank && p_la == rank) {
    return;
  }
  p4est_tree_ghost_neigh_add_r (p4est, from_tree, insul, insul_max,
                                which_tree, n->level, touch, p_fi, p_la, &fi,
                                &la, proc_hash, procs);
}

/* loop over the insuation neighbors of \a q , looking for neighboring
 * processes */
static void
p4est_tree_ghost_neighborhood_insert (p4est_t * p4est, p4est_topidx_t t,
                                      p4est_topidx_t flt, p4est_topidx_t llt,
                                      p4est_quadrant_t * q,
                                      p4est_tree_neigh_info_t * info,
                                      int *proc_hash,
                                      p4est_quadrant_t * sortquad,
                                      sc_array_t * procs)
{
  p4est_qcoord_t      h = P4EST_QUADRANT_LEN (q->level);
  int                 i, j, k, l, insulmax;
  int                 insort0 = 0;
  int                 insort1 = 0;
#ifdef P4_TO_P8
  int                 kstart = 0, kend = 3;
#else
  int                 kstart = 1, kend = 2;
#endif

  insulmax = info->offset[P4EST_INSUL];
  insort0 = (t == flt && p4est_quadrant_is_ancestor (&sortquad[0], q));
  insort1 = (t == llt && p4est_quadrant_is_ancestor (&sortquad[1], q));
  for (k = kstart, l = 0; k < kend; k++) {
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++, l++) {
        p4est_quadrant_t    n;
        p4est_quadrant_t    nn, nq;
        int                 rootinsul;
        size_t              s;

        if (l == P4EST_INSUL / 2) {
          continue;
        }

        n = *q;
        n.x += (i - 1) * h;
        n.y += (j - 1) * h;
#ifdef P4_TO_P8
        n.z += (k - 1) * h;
#endif

        rootinsul = p4est_root_insul (&n);
        if (rootinsul == P4EST_INSUL / 2) {
          if (insort0 && p4est_quadrant_is_ancestor (&sortquad[0], &n)) {
            /* do not search for ghost neighbors inside the region
             * where communication has been handled by the sort phase */
            continue;
          }
          if (insort1 && p4est_quadrant_is_ancestor (&sortquad[1], &n)) {
            /* do not search for ghost neighbors inside the region
             * where communication has been handled by the sort phase */
            continue;
          }
          p4est_tree_ghost_neighborhood_add (p4est, t,
                                             info->offset[rootinsul],
                                             insulmax, t, q, &n, proc_hash,
                                             procs);
        }
        else {
          size_t              offset = info->offset[rootinsul];
          size_t              stop = info->offset[rootinsul + 1];

          for (s = offset; s < stop; s++) {
            p4est_tree_neigh_t *tn =
              (p4est_tree_neigh_t *) sc_array_index (&info->tnarray, s);

            p4est_quadrant_utransform (q, &nq, &(tn->u[0]), 0);
            p4est_quadrant_utransform (&n, &nn, &(tn->u[0]), 0);
            p4est_tree_ghost_neighborhood_add (p4est, t, s, insulmax, tn->nt,
                                               &nq, &nn, proc_hash, procs);
          }
        }
      }
    }
  }
}

/* Compute the ghost neighbors of the process (ignoring neighbors who are only
 * neighbors that are "close" to the start / end of this process's range,
 * i.e. within the minlevel ancestor of first_desc / last_desc */
static void
p4est_tree_ghost_neighborhood (p4est_balance_obj_t * bobj, p4est_t * p4est,
                               p4est_topidx_t t, p4est_topidx_t flt,
                               p4est_topidx_t llt, int *proc_hash,
                               sc_array_t * procs, int *minlevel,
                               p4est_tree_neigh_info_t * info)
{
  p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
  p4est_quadrant_t    f = tree->first_desc;
  p4est_quadrant_t    l = tree->last_desc;
  p4est_quadrant_t    a, af, al, fstop, temp;
  p4est_quadrant_t    sortquad[2];
  int                 alevel, fid, lid, idxor, dir, fstopid;
  sc_flopinfo_t       snap;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);

  p4est_nearest_common_ancestor (&f, &l, &a);

  alevel = a.level;

  p4est_quadrant_first_descendant (&a, &af, P4EST_QMAXLEVEL);
  p4est_quadrant_last_descendant (&a, &al, P4EST_QMAXLEVEL);

  if (t == flt) {
    p4est_quadrant_ancestor (&f, minlevel[0], &sortquad[0]);
  }
  if (t == llt) {
    p4est_quadrant_ancestor (&l, minlevel[1], &sortquad[1]);
  }
  if (p4est_quadrant_is_equal (&f, &af) && p4est_quadrant_is_equal (&l, &al)) {
    p4est_tree_ghost_neighborhood_insert (p4est, t, flt, llt, &a, info,
                                          proc_hash, sortquad, procs);
    return;
  }
  /* loop over the implicit coarsest tree for this process's range,
   * inserting ghost neighbors based on the insulation layer neighbors of the
   * quads in the loop */
  fid = p4est_quadrant_ancestor_id (&f, alevel + 1);
  lid = p4est_quadrant_ancestor_id (&l, alevel + 1);
  P4EST_ASSERT (lid > fid);
  idxor = lid ^ fid;
  dir = SC_LOG2_8 (idxor);
  P4EST_ASSERT (dir >= 0);
  fstopid = (lid >> dir) << dir;
  p4est_quadrant_child (&a, &fstop, fstopid);
  while (f.level > alevel + 1 && p4est_quadrant_child_id (&f) == 0) {
    p4est_quadrant_parent (&f, &f);
  }
  while (l.level > alevel + 1 &&
         p4est_quadrant_child_id (&l) == P4EST_CHILDREN - 1) {
    p4est_quadrant_parent (&l, &l);
  }
  P4EST_ASSERT (!p4est_quadrant_overlaps (&f, &fstop)
                && p4est_quadrant_compare (&f, &fstop) < 0);
  P4EST_ASSERT (p4est_quadrant_compare (&fstop, &l) <= 0);
  for (;;) {
    p4est_tree_ghost_neighborhood_insert (p4est, t, flt, llt, &f, info,
                                          proc_hash, sortquad, procs);
    fid = p4est_quadrant_child_id (&f);
    while (fid == P4EST_CHILDREN - 1) {
      p4est_quadrant_parent (&f, &f);
      fid = p4est_quadrant_child_id (&f);
    }
    p4est_quadrant_sibling (&f, &temp, fid + 1);
    f = temp;
    if (p4est_quadrant_is_equal (&f, &fstop)) {
      break;
    }
  }
  for (;;) {
    while (p4est_quadrant_is_ancestor (&fstop, &l)) {
      fstop.level++;
    }
    p4est_tree_ghost_neighborhood_insert (p4est, t, flt, llt, &fstop, info,
                                          proc_hash, sortquad, procs);
    if (p4est_quadrant_is_equal (&fstop, &l)) {
      break;
    }
    fstopid = p4est_quadrant_child_id (&fstop);
    P4EST_ASSERT (fstopid < P4EST_CHILDREN - 1);
    p4est_quadrant_sibling (&fstop, &temp, fstopid + 1);
    fstop = temp;
  }
  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}

static int
p4est_balance_sort_min_insulation_level (p4est_balance_obj_t * bobj,
                                         p4est_t * p4est,
                                         p4est_tree_neigh_info_t * info,
                                         p4est_quadrant_t * q)
{
  int                 i, j, k, l, m, p;
#ifdef P4_TO_P8
  int                 kstart = 0, kend = 3;
#else
  int                 kstart = 1, kend = 2;
#endif
  p4est_quadrant_t    fd, ld;
  p4est_quadrant_t   *nd;

  if (bobj->use_root) {
    return 0;
  }
  for (l = 0; l < P4EST_QMAXLEVEL; l++) {
    p4est_quadrant_t    a, n;
    p4est_qcoord_t      h;

    p4est_quadrant_ancestor (q, l, &a);
    h = P4EST_QUADRANT_LEN (l);
    n.level = l;
    n.p.which_tree = q->p.which_tree;
    for (m = 0, k = kstart; k < kend; k++) {
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++, m++) {
          int                 rootinsul;

          n.x = a.x + (i - 1) * h;
          n.y = a.y + (j - 1) * h;
#ifdef P4_TO_P8
          n.z = a.z + (j - 1) * h;
#endif

          rootinsul = p4est_root_insul (&n);

          if (rootinsul == P4EST_INSUL / 2) {
            /* in root */
            p4est_quadrant_first_descendant (&n, &fd, P4EST_QMAXLEVEL);
            p4est_quadrant_last_descendant (&n, &ld, P4EST_QMAXLEVEL);
            ld.p.which_tree = n.p.which_tree;
            p =
              p4est_comm_find_owner (p4est, q->p.which_tree, &fd,
                                     p4est->mpirank);
            nd = &p4est->global_first_position[p + 1];
            if (p4est_quadrant_compare_piggy (&ld, nd) < 0) {
              /* n is entirely contained in one process */
              return l;
            }
          }
          else {
            size_t              s, s_start, s_end;

            s_start = info->offset[rootinsul];
            s_end = info->offset[rootinsul + 1];
            for (s = s_start; s < s_end; s++) {
              p4est_quadrant_t    nn;
              p4est_tree_neigh_t *tn =
                (p4est_tree_neigh_t *) sc_array_index (&(info->tnarray), s);

              p4est_quadrant_utransform (&n, &nn, &(tn->u[0]), 0);
              p4est_quadrant_first_descendant (&nn, &fd, P4EST_QMAXLEVEL);
              p4est_quadrant_last_descendant (&nn, &ld, P4EST_QMAXLEVEL);
              ld.p.which_tree = tn->nt;
              p = p4est_comm_find_owner (p4est, tn->nt, &fd, p4est->mpirank);
              nd = &p4est->global_first_position[p + 1];
              if (p4est_quadrant_compare_piggy (&ld, nd) < 0) {
                /* nn is entirely contained in one process */
                return l;
              }
            }
          }
        }
      }
    }
  }
  SC_ABORT_NOT_REACHED ();
  return -1;
}

static void
sc_array_crop (sc_array_t * array, size_t a_start, size_t a_end)
{
  P4EST_ASSERT (array);
  P4EST_ASSERT (a_end >= a_start);
  P4EST_ASSERT (a_end <= array->elem_count);
  if (a_end == a_start) {
    array->elem_count = 0;
    return;
  }
  if (a_start) {
    memmove (array->array, sc_array_index (array, a_start),
             (a_end - a_start) * array->elem_size);
  }
  array->elem_count = (a_end - a_start);
}

static void
sc_array_push_copy (sc_array_t * array, size_t a_start, size_t a_end,
                    sc_array_t * result)
{
  P4EST_ASSERT (result && array != result);
  P4EST_ASSERT (result->elem_size == array->elem_size);
  P4EST_ASSERT (a_end >= a_start);

  if (a_end > a_start) {
    char               *dest = sc_array_push_count (result, a_end - a_start);
    memcpy (dest, sc_array_index (array, a_start),
            (a_end - a_start) * array->elem_size);
  }
}

static void
sc_array_push_cut (sc_array_t * array, size_t a_start, size_t a_end,
                   sc_array_t * result)
{
  P4EST_ASSERT (result && array != result);
  P4EST_ASSERT (result->elem_size == array->elem_size);
  P4EST_ASSERT (a_end >= a_start);

  if (a_end > a_start) {
    char               *dest = sc_array_push_count (result, a_end - a_start);

    memcpy (dest, sc_array_index (array, a_start),
            (a_end - a_start) * array->elem_size);
    if (a_end < array->elem_count) {
      memmove (sc_array_index (array, a_start), sc_array_index (array, a_end),
               (array->elem_count - a_end) * array->elem_size);
    }
    array->elem_count -= (a_end - a_start);
  }
}

enum
{
  P4EST_BALSORT_ISEND,
  P4EST_BALSORT_IPROBE,
  P4EST_BALSORT_IRECV,
  P4EST_BALSORT_MERGE,
  P4EST_BALSORT_WAIT,
  P4EST_BALSORT_SYNC,
};

static void
p4est_quadrant_array_merge_reduce (sc_array_t * A, sc_array_t * B)
{
  size_t              nnew = B->elem_count;
  size_t              nold = A->elem_count;
  p4est_quadrant_t   *a, *b, *c;
  p4est_quadrant_t    cp;
  size_t              i, j, k;

  if (!nnew) {
    return;
  }
  P4EST_ASSERT (nold == 0 || p4est_quadrant_array_is_reduced (A));
  P4EST_ASSERT (p4est_quadrant_array_is_reduced (B));
  (void) sc_array_push_count (A, nnew);
  if (!nold) {
    sc_array_copy (A, B);
    return;
  }
  a = (p4est_quadrant_t *) sc_array_index (A, nnew);
  b = (p4est_quadrant_t *) B->array;
  c = (p4est_quadrant_t *) A->array;
  memmove (a, c, nold * A->elem_size);
  for (i = 0, j = 0, k = 0; i + j < nold + nnew;) {
    p4est_quadrant_t   *q, qp;
    int                 comp;

    if (i == nold) {
      comp = -1;
    }
    else if (j == nnew) {
      comp = 1;
    }
    else {
      comp = p4est_quadrant_compare (b, a);
    }
    if (comp < 0) {
      q = b++;
      j++;
    }
    else if (comp > 0) {
      q = a++;
      i++;
    }
    else {
      q = a++;
      b++;
      i++;
      j++;
    }
    if (!q->level) {
      continue;
    }
    if (!k) {
      p4est_quadrant_parent (q, &cp);
      *c = *q;
      k++;
    }
    else {
      P4EST_ASSERT (p4est_quadrant_compare (c, q) < 0);
      p4est_quadrant_parent (q, &qp);
      if (p4est_quadrant_overlaps (&cp, &qp)) {
        if (q->level > c->level) {
          *c = *q;
          cp = qp;
        }
      }
      else {
        *(++c) = *q;
        cp = qp;
        k++;
      }
    }
  }
  A->elem_count = k;
  P4EST_ASSERT (p4est_quadrant_array_is_reduced (A));
}

/* divide up a sorted list of quadrants from a given tree between the
 * neighbor processes */
static void
p4est_balance_sort_divide (p4est_t * p4est, int n_neigh, const int *neigh_procs,
                           sc_array_t ** send_arrays,
                           p4est_topidx_t which_tree,
                           sc_array_t * tquads, int skip_self)
{
  int                 neigh;
  int                 limit = skip_self ? n_neigh : n_neigh + 1;

  for (neigh = 0; neigh < limit; neigh++) {
    int                 p = neigh < n_neigh ? neigh_procs[neigh] : p4est->mpirank;
    p4est_quadrant_t    fp, np, lb, ub;
    p4est_quadrant_t   *nq;
    ssize_t             sz;
    size_t              qstart, qend;
    size_t              nbuf, s;

    fp = p4est->global_first_position[p];
    np = p4est->global_first_position[p + 1];
    if (np.p.which_tree < which_tree || fp.p.which_tree > which_tree) {
      continue;
    }
    if (fp.p.which_tree == which_tree) {
      lb = fp;

      while (lb.level && !p4est_quadrant_child_id (&lb)) {
        lb.level--;
      }
      sz = p4est_find_lower_bound (tquads, &lb, 0);
      if (sz < 0) {
        sz = tquads->elem_count;
      }
      qstart = sz;
    }
    else {
      qstart = 0;
    }

    if (np.p.which_tree > which_tree) {
      qend = tquads->elem_count;
    }
    else {
      ub = np;

      while (ub.level && !p4est_quadrant_child_id (&ub)) {
        ub.level--;
      }
      sz = p4est_find_lower_bound (tquads, &ub, 0);
      if (sz < 0) {
        sz = tquads->elem_count;
      }
      qend = sz;
    }
    P4EST_ASSERT (qend >= qstart);
    if (qend > qstart) {
      sc_array_t *buf = send_arrays[neigh];

      nbuf = buf->elem_count;
      sc_array_push_copy (tquads, qstart, qend, buf);
      nq = p4est_quadrant_array_index (buf, nbuf);
      for (s = 0; s < qend - qstart; s++) {
        nq[s].pad8 = 0;
        nq[s].pad16 = 0;
        nq[s].p.user_data = NULL;
        nq[s].p.which_tree = which_tree;
      }
    }
  }
}

static void
p4est_balance_sort_compute_pattern_sort (p4est_balance_obj_t * bobj,
                                         p4est_t * p4est, int flt, int llt,
                                         int num_trees,
                                         p4est_tree_neigh_info_t * tinfo,
                                         int (*procrange)[2], int *minlevel,
                                         int *num_sort,
                                         p4est_quadrant_t * desc)
{
  p4est_topidx_t      t;
  sc_flopinfo_t       snap;
  int                 rank = p4est->mpirank;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);

  /* block balance local */
  for (t = flt; t <= llt; t++) {
    p4est_tree_neigh_info_t *info;

    /* Get the connectivity information for this tree */
    info = &tinfo[t - flt];
    p4est_tree_get_conn_info (bobj, p4est, t, info);
  }

  /* compute the processes I communicate with in the sort phase */
  if (llt >= flt) {
    int                 i;
    int                 empty;

    desc[0] = p4est->global_first_position[rank];
    desc[1] = p4est->global_first_position[rank + 1];
    empty = p4est_quadrant_is_equal_piggy (&desc[0], &desc[1]);
    if (desc[1].p.which_tree > llt) {
      p4est_qcoord_t      last =
        P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
      desc[1].p.which_tree--;

      desc[1].x = last;
      desc[1].y = last;
#ifdef P4_TO_P8
      desc[1].z = last;
#endif
    }
    else {
      uint64_t            id =
        p4est_quadrant_linear_id (&desc[1], P4EST_QMAXLEVEL);

      p4est_quadrant_set_morton (&desc[1], P4EST_QMAXLEVEL, id - 1);
    }

    for (i = 0; i < 2; i++) {
      p4est_quadrant_t    a, fd, ld;
      int                 t = i ? llt : flt;

      minlevel[i] =
        p4est_balance_sort_min_insulation_level (bobj, p4est, &tinfo[t - flt],
                                                 &desc[i]);
      p4est_quadrant_ancestor (&desc[i], minlevel[i], &a);
      p4est_quadrant_first_descendant (&a, &fd, P4EST_QMAXLEVEL);
      p4est_quadrant_last_descendant (&a, &ld, P4EST_QMAXLEVEL);
      procrange[i][0] = p4est_comm_find_owner (p4est, t, &fd, rank);
      procrange[i][1] = p4est_comm_find_owner (p4est, t, &ld, rank);
    }
    if (empty) {
      *num_sort = 1;
      if (procrange[0][0] != procrange[1][0]) {
        /* communicate with no one */
        procrange[0][0] = procrange[1][0] = -1;
        procrange[0][1] = procrange[1][1] = -2;
        *num_sort = 0;
      }
    }
    else {
      if (procrange[0][0] == procrange[1][0] &&
          procrange[0][1] == procrange[1][1]) {
        *num_sort = 1;
      }
      else {
        *num_sort = 2;
      }
      P4EST_ASSERT (procrange[0][0] <= rank);
      P4EST_ASSERT (rank <= procrange[1][1]);
      if (procrange[0][1] > procrange[1][0]) {
        P4EST_ASSERT (procrange[0][0] == procrange[1][0]);
        P4EST_ASSERT (procrange[0][1] == procrange[1][1]);
      }
      else {
        P4EST_ASSERT (procrange[0][1] == rank);
        P4EST_ASSERT (procrange[1][0] == rank);
      }
    }
  }
  else {
    /* communicate with no one */
    procrange[0][0] = procrange[1][0] = -1;
    procrange[0][1] = procrange[1][1] = -2;
    num_sort = 0;
  }
#ifdef P4EST_ENABLE_DEBUG
  {
    int                *comm_out = P4EST_ALLOC_ZERO (int, p4est->mpisize);
    int                *comm_in = P4EST_ALLOC (int, p4est->mpisize);
    int                 i, j, mpiret;

    for (i = 0; i < 2; i++) {
      for (j = procrange[i][0]; j <= procrange[i][1]; j++) {
        comm_out[j] = 1;
      }
    }
    mpiret = sc_MPI_Alltoall (comm_out, 1, sc_MPI_INT, comm_in, 1, sc_MPI_INT,
                              p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    for (i = 0; i < p4est->mpisize; i++) {
      P4EST_ASSERT (comm_out[i] == comm_in[i]);
    }
    P4EST_FREE (comm_out);
    P4EST_FREE (comm_in);
  }
#endif
  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}

static void
p4est_balance_sort_compute_neigh (p4est_balance_obj_t * bobj,
                                  p4est_t * p4est,
                                  int flt, int llt, int num_trees,
                                  p4est_tree_neigh_info_t * tinfo,
                                  int *minlevel,
                                  p4est_quadrant_t * desc,
                                  p4est_neigh_t **neigh_p)
{
  p4est_topidx_t      t;
  sc_flopinfo_t       snap;
  size_t              s;
  int                *proc_hash;
  sc_array_t         *procs;

  if (bobj->neigh) {
    *neigh_p = bobj->neigh;
    return;
  }

  P4EST_BAL_FUNC_SNAP (bobj, &snap);

  procs = sc_array_new_count (sizeof (int), P4EST_INSUL);
  sc_array_truncate (procs);
  proc_hash = P4EST_ALLOC_ZERO (int, p4est->mpisize);

  /* block balance local */
  for (t = flt; t <= llt; t++) {
    p4est_tree_neigh_info_t *info;

    /* Get the connectivity information for this tree */
    info = &tinfo[t - flt];

    /* find ghost neighbors of this tree */
    if (num_trees) {
      p4est_tree_ghost_neighborhood (bobj, p4est, t, flt, llt, proc_hash,
                                     procs, minlevel, info);
    }
  }
  sc_array_sort (procs, sc_int_compare);
  for (s = 0; s < procs->elem_count; s++) {
    int                 p = *((int *) sc_array_index (procs, s));

    proc_hash[p] = s + 1;
  }

  *neigh_p = p4est_neigh_new (p4est, (int) procs->elem_count, (int *) procs->array);
  P4EST_FREE (proc_hash);
  sc_array_destroy (procs);
  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}

static void
p4est_balance_sort_local (p4est_balance_obj_t * bobj,
                          p4est_t * p4est,
                          p4est_connect_type_t btype,
                          int flt, int llt, int num_trees,
                          int num_sort, int *minlevel,
                          p4est_quadrant_t * desc,
                          sc_array_t * tree_bufs, sc_array_t * sort_bufs)
{
  p4est_topidx_t      t;
  const int8_t       *pre_adapt_flags = NULL;
  p4est_locidx_t      all_incount;
  int                 bound;
  sc_flopinfo_t       snap;
  sc_mempool_t       *qpool;
  sc_mempool_t       *list_alloc;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);
  pre_adapt_flags = p4est_balance_obj_get_adapt_flags (bobj);

  switch (p4est_connect_type_int (btype)) {
  case 0:
    bound = 1;
    break;
  case 1:
    bound = P4EST_DIM + 1;
    break;
  case P4EST_DIM:
    bound = (1 << P4EST_DIM);
    break;
#ifdef P4_TO_P8
  case 2:
    bound = 7;
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
  }

  qpool = p4est->quadrant_pool;
  list_alloc = sc_mempool_new (sizeof (sc_link_t));

  /* block balance local */
  for (t = flt, all_incount = 0; t < flt + num_trees; t++) {
    const int8_t       *this_pre_adapt_flags = NULL;
    p4est_tree_t       *tree;
    sc_array_t         *tquadrants;
    size_t              treecount;
    p4est_quadrant_t    root;
    size_t              ancstart, ancend, selfstart, selfend;
    sc_array_t         *inquads = &tree_bufs[t - flt];

    tree = p4est_tree_array_index (p4est->trees, t);
    tquadrants = &tree->quadrants;

    /* initial log message for this tree */
    P4EST_VERBOSEF ("Into balance tree %lld with %llu\n", (long long) t,
                    (unsigned long long) tquadrants->elem_count);

    p4est_nearest_common_ancestor (&tree->first_desc, &tree->last_desc,
                                   &root);
    if (t == flt && t == llt) {
      if (minlevel[0] <= root.level) {
        P4EST_ASSERT (minlevel[1] == minlevel[0]);
        p4est_quadrant_ancestor (&tree->first_desc, minlevel[0], &root);
      }
      /* is it worth it to split the local balance when it consists of two
       * regions that do not influence each other ? */
#if 0
      else if (minlevel[0] > root.level + 1 && minlevel[1] > root.level + 1) {
        int                 fid, lid, level, lolevel, hilevel;
        p4est_quadrant_t    lo, hi, ub, lb, temp;
        p4est_qcoord_t      H, diff[P4EST_DIM];

        fid = p4est_quadrant_ancestor_id (&tree->first_desc, root.level + 1);
        lid = p4est_quadrant_ancestor_id (&tree->last_desc, root.level + 1);
        P4EST_ASSERT (lid == fid + 1);

        p4est_quadrant_child (&root, &lo, fid);
        p4est_quadrant_last_descendant (&lo, &ub, P4EST_QMAXLEVEL);
        p4est_nearest_common_ancestor (&tree->first_desc, &ub, &lo);
        lolevel = SC_MIN (lo.level, minlevel[0]);

        p4est_quadrant_child (&root, &hi, lid);
        p4est_quadrant_first_descendant (&hi, &lb, P4EST_QMAXLEVEL);
        p4est_nearest_common_ancestor (&tree->last_desc, &lb, &hi);
        hilevel = SC_MIN (hi.level, minlevel[1]);
        level = SC_MIN (lolevel, hilevel);

        p4est_quadrant_ancestor (&tree->first_desc, level, &lo);
        p4est_quadrant_ancestor (&tree->last_desc, level, &hi);

        H = P4EST_QUADRANT_LEN (level);
        diff[0] = hi.x - lo.x;
        diff[1] = hi.y - lo.y;
#ifdef P4_TO_P8
        diff[2] = hi.z - lo.z;
#endif
        if (diff[0] < -H || diff[0] > H || diff[1] < -H || diff[1] > H ||
#ifdef P4_TO_P8
            diff[2] < -H || diff[2] > H ||
#endif
            0) {
          P4EST_DEBUGF ("Split root: [%d, %d], [%d, %d], [%d, %d, %d], %d\n",
                        lolevel, hilevel, minlevel[0], minlevel[1], diff[0],
                        diff[1],
#ifdef P4_TO_P8
                        diff[2],
#else
                        0,
#endif
                        H);
        }
      }
#endif
    }
    else if (t == flt) {
      p4est_quadrant_ancestor (&tree->first_desc,
                               SC_MIN (root.level, minlevel[0]), &root);
    }
    else if (t == llt) {
      p4est_quadrant_ancestor (&tree->first_desc,
                               SC_MIN (root.level, minlevel[1]), &root);
    }

    if (pre_adapt_flags) {
      this_pre_adapt_flags = &pre_adapt_flags[all_incount];
    }
    all_incount += tquadrants->elem_count;

    /* local balance first pass */
    p4est_quadrant_array_reduce (tquadrants, inquads, this_pre_adapt_flags);
    if (pre_adapt_flags) {
      p4est_quadrant_array_insert_endpoints (inquads, &root,
                                             &tree->first_desc,
                                             &tree->last_desc);
    }
    p4est_balance_kernel (inquads, NULL, root.level, bound, qpool,
                          list_alloc, NULL, NULL, NULL, NULL, NULL);

    treecount = inquads->elem_count;
    P4EST_VERBOSEF ("Balance tree %lld A %llu (reduced)\n",
                    (long long) t, (unsigned long long) treecount);

    if (t > flt && t < llt) {
      continue;
    }

    if (t == flt) {
      ssize_t             sz;
      p4est_quadrant_t    fdesc, anc;

      fdesc = tree->first_desc;
      while (fdesc.level > 0 && !p4est_quadrant_child_id (&fdesc)) {
        fdesc.level--;
      }
      sz = p4est_find_lower_bound (inquads, &fdesc, 0);
      if (sz < 0) {
        sz = inquads->elem_count;
      }
      selfstart = sz;

      p4est_quadrant_ancestor (&desc[0], minlevel[0], &anc);

      sz = p4est_find_lower_bound (inquads, &anc, 0);
      if (sz < 0) {
        sz = inquads->elem_count;
      }
      ancstart = sz;

      P4EST_ASSERT (ancstart <= selfstart);
      sc_array_push_cut (inquads, ancstart, selfstart, &sort_bufs[0]);
      selfstart = ancstart;
    }
    else {
      selfstart = 0;
    }
    if (t == llt) {
      ssize_t             sz;
      p4est_quadrant_t    anc, ub;

      sz = p4est_find_higher_bound (inquads, &tree->last_desc, 0);
      if (sz < 0) {
        sz = 0;
      }
      else {
        sz++;
      }
      selfend = sz;

      p4est_quadrant_ancestor (&desc[1], minlevel[1], &anc);
      p4est_quadrant_last_descendant (&anc, &ub, P4EST_QMAXLEVEL);
      sz = p4est_find_higher_bound (inquads, &ub, 0);
      if (sz < 0) {
        sz = 0;
      }
      else {
        sz++;
      }
      ancend = sz;
      sc_array_push_cut (inquads, selfend, ancend, &sort_bufs[1]);
    }
    else {
      selfend = inquads->elem_count;
    }
#if 0
    sc_array_crop (inquads, selfstart, selfend);
#endif
  }
  sc_mempool_destroy (list_alloc);

  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}

static void
p4est_balance_sort_sort (p4est_balance_obj_t * bobj,
                         p4est_t * p4est, int flt, int llt,
                         int num_trees, int num_sort,
                         int (*procrange)[2], sc_array_t * sort_bufs,
                         sc_statistics_t * stats)
{
  sc_flopinfo_t       snap;
  int                 rank = p4est->mpirank;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);
#ifdef P4EST_ENABLE_MPI
  if (stats) {
    if (!sc_statistics_has (stats, "balance sort send")) {
      sc_statistics_add_empty (stats, "balance sort send");
    }
    if (!sc_statistics_has (stats, "balance sort recv")) {
      sc_statistics_add_empty (stats, "balance sort recv");
    }
  }
  /* block balance parallel sort */
  if (num_sort > 0) {
    int                 i;
    MPI_Request         sendreq[2];
    MPI_Request         recvreq[2];
    int                 mpiret;
    int                 bit[2] = { 0, 0 };
    int                 proc[2] = { -1, -1 };
    int                 state[2];
    sc_array_t          recv_bufs[2], send_bufs[2];
    int                 num_active;

    sc_array_init (&recv_bufs[0], sizeof (p4est_quadrant_t));
    sc_array_init (&recv_bufs[1], sizeof (p4est_quadrant_t));
    sc_array_init (&send_bufs[0], sizeof (p4est_quadrant_t));
    sc_array_init (&send_bufs[1], sizeof (p4est_quadrant_t));

    bit[0] = bit[1] = -1;
    for (i = 0; i < 2; i++) {
      int                 shift =
        SC_LOG2_64 (procrange[i][1] - procrange[i][0]);
      if (shift >= 0) {
        bit[i] = 1 << shift;
      }
      else {
        bit[i] = 0;
      }
      state[i] = P4EST_BALSORT_ISEND;
    }

    num_active = 2;
    while (bit[0] > 0 || bit[1] > 0) {
      for (i = 0; i < 2; i++) {
        int                 flg;
        int                 rcount;
        MPI_Status          status;

        if (!bit[i]) {
          num_active = 1;
          continue;
        }

        switch (state[i]) {
        case P4EST_BALSORT_ISEND:
          /* all receives have been merged, all sends have completed:
           * initiate another send and receive */
          if (num_sort == 1 && bit[i ^ 1] > bit[i]) {
            num_active = 1;
            state[i] = P4EST_BALSORT_SYNC;
            break;
          }
          proc[i] = i ? rank + bit[i] : rank - bit[i];
          if (proc[i] >= procrange[i][0] && proc[i] <= procrange[i][1]) {
            /* prepare the send buffer */
            ssize_t             sz;
            p4est_quadrant_t    nextq =
              p4est->global_first_position[i ? proc[i] : proc[i] + 1];

            while (nextq.level > 0 && p4est_quadrant_child_id (&nextq) == 0) {
              nextq.level--;
            }
            sz = p4est_find_lower_bound (&sort_bufs[i], &nextq, 0);
            if (sz < 0) {
              sz = sort_bufs[i].elem_count;
            }
            if (!i) {
              sc_array_push_cut (&sort_bufs[i], 0, (size_t) sz,
                                 &send_bufs[i]);
            }
            else {
              sc_array_push_cut (&sort_bufs[i], (size_t) sz,
                                 sort_bufs[i].elem_count, &send_bufs[i]);
            }
            mpiret =
              MPI_Isend (send_bufs[i].array,
                         send_bufs[i].elem_count * sizeof (p4est_quadrant_t),
                         MPI_BYTE, proc[i], P4EST_COMM_BALANCE_SORT_SORT,
                         p4est->mpicomm, &sendreq[i]);
            SC_CHECK_MPI (mpiret);
            if (stats) {
              sc_statistics_accumulate (stats, "balance sort send",
                                        (double) send_bufs[i].elem_count);
            }
            state[i] = P4EST_BALSORT_IPROBE;
          }
          else {
            bit[i] >>= 1;
          }
          break;
        case P4EST_BALSORT_SYNC:
          if (bit[i ^ 1] <= bit[i]) {
            num_active = 2;
            state[i] = P4EST_BALSORT_ISEND;
          }
          break;
        case P4EST_BALSORT_IPROBE:
          /* a send has been initiated, probe for a matching receive */
          P4EST_ASSERT (proc[i] >= procrange[i][0]
                        && proc[i] <= procrange[i][1]);
          if (num_active == 1) {
            mpiret =
              MPI_Probe (proc[i], P4EST_COMM_BALANCE_SORT_SORT,
                         p4est->mpicomm, &status);
            SC_CHECK_MPI (mpiret);
            flg = 1;
          }
          else {
            mpiret =
              MPI_Iprobe (proc[i], P4EST_COMM_BALANCE_SORT_SORT,
                          p4est->mpicomm, &flg, &status);
            SC_CHECK_MPI (mpiret);
          }
          if (flg) {
            p4est_quadrant_t   *dest;

            mpiret = MPI_Get_count (&status, MPI_BYTE, &rcount);
            SC_CHECK_MPI (mpiret);
            dest =
              (p4est_quadrant_t *) sc_array_push_count (&recv_bufs[i],
                                                        (size_t) (rcount /
                                                                  sizeof
                                                                  (p4est_quadrant_t)));
            mpiret =
              MPI_Irecv (dest, rcount, MPI_BYTE, proc[i],
                         P4EST_COMM_BALANCE_SORT_SORT, p4est->mpicomm,
                         &recvreq[i]);
            SC_CHECK_MPI (mpiret);
            if (stats) {
              sc_statistics_accumulate (stats, "balance sort recv",
                                        (double) (rcount /
                                                  sizeof (p4est_quadrant_t)));
            }
            state[i] = P4EST_BALSORT_IRECV;
          }
          break;
        case P4EST_BALSORT_IRECV:
          /* a recv request has been initiated, wait for it to receive */
          if (num_active == 1) {
            mpiret = MPI_Wait (&recvreq[i], MPI_STATUS_IGNORE);
            SC_CHECK_MPI (mpiret);
            flg = 1;
          }
          else {
            mpiret = MPI_Test (&recvreq[i], &flg, MPI_STATUS_IGNORE);
            SC_CHECK_MPI (mpiret);
          }
          if (flg) {
            /* merge sort in the receive */
            if (num_sort == 2) {
              /* everything that I have received is for me */
              p4est_quadrant_array_merge_reduce (&sort_bufs[i],
                                                 &recv_bufs[i]);
            }
            else {
              /* what I have received goes to the opposite side */
              p4est_quadrant_array_merge_reduce (&sort_bufs[i ^ 1],
                                                 &recv_bufs[i]);
            }
            sc_array_truncate (&recv_bufs[i]);
            state[i] = P4EST_BALSORT_WAIT;
          }
          break;
        case P4EST_BALSORT_WAIT:
          /* wait for a send to finish */
          if (num_active == 1) {
            mpiret = MPI_Wait (&sendreq[i], MPI_STATUS_IGNORE);
            SC_CHECK_MPI (mpiret);
            flg = 1;
          }
          else {
            mpiret = MPI_Test (&sendreq[i], &flg, MPI_STATUS_IGNORE);
            SC_CHECK_MPI (mpiret);
          }
          if (flg) {
            bit[i] >>= 1;
            sc_array_truncate (&send_bufs[i]);
            state[i] = P4EST_BALSORT_ISEND;
          }
          break;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
    }

    if (num_sort == 1) {
      p4est_quadrant_array_merge_reduce (&sort_bufs[0], &sort_bufs[1]);
      sc_array_truncate (&sort_bufs[1]);
    }

    sc_array_reset (&recv_bufs[1]);
    sc_array_reset (&recv_bufs[0]);
    sc_array_reset (&send_bufs[1]);
    sc_array_reset (&send_bufs[0]);
  }
#endif
  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}

static void
p4est_balance_sort_merge_seeds (p4est_balance_obj_t * bobj,
                                p4est_t * p4est, p4est_connect_type_t btype,
                                int flt, int llt, int num_trees,
                                int num_sort, int *minlevel,
                                p4est_quadrant_t * desc,
                                sc_array_t * tree_bufs,
                                sc_array_t * sort_bufs)
{
  sc_flopinfo_t       snap;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);
#ifdef P4EST_ENABLE_MPI
  if (num_sort > 0) {
    int                 i;

    /* combine sort results back into the original trees */
    for (i = 0; i < num_sort; i++) {
      p4est_quadrant_t    anc;
      size_t              nq, s;
      sc_array_t          seeds;
      sc_array_t          newquads;

      p4est_quadrant_ancestor (&desc[i], minlevel[i], &anc);

      sc_array_init (&seeds, sizeof (p4est_quadrant_t));
      sc_array_init (&newquads, sizeof (p4est_quadrant_t));
      nq = sort_bufs[i].elem_count;
      for (s = 0; s < nq; s++) {
        p4est_quadrant_t   *q =
          (p4est_quadrant_t *) sc_array_index (&sort_bufs[i], s);
        p4est_quadrant_t   *prev =
          s ? (p4est_quadrant_t *) sc_array_index (&sort_bufs[i],
                                                   s - 1) : NULL;
        p4est_quadrant_t   *next =
          (s < nq - 1) ? (p4est_quadrant_t *) sc_array_index (&sort_bufs[i],
                                                              s + 1) : NULL;
        p4est_quadrant_t    lb, ub;
        int32_t             touch;
        int32_t             mask = 1;

        P4EST_ASSERT (p4est_quadrant_is_ancestor (&anc, q));
        if (prev) {
          p4est_quadrant_t    a, qstop;
          int                 qid, pid, idxor, dir, qstopid;

          p4est_nearest_common_ancestor (q, prev, &a);
          P4EST_ASSERT (a.level < q->level && a.level < prev->level);
          qid = p4est_quadrant_ancestor_id (q, a.level + 1);
          pid = p4est_quadrant_ancestor_id (prev, a.level + 1);
          P4EST_ASSERT (qid > pid);
          idxor = qid ^ pid;
          dir = SC_LOG2_8 (idxor);
          P4EST_ASSERT (dir >= 0);
          qstopid = (qid >> dir) << dir;
          p4est_quadrant_child (&a, &qstop, qstopid);
          p4est_quadrant_first_descendant (&qstop, &lb, P4EST_QMAXLEVEL);
        }
        if (next) {
          p4est_quadrant_t    a, qstop;
          int                 qid, nid, idxor, dir, qstopid;

          p4est_nearest_common_ancestor (q, next, &a);
          P4EST_ASSERT (a.level < q->level && a.level < next->level);
          qid = p4est_quadrant_ancestor_id (q, a.level + 1);
          nid = p4est_quadrant_ancestor_id (next, a.level + 1);
          P4EST_ASSERT (nid > qid);
          idxor = qid ^ nid;
          dir = SC_LOG2_8 (idxor);
          P4EST_ASSERT (dir >= 0);
          qstopid = ((nid >> dir) << dir) - 1;
          p4est_quadrant_child (&a, &qstop, qstopid);
          p4est_quadrant_last_descendant (&qstop, &ub, P4EST_QMAXLEVEL);
        }
        touch =
          p4est_find_range_boundaries (prev ? &lb : NULL, next ? &ub : NULL,
                                       anc.level, NULL,
#ifdef P4_TO_P8
                                       NULL,
#endif
                                       NULL);
        {
          int                 f;

          for (f = 0; f < P4EST_FACES; f++, mask <<= 1) {
            p4est_quadrant_t    n;
            int                 split;

            if (!(touch & mask)) {
              continue;
            }
            p4est_quadrant_face_neighbor (&anc, f, &n);
            split = p4est_balance_seeds (q, &n, btype, &seeds);
            if (split) {
              sc_array_push_cut (&seeds, 0, seeds.elem_count, &newquads);
            }
          }
        }
#ifdef P4_TO_P8
        {
          int                 e;

          for (e = 0; e < P8EST_EDGES; e++, mask <<= 1) {
            p4est_quadrant_t    n;
            int                 split;

            if (!(touch & mask)) {
              continue;
            }
            p8est_quadrant_edge_neighbor (&anc, e, &n);
            split = p4est_balance_seeds (q, &n, btype, &seeds);
            if (split) {
              sc_array_push_cut (&seeds, 0, seeds.elem_count, &newquads);
            }
          }
        }
#endif
        {
          int                 c;

          for (c = 0; c < P4EST_CHILDREN; c++, mask <<= 1) {
            p4est_quadrant_t    n;
            int                 split;

            if (!(touch & mask)) {
              continue;
            }
            p4est_quadrant_corner_neighbor (&anc, c, &n);
            split = p4est_balance_seeds (q, &n, btype, &seeds);
            if (split) {
              sc_array_push_cut (&seeds, 0, seeds.elem_count, &newquads);
            }
          }
        }
      }
      sc_array_sort (&newquads, p4est_quadrant_compare);
      p4est_quadrant_array_reduce (&newquads, NULL, NULL);
      p4est_quadrant_array_merge_reduce (&sort_bufs[i], &newquads);
      p4est_quadrant_array_merge_reduce (&tree_bufs
                                         [desc[i].p.which_tree - flt],
                                         &sort_bufs[i]);
      sc_array_reset (&seeds);
      sc_array_reset (&newquads);
    }

  }
#endif
  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}

static int8_t
p4est_process_level (p4est_t *p4est, int p)
{
  p4est_quadrant_t f, f_mid, l_mid, l;
  p4est_quadrant_t af, al, afl, alf, aff, all;
  int8_t ancestor_level, f_level, l_level;
  P4EST_ASSERT (!p4est_comm_is_empty(p4est, p));

  P4EST_QUADRANT_INIT(&f_mid);
  P4EST_QUADRANT_INIT(&l_mid);

  f = p4est->global_first_position[p];
  l = p4est->global_first_position[p+1];
  if (l.x == 0 && l.y == 0 &&
#ifdef P4_TO_P8
      l.z == 0 &&
#endif
      1) {
    l.p.which_tree--;
    l.x = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN(P4EST_QMAXLEVEL);
    l.y = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN(P4EST_QMAXLEVEL);
#ifdef P4_TO_P8
    l.z = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN(P4EST_QMAXLEVEL);
#endif
  } else {
    uint64_t id;

    id = p4est_quadrant_linear_id (&l, P4EST_QMAXLEVEL);
    p4est_quadrant_set_morton (&l, P4EST_QMAXLEVEL, id - 1);
  }
  P4EST_ASSERT (p4est_quadrant_compare_piggy (&f, &l) <= 0);
  if (l.p.which_tree > f.p.which_tree + 1) {
    /* the range contains a root */
    return 0;
  }
  if (l.p.which_tree == f.p.which_tree + 1) {
    ancestor_level = -1;
    f_mid.x = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN(P4EST_QMAXLEVEL);
    f_mid.y = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN(P4EST_QMAXLEVEL);
#ifdef P4_TO_P8
    f_mid.z = P4EST_ROOT_LEN - P4EST_QUADRANT_LEN(P4EST_QMAXLEVEL);
#endif
    f_mid.level = P4EST_QMAXLEVEL;
    l_mid.x = 0;
    l_mid.y = 0;
#ifdef P4_TO_P8
    l_mid.z = 0;
#endif
    l_mid.level = P4EST_QMAXLEVEL;
  }
  else {
    int8_t f_id, l_id;
    p4est_quadrant_t a;

    p4est_nearest_common_ancestor (&f, &l, &a);
    ancestor_level = a.level;

    p4est_quadrant_corner_descendant (&a, &af, 0, P4EST_QMAXLEVEL);
    p4est_quadrant_corner_descendant (&a, &al, P4EST_CHILDREN - 1, P4EST_QMAXLEVEL);
    if (p4est_quadrant_is_equal (&f, &af) &&
        p4est_quadrant_is_equal (&l, &al)) {
      /* the range is equivalent to the least common ancestor */
      return ancestor_level;
    }

    f_id = p4est_quadrant_ancestor_id (&f, ancestor_level + 1);
    l_id = p4est_quadrant_ancestor_id (&l, ancestor_level + 1);
    if (l_id > f_id + 1) {
      /* a full child of the ancestor is present in the range */
      return ancestor_level + 1;
    }
  }
  p4est_quadrant_ancestor (&f, ancestor_level + 1, &af);
  p4est_quadrant_corner_descendant (&af, &afl, P4EST_CHILDREN - 1, P4EST_QMAXLEVEL);
  p4est_nearest_common_ancestor(&f, &afl, &af);
  p4est_quadrant_corner_descendant (&af, &aff, 0, P4EST_QMAXLEVEL);
  if (p4est_quadrant_is_equal (&f, &aff)) {
    f_level = af.level;
  }
  else {
    f_level = af.level + 1;
  }

  p4est_quadrant_ancestor (&l, ancestor_level + 1, &al);
  p4est_quadrant_corner_descendant (&al, &alf, 0, P4EST_QMAXLEVEL);
  p4est_nearest_common_ancestor(&alf, &l, &al);
  p4est_quadrant_corner_descendant (&al, &all, P4EST_CHILDREN - 1, P4EST_QMAXLEVEL);
  if (p4est_quadrant_is_equal (&l, &all)) {
    l_level = al.level;
  }
  else {
    l_level = al.level + 1;
  }
  return SC_MIN(l_level, f_level);
}

enum {
  NEIGH_SAME = 0,
  NEIGH_LARGER,
  NEIGH_SMALLER,
};

typedef struct neigh_info
{
  int root_proc;
  int leaf_proc;
  size_t count;
}
neigh_info_t;

typedef struct neigh_entry
{
  union neigh_type {
    neigh_info_t info;
    p4est_quadrant_t quad;
  } p;
}
neigh_entry_t;

/* we will use tags to keep the multiple tree communications separate */
enum
{
  P4EST_COMM_BALANCE_SORT_NEIGH_TREE = P4EST_COMM_TAG_LAST,
};

static void
p4est_balance_sort_neigh (p4est_balance_obj_t * bobj,
                          p4est_t * p4est, int flt, int llt, int num_trees,
                          p4est_tree_neigh_info_t * tinfo,
                          sc_array_t * tree_bufs,
                          p4est_neigh_t *neigh,
                          sc_statistics_t * stats)
{
  p4est_topidx_t      t;
  sc_flopinfo_t       snap;
  int                 rank = p4est->mpirank;
  int                 n_neigh, nq;
  const int          *neigh_procs;
  sc_array_t        **send_arrays;
  sc_array_t         *recv_array;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);

  p4est_neigh_get_procs (neigh, &n_neigh, &neigh_procs);

  send_arrays = P4EST_ALLOC (sc_array_t *,n_neigh + 1);
  for (nq = 0; nq < n_neigh + 1; nq++) {
    send_arrays[nq] = sc_array_new (sizeof (p4est_quadrant_t));
  }
  recv_array = send_arrays[n_neigh];

  /* ghost neighbor balance */
  /* loop over trees: transform non-root quads into their rooted version in
   * neighbor trees */
  for (t = 0; t < num_trees; t++) {
    int                 i, j, k, l;
#ifdef P4_TO_P8
    int                 kstart = 0, kend = 3;
#else
    int                 kstart = 1, kend = 2;
#endif
    sc_array_t         *tquads = &tree_bufs[t];
    sc_array_t          orig_view;
    sc_array_t          tform_quads;
    p4est_qcoord_t      H = P4EST_ROOT_LEN;
    p4est_tree_neigh_info_t *info = &tinfo[t];
    p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t + flt);

    sc_array_init (&tform_quads, sizeof (p4est_quadrant_t));
    for (k = kstart, l = 0; k < kend; k++) {
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++, l++) {
          p4est_quadrant_t    insul, ub;
          size_t              sz;
          size_t              qstart, qend;

          insul.x = (i - 1) * H;
          insul.y = (j - 1) * H;
#ifdef P4_TO_P8
          insul.z = (k - 1) * H;
#endif
          insul.level = 0;
          p4est_quadrant_last_descendant (&insul, &ub, P4EST_QMAXLEVEL);
          sz = p4est_find_lower_bound (tquads, &insul, 0);
          if (sz < 0) {
            sz = tquads->elem_count;
          }
          qstart = (size_t) sz;
          sz = p4est_find_higher_bound (tquads, &ub, 0);
          if (sz < 0) {
            sz = 0;
          }
          else {
            sz++;
          }
          qend = (size_t) sz;
          if (qend <= qstart) {
            continue;
          }
          sc_array_init_view (&orig_view, tquads, qstart, qend - qstart);
          if (l == P4EST_INSUL / 2) {
            p4est_balance_sort_divide (p4est, n_neigh, neigh_procs,
                                       send_arrays,
                                       t + flt, &orig_view, 1);
          }
          else {
            size_t              z, sq;
            size_t              no = orig_view.elem_count;

            for (z = info->offset[l]; z < info->offset[l + 1]; z++) {
              p4est_tree_neigh_t *tn =
                (p4est_tree_neigh_t *) sc_array_index (&info->tnarray, z);
              p4est_quadrant_t   *oq =
                (p4est_quadrant_t *) sc_array_index (&orig_view, 0);
              p4est_quadrant_t   *trq;

              trq =
                (p4est_quadrant_t *) sc_array_push_count (&tform_quads, no);

              for (sq = 0; sq < no; sq++) {
                p4est_quadrant_t    temp;
                p4est_quadrant_utransform (&oq[sq], &temp, &(tn->u[0]), 0);

                p4est_quadrant_sibling (&temp, &trq[sq], 0);
                P4EST_ASSERT (p4est_quadrant_is_valid (&trq[sq]));
              }
              sc_array_sort (&tform_quads, p4est_quadrant_compare);
              p4est_balance_sort_divide (p4est, n_neigh, neigh_procs,
                                         send_arrays,
                                         tn->nt, &tform_quads, 0);
              sc_array_truncate (&tform_quads);
            }
          }
        }
      }
    }

    if (stats) {
      if (!sc_statistics_has (stats, "balance neigh")) {
        sc_statistics_add_empty (stats, "balance neigh");
      }
      sc_statistics_accumulate (stats, "balance neigh", (double) n_neigh);
      if (!sc_statistics_has (stats, "balance neigh send")) {
        sc_statistics_add_empty (stats, "balance neigh send");
      }
      if (!sc_statistics_has (stats, "balance neigh recv")) {
        sc_statistics_add_empty (stats, "balance neigh recv");
      }
    }

    {
      p4est_quadrant_t    fd = tree->first_desc;
      ssize_t             sz;
      size_t              qstart, qend;

      while (fd.level > 0 && !p4est_quadrant_child_id (&fd)) {
        fd.level--;
      }
      sz = p4est_find_lower_bound (tquads, &fd, 0);
      if (sz < 0) {
        sz = tquads->elem_count;
      }
      qstart = sz;
      sz = p4est_find_higher_bound (tquads, &tree->last_desc, 0);
      if (sz < 0) {
        sz = 0;
      }
      else {
        sz++;
      }
      qend = sz;
      if (qend > qstart) {
        sc_array_crop (tquads, qstart, qend);
      }
      else {
        sc_array_truncate (tquads);
      }
    }
    sc_array_reset (&tform_quads);
  }

  p4est_neigh_alltoallx (neigh, send_arrays, recv_array, NULL);

  /* exchange, merge and complete */
  if (num_trees) {
    sc_array_sort (recv_array, p4est_quadrant_compare_piggy);
#ifdef P4EST_ENABLE_DEBUG
    if (recv_array->elem_count) {
      p4est_quadrant_t   *first = p4est_quadrant_array_index (recv_array, 0);
      p4est_quadrant_t   *last =
        p4est_quadrant_array_index (recv_array, recv_array->elem_count - 1);
      p4est_quadrant_t   *gf = &p4est->global_first_position[rank];
      p4est_quadrant_t   *gn = &p4est->global_first_position[rank + 1];

      P4EST_ASSERT (first->p.which_tree >= gf->p.which_tree);
      P4EST_ASSERT (last->p.which_tree <= gn->p.which_tree);
      if (first->p.which_tree == gf->p.which_tree) {
        P4EST_ASSERT (p4est_quadrant_compare (gf, first) < 0
                      || (gf->x == first->x && gf->y == first->y &&
#ifdef P4_TO_P8
                          gf->z == first->z &&
#endif
                          1));

      }
      if (last->p.which_tree == gn->p.which_tree) {
        P4EST_ASSERT (p4est_quadrant_compare (last, gn) < 0);
      }
    }
#endif
    {
      size_t              s, n;

      n = recv_array->elem_count;
      for (s = 0; s < n;) {
        size_t              z;
        p4est_quadrant_t   *q = p4est_quadrant_array_index (recv_array, s);
        p4est_topidx_t      t = q->p.which_tree;
        sc_array_t          view;

        for (z = s + 1; z < n; z++) {
          p4est_quadrant_t   *r = p4est_quadrant_array_index (recv_array, z);

          if (r->p.which_tree != t) {
            break;
          }
        }
        sc_array_init_view (&view, recv_array, s, z - s);
        p4est_quadrant_array_reduce (&view, NULL, NULL);
        p4est_quadrant_array_merge_reduce (&tree_bufs[t - flt], &view);
        s = z;
      }
    }
  }

  for (nq = 0; nq < n_neigh + 1; nq++) {
    sc_array_destroy (send_arrays[nq]);
  }
  P4EST_FREE (send_arrays);
  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}

static void
p4est_balance_sort_complete (p4est_balance_obj_t * bobj,
                             p4est_t * p4est, int flt, int llt, int num_trees,
                             sc_array_t * tree_bufs,
                             p4est_init_t init_fn, p4est_replace_t replace_fn)
{
  p4est_topidx_t      t;
  sc_flopinfo_t       snap;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);

  if (num_trees) {
    sc_array_t          newcomplete;

    sc_array_init (&newcomplete, sizeof (p4est_quadrant_t));
    p4est->local_num_quadrants = 0;
    for (t = flt; t <= llt; t++) {
      p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, t);
      p4est_quadrant_t    dom;

      p4est_nearest_common_ancestor (&tree->first_desc, &tree->last_desc,
                                     &dom);
      (void) sc_array_push_count (&newcomplete,
                                  tree_bufs[t - flt].elem_count);
      sc_array_truncate (&newcomplete);
      p4est_complete_kernel (&tree_bufs[t - flt], &dom, &tree->first_desc,
                             &tree->last_desc, &newcomplete);
      p4est_subtree_replace (p4est, t, &newcomplete, init_fn, replace_fn);
      p4est->local_num_quadrants += newcomplete.elem_count;
      sc_array_truncate (&newcomplete);
    }
    sc_array_reset (&newcomplete);
  }

  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}

void
p4est_balance_sort (p4est_balance_obj_t * bobj, p4est_t * p4est)
{
  p4est_topidx_t      t, num_trees;
  p4est_topidx_t      flt = p4est->first_local_tree;
  p4est_topidx_t      llt = p4est->last_local_tree;
  p4est_tree_neigh_info_t *tinfo = NULL;
  int                 num_sort = -1;
  sc_flopinfo_t       snap;
  int                 minlevel[2];
  int                 procrange[2][2];
  p4est_quadrant_t    desc[2];
  int                 rank = p4est->mpirank;
  sc_array_t         *tree_bufs;
  sc_array_t          sort_bufs[2];
  p4est_connect_type_t btype;
  p4est_init_t        init_fn;
  p4est_replace_t     replace_fn;
  p4est_neigh_t      *neigh;
  sc_statistics_t    *stats = NULL;

  P4EST_BAL_FUNC_SNAP (bobj, &snap);

  btype = p4est_balance_obj_get_connect (bobj);
  init_fn = p4est_balance_obj_get_init (bobj);
  replace_fn = p4est_balance_obj_get_replace (bobj);
  stats = p4est_balance_obj_get_stats (bobj);

  /* first figure out the parallel neighborhood */
  num_trees = SC_MAX (0, llt + 1 - flt);
  if (!num_trees) {
    p4est_quadrant_t    q;

    q = p4est->global_first_position[rank];
    P4EST_ASSERT (p4est_quadrant_is_equal_piggy
                  (&q, &p4est->global_first_position[rank + 1]));
    if (q.x != 0 || q.y != 0 ||
#ifdef P4_TO_P8
        q.z != 0 ||
#endif
        0) {
      flt = llt = q.p.which_tree;
    }
  }

  {
    sc_statistics_t    *stats = NULL;

    stats = p4est_balance_obj_get_stats (bobj);
    if (stats) {
      if (!sc_statistics_has (stats, "balance neigh time")) {
        sc_statistics_add_empty (stats, "balance neigh time");
      }
      if (!sc_statistics_has (stats, "balance neigh recv")) {
        sc_statistics_add_empty (stats, "balance neigh recv");
      }
      if (!sc_statistics_has (stats, "balance neigh send")) {
        sc_statistics_add_empty (stats, "balance neigh send");
      }
      if (!sc_statistics_has (stats, "balance sort time")) {
        sc_statistics_add_empty (stats, "balance sort time");
      }
    }
  }

  tinfo = P4EST_ALLOC (p4est_tree_neigh_info_t, llt + 1 - flt);
  p4est_balance_sort_compute_pattern_sort (bobj, p4est, flt, llt, num_trees,
                                           tinfo,
                                           procrange, minlevel, &num_sort,
                                           desc);
  tree_bufs = P4EST_ALLOC (sc_array_t, num_trees);
  for (t = 0; t < num_trees; t++) {
    sc_array_init (&tree_bufs[t], sizeof (p4est_quadrant_t));
  }

  sc_array_init (&sort_bufs[0], sizeof (p4est_quadrant_t));
  sc_array_init (&sort_bufs[1], sizeof (p4est_quadrant_t));

  p4est_balance_sort_local (bobj, p4est, btype, flt, llt, num_trees,
                            num_sort, minlevel, desc, tree_bufs, sort_bufs);

  p4est_balance_sort_sort (bobj, p4est, flt, llt, num_trees, num_sort,
                           procrange, sort_bufs, stats);

  p4est_balance_sort_merge_seeds (bobj, p4est, btype, flt, llt, num_trees,
                                  num_sort, minlevel, desc,
                                  tree_bufs, sort_bufs);

  sc_array_reset (&sort_bufs[1]);
  sc_array_reset (&sort_bufs[0]);

  p4est_balance_sort_compute_neigh (bobj, p4est, flt, llt, num_trees,
                                    tinfo, minlevel, desc, &neigh);
  p4est_balance_sort_neigh (bobj, p4est, flt, llt, num_trees, tinfo, tree_bufs,
                            neigh, stats);
  if (neigh != bobj->neigh) {
    p4est_neigh_destroy (neigh);
  }

  p4est_balance_sort_complete (bobj, p4est, flt, llt, num_trees, tree_bufs,
                               init_fn, replace_fn);

  for (t = 0; t < num_trees; t++) {
    sc_array_reset (&tree_bufs[t]);
  }
  P4EST_FREE (tree_bufs);
  for (t = flt; t <= llt; t++) {
    sc_array_reset (&(tinfo[t - flt].tnarray));
  }
  P4EST_FREE (tinfo);
  P4EST_BAL_FUNC_SHOT (bobj, &snap);
}
