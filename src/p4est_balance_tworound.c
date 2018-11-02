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
#include <p8est_algorithms.h>
#include <p8est_communication.h>
#include <p8est_bits.h>
#else
#include <p4est_balance_obj.h>
#include <p4est_extended.h>
#include <p4est_algorithms.h>
#include <p4est_communication.h>
#include <p4est_bits.h>
#endif

static const int8_t fully_owned_flag = 0x01;
static const int8_t any_face_flag = 0x02;

typedef struct
{
  int8_t              have_first_count, have_first_load;
  int8_t              have_second_count, have_second_load;
  int                 recv_first_count, recv_second_count;
  int                 send_first_count, send_second_count;
  sc_array_t          send_first, send_second, recv_first, recv_second;
}
p4est_balance_peer_t;


static void
p4est_balance_response (p4est_t * p4est, p4est_balance_peer_t * peer,
                        p4est_connect_type_t balance, sc_array_t * borders,
						p4est_balance_obj_t * bobj)
{
  sc_array_t         *first_seeds = sc_array_new (sizeof (p4est_quadrant_t));

/* compute and uniqify overlap quadrants */
  p4est_tree_compute_overlap (p4est, &peer->recv_first,
                              &peer->send_second, balance, borders,
                              first_seeds);
  /* replace peer->recv_first with first_seeds */
  p4est_tree_uniqify_overlap (&peer->send_second);
  p4est_tree_uniqify_overlap (first_seeds);
  /* replace peer->recv_first with first_seeds */
  sc_array_resize (&peer->recv_first, first_seeds->elem_count);
  memcpy (peer->recv_first.array, first_seeds->array,
          first_seeds->elem_size * first_seeds->elem_count);
  sc_array_destroy (first_seeds);

  if (bobj->inspect) {
    bobj->inspect->balance_comm_sent += peer->send_second.elem_count;
    if (peer->send_second.elem_count) {
      bobj->inspect->balance_comm_nzpeers++;
    }
  }
}

/** Check if the insulation layer of a quadrant overlaps anybody.
 * If yes, the quadrant itself is scheduled for sending.
 * Both quadrants are in the receiving tree's coordinates.
 * \param [in]  qtree       Tree id of the receiving tree.
 * \param [in]  inter_tree  Boolean flag to specify inter-tree communication.
 * \param [in]  q           The quadrant to be sent if there is overlap.
 * \param [in]  insul       An insulation quadrant of \a q.
 * \param [in,out]  first_peer  Lowest peer, will be updated.
 * \param [in,out]  last_peer   Highest peer, will be updated.
 */
static void
p4est_balance_schedule (p4est_t * p4est, p4est_balance_peer_t * peers,
                        p4est_topidx_t qtree, int inter_tree,
                        const p4est_quadrant_t * q,
                        const p4est_quadrant_t * insul,
                        int *first_peer, int *last_peer)
{
  const int           rank = p4est->mpirank;
  int                 found;
  int                 back, pos;
  int                 owner, first_owner, last_owner;
  p4est_gloidx_t     *global_first_quadrant = p4est->global_first_quadrant;
  p4est_quadrant_t    ld, *s;
  p4est_balance_peer_t *peer;

  P4EST_QUADRANT_INIT (&ld);

  /* querying insul is equivalent to querying first descendant */
  first_owner = p4est_comm_find_owner (p4est, qtree, insul, rank);
  /* querying last descendant */
  p4est_quadrant_last_descendant (insul, &ld, P4EST_QMAXLEVEL);
  last_owner = p4est_comm_find_owner (p4est, qtree, &ld, rank);

  /* send to all processors possibly intersecting insulation */
  for (owner = first_owner; owner <= last_owner; ++owner) {
    if (owner == rank && !inter_tree) {
      /* do not send to self for the same tree */
      continue;
    }
    if (global_first_quadrant[owner] == global_first_quadrant[owner + 1]) {
      /* do not send to empty processors */
      continue;
    }
    peer = peers + owner;
    /* avoid duplicates in the send array */
    found = 0;
    for (back = 0; back < P4EST_INSUL - 1; ++back) {
      pos = (int) peer->send_first.elem_count - back - 1;
      if (pos < 0) {
        break;
      }
      s = (p4est_quadrant_t *) sc_array_index_int (&peer->send_first, pos);
      if (p4est_quadrant_is_equal (s, q) && s->p.piggy2.which_tree == qtree
          && s->p.piggy2.from_tree == q->p.piggy2.from_tree
          && s->pad16 == q->pad16) {
        found = 1;
        break;
      }
    }
    if (found) {
      continue;
    }

    /* copy quadrant into shipping list */
    s = p4est_quadrant_array_push (&peer->send_first);
    *s = *q;
    P4EST_ASSERT (p4est_quadrant_is_extended (q));
    s->p.piggy2.which_tree = qtree;     /* piggy back tree id */

    /* update lowest and highest peer */
    if (owner != rank) {
      *first_peer = SC_MIN (owner, *first_peer);
      *last_peer = SC_MAX (owner, *last_peer);
    }
  }
}

void
p4est_balance_tworound (p4est_balance_obj_t *bobj, p4est_t *p4est)
{
  const int           rank = p4est->mpirank;
  const int           num_procs = p4est->mpisize;
  int                 j, k, l, m, which;
#ifdef P4EST_ENABLE_DEBUG
  int                 facewhich;
#endif
  int                 face;
  int                 first_peer, last_peer;
  int                 quad_contact[P4EST_FACES];
  int                 any_boundary, tree_contact[P4EST_INSUL];
  int                 tree_fully_owned, full_tree[2];
  int8_t             *tree_flags;
  size_t              zz, treecount, ctree;
  size_t              localcount;
  size_t              qcount, qbytes;
  size_t              all_incount, all_outcount;
  p4est_qcoord_t      qh;
  const p4est_qcoord_t rh = P4EST_ROOT_LEN;
  p4est_topidx_t      qtree, nt;
  p4est_topidx_t      first_tree, last_tree;
  p4est_locidx_t      skipped;
  p4est_balance_peer_t *peers, *peer;
  p4est_tree_t       *tree;
  p4est_quadrant_t    mylow, nextlow;
  p4est_quadrant_t    tosend, insulq, tempq;
  p4est_quadrant_t   *q, *s;
  p4est_connectivity_t *conn = p4est->connectivity;
  sc_array_t         *qarray, *tquadrants;
  sc_array_t         *borders;
  const int8_t       *pre_adapt_flags = NULL;
#ifdef P4EST_ENABLE_DEBUG
  size_t              data_pool_size;
#endif
  int                 ftransform[P4EST_FTRANSFORM];
  int                 face_axis[3];     /* 3 not P4EST_DIM */
  int                 contact_face_only, contact_edge_only;
#ifdef P4_TO_P8
  int                 edge;
  size_t              etree;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta;
#endif
  int                 corner;
  p4est_corner_info_t ci;
  p4est_corner_transform_t *ct;
  sc_array_t         *cta;
#ifdef P4EST_ENABLE_MPI
#ifdef P4EST_ENABLE_DEBUG
  unsigned            checksum;
  sc_array_t          checkarray;
  p4est_gloidx_t      ltotal[2], gtotal[2];
#endif /* P4EST_ENABLE_DEBUG */
  int                 i;
  int                 mpiret, rcount;
  int                 request_first_send, request_first_recv, outcount;
  int                 request_second_send, total_send_count, total_recv_count;
  int                 send_zero[2], send_load[2];
  int                 recv_zero[2], recv_load[2];
  int                *wait_indices;
  int                *receiver_ranks, *sender_ranks;
  int                 num_receivers, num_senders;
  MPI_Request        *send_requests_first_load;
  MPI_Request        *send_requests_second_load;
  MPI_Request        *recv_requests_first_load;
  MPI_Request        *recv_requests_second_load;
  MPI_Status         *recv_statuses, *jstatus;
  sc_array_t         *receivers, *senders;
  sc_array_t         *in_counts, *out_counts;
  int                *icounts;
  sc_flopinfo_t       snap;
#endif /* P4EST_ENABLE_MPI */
  sc_notify_t        *notify = NULL;
  int                 own_notify = 0;
  p4est_connect_type_t btype;
  p4est_init_t init_fn;
  p4est_replace_t replace_fn;

  P4EST_FUNC_SNAP (p4est, &snap);

  btype = p4est_balance_obj_get_connect (bobj);
  init_fn = p4est_balance_obj_get_init (bobj);
  replace_fn = p4est_balance_obj_get_replace (bobj);

#ifndef P4_TO_P8
  P4EST_ASSERT (btype == P4EST_CONNECT_FACE || btype == P4EST_CONNECT_CORNER);
#else
  P4EST_ASSERT (btype == P8EST_CONNECT_FACE || btype == P8EST_CONNECT_EDGE ||
                btype == P8EST_CONNECT_CORNER);
#endif

#ifdef P4EST_ENABLE_DEBUG
  data_pool_size = 0;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }
#endif

  P4EST_QUADRANT_INIT (&mylow);
  P4EST_QUADRANT_INIT (&nextlow);
  P4EST_QUADRANT_INIT (&tosend);
  P4EST_QUADRANT_INIT (&insulq);
  P4EST_QUADRANT_INIT (&tempq);

  /* tree status flags (max 8 per tree) */
  tree_flags = P4EST_ALLOC (int8_t, conn->num_trees);
  for (nt = 0; nt < conn->num_trees; ++nt) {
    tree_flags[nt] = 0x00;
  }

  localcount = (size_t) (p4est->last_local_tree + 1 -
                         p4est->first_local_tree);
  borders = sc_array_new_size (sizeof (sc_array_t), localcount);
  for (zz = 0; zz < localcount; zz++) {
    qarray = (sc_array_t *) sc_array_index (borders, zz);
    sc_array_init (qarray, sizeof (p4est_quadrant_t));
  }

#ifdef P4EST_ENABLE_MPI
  send_requests_first_load = P4EST_ALLOC (MPI_Request, 4 * num_procs);
  recv_requests_first_load = send_requests_first_load + 1 * num_procs;
  send_requests_second_load = send_requests_first_load + 2 * num_procs;
  recv_requests_second_load = send_requests_first_load + 3 * num_procs;
  recv_statuses = P4EST_ALLOC (MPI_Status, num_procs);
  for (j = 0; j < num_procs; ++j) {
    send_requests_first_load[j] = MPI_REQUEST_NULL;
    recv_requests_first_load[j] = MPI_REQUEST_NULL;
    send_requests_second_load[j] = MPI_REQUEST_NULL;
    recv_requests_second_load[j] = MPI_REQUEST_NULL;
  }
  wait_indices = P4EST_ALLOC (int, num_procs);
#ifdef P4EST_ENABLE_DEBUG
  sc_array_init (&checkarray, 4);
#endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */

  /* allocate per peer storage and initialize requests */
  peers = P4EST_ALLOC (p4est_balance_peer_t, num_procs);
  for (j = 0; j < num_procs; ++j) {
    peer = peers + j;
    sc_array_init (&peer->send_first, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->send_second, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->recv_first, sizeof (p4est_quadrant_t));
    sc_array_init (&peer->recv_second, sizeof (p4est_quadrant_t));
    peer->send_first_count = peer->send_second_count = 0;
    peer->recv_first_count = peer->recv_second_count = 0;
    peer->have_first_count = peer->have_first_load = 0;
    peer->have_second_count = peer->have_second_load = 0;
  }
#ifdef P4_TO_P8
  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));

  /* compute first quadrant on finest level */
  mylow.x = p4est->global_first_position[rank].x;
  mylow.y = p4est->global_first_position[rank].y;
#ifdef P4_TO_P8
  mylow.z = p4est->global_first_position[rank].z;
#endif
  mylow.level = P4EST_QMAXLEVEL;

  /* and the first finest quadrant of the next processor */
  nextlow.x = p4est->global_first_position[rank + 1].x;
  nextlow.y = p4est->global_first_position[rank + 1].y;
#ifdef P4_TO_P8
  nextlow.z = p4est->global_first_position[rank + 1].z;
#endif
  nextlow.level = P4EST_QMAXLEVEL;

  /* start balance_A timing */
  if (bobj->inspect != NULL) {
    bobj->inspect->balance_A = -sc_MPI_Wtime ();
    bobj->inspect->balance_A_count_in = 0;
    bobj->inspect->balance_A_count_out = 0;
    bobj->inspect->use_B = 0;
    pre_adapt_flags = bobj->inspect->pre_adapt_flags;
  }

  /* loop over all local trees to assemble first send list */
  first_tree = p4est->first_local_tree;
  last_tree = p4est->last_local_tree;
  first_peer = rank;
  last_peer = rank;
  all_incount = 0;
  skipped = 0;
  for (nt = first_tree; nt <= last_tree; ++nt) {
    const int8_t       *this_pre_adapt_flags = NULL;
    p4est_comm_tree_info (p4est, nt, full_tree, tree_contact, NULL, NULL);
    tree_fully_owned = full_tree[0] && full_tree[1];
    any_boundary = 0;
    for (i = 0; i < P4EST_INSUL; ++i) {
      any_boundary |= tree_contact[i];
    }
    if (any_boundary) {
      tree_flags[nt] |= any_face_flag;
    }
    tree = p4est_tree_array_index (p4est->trees, nt);
    tquadrants = &tree->quadrants;
    if (pre_adapt_flags) {
      this_pre_adapt_flags = &pre_adapt_flags[all_incount];
    }
    all_incount += tquadrants->elem_count;

    /* initial log message for this tree */
    P4EST_VERBOSEF ("Into balance tree %lld with %llu\n", (long long) nt,
                    (unsigned long long) tquadrants->elem_count);

    /* local balance first pass */
    p4est_balance_subtree_ext (p4est, btype, nt, this_pre_adapt_flags,
                               init_fn, replace_fn);
    treecount = tquadrants->elem_count;
    P4EST_VERBOSEF ("Balance tree %lld A %llu\n",
                    (long long) nt, (unsigned long long) treecount);

    /* check if this tree is not shared with other processors */
    if (tree_fully_owned) {
      /* all quadrants in this tree are owned by me */
      tree_flags[nt] |= fully_owned_flag;
      if (!any_boundary) {
        /* this tree is isolated, no balance between trees */
        continue;
      }
    }

    if (borders != NULL) {
      qarray = (sc_array_t *) sc_array_index (borders,
                                              (size_t) (nt - first_tree));
    }
    else {
      qarray = NULL;
    }

    /* identify boundary quadrants and prepare them to be sent */
    for (zz = 0; zz < treecount; ++zz) {
      /* this quadrant may be on the boundary with a range of processors */
      q = p4est_quadrant_array_index (tquadrants, zz);
      qh = P4EST_QUADRANT_LEN (q->level);
      if (p4est_comm_neighborhood_owned (p4est, nt,
                                         full_tree, tree_contact, q)) {
        /* this quadrant's 3x3 neighborhood is owned by this processor */
        ++skipped;
        continue;
      }

      if (qarray != NULL) {
        s = (p4est_quadrant_t *) sc_array_push (qarray);
        *s = *q;
      }

#ifdef P4_TO_P8
      for (m = 0; m < 3; ++m) {
#if 0
      }
#endif
#else
      m = 0;
#endif
      for (k = 0; k < 3; ++k) {
        for (l = 0; l < 3; ++l) {
          which = m * 9 + k * 3 + l;    /* 2D: 0..8, 3D: 0..26 */
          /* exclude myself from the queries */
          if (which == P4EST_INSUL / 2) {
            continue;
          }
          /* may modify insulq below, never modify q itself! */
          insulq = *q;
          insulq.x += (l - 1) * qh;
          insulq.y += (k - 1) * qh;
#ifdef P4_TO_P8
          insulq.z += (m - 1) * qh;
#endif
          /* check boundary status of insulation quadrant */
          quad_contact[0] = (insulq.x < 0);
          quad_contact[1] = (insulq.x >= rh);
          face_axis[0] = quad_contact[0] || quad_contact[1];
          quad_contact[2] = (insulq.y < 0);
          quad_contact[3] = (insulq.y >= rh);
          face_axis[1] = quad_contact[2] || quad_contact[3];
#ifndef P4_TO_P8
          face_axis[2] = 0;
#else
          quad_contact[4] = (insulq.z < 0);
          quad_contact[5] = (insulq.z >= rh);
          face_axis[2] = quad_contact[4] || quad_contact[5];
          edge = -1;
#endif
          contact_edge_only = contact_face_only = 0;
          face = -1;
          if (face_axis[0] || face_axis[1] || face_axis[2]) {
            /* this quadrant is relevant for inter-tree balancing */
            if (!face_axis[1] && !face_axis[2]) {
              contact_face_only = 1;
              face = 0 + quad_contact[1];
#ifdef P4EST_ENABLE_DEBUG
              facewhich = P4EST_INSUL / 2 + (2 * quad_contact[1] - 1);
#endif
            }
            else if (!face_axis[0] && !face_axis[2]) {
              contact_face_only = 1;
              face = 2 + quad_contact[3];
#ifdef P4EST_ENABLE_DEBUG
              facewhich = P4EST_INSUL / 2 + 3 * (2 * quad_contact[3] - 1);
#endif
            }
#ifdef P4_TO_P8
            else if (!face_axis[0] && !face_axis[1]) {
              contact_face_only = 1;
              face = 4 + quad_contact[5];
#ifdef P4EST_ENABLE_DEBUG
              facewhich = P4EST_INSUL / 2 + 9 * (2 * quad_contact[5] - 1);
#endif
            }
            else if (!face_axis[0]) {
              contact_edge_only = 1;
              edge = 0 + 2 * quad_contact[5] + quad_contact[3];
            }
            else if (!face_axis[1]) {
              contact_edge_only = 1;
              edge = 4 + 2 * quad_contact[5] + quad_contact[1];
            }
            else if (!face_axis[2]) {
              contact_edge_only = 1;
              edge = 8 + 2 * quad_contact[3] + quad_contact[1];
            }
#endif
            if (contact_face_only) {
              /* square contact across a face */
              P4EST_ASSERT (!contact_edge_only);
              P4EST_ASSERT (face >= 0 && face < P4EST_FACES);
              P4EST_ASSERT (quad_contact[face]);
              qtree = p4est_find_face_transform (conn, nt, face, ftransform);
              if (qtree >= 0) {
                P4EST_ASSERT (tree_contact[facewhich]);
                p4est_quadrant_transform_face (q, &tosend, ftransform);
                tosend.p.piggy2.from_tree = nt;
                tosend.pad16 = face;
                p4est_quadrant_transform_face (&insulq, &tempq, ftransform);
                p4est_balance_schedule (p4est, peers, qtree, 1,
                                        &tosend, &tempq,
                                        &first_peer, &last_peer);
              }
              else {
                /* goes across a face with no neighbor */
                P4EST_ASSERT (!tree_contact[facewhich]);
              }
            }
#ifdef P4_TO_P8
            else if (contact_edge_only) {
              /* this quadrant crosses an edge */
              P4EST_ASSERT (!contact_face_only);
              P4EST_ASSERT (edge >= 0 && edge < P8EST_EDGES);
              p8est_find_edge_transform (conn, nt, edge, &ei);
              for (etree = 0; etree < eta->elem_count; ++etree) {
                et = p8est_edge_array_index (eta, etree);
                p8est_quadrant_transform_edge (q, &tosend, &ei, et, 0);
                tosend.p.piggy2.from_tree = nt;
                tosend.pad16 = edge;
                p8est_quadrant_transform_edge (&insulq, &tempq, &ei, et, 1);
                p4est_balance_schedule (p4est, peers, et->ntree, 1,
                                        &tosend, &tempq,
                                        &first_peer, &last_peer);
              }
            }
#endif
            else {
              /* this quadrant crosses a corner */
              P4EST_ASSERT (face_axis[0] && face_axis[1]);
              corner = quad_contact[1] + 2 * quad_contact[3];
#ifdef P4_TO_P8
              P4EST_ASSERT (face_axis[2]);
              corner += 4 * quad_contact[5];
#endif
              P4EST_ASSERT (p4est_quadrant_touches_corner (q, corner, 1));
              P4EST_ASSERT (p4est_quadrant_touches_corner
                            (&insulq, corner, 0));
              p4est_find_corner_transform (conn, nt, corner, &ci);
              for (ctree = 0; ctree < cta->elem_count; ++ctree) {
                ct = p4est_corner_array_index (cta, ctree);
                tosend = *q;
                p4est_quadrant_transform_corner (&tosend, (int) ct->ncorner,
                                                 0);
                tosend.p.piggy2.from_tree = nt;
                tosend.pad16 = corner;
                tempq = insulq;
                p4est_quadrant_transform_corner (&tempq, (int) ct->ncorner,
                                                 1);
                p4est_balance_schedule (p4est, peers, ct->ntree, 1,
                                        &tosend, &tempq, &first_peer,
                                        &last_peer);
              }
            }
          }
          else {
            /* no inter-tree contact */
            tosend = *q;
            tosend.p.piggy2.from_tree = nt;
            tosend.pad16 = -1;
            p4est_balance_schedule (p4est, peers, nt, 0,
                                    &tosend, &insulq, &first_peer,
                                    &last_peer);
          }
        }
      }
#ifdef P4_TO_P8
#if 0
      {
#endif
      }
#endif
    }
    tquadrants = NULL;          /* safeguard */
  }

  /* end balance_A, start balance_comm */
  if (bobj->inspect != NULL) {
    bobj->inspect->balance_A += sc_MPI_Wtime ();
    bobj->inspect->balance_comm = -sc_MPI_Wtime ();
    bobj->inspect->balance_comm_sent = 0;
    bobj->inspect->balance_comm_nzpeers = 0;
    for (k = 0; k < 2; ++k) {
      bobj->inspect->balance_load_sends[k] = 0;
      bobj->inspect->balance_load_receives[k] = 0;
      bobj->inspect->balance_zero_sends[k] = 0;
      bobj->inspect->balance_zero_receives[k] = 0;
    }
    notify = bobj->inspect->notify;
  }
  if (!notify) {
    notify = sc_notify_new (bobj->mpicomm);
    own_notify = 1;
  }

#ifdef P4EST_ENABLE_MPI
  /* encode and distribute the asymmetric communication pattern */
  receivers = sc_array_new (sizeof (int));
  senders = sc_array_new (sizeof (int));
  in_counts = sc_array_new (sizeof (int));
  out_counts = sc_array_new (sizeof (int));
#ifdef P4EST_ENABLE_DEBUG
  for (j = 0; j < first_peer; j++) {
    P4EST_ASSERT (peers[j].send_first.elem_count == 0);
  }
#endif
  for (num_receivers = 0, j = first_peer; j <= last_peer; j++) {
    if (j == rank) {
      continue;
    }
    if (peers[j].send_first.elem_count) {
      num_receivers++;
    }
  }
#ifdef P4EST_ENABLE_DEBUG
  for (j = last_peer + 1; j < num_procs; j++) {
    P4EST_ASSERT (peers[j].send_first.elem_count == 0);
  }
#endif
  sc_array_resize (receivers, (size_t) num_receivers);
  sc_array_resize (in_counts, (size_t) num_receivers);
  receiver_ranks = (int *) receivers->array;
  icounts = (int *) in_counts->array;
  for (num_receivers = 0, j = first_peer; j <= last_peer; j++) {
    if (j == rank) {
      continue;
    }
    if (peers[j].send_first.elem_count) {
      receiver_ranks[num_receivers] = j;
      icounts[num_receivers++] = (int) peers[j].send_first.elem_count;
    }
  }

  if (bobj->inspect && bobj->inspect->stats) {
    sc_notify_set_stats (notify, bobj->inspect->stats);
  }
  sc_notify_payload (receivers, senders, in_counts, out_counts, 0, notify);
  sc_array_destroy (in_counts);
  sender_ranks = (int *) senders->array;
  num_senders = (int) senders->elem_count;

  /*
   * loop over all peers and send first round of quadrants
   * for intra-tree balancing, each load is contained in one tree
   */
  total_send_count = total_recv_count = 0;
  request_first_send = 0;
  request_first_recv = 0;
  request_second_send = 0;
  send_zero[0] = send_load[0] = recv_zero[0] = recv_load[0] = 0;
  send_zero[1] = send_load[1] = recv_zero[1] = recv_load[1] = 0;

  /* Use sender_ranks array to recv from them */
  icounts = (int *) out_counts->array;
  for (k = 0; k < num_senders; k++) {
    j = sender_ranks[k];
    qcount = icounts[k];
    P4EST_ASSERT (j >= 0 && j < num_procs && j != rank);
    peer = peers + j;
    sc_array_resize (&(peer->recv_first), (size_t) qcount);
    /* first send number of quadrants to be expected */
    if (qcount > 0) {
      P4EST_LDEBUGF ("Balance A recv %llu quadrants from %d\n",
                     (unsigned long long) qcount, j);
      total_recv_count += qcount;
      qbytes = qcount * sizeof (p4est_quadrant_t);
      peer->recv_first_count = qcount;
      mpiret = MPI_Irecv (peer->recv_first.array, (int) qbytes, MPI_BYTE,
                          j, P4EST_COMM_BALANCE_FIRST_LOAD,
                          p4est->mpicomm, &recv_requests_first_load[j]);
      SC_CHECK_MPI (mpiret);
      ++recv_load[0];
      request_first_recv++;
    }
    else {
      ++recv_zero[0];
    }
  }
  sc_array_destroy (out_counts);
  /* Use receiver_ranks array to send to them */
  for (k = 0; k < num_receivers; ++k) {
    j = receiver_ranks[k];
    P4EST_ASSERT (j >= first_peer && j <= last_peer && j != rank);
    peer = peers + j;
    qcount = peer->send_first.elem_count;

    /* first send number of quadrants to be expected */
    if (qcount > 0) {
      P4EST_LDEBUGF ("Balance A send %llu quadrants to %d\n",
                     (unsigned long long) qcount, j);
      ++send_load[0];
    }
    else {
      ++send_zero[0];
    }
    peer->send_first_count = (int) qcount;
    /* sort and send the actual quadrants and post receive for reply */
    if (qcount > 0) {
      sc_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);

#ifdef P4EST_ENABLE_DEBUG
      checksum = p4est_quadrant_checksum (&peer->send_first, &checkarray, 0);
      P4EST_LDEBUGF ("Balance A send checksum 0x%08x to %d\n", checksum, j);
#endif /* P4EST_ENABLE_DEBUG */

      total_send_count += qcount;
      qbytes = qcount * sizeof (p4est_quadrant_t);
      mpiret = MPI_Isend (peer->send_first.array, (int) qbytes, MPI_BYTE,
                          j, P4EST_COMM_BALANCE_FIRST_LOAD,
                          p4est->mpicomm, &send_requests_first_load[j]);
      SC_CHECK_MPI (mpiret);
      ++request_first_send;
    }
  }
  peer = NULL;

  /* wait for quadrant counts and post receive and send for quadrants */
  while (request_first_recv > 0) {
    mpiret = MPI_Waitsome (num_procs, recv_requests_first_load,
                           &outcount, wait_indices, recv_statuses);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount != MPI_UNDEFINED);
    P4EST_ASSERT (outcount > 0);
    for (i = 0; i < outcount; ++i) {
      /* retrieve sender's rank */
      j = wait_indices[i];
      jstatus = &recv_statuses[i];
      wait_indices[i] = -1;
      P4EST_ASSERT (j != rank && 0 <= j && j < num_procs);
      P4EST_ASSERT (jstatus->MPI_SOURCE == j);

      /* check if we are in receiving count or load */
      peer = peers + j;
      /* verify received size */
      P4EST_ASSERT (jstatus->MPI_TAG == P4EST_COMM_BALANCE_FIRST_LOAD);
      P4EST_ASSERT (peer->recv_first_count > 0);
      mpiret = MPI_Get_count (jstatus, MPI_BYTE, &rcount);
      SC_CHECK_MPI (mpiret);
      SC_CHECK_ABORTF (rcount ==
                       peer->recv_first_count *
                       (int) sizeof (p4est_quadrant_t),
                       "Receive load mismatch A %d %dx%llu", rcount,
                       peer->recv_first_count,
                       (unsigned long long) sizeof (p4est_quadrant_t));

      /* received load, close this request */
      peer->have_first_load = 1;
      --request_first_recv;

#ifdef P4EST_ENABLE_DEBUG
      checksum = p4est_quadrant_checksum (&peer->recv_first, &checkarray, 0);
      P4EST_LDEBUGF ("Balance A recv checksum 0x%08x from %d\n", checksum, j);
#endif /* P4EST_ENABLE_DEBUG */

      /* process incoming quadrants to interleave with communication */
      p4est_balance_response (p4est, peer, btype, borders, bobj);
      qcount = peer->send_second.elem_count;
      if (qcount > 0) {
        P4EST_LDEBUGF ("Balance B send %llu quadrants to %d\n",
                       (unsigned long long) qcount, j);
        ++send_load[1];
      }
      else {
        ++send_zero[1];
      }
      peer->send_second_count = (int) qcount;
      if (qcount > 0) {
#ifdef P4EST_ENABLE_DEBUG
        checksum =
          p4est_quadrant_checksum (&peer->send_second, &checkarray, 0);
        P4EST_LDEBUGF ("Balance B send checksum 0x%08x to %d\n", checksum, j);
#endif /* P4EST_ENABLE_DEBUG */
      }
      total_send_count += qcount;
      qbytes = qcount * sizeof (p4est_quadrant_t);
      mpiret = MPI_Isend (peer->send_second.array, (int) qbytes, MPI_BYTE,
                          j, P4EST_COMM_BALANCE_SECOND_LOAD,
                          p4est->mpicomm, &send_requests_second_load[j]);
      SC_CHECK_MPI (mpiret);
      ++request_second_send;
    }
  }
  for (j = 0; j < num_procs; ++j) {
    /* we cannot have exited this loop unless the first recv is complete */
    P4EST_ASSERT (recv_requests_first_load[j] == MPI_REQUEST_NULL);
  }
#endif /* P4EST_ENABLE_MPI */

  /* simulate send and receive with myself across tree boundaries */
  peer = peers + rank;
  sc_array_sort (&peer->send_first, p4est_quadrant_compare_piggy);
  qcount = peer->send_first.elem_count;
  peer->recv_first_count = peer->send_first_count = (int) qcount;
  qbytes = qcount * sizeof (p4est_quadrant_t);
  qarray = &peer->recv_first;
  sc_array_resize (qarray, qcount);
  memcpy (qarray->array, peer->send_first.array, qbytes);
  p4est_balance_response (p4est, peer, btype, borders, bobj);
  qcount = peer->send_second.elem_count;
  peer->recv_second_count = peer->send_second_count = (int) qcount;
  qbytes = qcount * sizeof (p4est_quadrant_t);
  qarray = &peer->recv_second;
  sc_array_resize (qarray, qcount);
  memcpy (qarray->array, peer->send_second.array, qbytes);

#ifdef P4EST_ENABLE_MPI
  /* receive second round appending to the same receive buffer */
  for (k = 0; k < request_first_send; k++) {
    MPI_Status          status;

    mpiret =
      MPI_Probe (MPI_ANY_SOURCE, P4EST_COMM_BALANCE_SECOND_LOAD,
                 p4est->mpicomm, &status);
    SC_CHECK_MPI (mpiret);
    j = status.MPI_SOURCE;
    peer = peers + j;
    mpiret = MPI_Get_count (&status, MPI_BYTE, &rcount);
    qcount = rcount / sizeof (p4est_quadrant_t);
    sc_array_resize (&(peer->recv_second), (size_t) qcount);
    total_recv_count += qcount;
    if (qcount) {
      recv_load[1]++;
    }
    else {
      recv_zero[1]++;
    }
    peer->recv_second_count = qcount;
    mpiret = MPI_Irecv (peer->recv_second.array, rcount, MPI_BYTE,
                        j, P4EST_COMM_BALANCE_SECOND_LOAD, p4est->mpicomm,
                        &recv_requests_second_load[j]);
    SC_CHECK_MPI (mpiret);
  }
  /* we cannot have initiated the second receive unless the first send is
   * complete */
  /* wait for the second recv to finish */
  mpiret = MPI_Waitall (num_procs, send_requests_first_load,
                        MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  /* wait for the second recv to finish */
  mpiret = MPI_Waitall (num_procs, recv_requests_second_load,
                        MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  for (j = 0; j < num_procs; ++j) {
    P4EST_ASSERT (recv_requests_second_load[j] == MPI_REQUEST_NULL);
#ifdef P4EST_ENABLE_DEBUG
    peer = peers + j;
    if (peer->recv_second.elem_count) {
      checksum = p4est_quadrant_checksum (&peer->recv_second, &checkarray, 0);
      P4EST_LDEBUGF ("Balance B recv checksum 0x%08x from %d\n", checksum, j);
    }
#endif
  }

  /* print buffer statistics */
  P4EST_VERBOSEF ("first send Z %d L %d recv Z %d L %d\n",
                  send_zero[0], send_load[0], recv_zero[0], recv_load[0]);
  P4EST_VERBOSEF ("second send Z %d L %d recv Z %d L %d\n",
                  send_zero[1], send_load[1], recv_zero[1], recv_load[1]);
  P4EST_VERBOSEF ("total send %d recv %d\n", total_send_count,
                  total_recv_count);
  for (j = 0; j < num_procs; ++j) {
    peer = peers + j;
    if (peer->send_first_count > 0 || peer->recv_first_count > 0 ||
        peer->send_second_count > 0 || peer->recv_second_count > 0) {
      P4EST_VERBOSEF ("peer %d first S %llu R %d second S %llu R %d\n",
                      j, (unsigned long long) peer->send_first_count,
                      peer->recv_first_count,
                      (unsigned long long) peer->send_second_count,
                      peer->recv_second_count);
    }
  }
#endif /* P4EST_ENABLE_MPI */

  /* end balance_comm, start balance_B */
  if (bobj->inspect != NULL) {
    bobj->inspect->balance_comm += sc_MPI_Wtime ();
    bobj->inspect->balance_B = -sc_MPI_Wtime ();
    bobj->inspect->balance_B_count_in = 0;
    bobj->inspect->balance_B_count_out = 0;
    bobj->inspect->use_B = 1;
#ifdef P4EST_ENABLE_MPI
    for (k = 0; k < 2; ++k) {
      bobj->inspect->balance_load_sends[k] = send_load[k];
      bobj->inspect->balance_load_receives[k] = recv_load[k];
      bobj->inspect->balance_zero_sends[k] = send_zero[k];
      bobj->inspect->balance_zero_receives[k] = recv_zero[k];
    }
#endif
  }
  if (own_notify) {
    sc_notify_destroy (notify);
  }

  /* merge received quadrants */
  for (j = 0; j < num_procs; ++j) {
    size_t              fcount;

    /* access peer information */
    peer = peers + j;
    fcount = peer->recv_first.elem_count;
    qcount = fcount + peer->recv_second.elem_count;
    P4EST_ASSERT (peer->send_first_count ==
                  (int) peer->send_first.elem_count);
    P4EST_ASSERT (peer->send_second_count ==
                  (int) peer->send_second.elem_count);
    P4EST_ASSERT (peer->recv_second_count ==
                  (int) peer->recv_second.elem_count);
    if (qcount == 0) {
      continue;
    }

    /* merge received quadrants into correct tree */
    for (zz = 0; zz < qcount; ++zz) {
      s = zz < fcount ? p4est_quadrant_array_index (&peer->recv_first, zz) :
        p4est_quadrant_array_index (&peer->recv_second, zz - fcount);
      P4EST_ASSERT (p4est_quadrant_is_extended (s));
      qtree = s->p.piggy2.which_tree;
      if (qtree < first_tree || qtree > last_tree) {
        /* this is a corner/edge quadrant from the second pass of balance */
        continue;
      }
      if (borders == NULL) {
        tree = p4est_tree_array_index (p4est->trees, qtree);
        q = p4est_quadrant_array_push (&tree->quadrants);
        *q = *s;
        ++tree->quadrants_per_level[q->level];
        tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q->level);
        ++p4est->local_num_quadrants;
        p4est_quadrant_init_data (p4est, qtree, q, init_fn);
      }
      else {
        qarray = (sc_array_t *) sc_array_index (borders,
                                                (int) (qtree - first_tree));
        q = p4est_quadrant_array_push (qarray);
        *q = *s;
      }
    }
  }

  /* rebalance and clamp result back to original tree boundaries */
  p4est->local_num_quadrants = 0;
  for (nt = first_tree; nt <= last_tree; ++nt) {
    /* check if we are the only processor in an isolated tree */
    tree = p4est_tree_array_index (p4est->trees, nt);
    tree->quadrants_offset = p4est->local_num_quadrants;
    tquadrants = &tree->quadrants;
    treecount = tquadrants->elem_count;
    if (!(tree_flags[nt] & fully_owned_flag) ||
        (tree_flags[nt] & any_face_flag)) {
      /* we have most probably received quadrants, run sort and balance */
      /* balance the border, add it back into the tree, and linearize */
      p4est_balance_border (p4est, btype, nt, init_fn, replace_fn, borders);
      P4EST_VERBOSEF ("Balance tree %lld B %llu to %llu\n",
                      (long long) nt,
                      (unsigned long long) treecount,
                      (unsigned long long) tquadrants->elem_count);
    }
    p4est->local_num_quadrants += tquadrants->elem_count;
    tquadrants = NULL;          /* safeguard */
  }
  if (last_tree >= 0) {
    for (; nt < conn->num_trees; ++nt) {
      tree = p4est_tree_array_index (p4est->trees, nt);
      tree->quadrants_offset = p4est->local_num_quadrants;
    }
  }

  /* end balance_B */
  if (bobj->inspect != NULL) {
    bobj->inspect->balance_B += sc_MPI_Wtime ();
  }

#ifdef P4EST_ENABLE_MPI
  /* wait for all send operations */
  if (request_second_send > 0) {
    mpiret = MPI_Waitall (num_procs,
                          send_requests_second_load, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }

  /* compute global sum of send and receive counts */
#ifdef P4EST_ENABLE_DEBUG
  gtotal[0] = gtotal[1] = 0;
  ltotal[0] = (p4est_gloidx_t) total_send_count;
  ltotal[1] = (p4est_gloidx_t) total_recv_count;
  mpiret = MPI_Reduce (ltotal, gtotal, 2, P4EST_MPI_GLOIDX,
                       MPI_SUM, 0, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_STATISTICSF ("Global number of shipped quadrants %lld\n",
                            (long long) gtotal[0]);
  P4EST_ASSERT (rank != 0 || gtotal[0] == gtotal[1]);
#endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */

  /* loop over all local trees to finalize balance */
  all_outcount = 0;
  for (nt = first_tree; nt <= last_tree; ++nt) {
    tree = p4est_tree_array_index (p4est->trees, nt);
    all_outcount += tree->quadrants.elem_count;

    /* final log message for this tree */
    P4EST_VERBOSEF ("Done balance tree %lld now %llu\n", (long long) nt,
                    (unsigned long long) tree->quadrants.elem_count);
  }

  /* cleanup temporary storage */
  P4EST_FREE (tree_flags);
  for (j = 0; j < num_procs; ++j) {
    peer = peers + j;
    sc_array_reset (&peer->send_first);
    sc_array_reset (&peer->send_second);
    sc_array_reset (&peer->recv_first);
    sc_array_reset (&peer->recv_second);
  }
  P4EST_FREE (peers);

  if (borders != NULL) {
    for (zz = 0; zz < localcount; zz++) {
      qarray = (sc_array_t *) sc_array_index (borders, zz);
      sc_array_reset (qarray);
    }
    sc_array_destroy (borders);
  }

#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
  sc_array_reset (cta);

#ifdef P4EST_ENABLE_MPI
  sc_array_destroy (receivers);
  sc_array_destroy (senders);
  P4EST_FREE (send_requests_first_load);
  P4EST_FREE (recv_statuses);
  P4EST_FREE (wait_indices);
#ifdef P4EST_ENABLE_DEBUG
  sc_array_reset (&checkarray);
#endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */

  /* some sanity checks */
  P4EST_ASSERT ((p4est_locidx_t) all_outcount == p4est->local_num_quadrants);
  P4EST_ASSERT (pre_adapt_flags || all_outcount >= all_incount);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + all_outcount - all_incount ==
                  p4est->user_data_pool->elem_count);
  }
  P4EST_VERBOSEF ("Balance skipped %lld\n", (long long) skipped);

  P4EST_FUNC_SHOT (p4est, &snap);
}


