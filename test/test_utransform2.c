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
#include <p4est_bits.h>
#include <p4est_connectivity.h>
#else
#include <p8est_bits.h>
#include <p8est_connectivity.h>
#endif

int
main (int argc, char **argv)
{
  const char         *names[] = {
#ifndef P4_TO_P8
    "brick23",
    "corner",
    "cubed",
    "disk",
    "moebius",
    "periodic",
    "pillow",
    "rotwrap",
    "star",
    "unit",
#else
    "brick235",
    "periodic",
    "rotcubes",
    "rotwrap",
    "shell",
    "sphere",
    "twocubes",
    "twowrap",
    "unit",
    "edge",
#endif
  };
  int                 i, num_names;

  sc_init (sc_MPI_COMM_NULL, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  num_names = sizeof (names) / sizeof (names[0]);

  for (i = 0; i < num_names; i++) {
    p4est_connectivity_t *conn;
    p4est_topidx_t      j, num_trees;
    p4est_quadrant_t    root, desc, neigh, tran, utran, utran2;

    conn = p4est_connectivity_new_byname (names[i]);
    num_trees = conn->num_trees;

    P4EST_QUADRANT_INIT (&root);
    P4EST_QUADRANT_INIT (&desc);
    P4EST_QUADRANT_INIT (&neigh);
    P4EST_QUADRANT_INIT (&tran);
    P4EST_QUADRANT_INIT (&utran);
    P4EST_QUADRANT_INIT (&utran2);
    root.x = 0;
    root.y = 0;
#ifdef P4_TO_P8
    root.z = 0;
#endif
    root.level = 0;

    for (j = 0; j < num_trees; j++) {
      int                 k, l;
      int                 utransform[P4EST_UTRANSFORM];

      for (k = 0; k < P4EST_FACES; k++) {
        int                 ftransform[P4EST_FTRANSFORM];
        p4est_topidx_t      nt;

        nt = p4est_find_face_transform (conn, j, k, ftransform);
        if (nt < 0) {
          continue;
        }

        p4est_face_transform_to_utransform (ftransform, utransform);

        for (l = 0; l < P4EST_HALF; l++) {
          int                 c = p4est_face_corners[k][l];

          p4est_quadrant_corner_descendant (&root, &desc, c, P4EST_QMAXLEVEL);
          p4est_quadrant_transform_face (&desc, &tran, ftransform);
          p4est_quadrant_utransform (&desc, &utran, utransform, 0);
          SC_CHECK_ABORTF (p4est_quadrant_is_equal (&tran, &utran),
                           "Fail conn %s, tree %d, face %d, corner %d\n",
                           names[i], (int) j, k, l);
          p4est_quadrant_utransform (&utran, &utran2, utransform, 1);
          SC_CHECK_ABORTF (p4est_quadrant_is_equal (&utran2, &desc),
                           "Fail conn %s, tree %d, face %d, corner %d\n",
                           names[i], (int) j, k, l);

          p4est_quadrant_face_neighbor (&desc, k, &neigh);
          p4est_quadrant_transform_face (&neigh, &tran, ftransform);
          p4est_quadrant_utransform (&neigh, &utran, utransform, 0);
          SC_CHECK_ABORTF (p4est_quadrant_is_equal (&tran, &utran),
                           "Fail conn %s, tree %d, face %d, corner %d\n",
                           names[i], (int) j, k, l);
          p4est_quadrant_utransform (&utran, &utran2, utransform, 1);
          SC_CHECK_ABORTF (p4est_quadrant_is_equal (&utran2, &neigh),
                           "Fail conn %s, tree %d, face %d, corner %d\n",
                           names[i], (int) j, k, l);
        }
      }
#ifdef P4_TO_P8
      for (k = 0; k < P8EST_EDGES; k++) {
        p8est_edge_info_t   ei;
        sc_array_t         *eta = &ei.edge_transforms;
        p8est_edge_transform_t *et;
        size_t              n, nneigh;

        sc_array_init (eta, sizeof (p8est_edge_transform_t));

        p8est_find_edge_transform (conn, j, k, &ei);

        nneigh = eta->elem_count;
        et = (p8est_edge_transform_t *) eta->array;

        for (n = 0; n < nneigh; n++) {
          p8est_edge_transform_to_utransform (&et[n], k, utransform);
          for (l = 0; l < 2; l++) {
            int                 c = p8est_edge_corners[k][l];

            p4est_quadrant_corner_descendant (&root, &desc, c,
                                              P4EST_QMAXLEVEL);
            p8est_quadrant_transform_edge (&desc, &tran, &ei, &et[n], 0);
            p4est_quadrant_utransform (&desc, &utran, utransform, 0);
            SC_CHECK_ABORTF (p4est_quadrant_is_equal (&tran, &utran),
                             "Fail conn %s, tree %d, edge %d, corner %d\n",
                             names[i], (int) j, k, l);
            p4est_quadrant_utransform (&utran, &utran2, utransform, 1);
            SC_CHECK_ABORTF (p4est_quadrant_is_equal (&utran2, &desc),
                             "Fail conn %s, tree %d, edge %d, corner %d\n",
                             names[i], (int) j, k, l);

            p8est_quadrant_edge_neighbor (&desc, k, &neigh);
            p8est_quadrant_transform_edge (&neigh, &tran, &ei, &et[n], 1);
            p4est_quadrant_utransform (&neigh, &utran, utransform, 0);
            SC_CHECK_ABORTF (p4est_quadrant_is_equal (&tran, &utran),
                             "Fail conn %s, tree %d, edge %d, corner %d\n",
                             names[i], (int) j, k, l);
            p4est_quadrant_utransform (&utran, &utran2, utransform, 1);
            SC_CHECK_ABORTF (p4est_quadrant_is_equal (&utran2, &neigh),
                             "Fail conn %s, tree %d, edge %d, corner %d\n",
                             names[i], (int) j, k, l);
          }
        }

        sc_array_reset (eta);
      }
#endif
      for (k = 0; k < P4EST_CHILDREN; k++) {
        p4est_corner_info_t ci;
        sc_array_t         *cta = &ci.corner_transforms;
        p4est_corner_transform_t *ct;
        size_t              n, nneigh;

        sc_array_init (cta, sizeof (p4est_corner_transform_t));

        p4est_find_corner_transform (conn, j, k, &ci);

        nneigh = cta->elem_count;
        ct = (p4est_corner_transform_t *) cta->array;

        for (n = 0; n < nneigh; n++) {
          p4est_corner_transform_to_utransform (&ct[n], k, utransform);

          p4est_quadrant_corner_descendant (&root, &desc, k, P4EST_QMAXLEVEL);
          tran = desc;
          p4est_quadrant_transform_corner (&tran, ct[n].ncorner, 0);
          p4est_quadrant_utransform (&desc, &utran, utransform, 0);
          SC_CHECK_ABORTF (p4est_quadrant_is_equal (&tran, &utran),
                           "Fail conn %s, tree %d, corner %d, corner %d\n",
                           names[i], (int) j, k, l);
          p4est_quadrant_utransform (&utran, &utran2, utransform, 1);
          SC_CHECK_ABORTF (p4est_quadrant_is_equal (&utran2, &desc),
                           "Fail conn %s, tree %d, corner %d, corner %d\n",
                           names[i], (int) j, k, l);

          p4est_quadrant_corner_neighbor (&desc, k, &neigh);
          tran = neigh;
          p4est_quadrant_transform_corner (&tran, ct[n].ncorner, 1);
          p4est_quadrant_utransform (&neigh, &utran, utransform, 0);
          SC_CHECK_ABORTF (p4est_quadrant_is_equal (&tran, &utran),
                           "Fail conn %s, tree %d, corner %d\n", names[i],
                           (int) j, k);
          p4est_quadrant_utransform (&utran, &utran2, utransform, 1);
          SC_CHECK_ABORTF (p4est_quadrant_is_equal (&utran2, &neigh),
                           "Fail conn %s, tree %d, corner %d\n", names[i],
                           (int) j, k);
        }

        sc_array_reset (cta);
      }
    }

    p4est_connectivity_destroy (conn);
  }

  sc_finalize ();

  return 0;
}
