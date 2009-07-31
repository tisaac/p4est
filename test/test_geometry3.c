/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007-2009 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <p8est_geometry.h>

static void
test_tree (p8est_geometry_t * geom, p4est_topidx_t which_tree,
           const double bounds[6], int N)
{
  const double        epsilon = 1e-8;
  int                 i, j, k, l, m;
  double              h[3];
  double              xyz[4][3], XYZ[4][3];
  double              Jgeom[3][3], Jdisc[3][3];
  double              detD, detJ, fd, fg;
  double              diffD, diffJ, maxJ;

  for (l = 0; l < 3; ++l) {
    h[l] = bounds[2 * l + 1] - bounds[2 * l];
  }

  for (i = 0; i < N; ++i) {
    xyz[3][2] = ((N - (i + .5)) * bounds[4] + (i + .5) * bounds[5]) / N;
    for (j = 0; j < N; ++j) {
      xyz[3][1] = ((N - (j + .5)) * bounds[2] + (j + .5) * bounds[3]) / N;
      for (k = 0; k < N; ++k) {
        xyz[3][0] = ((N - (k + .5)) * bounds[0] + (k + .5) * bounds[1]) / N;

        /* compute transformed point */
        geom->X (geom, which_tree, xyz[3], XYZ[3]);

        /* offset point three times in the coordinate directions */
        for (l = 0; l < 3; ++l) {       /* l runs in domain */
          memcpy (xyz[l], xyz[3], 3 * sizeof (double));
          xyz[l][l] += epsilon * h[l];
          geom->X (geom, which_tree, xyz[l], XYZ[l]);

          for (m = 0; m < 3; ++m) {     /* m runs in image domain */
            Jdisc[m][l] = (XYZ[l][m] - XYZ[3][m]) / (epsilon * h[l]);
          }
        }
        detJ = geom->J (geom, which_tree, xyz[3], Jgeom);
        detD = geom->D (geom, which_tree, xyz[3]);

        /* compare results */
        diffD = fabs (detD - detJ) / SC_MAX (detD, detJ);
        SC_CHECK_ABORTF (diffD < 1e-8, "Determinant mismatch %lld %g %g %g",
                         (long long) which_tree, xyz[3][0], xyz[3][1],
                         xyz[3][2]);

        diffJ = maxJ = 0.;
        for (l = 0; l < 3; ++l) {
          for (m = 0; m < 3; ++m) {
            fd = fabs (Jdisc[m][l]);
            fg = fabs (Jgeom[m][l]);
            maxJ = SC_MAX (maxJ, fd);
            maxJ = SC_MAX (maxJ, fg);
            diffJ += fabs (Jdisc[m][l] - Jgeom[m][l]);
          }
        }
        diffJ /= 3 * 3 * maxJ;
        SC_CHECK_ABORTF (diffJ < 100. * epsilon,
                         "Jacobian mismatch %lld %g %g %g",
                         (long long) which_tree, xyz[3][0], xyz[3][1],
                         xyz[3][2]);
      }
    }
  }
}

static void
test_identity (int N)
{
  const double        bounds[6] = { -1, 1, -1, 1, -1, 1 };
  p8est_geometry_t   *geom;

  P4EST_STATISTICS ("Test identity\n");

  geom = p8est_geometry_new_identity ();
  test_tree (geom, 0, bounds, N);

  P4EST_FREE (geom);
}

static void
test_shell (int N)
{
  const double        bounds[6] = { -1, 1, -1, 1, 1, 2 };
  p4est_topidx_t      jt;
  p8est_geometry_t   *geom;

  P4EST_STATISTICS ("Test shell\n");

  geom = p8est_geometry_new_shell (6371, 3480);
  for (jt = 0; jt < 24; ++jt) {
    test_tree (geom, jt, bounds, N);
  }

  P4EST_FREE (geom);
}

static void
test_sphere (int N)
{
  /* *INDENT-OFF* */
  const double        boundsS[6] = { -1,  1, -1,  1,  1,  2 };
  const double        boundsC[6] = { -1,  1, -1,  1, -1,  1 };
  /* *INDENT-ON* */
  p4est_topidx_t      jt;
  p8est_geometry_t   *geom;

  P4EST_STATISTICS ("Test sphere\n");

  geom = p8est_geometry_new_sphere (3, 2, 1);
  for (jt = 0; jt < 13; ++jt) {
    test_tree (geom, jt, jt < 12 ? boundsS : boundsC, N);
  }

  P4EST_FREE (geom);
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 N;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;

  sc_init (mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  N = 43;
  test_identity (N);
  test_shell (N);
  test_sphere (N);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}