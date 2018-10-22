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
#include <p4est_balance_obj.h>
#else
#include <p8est_balance_obj.h>
#endif

struct p4est_balance_obj_s
{
  sc_MPI_Comm            comm;
  p4est_balance_method_t method;
  sc_statistics_t        *stats;
  p4est_connect_type_t   connect;
  p4est_init_t           init_fn;
  p4est_replace_t        replace_fn;
};

const char *p4est_balance_method_strings[] = {
  P4EST_BALANCE_STR_SORT,
  P4EST_BALANCE_STR_TWOROUND,
};

p4est_balance_obj_t *
p4est_balance_obj_new (sc_MPI_Comm mpicomm)
{
  p4est_balance_obj_t *bal;

  bal = SC_ALLOC_ZERO (p4est_balance_obj_t, 1);

  bal->comm = mpicomm;
  bal->method = P4EST_BALANCE_DEFAULT;
  bal->connect = P4EST_CONNECT_FULL;
  return bal;
}

void
p4est_balance_obj_destroy (p4est_balance_obj_t *obj)
{
  SC_FREE (obj);
}

sc_MPI_Comm
p4est_balance_obj_get_comm (p4est_balance_obj_t *obj)
{
  return obj->comm;
}

void
p4est_balance_set_stats (p4est_balance_obj_t *bobj,
                         sc_statistics_t *stats)
{
  bobj->stats = stats;
}

sc_statistics_t *
p4est_balance_get_stats (p4est_balance_obj_t *bobj)
{
  return bobj->stats;
}


