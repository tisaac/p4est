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
#include <p4est_algorithms.h>
#include <p4est_communication.h>
#else
#include <p8est_balance_obj.h>
#include <p8est_algorithms.h>
#include <p8est_communication.h>
#endif

struct p4est_balance_obj_s
{
  sc_MPI_Comm         comm;
  p4est_balance_method_t method;
  sc_statistics_t    *stats;
  p4est_connect_type_t connect;
  p4est_init_t        init_fn;
  p4est_replace_t     replace_fn;
};

const char         *p4est_balance_method_strings[] = {
  P4EST_BALANCE_STR_SORT,
  P4EST_BALANCE_STR_TWOROUND,
};

p4est_balance_method_t p4est_balance_method_default = P4EST_BALANCE_TWOROUND;

p4est_balance_obj_t *
p4est_balance_obj_new (sc_MPI_Comm mpicomm)
{
  p4est_balance_obj_t *bal;

  bal = SC_ALLOC_ZERO (p4est_balance_obj_t, 1);

  bal->comm = mpicomm;
  bal->connect = P4EST_CONNECT_FULL;
  p4est_balance_obj_set_method (bal, P4EST_BALANCE_DEFAULT);
  return bal;
}

void
p4est_balance_obj_destroy (p4est_balance_obj_t * obj)
{
  SC_FREE (obj);
}

sc_MPI_Comm
p4est_balance_obj_get_comm (p4est_balance_obj_t * obj)
{
  return obj->comm;
}

void
p4est_balance_obj_set_stats (p4est_balance_obj_t * bobj, sc_statistics_t * stats)
{
  bobj->stats = stats;
}

sc_statistics_t    *
p4est_balance_obj_get_stats (p4est_balance_obj_t * bobj)
{
  return bobj->stats;
}

void
p4est_balance_obj_set_method (p4est_balance_obj_t * bobj,
                              p4est_balance_method_t in_method)
{
  p4est_balance_method_t current_method;

  current_method = p4est_balance_obj_get_method (bobj);
  if (in_method == P4EST_BALANCE_DEFAULT) {
    in_method = p4est_balance_method_default;
  }
  P4EST_ASSERT (in_method >= 0 && in_method < P4EST_BALANCE_NUM_METHODS);
  if (current_method != in_method) {
    bobj->method = in_method;
    /* initialize data */
    switch (in_method) {
    default:
      break;
    }
  }
}

p4est_balance_method_t
p4est_balance_obj_get_method (p4est_balance_obj_t * bobj)
{
  return bobj->method;
}

void
p4est_balance_obj_set_connect (p4est_balance_obj_t * bobj,
                               p4est_connect_type_t connect)
{
  bobj->connect = connect;      //TODO checks here?
}

p4est_connect_type_t
p4est_balance_obj_get_connect (p4est_balance_obj_t * bobj)
{
  return bobj->connect;
}

void
p4est_balance_obj_set_init (p4est_balance_obj_t * bobj, p4est_init_t init_fn)
{
  bobj->init_fn = init_fn;
}

p4est_init_t
p4est_balance_obj_get_init (p4est_balance_obj_t * bobj)
{
  return bobj->init_fn;
}

void
p4est_balance_obj_set_replace (p4est_balance_obj_t * bobj,
                               p4est_replace_t replace_fn)
{
  bobj->replace_fn = replace_fn;
}

p4est_replace_t
p4est_balance_obj_get_replace (p4est_balance_obj_t * bobj)
{
  return bobj->replace_fn;
}

extern void p4est_balance_sort (p4est_balance_obj_t * bobj, p4est_t *p4est);
extern void p4est_balance_tworound (p4est_balance_obj_t * bobj, p4est_t *p4est);

void
p4est_balance_obj (p4est_balance_obj_t * bobj, p4est_t * p4est)
{
  p4est_balance_method_t method = p4est_balance_obj_get_method (bobj);
  p4est_connect_type_t ctype = p4est_balance_obj_get_connect (bobj);
  p4est_gloidx_t      old_gnq;
  const int8_t       *pre_adapt_flags = NULL;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING
                            "_balance %s with %lld total quadrants\n",
                            p4est_connect_type_string (ctype),
                            (long long) p4est->global_num_quadrants);
  p4est_log_indent_push ();
  P4EST_GLOBAL_INFOF ("Balance algorithm: %s\n", p4est_balance_method_strings[method]);
  P4EST_ASSERT (p4est_is_valid (p4est));
  /* remember input quadrant count; it will not decrease */
  old_gnq = p4est->global_num_quadrants;
  switch (method) {
  case P4EST_BALANCE_SORT:
    p4est_balance_sort (bobj, p4est);
    break;
  case P4EST_BALANCE_TWOROUND:
    p4est_balance_tworound (bobj, p4est);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  /* compute global number of quadrants */
  p4est_comm_count_quadrants (p4est);
  if (p4est->inspect) {
    pre_adapt_flags = p4est->inspect->pre_adapt_flags;
  }
  P4EST_ASSERT (pre_adapt_flags || p4est->global_num_quadrants >= old_gnq);
  if (pre_adapt_flags || old_gnq != p4est->global_num_quadrants) {
    ++p4est->revision;
  }
  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (p4est_is_balanced (p4est, ctype));
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P4EST_STRING
                            "_balance with %lld total quadrants\n",
                            (long long) p4est->global_num_quadrants);
}

