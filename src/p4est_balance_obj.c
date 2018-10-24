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

p4est_balance_method_t   p4est_balance_method_default = P4EST_BALANCE_SORT; //TODO: sort as default?

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

p4est_balance_method_t
p4est_balance_get_method (p4est_balance_obj_t * bobj)
{
  return bobj->method;
}

void      p4est_balance_obj_set_method (p4est_balance_obj_t *bobj,
                                        p4est_balance_method_t in_method)
{
	p4est_balance_method_t    current_method;

	current_method = p4est_balance_get_method (bobj);
	if (in_method == P4EST_BALANCE_DEFAULT) {
	  in_method = p4est_balance_method_default;
	}
	P4EST_ASSERT (in_method >= 0 && in_type < P4EST_BALANCE_NUM_METHODS);
	if (current_method != in_method){
		bobj->method = in_method;
		/* initialize data */
		switch (in_method) {
			case P4EST_BALANCE_SORT:
				p4est_balance_sort_init (bobj); //????
				break;
			case P4EST_BALANCE TWOROUND:
				p4est_balance_tworound_init (bobj); //????
				break;
			default:
				SC_ABORT_NOT_REACHED();
		}
	}
}

p4est_balance_method_t p4est_balance_obj_get_method (p4est_balance_obj_t *bobj)
{
  return bobj->method; 
}

void      p4est_balance_obj_set_connect_type (p4est_balance_obj_t *bobj,
                                              p4est_connect_type_t connect)
{
	bobj->connect = connect;//TODO checks here?
}

p4est_connect_type_t p4est_balance_obj_get_connect (p4est_balance_obj_t *bobj)
{
  return bobj->connect;
}

void      p4est_balance_obj_set_init (p4est_balance_obj_t *bobj,
                                      p4est_init_t init_fn)
{
 //TODO
}

p4est_init_t p4est_balance_obj_get_init (p4est_balance_obj_t *bobj)
{
  return bobj->init_fn;
}

void      p4est_balance_obj_set_replace (p4est_balance_obj_t *bobj,
                                         p4est_replace_t replace_fn)
{
 //TODO 
}

p4est_replace_t p4est_balance_obj_get_replace (p4est_balance_obj_t *bobj)
{
  return bobj->replace_fn;
}

void      p4est_balance_obj (p4est_balance_obj_t * bobj, p4est_t *p4est)
{
  p4est_balance_method_t	method = p4est_balance_obj_get_method (bobj);
  p4est_connect_type_t		ctype = p4est_balance_obj_get_connect (bobj);
  p4est_init_t				init_fn = p4est_balance_obj_get_init (bobj);
  p4est_replace_t			replace_fn = p4est_balance_obj_get_replace (bobj);

  switch (method) {
	  case P4EST_BALANCE_SORT:
		  p4est_balance_sort (p4est, ctype, init_fn, replace_fn);
		  break;
	  case P4EST_BALANCE_TWOROUND:
		  p4est_balance_tworound (p4est, ctype, init_fn, replace_fn);
		  break;
	  default:
		 SC_ABORT_NOT_REACHED (); 
  }
}

static void
p4est_balance_init (p4est_balance_obj_t * bobj)
{
  sc_MPI_Comm			comm;
  int					success, mpisize, mpirank;
 
  success = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (success);
  bobj->method
  notify->data.nary.mpisize = mpisize;
  success = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (success);
  notify->data.nary.mpirank = mpirank;
  sc_notify_nary_set_widths (notify, sc_notify_nary_ntop_default,
                             sc_notify_nary_nint_default,
                             sc_notify_nary_nbot_default);
}


							
