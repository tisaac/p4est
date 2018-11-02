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

/********************************************************************
 *                          IMPORTANT NOTE                          *
 *                                                                  *
 * These interfaces are intended for those who like finer control.  *
 * The API offers extended versions of some basic p4est functions.  *
 * The API may change without notice.                               *
 ********************************************************************/

/** \file p4est_balance_obj.h
 *
 *
 * \ingroup p4est
 */

#ifndef P4EST_BALANCE_OBJ_H
#define P4EST_BALANCE_OBJ_H

#include <p4est.h>
#include <p4est_extended.h>
#include <sc_statistics.h>

SC_EXTERN_C_BEGIN;

typedef struct p4est_balance_obj_s p4est_balance_obj_t;

typedef enum
{
  P4EST_BALANCE_DEFAULT = -1,
  P4EST_BALANCE_SORT = 0,
  P4EST_BALANCE_TWOROUND,
  P4EST_BALANCE_NUM_METHODS
}
p4est_balance_method_t;

#define P4EST_BALANCE_STR_SORT     "sort"
#define P4EST_BALANCE_STR_TWOROUND "tworound"

/** The method that will be used if the user chooses P4EST_BALANCE_DEFAULT */
extern p4est_balance_method_t p4est_balance_method_default;

/** Names for each notify method */
extern const char  *p4est_balance_method_strings[P4EST_BALANCE_NUM_METHODS];

p4est_balance_obj_t *p4est_balance_obj_new (sc_MPI_Comm mpicomm);
void                p4est_balance_obj_destroy (p4est_balance_obj_t * bobj);

sc_MPI_Comm         p4est_balance_obj_get_comm (p4est_balance_obj_t * notify);

void                p4est_balance_obj_set_stats (p4est_balance_obj_t * bobj,
                                                 sc_statistics_t * stats);
sc_statistics_t    *p4est_balance_obj_get_stats (p4est_balance_obj_t * bobj);

void                p4est_balance_obj_set_method (p4est_balance_obj_t * bobj,
                                                  p4est_balance_method_t
                                                  method);
p4est_balance_method_t p4est_balance_obj_get_method (p4est_balance_obj_t *
                                                     bobj);

void                p4est_balance_obj_set_connect (p4est_balance_obj_t *
                                                   bobj,
                                                   p4est_connect_type_t
                                                   connect);
p4est_connect_type_t p4est_balance_obj_get_connect (p4est_balance_obj_t *
                                                    bobj);

void                p4est_balance_obj_set_init (p4est_balance_obj_t * bobj,
                                                p4est_init_t init_fn);
p4est_init_t        p4est_balance_obj_get_init (p4est_balance_obj_t * bobj);

void                p4est_balance_obj_set_replace (p4est_balance_obj_t * bobj,
                                                   p4est_replace_t
                                                   replace_fn);
p4est_replace_t     p4est_balance_obj_get_replace (p4est_balance_obj_t *
                                                   bobj);

/* TODO: make balance_sort use these flags */
void                p4est_balance_obj_set_adapt_flags (p4est_balance_obj_t *bobj,
                                                       const int8_t *flags);
const int8_t *      p4est_balance_obj_get_adapt_flags (p4est_balance_obj_t *bobj);

void                p4est_balance_obj_set_inspect (p4est_balance_obj_t * bobj, p4est_inspect_t * inspect);

p4est_inspect_t *   p4est_balance_obj_get_inspect (p4est_balance_obj_t * bobj);

void                p4est_balance_obj_set_notify (p4est_balance_obj_t * bobj, sc_notify_t * notify);

sc_notify_t *       p4est_balance_obj_get_notify (p4est_balance_obj_t * bobj);

void                p4est_balance_obj (p4est_balance_obj_t * bobj,
                                       p4est_t * p4est);

#endif /* P4EST_BALANCE_OBJ_H */
