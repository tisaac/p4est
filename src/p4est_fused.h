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

/** \file p4est_fused.h
 *
 * Routines for fused adaptivity operations
 *
 * \ingroup p4est
 */

#ifndef P4EST_FUSED_H
#define P4EST_FUSED_H

#include <p4est.h>
#include <p4est_ghost.h>
#include <p4est_extended.h>
#include <p4est_balance_obj.h>

SC_EXTERN_C_BEGIN;

/** Coarsen, refine, balance, repartition, and reconstruct the ghost layer of
 * a forest.
 *
 */

enum
{
  P4EST_FUSED_KEEP = 0,
  P4EST_FUSED_REFINE,
  P4EST_FUSED_COARSEN
};

/** This is a fused adaptation loop equivalent to the following code:
 *
 * if (*p4est_out != p4est) {
 *   *p4est_out = p4est_copy (p4est, copy_data);
 * }
 * // coarsen families where every child's adapt_flag is P4EST_FUSED_COARSEN
 * p4est_coarsen_ext (*p4est_out, 0, 0, ..., init_fn, replace_fn);
 * // refine all quadrants marked P4EST_FUSED_REFINE
 * p4est_refine_ext (*p4est_out, 0, P4EST_QMAXLEVEL, ..., init_fn, *replace_fn);
 * p4est_balance_ext (*p4est_out, balance_type, init_fn, replace_fn);
 * if (repartition) {
 *   p4est_partition (*p4est_out, partition_for_coarsening, weight_fn);
 * }
 * if (ghost_layer_width > 0) {
 *   int i;
 *
 *   *ghost_out = p4est_ghost_new (*p4est_out, ghost_type);
 *   for (i = 1; i < ghost_layer_width; i++) {
 *     p4est_ghost_expand (*p4est_out, *ghost_out);
 *   }
 * }
 *
 */
void                p4est_adapt_fused (p4est_t * p4est,
                                       const int8_t * adapt_flag,
                                       int copy_data,
                                       p4est_balance_obj_t * bobj,
                                       int repartition,
                                       int partition_for_coarsening,
                                       int ghost_layer_width,
                                       p4est_connect_type_t ghost_type,
                                       p4est_weight_t weight_fn,
                                       p4est_init_t init_fn,
                                       p4est_replace_t replace_fn,
                                       p4est_t ** p4est_out,
                                       p4est_ghost_t ** ghost_out);

void                p4est_adapt_fused_reference (p4est_t * p4est,
                                                 const int8_t * adapt_flag,
                                                 int copy_data,
                                                 p4est_balance_obj_t * bobj,
                                                 int repartition,
                                                 int partition_for_coarsening,
                                                 int ghost_layer_width,
                                                 p4est_connect_type_t
                                                 ghost_type,
                                                 p4est_weight_t weight_fn,
                                                 p4est_init_t init_fn,
                                                 p4est_replace_t replace_fn,
                                                 p4est_t ** p4est_out,
                                                 p4est_ghost_t ** ghost_out);

SC_EXTERN_C_END;

#endif /* !P4EST_FUSED_H */
