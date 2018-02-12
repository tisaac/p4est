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

/** \file p8est_fused.h
 *
 * Routines for fused adaptivity operations
 *
 * \ingroup p8est
 */

#ifndef P8EST_FUSED_H
#define P8EST_FUSED_H

#include <p8est.h>
#include <p8est_ghost.h>
#include <p8est_extended.h>

SC_EXTERN_C_BEGIN;

/** Coarsen, refine, balance, repartition, and reconstruct the ghost layer of
 * a forest.
 *
 */

enum
{
  P8EST_FUSED_KEEP = 0,
  P8EST_FUSED_REFINE,
  P8EST_FUSED_COARSEN
};

/** This is a fused adaptation loop equivalent to the following code:
 *
 * if (*p8est_out != p8est) {
 *   *p8est_out = p8est_copy (p8est, copy_data);
 * }
 * // coarsen families where every child's adapt_flag is P8EST_FUSED_COARSEN
 * p8est_coarsen_ext (*p8est_out, 0, 0, ..., init_fn, replace_fn);
 * // refine all quadrants marked P8EST_FUSED_REFINE
 * p8est_refine_ext (*p8est_out, 0, P8EST_QMAXLEVEL, ..., init_fn, *replace_fn);
 * p8est_balance_ext (*p8est_out, balance_type, init_fn, replace_fn);
 * if (repartition) {
 *   p8est_partition (*p8est_out, partition_for_coarsening, weight_fn);
 * }
 * if (ghost_layer_width > 0) {
 *   int i;
 *
 *   *ghost_out = p8est_ghost_new (*p8est_out, ghost_type);
 *   for (i = 1; i < ghost_layer_width; i++) {
 *     p8est_ghost_expand (*p8est_out, *ghost_out);
 *   }
 * }
 *
 */
void                p8est_adapt_fused (p8est_t * p8est,
                                       const int8_t * adapt_flag,
                                       int copy_data,
                                       p8est_connect_type_t balance_type,
                                       int repartition,
                                       int partition_for_coarsening,
                                       int ghost_layer_width,
                                       p8est_connect_type_t ghost_type,
                                       p8est_weight_t weight_fn,
                                       p8est_init_t init_fn,
                                       p8est_replace_t replace_fn,
                                       p8est_t ** p8est_out,
                                       p8est_ghost_t ** ghost_out);

void                p8est_adapt_fused_reference (p8est_t * p8est,
                                                 const int8_t * adapt_flag,
                                                 int copy_data,
                                                 p8est_connect_type_t
                                                 balance_type,
                                                 int repartition,
                                                 int partition_for_coarsening,
                                                 int ghost_layer_width,
                                                 p8est_connect_type_t
                                                 ghost_type,
                                                 p8est_weight_t weight_fn,
                                                 p8est_init_t init_fn,
                                                 p8est_replace_t replace_fn,
                                                 p8est_t ** p8est_out,
                                                 p8est_ghost_t ** ghost_out);

SC_EXTERN_C_END;

#endif /* !P8EST_FUSED_H */
