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

#ifndef P4EST_BALANCE_OBJ_IMPL_H
#define P4EST_BALANCE_OBJ_IMPL_H

#ifndef P4_TO_P8
#include <p4est.h>
#else
#include <p8est.h>
#endif

#include <sc_notify.h>
#include <sc_flops.h>

struct p4est_balance_obj_s
{
  sc_MPI_Comm         comm;
  p4est_balance_method_t method;
  sc_statistics_t    *stats;
  p4est_connect_type_t connect;
  p4est_init_t        init_fn;
  p4est_replace_t     replace_fn;
  p4est_inspect_t    *inspect;
  sc_notify_t        *notify;
  int                 flop_started;
  sc_flopinfo_t       flop;
  const int8_t       *adapt_flags;
};

#define P4EST_BAL_FUNC_SNAP(bobj,snap)                       \
  do {                                                       \
    if ((bobj)->stats) {                                     \
      if (!(bobj)->inspect->flop_started) {                   \
        (bobj)->flop_started = 1;                            \
        sc_flops_start (&((bobj)->flop));                    \
      }                                                      \
      SC_FUNC_SNAP ((bobj)->stats, &((bobj)->flop), (snap)); \
    }                                                        \
  } while (0)

#define P4EST_BAL_FUNC_SHOT(bobj,snap)                       \
  do {                                                       \
    if ((bobj)->stats) {                                     \
      SC_FUNC_SHOT ((bobj)->stats, &((bobj)->flop), (snap)); \
    }                                                        \
  } while (0)

#endif
