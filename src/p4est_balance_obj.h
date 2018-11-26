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

/** Opaque object for controlling the behavior of 2:1 balance algorithms */
typedef struct p4est_balance_obj_s p4est_balance_obj_t;

/** A list of methods (algorithms) for applying 2:1 balance */
typedef enum
{
  P4EST_BALANCE_DEFAULT = -1, /**< The default is currently to use the two-round algorithm, but that may change in the future: explicitly set a particular algorithm if it is needed */
  P4EST_BALANCE_SORT = 0, /**< Apply by having each process balance the full forest (in compressed form) and sorting the result leaves */
  P4EST_BALANCE_TWOROUND, /**< Each process identifies other processes that it must communicate with and initiates a send/response communication pattern */
  P4EST_BALANCE_NUM_METHODS /**< The numbr of methods available (invalid as an input to #p4est_balance_obj_set_method()) */
}
p4est_balance_method_t;

#define P4EST_BALANCE_STR_SORT     "sort"
#define P4EST_BALANCE_STR_TWOROUND "tworound"

/** The method that will be used if the user chooses P4EST_BALANCE_DEFAULT.
 * The user is free  to change this to any #p4est_balance_method_t that
 * is not #P4EST_BALANCE_DEFAULT or #P4EST_BALANCE_NUM_METHODS */
extern p4est_balance_method_t p4est_balance_method_default;

/** Names for each notify method, which can be used (for example)
 * to specify a method from the command line */
extern const char  *p4est_balance_method_strings[P4EST_BALANCE_NUM_METHODS];

/** Create a balance object that can be customized and used to balance in
 * #p4est_balance_obj().
 *
 * \param[in] mpicomm  The communicator on which all processes are creatin
 *                     the same object.  (Should be the same as for the
 *                     #p4est_t that it will act on).
 * \return             A #p4est_balance_obj_t pointer that can be destroyed
 *                     with #p4est_balance_obj_destroy().
 */
p4est_balance_obj_t *p4est_balance_obj_new (sc_MPI_Comm mpicomm);

/** Destroy a balance object.
 *
 * \param[in/out] bobj  The balance object to be destroyed.
 */
void                p4est_balance_obj_destroy (p4est_balance_obj_t * bobj);

/** Get the sc_MPI_Comm * object.
 *
 * \param [in/out] bobj       Balance object from which communicator will be
 *                            fetched.
 * \return                    A sc_MPI_Comm communicator.
 */
sc_MPI_Comm         p4est_balance_obj_get_comm (p4est_balance_obj_t * notify);

/** Set a sc_statistics_t * object from balance object.  If set,
 * statistics about runtime of balance routines and measures of balance
 * communications will be gathered.
 *
 * \param[in/out]     bobj        The balance object.
 * \param[in]         stats       The sc_statistics_t * object.
 *                                \a bobj does not take ownership of
 *                                \a stats and will not destroy it.
 */
void                p4est_balance_obj_set_stats (p4est_balance_obj_t * bobj,
                                                 sc_statistics_t * stats);

/** Get the sc_statistics_t * object from balance object (see
 * #p4est_balance_obj_set_stats()).
 *
 * \param [in] bobj       The balance object
 * \return                A sc_statistics_t * pointer.
 *                        The user should not destroy the pointer
 *                        while \a bobj is still using it.
 */
sc_statistics_t    *p4est_balance_obj_get_stats (p4est_balance_obj_t * bobj);

/** Set a #p4est_balance_method_t for a balance object.
 * The methods will all produce the same 2:1 balanced octree, but
 * using different algorithms that maybe optimized for different forests
 * and/or hardware.
 *
 * \param[in/out] bobj        The balance object.
 * \param[in]     method      A #p4est_balance_method_t.
 *                            If \a method is #P4EST_BALANCE_DEFAULT,
 *                            the value of #p4est_balance_method_default
 *                            will be used.
 */
void                p4est_balance_obj_set_method (p4est_balance_obj_t * bobj,
                                                  p4est_balance_method_t
                                                  method);

/** Get the #p4est_balance_method_t of a balance object (see
 * #p4est_balance_obj_set_method())
 *
 * \param [in] bobj       The balance object
 * \return                The #p4est_balance_method_t used by \a bobj.
 */
p4est_balance_method_t p4est_balance_obj_get_method (p4est_balance_obj_t *
                                                     bobj);

/** Set the #p4est_connect_type_t of balance object.
 * 2:1 balance may be defined across different types of interfaces:
 *
 * #P4EST_CONNECT_FULL or #P4EST_CONNECT_FACE means that the result
 * of #p4est_balance_obj() will have 2:1 balance (a difference of
 * refinement level of at most one) between quadrants that neighbor across faces
 * or corners.
 *
 * #P4EST_CONNECT_FACE means that the result will have 2:1 balance
 * only across faces.
 *
 * \param[in]     bobj        The balance object.
 * \param[in]     connect     The #p4est_connect_type_t.
 */
void                p4est_balance_obj_set_connect (p4est_balance_obj_t *
                                                   bobj,
                                                   p4est_connect_type_t
                                                   connect);

/** Get the #p4est_connect_type_t from a balance object (see
 * #p4est_balance_obj_set_connect()).
 *
 * \param [in] bobj       The balance object
 * \return                A #p4est_connect_type_t.
 */
p4est_connect_type_t p4est_balance_obj_get_connect (p4est_balance_obj_t *
                                                    bobj);

/** Set a #p4est_init_t callback function from balance object.
 *
 * If specified, this will be called once on each new quadrant
 * created in #p4est_balance_obj().
 *
 * \param[in]     bobj        The balance object.
 * \param[in]     init_fn     The #p4est_init_t.
 */
void                p4est_balance_obj_set_init (p4est_balance_obj_t * bobj,
                                                p4est_init_t init_fn);

/** Get the #p4est_init_t callback function from balance object (see
 * #p4est_balance_obj_set_init())
 *
 * \param [in] bobj       The balance object
 * \return                A #p4est_init_t.
 */
p4est_init_t        p4est_balance_obj_get_init (p4est_balance_obj_t * bobj);

/** Set a p4est_replace_t * callback function from balance object.
 *
 * If specified, this will be called once (recursively) for each
 * parent/children set that is updated in #p4est_balance_obj().
 *
 * \param[in]     bobj        The balance object.
 * \param[in]     replace_fn     The p4est_replace_t * pointer.
 */
void                p4est_balance_obj_set_replace (p4est_balance_obj_t * bobj,
                                                   p4est_replace_t
                                                   replace_fn);

/** Get the p4est_replace_t * callback function from balance object.
 *
 * \param [in] bobj       The balance object
 * \return           A p4est_replace_t * pointer.
 */
p4est_replace_t     p4est_balance_obj_get_replace (p4est_balance_obj_t *
                                                   bobj);

/** Set a adapt flags * value from balance object.
 *
 * If non-recursive adaptivity is required,
 * rather than using the callback-based interface
 * (#p4est_coarsen() / #p4est_refine()), a
 * list of flags of type #p4est_adapt_flag_t can
 * be passed to a balance object, and the resulting
 * coarsening/refinement will be done during a pre-processing step in
 * #p4est_balance_obj().
 *
 * \param[in]     bobj            The balance object.
 * \param[in]     const int8_t    A const int8_t pointer to flags array.
 */
void                p4est_balance_obj_set_adapt_flags (p4est_balance_obj_t *
                                                       bobj,
                                                       const int8_t * flags);

/** Get the adapt flags array * from balance object.
 *
 * \param [in] bobj       The balance object
 * \return           A const_int8_t * array pointer of adapt flags.
 */
const int8_t       *p4est_balance_obj_get_adapt_flags (p4est_balance_obj_t *
                                                       bobj);

/** Set a p4est_inspect_t * object from balance object.
 *
 * \param[in]     bobj            The balance object.
 * \param[in]     p4est_inspect_t The p4est_inspect_t * pointer.
 */
void                p4est_balance_obj_set_inspect (p4est_balance_obj_t * bobj,
                                                   p4est_inspect_t * inspect);

/** Get the p4est_inspect_t * from balance object.
 *
 * \param [in] bobj       The balance object
 * \return           A p4est_inspect_t * pointer.
 */
p4est_inspect_t    *p4est_balance_obj_get_inspect (p4est_balance_obj_t *
                                                   bobj);

/** Set a sc_notify_t * object from balance object.
 *
 * \param[in]     bobj            The balance object.
 * \param[in]     sc_notify_t     The sc_notify_t * pointer.
 */
void                p4est_balance_obj_set_notify (p4est_balance_obj_t * bobj,
                                                  sc_notify_t * notify);

/** Get the sc_notify_t * from balance object.
 *
 * \param [in] bobj       The balance object
 * \return           The sc_notify_t * pointer.
 */
sc_notify_t        *p4est_balance_obj_get_notify (p4est_balance_obj_t * bobj);

/** Set a use_root int from balance object.
 *
 * \param[in]     bobj            The balance object.
 * \param[in]     int          The use_root integer.
 */
void                p4est_balance_obj_sort_set_use_root (p4est_balance_obj_t *
                                                         bobj, int use_root);
/** Balance function that takes p4est and balance object as inputs.
 * \ It will perform fusion balancing operations based on methods
 * \ stored in balance object.
 *
 * \param[in]     bobj      The balance object.
 * \param[in]     p4est_t    The p4est object.
 */
void                p4est_balance_obj (p4est_balance_obj_t * bobj,
                                       p4est_t * p4est);

#endif /* P4EST_BALANCE_OBJ_H */
