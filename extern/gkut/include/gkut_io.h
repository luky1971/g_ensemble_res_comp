/*
 * Copyright 2016 Ahnaf Siddiqui
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.
 */

#ifndef GKUT_IO_H
#define GKUT_IO_H

#include "vec.h"
#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif

#define FRAMESTEP 1000 // The number of new frames by which to reallocate an array of length # trajectory frames

void gk_read_traj_t(const char *traj_fname, real **t, rvec ***x, matrix **box, int *nframes, int *natoms, output_env_t *oenv);
/* Reads a trajectory file.
 * real *t is each frame's time indexed [frame #].
 * rvec **x is position coordinates indexed x[frame #][atom #].
 * matrix *box is a 1D array of matrices indexed [frame #].
 * 2D memory is allocated for x and 1D memory is allocated for box.
 * Use free_traj(rvec **x) to free x, and sfree(matrix *box) to free box.
 */

void gk_read_traj(const char *traj_fname, rvec ***x, matrix **box, int *nframes, int *natoms, output_env_t *oenv);
/* Same as read_traj function above but does not return time information.
 */

void gk_free_traj(rvec **x, int nframes, int natoms);
/* Frees the memory associated with a trajectory of vectors.
 */

void gk_print_traj(rvec **x, int nframes, int natoms, const char *fname);
/* Prints the given 2D array of vectors to a text file with the given name.
 */

void gk_ndx_get_indx(const char *ndx_fname, int numgroups, atom_id ***indx, int **isize);
/* Get indexes for a given number of groups from a given index file.
 * Indexes are put in 2D array **indx and number of indexes per group are put in 1D array *isize.
 * 2D memory is allocated for indx and 1D memory is allocated for isize.
 */

void gk_filter_vecs(atom_id *indx, int isize, rvec *pre_x, rvec **new_x);
/* Creates a new 1D array of vectors new_x with only the coordinates from pre_x that are specified by the indexes in indx.
 * 1D memory is allocated for new_x.
 */

void gk_ndx_filter_traj(const char *ndx_fname, rvec **pre_x, rvec ***new_x, int nframes, int *natoms);
/* Creates a new trajectory with only the coordinates from the old trajectory (pre_x) that are specified by the index file.
 * 2D memory is allocated for new_x.
 */

int gk_read_topology(const char *top_fname, t_topology *top);
/* Reads topology information from a gro or tpr file into a t_topology struct.
 * Other formats are not currently supporeted.
 */

void gk_read_top_gro(const char *gro_fname, t_topology *top);
/* Reads topology information from a gro file into a t_topology struct.
 * Allocates memory for data within t_topology *top.
 * Use free_topology to free.
 */

void gk_free_topology(t_topology *top);
/* Frees the dynamic memory in a t_topology struct.
 */

void gk_read_top_tpr(const char *tpr_fname, gmx_mtop_t *mtop);
/* Reads topology information from a tpr file into a gmx_mtop_t struct.
 * Allocates memory for data within gmx_mtop_t *mtop.
 * Use free_mtop to free memory.
 */

void gk_free_mtop(gmx_mtop_t *mtop);
/* Frees the dynamic memory in a gmx_mtop_t struct.
 */

#endif // GKUT_IO_H