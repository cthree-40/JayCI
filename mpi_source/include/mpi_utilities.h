// File: mpi_utilities.h
/*
 * Utilities for MPI implementations of pjayci. The routines contained
 * herein modify global pjayci.x variables, and are responsible for
 * inter-node communication. Mostly interfaces for MPI Library routines.
 */

#ifndef mpi_utilities_h
#define mpi_utilities_h

/*
 * set_mpi_process_number_and_rank: set global variables $mpi_num_procs
 * and $mpi_proc_rank.
 */
void set_mpi_process_number_and_rank ();

#endif
