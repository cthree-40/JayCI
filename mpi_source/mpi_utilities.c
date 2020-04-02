// File: mpi_utilities.c
/*
 * Utilities for MPI implementations of pjayci. The routines contained
 * herein modify global pjayci.x variables, and are responsible for
 * inter-node communication. Mostly interfaces for MPI Library routines.
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "mpi_utilities.h"
#include "pjayci_global.h"
#include <mpi.h>
#include <ga.h>
#include <macdecls.h>

/*
 * mpi_error_check_msg: check for error. Print message if error has
 * occured.
 * Input:
 *  error    = error flag
 *  fcn_name = calling function name
 *  message  = error message to print
 */
void mpi_error_check_msg (int error, char *fcn_name, char *message)
{
        int eflag;
        MPI_Allreduce(&error, &eflag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (eflag != 0) {
                MPI_Allreduce(&error, &eflag, 1, MPI_INT, MPI_MINLOC,
                              MPI_COMM_WORLD);
                if (mpi_proc_rank == mpi_root) {
                        error_message(eflag, message, fcn_name);
                }
                MPI_Finalize();
                exit(-1);
        }
}

/*
 * mpi_split_work_array_1d: get first and last elements for partitioning
 * a 1d array amongst work processes.
 * Input:
 *  len  = length of vector (total)
 * Output:
 *  chunk = chunk size
 *  lo    = first element
 *  hi    = last  element
 */
void mpi_split_work_array_1d (int len, int *chunk, int *lo, int *hi)
{
        *chunk = 0;
        *chunk = len / mpi_num_procs;
        *lo = mpi_proc_rank * (*chunk);
        *hi = (mpi_proc_rank + 1) * (*chunk) - 1;
        if ((*chunk * mpi_num_procs) != len) {
                if (mpi_proc_rank == (mpi_num_procs - 1)) {
                        *hi = len - 1; // Last element of array
                }
        }
        return;
}

/*
 * print_gavectors2file: print a set of GA vectors to a file.
 * Input:
 *  hndl  = global arrays handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void print_gavectors2file(int hndl, int len, int dim, char *fname)
{
        FILE *fptr = NULL;

        GA_Sync();
        
        return;
}

/*
 * set_mpi_process_number_and_rank: set global variables $mpi_num_procs
 * and $mpi_proc_rank.
 */
void set_mpi_process_number_and_rank ()
{
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_proc_rank);
}

/*
 * set_ga_process_number_and_rank: set global variables $mpi_num_procs
 * and $mpi_proc_rank with GA wrappers.
 */
void set_ga_process_number_and_rank ()
{
        mpi_proc_rank = GA_Nodeid();
        mpi_num_procs = GA_Nnodes();
}
