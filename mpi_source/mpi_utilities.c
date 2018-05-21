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
 * set_mpi_process_number_and_rank: set global variables $mpi_num_procs
 * and $mpi_proc_rank.
 */
void set_mpi_process_number_and_rank ()
{
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_proc_rank);
}
