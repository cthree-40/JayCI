// File: pjayci.c
/*
 * Parallel implementation of JayCI determinant CI algorithm.
 *
 * By Christopher L Malbon
 * Yarkony Group, Dept of Chemistry
 * The Johns Hopkins University
 */
#include <stdio.h>
#include <stdlib.h>
#include "pjayci.h"
#include <mpi.h>

/*
 * Main driver.
 */
int main ()
{
        int error = 0;

        int mpi_num_procs  = 0; /* MPI: Number of mpi processes */
        int mpi_proc_rank  = 0; /* MPI: Process rank */
        const int mpi_root = 0; /* Root process always 0 */
        
        MPI_Init(&argc, &argv);



        MPI_Finalize();
        return error;
}
