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
#include "pjayci_global.h"
#include "mpi_utilities.h"
#include "execute_pjayci.h"
#include <mpi.h>
#include <ga.h>
#include <macdecls.h>

const int mpi_root = 0; /* MPI: Root rank is always 0 */

/*
 * Main driver.
 */
int main (int argc, char **argv)
{
        int error = 0;

        MPI_Init(&argc, &argv);
        GA_Initialize();
        
        set_mpi_process_number_and_rank();

        error = execute_pjayci();

        GA_Terminate();
        MPI_Finalize();

        return error;
}
