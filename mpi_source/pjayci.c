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
#include "binarystr.h"
#include "citruncate.h"
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
#ifdef BIGMEM
	long long int stack = 800000000, heap = 800000000;
#else
        long long int stack = 80000000,  heap = 80000000;
#endif
        
        MPI_Init(&argc, &argv);
        GA_Initialize();

        /* Initialize memory allocator (MA) */
        if (! MA_init(C_DBL, stack, heap)) GA_Error("MA_init failed",stack+heap);

        set_ga_process_number_and_rank();
	printf("Greetings from process: %d\n", mpi_proc_rank);
	fflush(stdout);
        error = execute_pjayci();

        GA_Terminate();
        MPI_Finalize();

        return error;
}

