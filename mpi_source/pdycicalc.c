// File: pdycicalc.c
/*
 * Calculate the dyson orbital for JayCI wavefunctions.
 *
 * By Christopher L Malbon
 * Yarkony Group, Department of Chemistry
 * The Johns Hopkins University
 */
#include <stdio.h>
#include <stdlib.h>
#include "mpi_utilities.h"
#include "run_pdycicalc.h"
#include <mpi.h>
#include <ga.h>
#include <macdecls.h>

const int mpi_root = 0; /* MPI: Root rank is always 0 */

/*
 * main: driver for program.
 */
int main (int argc, char **argv)
{
        int error = 0;
#ifdef BIGMEM
        long long int stack = 800000000, heap = 800000000;
#else
        long long int stack = 80000000, heap = 80000000;
#endif
        

        MPI_Init(&argc, &argv);
        GA_Initialize();

        set_ga_process_number_and_rank();

        /* Initialize memory allocator (MA) */
        if (! MA_init(C_DBL, stack, heap)) GA_Error("MA_init failed",stack+heap);
        
        error = run_pdycicalc();
        mpi_error_check_msg(error, "run_pdycicalc", "Error occured!");
        
        GA_Terminate();
        MPI_Finalize();

        return error;
}

