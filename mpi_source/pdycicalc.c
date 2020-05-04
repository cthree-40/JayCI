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
#include "pdyson_global.h"
#include "mpi_utilities.h"
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
        int stack = 8000000, heap = 8000000;

        MPI_Init(&argc, &argv);
        GA_Initialize();

        /* Initialize memory allocator (MA) */
        if (! MA_init(C_DBL, stack, heap)) GA_Error("MA_init failed",stack+heap);

        set_ga_process_number_and_rank();

        error = 0;

        GA_Terminate();
        MPI_Finalize();

        return error;
}
