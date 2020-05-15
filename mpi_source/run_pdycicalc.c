// File: run_pdycical.c
#include <stdio.h>
#include <stdlib.h>
#include "pdyson_global.h"
#include "mpi_utilities.h"
#include <mpi.h>
#include <ga.h>
#include <macdecls.h>

/*
 * run_pdycicalc: Execute CI dyson orbital calculation.
 */
int run_pdycicalc ()
{
        int error = 0; /* Error flag */

        int w0_hndl = 0, w1_hndl = 0; /* GA wavefunctions 0 and 1 */
        int w_dims[2] = {0, 0};       /* GA wavefunction dimension */
        int w_chunk[2]= {0, 0};       /* GA wavefunction chunk sizes */

        double **dyorb = NULL;        /* LOCAL dyson orbitals */
        double *dyorb_data = NULL;    /* LOCAL dyson orbital memory block */
        
        /* Read wavefunction input. */
        
        /* Allocate GLOBAL arrays: W0, W1, V0, V1 */

        /* Allocate LOCAL arrays: dyorb */
        
        /* Read civectors. */

        /* Compute dyson orbital */
        
        return error;
}
