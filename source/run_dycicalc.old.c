/* File: run_dycicalc.c
 *
 * run_dycicalc: Execute CI dyson orbital calculation.
 */
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "ioutil.h"
#include "execute_dycicalc.h"
#include "run_dycicalc.h"

/*
 * run_dycicalc: Execute CI dyson orbital calculation.
 */
int run_dycicalc()
{
        char wvfcn_file0[FLNMSIZE]; /* N+1 electron wavefunction file */
        char wvfcn_file1[FLNMSIZE]; /* N   electron wavefunction file */
        int nstates0 = 0; /* N+1 electron wavefunction states */
        int nelecs0  = 0; /* N+1 electron wavefunction electrons */
        int ndets0   = 0; /* N+1 electron wavefunction determinants */
        int norbs0   = 0; /* N+1 electron wavefunction orbitals */
        int ninto0   = 0; /* N+1 electron internal orbitals */
        int nstates1 = 0; /* N electron wavefunction states to read */
        int nelecs1  = 0; /* N electron wavefunction electrons */
        int ndets1   = 0; /* N electron wavefunction determinants */
        int norbs1   = 0; /* N electron wavefunction orbitals */
        int ninto1   = 0; /* N electron internal orbitals */
        int error = 0; /* Error flag */

        /* Check for namelist input file. Read namelist input file. Verify
         * that the wavefunction files are present in the work directory. */
        error = check_for_input_files();
        if (error != 0) {
                error_flag(error, "run_dycicalc: missing input files!");
                return error;
        }
        read_dysonorb_input(wvfcn_file0, wvfcn_file1, &nstates0, &nstates1,
                            &nelecs0, &nelecs1, &norbs0, &norbs1, &ndets0,
                            &ndets1, &ninto0, &ninto1, &error);
        if (error != 0) {
                error_flag(error, "run_dycicalc: error reading input!");
                return error;
        }
        error = check_for_file(wvfcn_file0, "r");
        if (error != 0) {
                error_flag(error, "run_dycicalc");
                return error;
        }
        error = check_for_file(wvfcn_file1, "r");
        if (error != 0) {
                error_flag(error, "run_dycicalc");
                return error;
        }

        /* Both wavefunction files exist, and we have read the main input
         * file. Perform dyson orbital calulcation between */
        // Ensure orbital number is the same for both wavefunctions.
        if (norbs1 != norbs0) {
                error_message("norbs1 != norbs0", "run_dycicalc");
                error = 1;
                return error;
        }
        error = execute_dycicalc(wvfcn_file0, wvfcn_file1, nstates0, nstates1,
                                 nelecs0, nelecs1, norbs1, ndets0,
                                 ndets1);
        
        return error;
}

/* 
 * check_for_input_files: Check for the proper input and output files. 
 * Proper files: dycicalc.in.
 */
int check_for_input_files()
{
        int error = 0;
        error = check_for_file("dycicalc.in", "r");
        return error;
}
