/* File: run_dycicalc.c
 *
 * run_dycicalc: Execute CI dyson orbital calculation.
 */
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "ioutil.h"
#include "run_dycicalc.h"

/*
 * run_dycicalc: Execute CI dyson orbital calculation.
 */
int run_dycicalc()
{
        char wvfcn_file0[FLNMSIZE]; /* N+1 electron wavefunction file */
        char wvfcn_file1[FLNMSIZE]; /* N   electron wavefunction file */
        int error = 0; /* Error flag */

        error = check_for_input_files();
        if (error != 0) {
                error_flag(error, "run_dycicalc: missing input files!");
                return error;
        }
        read_dysonorb_input(wvfcn_file0, wvfcn_file1, &error);
        if (error != 0) {
                error_flag(error, "run_dycicalc: error reading input!");
                return error;
        }
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
