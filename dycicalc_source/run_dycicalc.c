/* File: run_dycicalc.c
 *
 * run_dycicalc: Execute CI dyson orbital calculation.
 */
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "run_dycicalc.h"

/*
 * run_dycicalc: Execute CI dyson orbital calculation.
 */
int run_dycicalc()
{
        int error = 0; /* Error flag */

        error = check_for_input_files();
        if (error != 0) {
                error_flag(error, "run_dycicalc");
                return error;
        }

        return error;
}

/* 
 * check_for_input_files: Check for the proper input and output files. 
 */
int check_for_input_files();
{
        int error = 0;
        return error;
}
