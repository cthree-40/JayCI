// File: pjayci.c

/*
 * Main program file for pjayci.x.
 *
 *
 * By Christopher L Malbon (cthree-40)
 * Yarkony Group, Dept. of Chemistry, The Johns Hopkins University
 * 2016-04-29
 */

#include <stdio.h>
#include <stdlib.h>
#include "progheader.h"
#include "cmdline.h"
#include "run_jayci.h"

/*
 * main: Reads in command line arguments, then calls run_jayci().
 */
int main(int argc, char *argv[])
{
	int error = 0; /* Error handling */
	print_jayciheader();
	error = mpi_runjayci();
	return error;
}
    

    
    
