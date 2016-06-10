// File: jayci.c

/*
 * Main program file for jayci.x.
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
int main()
{
	int error = 0; /* Error handling */
	
	/* Print header. */
	print_jayciheader();
	
	/* Process command line arguments */
	//error = processcmdlineargs(argc, argv, &memory);
	//if (error != 0) exit(error);
	
	/* Call run_jayci. */
	error = run_jayci();
	if (error != 0) exit(error);
	
	return error;
}
    

    
    
