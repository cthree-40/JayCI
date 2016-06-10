// File: cmdline.c

/*
 * Process command line arguments.
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "cmdline.h"

/*
 * processcmdlineargs: processes command line arguments. Output of function
 * is error flag. Returns amoutn of memory to be allocated in bytes.
 */
int processcmdlineargs(int argc, char *argv[], long long int *memory)
{
    int i, c, error;
    
    if (argc == 1) {
	fprintf(stderr,"Usage: jayci.x -m [memory]\n");
	error = -1;
	error_flag(error, "processcmdlineargs");
	return error;
    }
    
    for (i = 1; i < argc; i++) {
	if (*argv[i] == '-') {
	    while ((c = *(++argv[i]))) {
		switch (c) {
		case 'm':
		    *memory = atoi(argv[i+1]);
		    break;
		default:
		    fprintf(stderr,"Illegal option: %c\n", c);
		    i = 100;
		    break;
		}
	    }
	}
    }

    if (argc != 3) {
	fprintf(stderr,"Usage: jayci.x -m [memory]\n");
	error = -1;
	error_flag(error, "processcmdlineargs");
	return error;
    } else {
	fprintf(stdout,"Allocated %lld bytes\n", *memory);
	error = 0;
	return error;
    }
}
    
