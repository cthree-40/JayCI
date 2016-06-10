// File: allocate_mem.c

/*
 * Allocate/Deallocate 2d array.
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "allocate_mem.h"

/*
 * allocate_mem_double: allocate a 2d array.
 */
int allocate_mem_double(double ***array, int n, int m)
{
	int i;
	*array = (double **) malloc(m * sizeof(double *));
	for (i = 0; i < m; i++) {
		if (((*array)[i] = (double *) malloc(n * sizeof(double))) == NULL)
		{
			error_flag(i, "allocate_mem_double");
			return i;
		}
	}
	i = 0;
	return i;
}

/*
 * deallocate_mem: deallocate 2d array
 */
int deallocate_mem(double ***array, int n, int m)
{
	int i = 0;
	for (i = 0; i < m; i++) {
		free((*array)[i]);
	}
	free(*array);
	return i;
}
	
