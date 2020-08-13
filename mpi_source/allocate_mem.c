// File: allocate_mem.c

/*
 * Allocate/Deallocate 2d array.
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "allocate_mem.h"
#include "mpi_utilities.h"
#include <mpi.h>

/*
 * allocate_mem_char: allocated a 2d character array.
 */
int allocate_mem_char(char ***array, int n, int m)
{
        int i = 0;
        *array = (char **) malloc(sizeof(char *) * m);
        for (i = 0; i < m; i++) {
                if (((*array)[i] = (char *) malloc(sizeof(char) * n)) == NULL)
                {
                        error_flag(mpi_proc_rank, i, "allocate_mem_double");
                        return i;
                }
        }
        i = 0;
        return i;
}

/*
 * allocate_mem_double: allocate a 2d array.
 * Array = M[n x m]
 */
int allocate_mem_double(double ***array, int n, int m)
{
	int i = 0;
	*array = (double **) malloc(m * sizeof(double *));
	for (i = 0; i < m; i++) {
		if (((*array)[i] = (double *) malloc(n * sizeof(double))) == NULL)
		{
			error_flag(mpi_proc_rank, i, "allocate_mem_double");
			return i;
		}
	}
	i = 0;
	return i;
}

/*
 * allocate_mem_double_cont: allocate a *contiguous* 2d array.
 * Array = M[n x m]
 */
double *allocate_mem_double_cont(double ***array, int n, int m)
{
        double *ptr = NULL; /* Data */
        int i = 0;
        if ((ptr = (double *) malloc(n * m * sizeof(double))) == NULL) {
                error_flag(mpi_proc_rank, (n * m), "allocate_mem_double_cont");
                return ptr;
        }
        *array = (double **) malloc(m * sizeof(double *));
        for (i = 0; i < m; i++) {
                (*array)[i] = ptr + i * n;
        }
        return ptr;
}

/*
 * allocate_mem_double_cont2: allocate a *contiguous* 2d array. Returns
 * 1-d array. Input **array is pointers to columns.
 * Array = M[n x m], n = rows, m = columns
 */
double *allocate_mem_double_cont2(double **array, int n, int m)
{
        double *ptr = NULL; /* Data */
        int i = 0;
        if ((ptr = (double *) malloc(n * m * sizeof(double))) == NULL) {
                error_flag(mpi_proc_rank, (n * m), "allocate_mem_double_cont");
                return ptr;
        }
        for (i = 0; i < m; i++) {
                array[i] = ptr + (i * n);
        }
        return ptr;
}

/*
 * allocate_mem_int: allocate a 2d integer array
 * Array = M[n x m]
 */
int allocate_mem_int(int ***array, int n, int m)
{
        int i = 0;
        *array = (int **) malloc(sizeof(int *) * m);
        for (i = 0; i < m; i++) {
                if (((*array)[i] = (int *) malloc(sizeof(int) * n)) == NULL)
                {
                        error_flag(mpi_proc_rank, i, "allocate_mem_int");
                        return i;
                }
        }
        i = 0;
        return i;
}

/*
 * allocate_mem_int_cont: allocate a *contiguous* 2d array.
 * Array = M[n x m]
 */
int *allocate_mem_int_cont(int ***array, int n, int m)
{
        int *ptr = NULL; /* Data */
        int i = 0;
        if ((ptr = (int *) malloc(n * m * sizeof(int))) == NULL) {
                error_flag(mpi_proc_rank, (n * m), "allocate_mem_double_cont");
                return ptr;
        }
        *array = (int **) malloc(m * sizeof(int *));
        for (i = 0; i < m; i++) {
                (*array)[i] = ptr + i * n;
        }
        return ptr;
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

/*
 * deallocate_mem_int: deallocate 2d int array.
 */
void deallocate_mem_int(int ***array, int m)
{
    for (int i = 0; i < m; i++) {
        free((*array)[i]);
    }
    free(*array);
    return;
}

/*
 * deallocate_mem_cont: deallocate *contiguous* 2d array
 */
void deallocate_mem_cont(double ***array, double *ptr)
{
        free(ptr);
        free(*array);
}
/*
 * deallocate_mem_cont_int: deallocate *contiguous* 2d array (C_INT)
 */
void deallocate_mem_cont_int(int ***array, int *ptr)
{
        free(ptr);
        free(*array);
}

