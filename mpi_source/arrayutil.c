// File: arrayutil.c
/*
 * arrayutil.c
 * -----------
 * Subfunctions for performing operations on arrays.
 *
 * cparray_1d1d: copy contents of 1d array into 1d array.
 * cparray_2d1d: copy contents of 2d array into 1d array.
 * cparray_2d2d: copy contents of 2d array into 2d array.
 * cparray_1d2d: copy contents of 1d array into 2d array.
 * find_pos_in_array_lnsrch: find position of integer in array. linear search.
 * init_dbl_array_0: initialize double  array to zero
 * init_dbl_2darray_0: initialize double 2d array to zero
 * init_int_array_0: initialize integer array to zero
 * ndiffs_array: find number of differences between two arrays
 * sort_array_fast_onesub: sort array with one index out of order
 *
 * NOTE: Concerning 2-d arrays: rows and columns correspond to the array in
 * "true" matrix form. This is the OPPOSITE of the arrays ACTUAL form in a
 * C program. This is done to keep array references consistent with how they
 * are used in the master program and in f90 subroutines.
 *
 * By Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 */

#include <stdio.h>
#include "iminmax.h"
#include "arrayutil.h"

/*
 * cparray_1d1d: copy contents of 1-d array to 1-d array.
 */
void cparray_1d1d(double *array1, int len1, double *array2, int len2)
{
	int i;
	/* Test bounds */
	if (len1 > len2) {
		fprintf(stderr, "WARNING ");
		fprintf(stderr, "cparray_1d1d: %d > %d!\n", len1, len2);
	}
	for (i = 0; i < len2; i++) {
		array2[i] = array1[i];
	}
	return;
}

/*
 * cparray_1d2d: copy contents of 1-d array to 2-d array.
 */
void cparray_1d2d(double *array_1d, double **array_2d, int rows, int cols)
{
	int i, j;
	int ptr;
	for (i = 0; i < cols; i++) {
		for (j = 0; j < rows; j++) {
			ptr = (i * rows) + j;
			array_2d[i][j] = array_1d[ptr];
		}
	}
	return;
}

/*
 * cparray_2d1d: copy contents of 2-d array to 1-d array.
 */
void cparray_2d1d(double **array_2d, int rows, int cols, double *array_1d)
{
	int i, j;
	int ptr;
	for (i = 0; i < cols; i++) {
		for (j = 0; j < rows; j++) {
			ptr = (i * rows) + j;
			array_1d[ptr] = array_2d[i][j];
		}
	}
	return;
}

/*
 * cparray_2d2d: copy contents of 2-d array into 2-d array.
 */
void cparray_2d2d(double **array1, int rows1, int cols1, double **array2,
		  int rows2, int cols2)
{
	int i, j;
	/* Check bounds */
	if (rows1 > rows2) {
		fprintf(stderr, "WARNING: ");
		fprintf(stderr, "cparray_2d2d: (rows) %d > %d\n",
			rows1, rows2);
	}
	if (cols1 > cols2) {
		fprintf(stderr, "WARNING: ");
		fprintf(stderr, "cparray_2d2d: (cols) %d > %d\n",
			cols1, cols2);
	}
	/* copy array */
	for (i = 0; i < cols2; i++) {
		for (j = 0; j < rows2; j++) {
			array2[i][j] = array1[i][j];
		}
	}
	return;
}


/* 
 * find_pos_in_array_lnsrch: find position of integer in array. linear search.
 */
int find_pos_in_array_lnsrch(int look, int *array, int len)
{
	int i = 0;
	for (i = 0; i < len; i++) {
		if (look == array[i]) {
			return i;
		}
	}
	return (-1);
}

/* 
 * init_dbl_array_0: initialize double  array to 0
 */
void init_dbl_array_0(double *array, int len)
{
	int i;
	for (i = 0; i < len; i++) {
		array[i] = 0.0;
	}
	return;
}

/*
 * init_dbl_array_1: initialize doulbe array to 1
 */
void init_dbl_array_1(double *array, int len)
{
        int i;
        for (i = 0; i < len; i++) {
                array[i] = 1.0;
        }
        return;
}

/* 
 * init_dbl_2darray_0: initialize double  array to 0
 */
void init_dbl_2darray_0(double **array, int len, int ndim)
{
	int i;
	for (i = 0; i < ndim; i++) {
	    init_dbl_array_0(array[i], len);
	}
	return;
}

/* 
 * init_int_array_0: initialize integer array to 0
 */
void init_int_array_0(int *array, int len)
{
	int i;
	for (i = 0; i < len; i++) {
		array[i] = 0;
	}
	return;
}

/* 
 * ndiffs_array: find number of similarities between two arrays
 * 
 * Input:
 *  ar1 = array 1
 *  ar2 = array 2
 *  ard1 = array 1 dimension
 *  ard2 = array 2 dimension
 * Returns:
 *  ndiffs = number of differences 
 */
int ndiffs_array(int *ar1, int *ar2, int ard1, int ard2)
{
	int ndiffs = 0;
	int test;
	int i, j;
	
	for (i = 0; i < ard1; i++) {
		test = 0;
		for (j = 0; j < ard2; j++) {
			if (ar1[i] != ar2[j]) {
				test++;
			} else {
				test = 0;
				break;
			}
		}
		if (test != 0) ndiffs++;
	}
	return ndiffs;
}

/*
 * remove_leq_int: remove values from array less than or equal to input val
 * On output, *nelem = new number of values.
 * Input:
 *   val   = value to compare
 *   list  = list of (int)
 *   nelem = number of elemetns in list (changed on output)
 *   scr   = scratch array to hold values
 */
void remove_leq_int(int val, int *list, int *nelem, int *scr)
{
    int len = *nelem;
    int num = 0; /* Number of elements > val */
    int i;
    for (i = 0; i < len; i++) {
        /* If greater than add to scr */
        if (list[i] > val) {
            scr[num] = list[i];
            num++;
        }
        /* Zero out list as you go */
        list[i] = 0;
    }
    /* Copy over elements */
    for (i = 0; i < num; i++) {
        list[i] = scr[i];
    }
    /* Reset nelem to reflect new list */
    *nelem = num;
    return;
}

/*
 * sort_array: sort ordered (except for first or last element) array.
 */
int sort_array(int updown, int *list, int list_len, int *new_loc)
{
	int i;
	int tmp;
	int pindx = 1; /* permutational index */
	if (updown > 0) {
		/* first element is out of order */
		*new_loc = 0;
		i = *new_loc + 1;
		while (list[*new_loc] > list[i] && i < list_len) {
			/* swap elements */
			tmp = list[i];
			list[i] = list[*new_loc];
			list[*new_loc] = tmp;
			/* iterate i and new_loc */
			(*new_loc)++;
			i++;
			/* multiply parity by -1 */
			pindx = pindx * (-1);
		}
	} else if (updown < 0) {
		/* last element is out of order */
		*new_loc = list_len - 1;
		i = *new_loc - 1;
		while(list[*new_loc] < list[i] && i >= 0) {
			/* swap elements */
			tmp = list[i];
			list[i] = list[*new_loc];
			list[*new_loc] = tmp;
			/* iteratate i and new_loc */
			(*new_loc)--;
			i--;
			/* multiply parity by -1 */
			pindx = pindx * (-1);
		}
	}
	return pindx;
}

/*
 * sort_array_fast_onesub: sort array with one index out of order. Returns
 * parity of sort and updates index to new location in list.
 */
int sort_array_fast_onesub(int *list, int list_len, int *index)
{
	int htest = 0, ltest = 0; /* bounds limits for list */
	int new_loc = 0; /* new location of list[index] */
	int pindx = 1; /* permutational index */

	/* Get bounds limits */
	htest = int_min((list_len - 1), (*index + 1));
	ltest = int_max(0, (*index - 1));
	
	/* Check sort direction */
	if (list[*index] > list[htest]) {
		/* Sorting up */
		pindx = sort_array(1, (list+(*index)), (list_len - *index),
				   &new_loc);
		*index = *index + new_loc;
	} else if (list[*index] < list[ltest]) {
		/* Sorting down */
		pindx = sort_array(-1, list, (*index + 1), &new_loc); 
		*index = new_loc;
	}
    
	return pindx;
}
