// File: arrayutil.c
/*
 * arrayutil.c
 * -----------
 * Subfunctions for performing operations on arrays.
 *
 * find_pos_in_array_lnsrch: find position of integer in array. linear search.
 * init_dbl_array_0: initialize double  array to zero
 * init_int_array_0: initialize integer array to zero
 * ndiffs_array: find number of differences between two arrays
 *
 * By Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 */

#include <stdio.h>
#include "arrayutil.h"

/* 
 * find_pos_in_array_lnsrch: find position of integer in array. linear search.
 */
int find_pos_in_array_lnsrch(int look, int *array, int len)
{
	int loc = 0; /* location */
	int i;
	for (i = 0; i < len; i++) {
		if (look == array[i]) {
			loc = i;
			return loc;
		}
	}
	return loc;
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
