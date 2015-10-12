// File: arrayutil.c
/*********************************************************************
 * arrayutil.c
 * -----------
 * Subfunctions for performing operations on arrays.
 *
 * ndiffs_array: find number of differences between two arrays
 *
 * By Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include "arrayutil.h"

/* ndiffs_array: find number of similarities between two arrays
 * -------------------------------------------------------------------
 * Input:
 *  ar1 = array 1
 *  ar2 = array 2
 *  ard1 = array 1 dimension
 *  ard2 = array 2 dimension
 * Returns:
 *  ndiffs = number of differences */
int ndiffs_array(int *ar1, int *ar2, int ard1, int ard2)
{
	int ndiffs = 0;
	int test;
	int i, j;
	
	for (i = 0; i < ard1; i++) {
		test = 0;
		for (j = 0; j < ard2; j++) {
			if (ar1[i] != ar2[j]) test++;
		}
		if (test != 0) ndiffs++;
	}
	return ndiffs;
}
