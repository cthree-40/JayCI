// File: arrayutil.h
#ifndef arrayutil_h
#define arrayutil_h

/* 
 * find_pos_in_array_lnsrch: find position of integer in array. linear search.
 */
int find_pos_in_array_lnsrch(
	int look,   /* integer to look for */ 
	int *array, /* array to look in */
	int len);   /* length of array */

/* 
 * init_dbl_array_0: initialize a double array to 0.0
 */
void init_dbl_array_0(
	double *array,  /* array to initialize */
	int len);       /* dimension of array  */

void init_int_array_0(
	int *array,
	int len);

int ndiffs_array(
	int *ar1,
	int *ar2,
	int ard1,
	int ard2);

#endif
