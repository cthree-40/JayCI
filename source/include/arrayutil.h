// File: arrayutil.h
#ifndef arrayutil_h
#define arrayutil_h

/*
 * cparray_1d1d: copy contents of 1-d array to 1-d array.
 */
void cparray_1d1d(
	double *array1, /* Array to copy from */
	int len1, /* Length of array 1 */
	double *array2, /* Array to copy TO */
	int len2 /* Length of array 2 */
	);

/*
 * cparray_1d2d: copy contents of 1-d array to 2-d array.
 */
void cparray_1d2d(
    double *array_1d,
    double **array_2d,
    int rows,
    int cols
    );

/*
 * cparray_2d1d: copy contents of 2-d array to 1-d array.
 */
void cparray_2d1d(
    double **array_2d,
    int rows,
    int cols,
    double *array_1d
    );

/*
 * cparray_2d2d: copy contents of 2-d array into 2-d array.
 */
void cparray_2d2d(
	double **array1, /* Array to copy FROM */
	int rows1, /* Rows of array 1 */
	int cols1, /* Columns of array 1 */
	double **array2, /* Array to copy TO */
	int rows2, /* Rows of array 2 */
	int cols2 /* Columns of array 2 */
	);

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

/*
 * init_dbl_array_1: initialize doulbe array to 1
 */
void init_dbl_array_1(
        double *array,
        int len
        );

/* 
 * init_dbl_2darray_0: initialize double  array to 0
 */
void init_dbl_2darray_0(
    double **array,
    int len,
    int ndim
    );

void init_int_array_0(
	int *array,
	int len);

int ndiffs_array(
	int *ar1,
	int *ar2,
	int ard1,
	int ard2);

/*
 * sort_array: sort ordered (except for first or last element) array.
 */
int sort_array(
	int updown,
	int *list,
	int list_len,
	int *new_loc
	);

/*
 * sort_array_fast_onesub: sort array with one index out of order. Returns
 * parity of sort and updates index to new location in list.
 */
int sort_array_fast_onesub(
	int *list,
	int list_len,
	int *index
	);

#endif
