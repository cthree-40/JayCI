// File: initarray.c
#include <stdio.h>
#include <stdlib.h>
/* initarray: initializes an array, seting all elements to 0.0 */
void initarray(double *array, int dim)
{
    int i;
    for (i = 0; i < dim; i++) {
	array[i] = 0.0;
    }
    return;
}
