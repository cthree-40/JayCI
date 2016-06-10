/* test.c */
/*
 * Test driver for subfunctions of jayci.x
 */

#include <stdio.h>
#include <stdlib.h>
#include "arrayutil.h"
#include "allocate_mem.h"
#include "mathutil.h"
#include "ioutil.h"

int main()
{
	int *list = NULL;
	int i = 0;
	int indx1;
	int indx2;

	list = malloc(10 * sizeof(int));
	indx = (int) 4;
	
	for (i = 0; i < 10; i++) {
		list[i] = i + 3;
	}

	list[indx] = 0;

	for (i = 0; i < 10; i++) {
		printf(" %d", list[i]);
	}
	printf("\n");

	i = sort_array_fast_onesub(list, 10, &indx);
	printf("pindx = %d\n", i);
	for (i = 0; i < 10; i++) {
		printf(" %d", list[i]);
	}
	printf("\n");

	return 0;
}
