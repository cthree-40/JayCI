// FILE: moindex.c
/*********************************************************************
 * moindex.c
 * ---------
 * Contains routines for index MO integrals.
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include "moindex.h"

int index1e(int i, int j)
/* index1e
 * -------
 * Index 1-e integral. Lower-triangle packing.
 */
{
     int result;
     
     if (i > j)
	  result = (((i - 1) * i) / 2) + j;
     else
	  result = (((j - 1) * j) / 2) + i;

     return result;
}

int index2e(int i, int j, int k, int l)
/* index2e
 * -------
 * Index 2-e integral. Lower triangle packing.
 */
{
     int result;

     result = index1e(index1e(i, j), index1e(k, l));

     return result;
}
