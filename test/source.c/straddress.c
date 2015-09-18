// File: straddress.c
/********************************************************************
 * straddress
 * ----------
 * Contains string addressing utilities
 *
 * Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "combinatorial.h"
#include "straddress.h"
#include "binary.h"

/* subfunctions
 * ------------
 * int  str2adr(int *ostr, int elec, int orbitals)
 */

void str2det(int *astr, int *bstr, int *dstr, int orb)
{
     int i;

     /* build high order bit */
     for (i = 0; i < orb; i++) {
	  dstr[i] = astr[i];
     }
     /* build low order bit */
     for (i = 0; i < orb; i++) {
	  dstr[i + orb] = bstr[i];
     }

     return;
}
