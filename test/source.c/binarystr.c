// File: binarystr.c
/*********************************************************************
 * binarystr.c
 * -----------
 * Subfunctions for binary digit implementation of orbital strings
 *
 * str2occstr: convert orbital index string to occstr type
 * 
 * By Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "binarystr.h"

/* str2occstr: convert orbital index string -> occstr type
 * -------------------------------------------------------------------
 * Input:
 *  istr = orbital index string
 *  elec = electrons
 * ndocc = number of doubly-occupied orbitals
 * nactv = number of active orbitals
 *  ostr = occupation string  */
struct occstr str2occstr(int *istr, int elec, int ndocc, int nactv)
{
     struct occstr ostr;
     int i;
     int vptr = 0;

     ostr.byte1 = 0;
     ostr.byte2 = 0;
     ostr.virtx[0] = 0;
     ostr.virtx[1] = 0;
     
     /* test if active space is larger than 64 orbitals.
      * if less, occupation information can be stored entirely in
      * ostr.byte1 */
     if (ndocc + nactv <= 64) {

	  for (i = 0; i < elec; i++) {
	       if (istr[i] <= (ndocc + nactv)) {
		    ostr.byte1 = ostr.byte1 + pow(2, (istr[i] - 1));
	       } else {
		    ostr.virtx[vptr] = istr[i];
		    vptr++;
	       }
	  }

	  return ostr;

     } else if (ndocc + nactv > 64) {

	  for (i = 0; i < elec; i++) {
	       if (istr[i] <= 64) {
		    ostr.byte1 = ostr.byte1 + pow(2, (istr[i] - 1));
	       } else if (istr[i] > 64 && istr[i] <= (ndocc + nactv)) {
		    ostr.byte2 = ostr.byte2 + pow(2, (istr[i] - 1));
	       } else {
		    ostr.virtx[vptr] = istr[i];
		    vptr++;
	       }
	  }

	  return ostr;

     } else {
	  fprintf(stderr,"*** ERROR: Active space > 128 orbitals! ***\n");
	  exit(1);
     }
	  
}

     
