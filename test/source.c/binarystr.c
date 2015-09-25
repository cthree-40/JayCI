// File: binarystr.c
/*********************************************************************
 * binarystr.c
 * -----------
 * Subfunctions for binary digit implementation of orbital strings
 *
 * str2occstr: convert orbital index string to occstr type
 * comparedets: compare two determinants 
 *
 * By Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bitutil.h"
#include "binarystr.h"

/* comparedets: compare two determinants, returning differences
 * -------------------------------------------------------------------
 * Input:
 *  deti = determinant 1
 *  detj = determinant j
 * Output:
 *  axi = alpha excitations initial orbitals
 *  axf = alpha excitations final orbitals
 *  bxi = beta  excitations initial orbitals
 *  bxf = beta  excitations final orbitals
 * Returns:
 *  numx = number of excitations */
int comparedets(struct det deti, struct det detj,
		long long int *axi, long long int *axf,
		long long int *bxi, long long int *bxf)
{
    int numx;

    /* .. local scalars ..
     * same  = output value of nsamebytes() calls
     * diffs = output of ndiffbytes calls 
     * sames = output of nsamebytes calls */
    int same; 
    long long int diffs, sames;

    numx = 0;
    /* compare alpha strings */
    numx = ndiffbytes(deti.astr.byte1, detj.astr.byte1, &diffs);
    same = nsamebytes(deti.astr.byte1, diffs, axi);
    same = nsamebytes(detj.astr.byte1, diffs, axf);
    
    /* compare beta strings */
    numx = numx + ndiffbytes(deti.bstr.byte1, detj.bstr.byte1, &diffs);
    same = nsamebytes(deti.bstr.byte1, diffs, bxi);
    same = nsamebytes(detj.bstr.byte1, diffs, bxf);
    
#ifdef BIGCAS
    if ((ndocc + nactv) > 64) {
	numx = numx + ndiffbytes(deti.astr.byte2, detj.astr.byte2, &diffs);

	numx = numx + ndiffbytes(deti.bstr.byte2, detj.bstr.byte2, &diffs);
    }
#endif

    numx = numx / 2;

    return numx;
}
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
#ifndef BIGCAS
     for (i = 0; i < elec; i++) {
	 if (istr[i] <= (ndocc + nactv)) {
	     ostr.byte1 = ostr.byte1 + pow(2, (istr[i] - 1));
	 } else {
	     ostr.virtx[vptr] = istr[i];
	     vptr++;
	 }
     }

     return ostr;
#endif

#ifdef BIGCAS
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
#endif
}

     
