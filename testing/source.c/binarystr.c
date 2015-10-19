// File: binarystr.c
/*********************************************************************
 * binarystr.c
 * -----------
 * Subfunctions for binary digit implementation of orbital strings
 *
 * comparedets_cas: compare two cas-flagged determinants 
 * comparedets_ncas: compare two non-cas-flagged determinants 
 * str2occstr: convert orbital index string to occstr type
 * 
 *
 * By Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrayutil.h"
#include "bitutil.h"
#include "binarystr.h"

/* 
 * comparedets_cas: compare CAS bytes of two determinants, 
 *                  returning differences
 * 
 * Input:
 *  deti  = determinant 1
 *  detj  = determinant j
 * Output:
 *  numax = number of alpha excitations
 *  numbx = number of beta  excitations
 *  axi   = alpha excitations initial orbitals
 *  axf   = alpha excitations final orbitals
 *  bxi   = beta  excitations initial orbitals
 *  bxf   = beta  excitations final orbitals
 * Returns:
 *  numx  = number of excitations 
 */
int comparedets_cas(struct det deti, struct det detj,
		    int *numax, int *numbx, long long int *axi, 
		    long long int *axf, long long int *bxi, 
		    long long int *bxf)
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
	*numax = ndiffbytes(deti.astr.byte1, detj.astr.byte1, &diffs);
	same = nsamebytes(deti.astr.byte1, diffs, axi);
	same = nsamebytes(detj.astr.byte1, diffs, axf);
	numx = numx + *numax;
	*numax = *numax / 2;
	
	/* compare beta strings */
	*numbx = ndiffbytes(deti.bstr.byte1, detj.bstr.byte1, &diffs);
	same = nsamebytes(deti.bstr.byte1, diffs, bxi);
	same = nsamebytes(detj.bstr.byte1, diffs, bxf);
	numx = numx + *numbx;
	*numbx = *numbx / 2;
	
#ifdef BIGCAS
	if ((ndocc + nactv) > 64) {
		*numax = *numax + 
			ndiffbytes(deti.astr.byte2, detj.astr.byte2, &diffs);
		*numbx = *numbx + 
			ndiffbytes(deti.bstr.byte2, detj.bstr.byte2, &diffs);
		numx = numx + *numax;
		numx = numx + *numbx;
	}
#endif
	/* Note: numx is NOT excitation number. It is the differences in 
	 * binary representations of the two bytes. You must divide the
	 * final numx, including virtual excitations , to get excitation 
	 * number. */
	return numx;
}

/* 
 * comparedets_ncas: compare determinants that are not cas-flaged
 * 
 * Input:
 *  deti  = determinant 1
 *  detj  = determinant j
 * Output:
 *  numaxc = number of CAS alpha excitations
 *  numbxc = number of CAS beta  excitations
 *  numaxv = number of virtual alpha excitations
 *  numbxv = number of virtual beta  excitations
 *  numaxcv = number of CAS<->virutal orbitals
 *  numbxcv = number of virtual<->CAS orbitals
 *  axi   = alpha excitations initial orbitals
 *  axf   = alpha excitations final orbitals
 *  bxi   = beta  excitations initial orbitals
 *  bxf   = beta  excitations final orbitals
 * Returns:
 *  numx  = number of excitations 
 */
int comparedets_ncas(struct det deti, struct det detj, 
		     int *numaxc,  int *numbxc,
		     int *numaxv,  int *numbxv,
		     int *numaxcv, int *numbxcv,
		     long long int *axi, long long int *axf, 
		     long long int *bxi, long long int *bxf)
{
	int numx = 0;
	/* .. local scalars ..
	 * samei = electrons in CAS of string of determinant i
	 * samej = electrons in CAS of string of determinant j
	 * diffs = output of ndiffbytes calls
	 * sames = output of nsamebytes calls */
	int samei = 0, samej = 0;
	long long int diffs = 0, sames = 0;
	
	*axi = 0;
	*bxi = 0;
	*axf = 0;
	*bxf = 0;

	*numaxv = compute_virt_diffs(deti.astr, detj.astr);
	*numbxv = compute_virt_diffs(deti.bstr, detj.bstr);
	numx = numx + *numaxv + *numbxv;
	if (numx > 2) return numx;
	
	/* Compare CAS byte. 
	 * After comparing check for CAS<->Virtual excitations. 
	 *
	 * This is done by comparing the number of electrons in the CAS 
	 * spaces of determinant i and determinant j. 
	 * 
	 * When we divide num{a,b}xc by 2 the result is rounded down
	 * giving us only CAS->CAS excitations, which is what we want. 
	 *
	 * We then subtract the CAS->Virt excitation number from 
	 * the differences in the virtual orbitals.*/
	*numaxc = ndiffbytes(deti.astr.byte1, detj.astr.byte1, &diffs);
	samei = nsamebytes(deti.astr.byte1, diffs, *(&axi));
	samej = nsamebytes(detj.astr.byte1, diffs, *(&axf));
	if (*axi == 0 || *axf == 0) {
		*numaxc = 0;
	}
	*numaxcv = abs(samei - samej);
	*numaxc = *numaxc / 2;
	*numaxv = *numaxv - *numaxcv;
	
	*numbxc = ndiffbytes(deti.bstr.byte1, detj.bstr.byte1, &diffs);
	samei = nsamebytes(deti.bstr.byte1, diffs, *(&bxi));
	samej = nsamebytes(detj.bstr.byte1, diffs, *(&bxf));
	if (*bxi == 0 || *bxf == 0) {
		*numbxc = 0;
	}
	*numbxcv = abs(samei - samej);
	*numbxc = *numbxc / 2;
	*numbxv = *numbxv - *numbxcv;
	
	numx = *numaxc + *numaxcv + *numaxv + *numbxc + *numbxcv + *numbxv;
	return numx;
}

/*
 * compute_virt_diffs: compute virtual orbital differences 
 * ------------------------------------------------------------------
 * Input:
 *  ostri = occupation string i
 *  ostrj = occupation string j
 * Returns:
 *  numxv = number of virtual orbital differences
 */
int compute_virt_diffs(struct occstr ostri, struct occstr ostrj)
{
	int numxv;
	int i;

	if (ostri.nvrtx >= ostrj.nvrtx ) {
		i = ostri.nvrtx;
	} else {
		i = ostrj.nvrtx;
	}
	
	if (ostri.nvrtx + ostrj.nvrtx == 0) {
		numxv = 0;
		return numxv;
	} else if (ostri.nvrtx == 0) {
		numxv = ostrj.nvrtx;
		return numxv;
	} else if (ostrj.nvrtx == 0) {
		numxv = ostri.nvrtx;
		return numxv;
	} else {
		/* both strings have virtual orbitals */
		numxv = ndiffs_array(
			ostri.virtx, ostrj.virtx, i, i);
		return numxv;
	}
}
/* 
 * str2occstr: convert orbital index string -> occstr type
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
     ostr.nvrtx = 0;
     
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
	     ostr.nvrtx = vptr;
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
	     ostr.nvrtx = vptr;
	 }
     }
     
     return ostr;
#endif
}

     
