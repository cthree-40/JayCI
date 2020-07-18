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
#include "iminmax.h"

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
		    long long int *bxf, int nactv)
{
	int numx = 0;
	int samei = 0, samej = 0; 
	long long int diffs = 0;
	
	numx = 0;
        *axi = 0;
        *axf = 0;
        *bxi = 0;
        *bxf = 0;
	/* compare alpha strings */
	*numax = ndiffbytes(deti.astr.byte1, detj.astr.byte1, 
			    nactv, &diffs);
	samei = nsamebytes(deti.astr.byte1, diffs, nactv, axi);
	samej = nsamebytes(detj.astr.byte1, diffs, nactv, axf);
	*numax = int_min(samei, samej);
	numx = numx + *numax;
        if (numx > 2) return numx;
        
	/* compare beta strings */
	*numbx = ndiffbytes(deti.bstr.byte1, detj.bstr.byte1, 
			    nactv, &diffs);
	samei = nsamebytes(deti.bstr.byte1, diffs, nactv, bxi);
	samej = nsamebytes(detj.bstr.byte1, diffs, nactv, bxf);
	*numbx = int_min(samei, samej);
	numx = numx + *numbx;
	
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
		     long long int *bxi, long long int *bxf,
		     int nactv)
{
	int numx = 0;
	/* .. local scalars ..
	 * samei = electrons in CAS of string of determinant i
	 * samej = electrons in CAS of string of determinant j
	 * diffs = output of ndiffbytes calls
	 * sames = output of nsamebytes calls */
	int samei = 0, samej = 0;
	long long int diffs = 0;
	
	*axi = 0;
	*bxi = 0;
	*axf = 0;
	*bxf = 0;
	
	*numaxv = compute_virt_diffs(deti.astr, detj.astr);
	*numbxv = compute_virt_diffs(deti.bstr, detj.bstr);
	numx = numx + *numaxv + *numbxv;
	if (numx > 2) {
		return numx;
	} 
	
	/* get CAS->Virt excitations */
	*numaxcv = abs(deti.astr.nvrtx - detj.astr.nvrtx);
	*numbxcv = abs(deti.bstr.nvrtx - detj.bstr.nvrtx);
	
	if (*numaxcv + *numbxcv == 0 && numx == 2) {
		if ((deti.astr.byte1 ^ detj.astr.byte1) != 0x00) {
			numx+=10;
			return numx;
		}
		if ((deti.bstr.byte1 ^ detj.bstr.byte1) != 0x00) {
			numx+=10;
			return numx;
		}
		*numaxc = 0;
		*numbxc = 0;
		return numx;
	}
	
	/* Compare CAS byte. */
	*numaxc = ndiffbytes(deti.astr.byte1, detj.astr.byte1, 
			     nactv, &diffs);
	samei = nsamebytes(deti.astr.byte1, diffs, nactv, *(&axi));
	samej = nsamebytes(detj.astr.byte1, diffs, nactv, *(&axf));
	*numaxc = int_min(samei, samej);
	*numaxv = *numaxv - *numaxcv;
	
	numx = *numaxc + *numaxcv + *numaxv;
	if (numx > 2) return numx;
	
	*numbxc = ndiffbytes(deti.bstr.byte1, detj.bstr.byte1, 
			     nactv, &diffs);
	samei = nsamebytes(deti.bstr.byte1, diffs, nactv, *(&bxi));
	samej = nsamebytes(detj.bstr.byte1, diffs, nactv, *(&bxf));
	*numbxc = int_min(samei, samej);
	*numbxv = *numbxv - *numbxcv;
	
	numx = numx + *numbxc + *numbxcv + *numbxv;
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
	int numxv = 0;
	int i = 0;
	
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
 * get_string_eospace_info: get electon occupation space information given
 * a binary electron string.
 * -------------------------------------------------------------------
 * Input:
 *  str      = electron occupation string
 *  ndocc    = number of docc orbitals
 *  nactv    = number of CAS  orbtials
 * Output:
 *  nde = number of docc electrons
 *  nce = number of CAS electrons
 *  nve = number of virtual electrons
 */
void get_string_eospace_info(struct occstr str, int ndocc, int nactv, int *nde,
                             int *nce, int *nve)
{
        /* 64-bit byte to check docc and cas occupation */
        long long int intorb_check = 0x00;
        /* Scratch 64-bit byte */
        long long int scr = 0x00;
        int i;
        char tmp[65];
        /* Turn "on" DOCC and ACTV bits in intorb_check */
        for (i = 0; i < (ndocc + nactv); i++) {
                intorb_check = intorb_check + pow(2, i);
        }
        scr = intorb_check & str.byte1;
        /* Check DOCC orbitals */
        *nde = 0;
        for (i = 0; i < ndocc; i++) {
                if (scr & 0x01) {
                        *nde = *nde + 1;
                }
                scr = scr >> 1;
        }
        *nce = 0;
        /* Check ACTV orbitals */
        for (i = 0; i < nactv; i++) {
                if (scr &0x01) {
                        *nce = *nce + 1;
                }
                scr = scr >> 1;
        }
        *nve = str.nvrtx;
        return;
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
	struct occstr ostr = {.byte1 = 0x0, .virtx = {0},
			      .nvrtx = 0};
	int i;
	int vptr = 0;
	
	/* test if active space is larger than 64 orbitals.
	 * if less, occupation information can be stored entirely in
	 * ostr.byte1 */
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
}

/*
 * print_determinant: print a determinant
 */
void print_determinant(struct det d, int aelec, int belec)
{
        int stra[aelec];
        int strb[belec];
        int acnt = 0, bcnt = 0; /* Orbital occupation count */
        int i, j;
        
        init_int_array_0(stra, aelec);
        init_int_array_0(strb, belec);
        nonzerobits(d.astr.byte1, 64, stra);
        for (i = 0; i < aelec; i++) {
                if (stra[i] > 0) acnt++;
        }
        nonzerobits(d.bstr.byte1, 64, strb);
        for (i = 0; i < belec; i++) {
                if (strb[i] > 0) bcnt++;
        }

        j = 0;
        for (i = acnt; i < aelec; i++) {
                stra[i] = d.astr.virtx[j];
                j++;
        }
        j = 0;
        for (i = bcnt; i < belec; i++) {
                strb[i] = d.bstr.virtx[j];
                j++;
        }

        printf("a: ");
        for (i = 0; i < aelec; i++) {
                printf("%d ", stra[i]);
        }
        printf(" b: ");
        for (i = 0; i < belec; i++) {
                printf("%d ", strb[i]);
        }
        printf("\n");
}

/*
 * print_occstring: print occupation string
 *  byte1, virtx, nvrtx
 */
void print_occstring(struct occstr ostr, int nelec, int ndocc, int nactv)
{
        int *string = NULL;
        int i = 0, j = 0;
        int ecnt = 0;
        string = malloc(sizeof(int) * nelec);
        init_int_array_0(string, nelec);
        nonzerobits(ostr.byte1, (ndocc + nactv), string);
        for (i = 0; i < nelec; i++) {
                if (string[i] > 0) ecnt++;
        }
        for (i = ecnt; i < nelec; i++) {
                string[i] = ostr.virtx[j];
                j++;
        }
        
        printf("String:  ");
        for (i = 0; i < nelec; i++) {
                printf(" %d", string[i]);
        }
        printf("\n");
        free(string);
}

/*
 * init_detlist: initialize determinant list
 */
void init_detlist(struct det *dlist, int ndets)
{
	int i = 0;
	for (i = 0; i < ndets; i++) {
		init_occstr(dlist[i].astr);
		init_occstr(dlist[i].bstr);
		dlist[i].cas = 0;
	}
	return;
}

/*
 * init_occstr: initialize occupation string
 */
void init_occstr(struct occstr ostr)
{
	ostr.byte1 = 0x0;
	init_int_array_0(ostr.virtx, 2);
	ostr.nvrtx = 0;
	return;
}

     
