// File: action_util.c
/*
 * action_util
 * -----------
 * Utilities for performing Hv=c
 *
 * Subfunctions:
 *  
 *  hmatels: compute matrix element <i|H|j>
 *  make_orbital_strings_virt: make orbital strings with virtual occupations
 *  virtdiffs_1: find location of virtual orbital replacements
 * 
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 */

#include <stdio.h>
#include <math.h>
#include "bitutil.h"
#include "binarystr.h"
#include "action_util.h"

/* 
 * hmatels: compute matrix element <i|H|j>
 * ---------------------------------------
 * Input:
 *  deti = determinant i
 *  detj = determinant j
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 *  m1len   = length of 1-e integrals
 *  m2len   = length of 2-e integrals
 *  aelec   = alpha electrons
 *  belec   = beta  electrons
 * Returns:
 *  val = <i|H|j> 
 */
double hmatels(struct det deti, struct det detj, double *moints1,
	       double *moints2, int m1len, int m2len, int aelec,
	       int belec)
{
	double val;
    
	/* .. local scalars ..
	 * detdiff = differences between determinants
	 * numaxc = number of CAS alpha excitations
	 * numbxc = number of CAS beta excitations
	 * numaxv = number of virtual alpha excitations
	 * numbxv = number of virtual beta excitations
	 * axi = alpha excitations initial orbitals
	 * axf = alpha excitations final orbitals
	 * bxi = beta excitaitons initial orbitals
	 * bxf = beta excitations final orbitals */
    
	int detdiff; 
	int numaxc, numbxc, numaxv, numbxv, numaxcv, numbxcv;
	long long int axi, axf, bxi, bxf;
    
	val = 0.0;
	/* test if determinants are CAS-flagged. */ 
	if (deti.cas + detj.cas < 2) {
		detdiff = comparedets_ncas(
			deti, detj, &numaxc, &numbxc, &numaxv, &numbxv,
			&numaxcv, &numbxcv, &axi, &axf, &bxi, &bxf);
		printf(" detdiff = %d, numaxc = %d, numbxc = %d ", 
		       detdiff, numaxc, numbxc);
		printf(" numaxv = %d, numbxv = %d ", 
		       numaxv, numbxv);
		printf(" numaxcv = %d, numbxcv = %d ",
		       numaxcv, numbxcv);
		if (detdiff > 2) return val;
		val = 2.0;
	} else {
		detdiff = comparedets_cas(
			deti, detj, &numaxc, &numbxc, &axi, &axf, &bxi, &bxf);
		detdiff = detdiff / 2;
		printf(" detdiff = %d, numaxc = %d, numbxc = %d ", 
		       detdiff, numaxc, numbxc);
		if (detdiff > 2) return val;
		val = evaluate_dets_cas(
			detdiff, deti, detj, numaxc, numbxc, axi, axf,
			bxi, bxf, aelec, belec, moints1, moints2);
	}

	return val;
	
}

/* 
 * evaluate_dets_cas: evaluate matrix element <i|H|j>.
 * ---------------------------------------------------
 * Input:
 *  ndiff = number of differences (excitation level)
 *  deti  = determinant i
 *  detj  = determinant j
 *  numax = number of alpha excitations
 *  numbx = number of beta  excitations
 *  axi   = alpha excitation initial orbitals
 *  axf   = alpha excitation final orbitals
 *  bxi   = beta  excitation initial orbitals
 *  bxf   = beta  excitation final orbitals
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 * Returns:
 *  value = value of matrix element 
 */
double evaluate_dets_cas(int ndiff, struct det deti, struct det detj, 
			 int numax, int numbx, long long int axi, 
			 long long int axf, long long int bxi, long long int bxf, 
			 int aelec, int belec, double *moints1, double *moints2)
{
	double value;
	
	/* check number of differences and call appropriate subfunction */
	if (ndiff == 2) {
		/* 1,1 or 2,0 or 0,2 */
		if (numax == 1) {
			value = eval2_11_cas(axi, axf, bxi, bxf, moints2);
		} else if (numax == 2) {
			value = eval2_20_cas(axi, axf, moints2);
		} else {
			value = eval2_20_cas(bxi, bxf, moints2);
		}
	} else if (ndiff == 1) {
		if (numax == 1) {
			value = eval1_10_cas(
				deti.astr, axi, axf, deti.bstr, aelec, belec, 
				moints1, moints2);
		} else {
			value = eval1_10_cas(
				deti.bstr, bxi, bxf, deti.astr, belec, aelec,
				moints1, moints2);
		}
	} else {
		value = eval0_cas(deti, aelec, belec, moints1, moints2);
	}
	
	return value;
	
}

/* 
 * eval0_cas: evaluate diagonal elements that are CAS-flagged
 * ----------------------------------------------------------
 * Input:
 *  deti = determinant i
 *  detj = determinant j
 *  aelec = alpha electrons
 *  belec = beta electrons
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 * Returns:
 *  val = <i|H|i> 
 */
double eval0_cas(struct det deti, int aelec, int belec, 
		 double *moints1, double *moints2)
{
	double val;
	
	/* .. local scalars ..
	 * i1, i2 = integral indexes */
	int i1, i2;
	int i, j;
	
	/* .. local arrays ..
	 * eostr1 = electron orbital occupation string 1
	 * eostr2 = electron orbital occupation string 2 */
	int eostr1[aelec], eostr2[belec];
	
	/* form eostr1 and eostr2 */
#ifndef BIGCAS
	nonzerobits(deti.astr.byte1, eostr1);
	nonzerobits(deti.bstr.byte1, eostr2);
#endif

	val = 0.0;
	/* compute alpha contribution */
	for (i = 0; i < aelec; i++) {
		/* one electron integral contribution */
		i1 = index1e(eostr1[i], eostr1[i]);
		val = val + moints1[i1 - 1];
		for (j = 0; j < i; j++) {
			/* two electron integral contribution */
			i1 = index2e(eostr1[i], eostr1[i], eostr1[j], eostr1[j]);
			i2 = index2e(eostr1[i], eostr1[j], eostr1[i], eostr1[j]);
			val = val + (moints2[i1 - 1] - moints2[i2 - 1]);

		}
	}
	/* compute beta contribution */
	for (i = 0; i < belec; i++) {
		/* one electron integral contribution */
		i1 = index1e(eostr2[i], eostr2[i]);
		val = val + moints1[i1 - 1];
		for (j = 0; j < i; j++) {
			/* two electron integral contribution */
			i1 = index2e(eostr2[i], eostr2[i], eostr2[j], eostr2[j]);
			i2 = index2e(eostr2[i], eostr2[j], eostr2[i], eostr2[j]);
			val = val + (moints2[i1 - 1] - moints2[i2 - 1]);
		}
	}
	/* compute alpha and beta contribution */
	for (i = 0; i < aelec; i++) {
		for (j = 0; j < belec; j++) {
			i1 = index2e(eostr1[i], eostr1[i], eostr2[j], eostr2[j]);
			val = val + moints2[i1 - 1];
		}
	}
	return val;
}

/* 
 * eval1_10_cas: evaluate the matrix element of a single replacement
 * -----------------------------------------------------------------
 * Input:
 *  ostr1 = orbital occupation string of determinant i of (alpha/beta) electrons
 *  xi    = initial excitation orbital of (alpha/beta) electron
 *  xf    = final excitation orbital of (alpha/beta) electrons
 *  ostr2 = orbital occupation string of determinant i of (beta/alpha) electrons
 *  ne1   = number of (alpha/beta) electrons
 *  ne2   = number of (beta/alpha) electrons
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 * Returns:
 *  val = value of matrix element 
 */
double eval1_10_cas(struct occstr ostr1, long long int xi, long long int xf,
		    struct occstr ostr2, int ne1, int ne2, double *moints1,
		    double *moints2)
{
	double val;
	
	/* .. local scalars ..
	 * pindx = permuational index 
	 * io = initial orbital
	 * fo = final orbital 
	 * i1, i2 = integral indexes */
	int pindx;
	int io, fo;
	int i1, i2;
	
	/* .. local arrays ..
	 * eostr1 = electron orbital occupation string 1
	 * eostr2 = electron orbital occupation string 2 */
	int eostr1[ne1];
	int eostr2[ne2];
	
	int i;
	
	/* locate nonzero bits in xi and xf, and ostr1 and ostr2*/
	nonzerobits(xi, &io);
	nonzerobits(xf, &fo);
#ifndef BIGCAS
	nonzerobits(ostr1.byte1, eostr1);
	nonzerobits(ostr2.byte1, eostr2);
#endif
	
	/* compute permuation index */
	pindx = pow(-1, abs(fo - io));
	
	/* compute 1-e contribution */
	i1 = index1e(io, fo);
	val = pindx * moints1[i1 - 1];
	
	/* compute 2-e contributions */
	for (i = 0; i < ne1; i++) {
		if (eostr1[i] != io) {
			i1 = index2e(eostr1[i], eostr1[i], io, fo);
			i2 = index2e(eostr1[i], io, eostr1[i], fo);
			val = val + pindx * (moints2[i1 - 1] - moints2[i2 - 1]);
		}
	}
	for (i = 0; i < ne2; i++) {
		i1 = index2e(eostr2[i], eostr2[i], io, fo);
		val = val + pindx * moints2[i1 - 1];
	}
	
	return val;
}
    
/* 
 * eval2_11_cas: evaluate the matrix element of one replacement in two strings
 * ---------------------------------------------------------------------------
 * Input:
 *  axi = alpha initial orbitals
 *  axf = alpha final orbitals
 *  bxi = beta initial orbitals
 *  bxf = beta final orbitals
 *  moints2 = 2-e integrals
 * Returns:
 *  val = <i|H|j> = (axi(1),bxi(1)|axf(1),bxf(1)) 
 */
double eval2_11_cas(long long int axi, long long int axf, long long int bxi,
		    long long int bxf, double *moints2)
{
	double val;
    
	/* .. local scalars ..
	 * aio = alpha initial orbital
	 * afo = alpha final orbital
	 * bio = beta intital orbital
	 * bfo = beta final orbital
	 * pindx = alpha permuational index
	 * i1,i2 = integral indexes */
	int aio, afo, bio, bfo;
	int pindx;
	int i1, i2;
	
	/* locate nonzero bits in axi, axf, bxi, bxf */
	nonzerobits(axi, &aio);
	nonzerobits(axf, &afo);
	nonzerobits(bxi, &bio);
	nonzerobits(bxf, &bfo);
	
	/* compute permuation index */
	pindx = pow(-1, (abs(afo - aio) + abs(bfo - bio)));
	
	/* compute matrix element */
	i1 = index2e(aio, bio, afo, bfo);
	i2 = index2e(aio, afo, bio, bfo);
	val = pindx * (moints2[i1 - 1] - moints2[i2 - 1]);
	
	return val;
}

/* 
 * eval2_20_cas: evaluate the matrix element of two replacements in one string
 * ---------------------------------------------------------------------------
 * Input:
 *  xi = initial orbitals
 *  xf = final  orbitals
 *  moints2 = 2-e integrals **bounds not passed!!
 * Returns:
 *  val = <i|H|j> = (xi(1),xi(2)|xf(1),xf(2)) 
 */
double eval2_20_cas(long long int xi, long long int xf, double *moints2)
{
	double val;
	
	/* .. local scalars ..
	 * pindx = permuation index 
	 * i1,i2 = integral indexes */
	int pindx;
	int i1, i2;
	
	/* .. local arrays ..
	 * init_orbs = initial orbitals
	 * finl_orbs = final orbitals */
	int init_orbs[2], finl_orbs[2];
	
	/* locate nonzero bits in xi and xf to get initial and final orbitals */
	nonzerobits(xi, init_orbs);
	nonzerobits(xf, finl_orbs);
	
	/* compute permuation index */
	pindx = pow(-1, (abs(init_orbs[0] - finl_orbs[0]) 
			 + abs(init_orbs[1] - finl_orbs[1])));
	
	/* compute matrix element */
	i1 = index2e(init_orbs[0], init_orbs[1], finl_orbs[0], finl_orbs[1]);
	i2 = index2e(init_orbs[0], finl_orbs[0], init_orbs[1], finl_orbs[1]);
	val = pindx * (moints2[i1 - 1] - moints2[i2 - 1]);
	
	return val;
}

/* 
 * make_orbital_strings_virt: make orbital strings with virtual occupations
 * -----------------------------------------------------------------------------
 * Input:
 *  ostr1i = orbital occupation string
 * Output:
 *  eostr1 = electron orbital occupation string
 */
void make_orbital_strings_virt(struct occstr ostr1i, int *eostr1, int nelec1)
{
#ifndef BIGCAS
	nonzerobits(ostr1i.byte1, eostr1);
#endif
	if (eostr1[nelec1 - 2] == 0) {
		eostr1[nelec1 - 2] = ostr1i.virtx[0];
		eostr1[nelec1 - 1] = ostr1i.virtx[1];
	} else {
		eostr1[nelec1 - 1] = ostr1i.virtx[0];
	}
	return;
}

/* virtdiffs_1: find location of virtual orbital replacements
 * -------------------------------------------------------------------
 * Input:
 *  vxi = deti virtual orbitals
 *  vxj = detj virtual orbitals
 * Output:
 *  ifo = differences */
void virtdiffs_1(int *vxi, int *vxj, int *ifo)
{
	/* locate where differences are 
	 * possible differences:
	 * b c  b a  a a
	 * a a  a c  b c */
	if (vxi[0] == vxj[0]) {
		ifo[0] = vxi[1];
		ifo[1] = vxj[1];
	} else {
		if (vxi[0] == vxj[1]) {
			ifo[0] = vxi[1];
			ifo[1] = vxj[0];
		} else {
			ifo[0] = vxi[0];
			ifo[1] = vxj[0];
		}
	}
}
