/* File: action_util.c */
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
#include <stdlib.h>
#include <math.h>
#include "arrayutil.h"
#include "bitutil.h"
#include "binarystr.h"
#include "moindex.h"
#include "action_util.h"
#include "binary.h"

/* 
 * cas_to_virt_replacements: compute excitations for cas<->virt replacements
 */
void cas_to_virt_replacements(int ncreps, int ncr, int nvr, long long int xi, 
			      long long int xf, int *restrict vxi, int *restrict vxj, 
			      int *restrict reps, int ninto)
{
	if (ncreps == 2) {
		if (xi == 0) {
			reps[0] = vxi[0];
			reps[1] = vxi[1];
			nonzerobits(xf, ninto, &(reps[2]));
		} else {
			nonzerobits(xi, ninto, &(reps[0]));
			reps[2] = vxj[0];
			reps[3] = vxj[1];
		}
	} else {
		/* test if cas->virt is the only excitation */
		if (nvr + ncr == 0) {
			if (xi != 0) {
				nonzerobits(xi, ninto, &(reps[0]));
				virtdiffs_single_cas_to_virt(vxi, vxj, &(reps[2]));
			} else {
				nonzerobits(xf, ninto, &(reps[2]));
				virtdiffs_single_cas_to_virt(vxj, vxi, &(reps[0]));
			}
		} else if (ncr == 1) {
			nonzerobits(xi, ninto, &(reps[0]));
			if (reps[1] == 0) {
				nonzerobits(xf, ninto, &(reps[2]));
				virtdiffs_single_cas_to_virt(vxj, vxi, &(reps[1]));
			} else {
				nonzerobits(xf, ninto, &(reps[2]));
				virtdiffs_single_cas_to_virt(vxi, vxj, &(reps[3]));
			}
		} else {
			if (vxi[1] == 0) {
				nonzerobits(xi, ninto, &(reps[0]));
				reps[1] = vxi[0];
				reps[2] = vxj[0];
				reps[3] = vxj[1];
			} else {
				reps[0] = vxi[0];
				reps[1] = vxi[1];
				nonzerobits(xf, ninto, &(reps[2]));
				reps[3] = vxj[0];
			}
		}
	}
	return;
}

/*
 * compute_hv: perform Hv=c
 */
void compute_hv(struct det *dlist, int ndets, double *moints1, double *moints2,
		int aelec, int belec, double *restrict v, double *restrict c,
		int ninto)
{
	int i, j;
	double valij;

	for (i = 0; i < ndets; i++) {
		c[i] = c[i] + hmatels(dlist[i],dlist[i], moints1, moints2, aelec,
				      belec, ninto) * v[i];
		for (j = i+1; j < ndets; j++) {
			valij = hmatels(dlist[i], dlist[j], moints1, moints2,
					aelec, belec, ninto);
			c[i] = c[i] + valij * v[j];
			c[j] = c[j] + valij * v[i];
		}
	}
	return;
}

/* 
 * hmatels: compute matrix element <i|H|j>
 * ---------------------------------------
 * Input:
 *  deti = determinant i
 *  detj = determinant j
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 *  aelec   = alpha electrons
 *  belec   = beta  electrons
 * Returns:
 *  val = <i|H|j> 
 */
double hmatels(struct det deti, struct det detj, double *moints1,
	       double *moints2, int aelec, int belec, int nactv)
{
	double val = 0.0;
    
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
    
	int detdiff = 0; 
	int numaxc = 0, numbxc  = 0, numaxv  = 0, 
	    numbxv = 0, numaxcv = 0, numbxcv = 0;
	long long int axi = 0, axf = 0, bxi = 0, bxf = 0;
    
	/* test if determinants are CAS-flagged. */ 
	if (deti.cas + detj.cas < 2) {
		detdiff = comparedets_ncas(
			deti, detj, &numaxc, &numbxc, &numaxv, &numbxv,
			&numaxcv, &numbxcv, &axi, &axf, &bxi, &bxf, nactv);
		if (detdiff > 2) return val;
		val = evaluate_dets_ncas(
			detdiff, deti, detj, numaxc, numbxc, numaxcv,
			numbxcv, numaxv, numbxv, axi, axf, bxi, bxf,
			aelec, belec, moints1, moints2, nactv);
	} else {
		detdiff = comparedets_cas(
			deti, detj, &numaxc, &numbxc, &axi, &axf, &bxi, &bxf,
			nactv);
		if (detdiff > 2) return val;
		val = evaluate_dets_cas(
			detdiff, deti, detj, numaxc, numbxc, axi, axf,
			bxi, bxf, aelec, belec, moints1, moints2, nactv);
	}

	return val;
	
}

/* 
 * evaluate_dets_cas: evaluate matrix element <i|H|j>.
 *
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
			 int aelec, int belec, double *moints1, double *moints2,
			 int ninto)
{
	double value;
	
	/* check number of differences and call appropriate subfunction */
	if (ndiff == 2) {
		/* 1,1 or 2,0 or 0,2 */
		if (numax == 1) {
			value = eval2_11_cas(axi, axf, bxi, bxf, moints2, 
                                        ninto, deti.astr.byte1, 
                                        deti.bstr.byte1);
		} else if (numax == 2) {
			value = eval2_20_cas(axi, axf, moints2, ninto,
                                        deti.astr.byte1);
		} else {
			value = eval2_20_cas(bxi, bxf, moints2, ninto,
                                        deti.bstr.byte1);
		}
	} else if (ndiff == 1) {
		if (numax == 1) {
			value = eval1_10_cas(
				deti.astr, axi, axf, deti.bstr, aelec, belec, 
				moints1, moints2, ninto);
		} else {
			value = eval1_10_cas(
				deti.bstr, bxi, bxf, deti.astr, belec, aelec,
				moints1, moints2, ninto);
		}
	} else {
		value = eval0_cas(deti, aelec, belec, moints1, moints2, ninto);
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
		 double *moints1, double *moints2, int ninto)
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
	nonzerobits(deti.astr.byte1, ninto, eostr1);
	nonzerobits(deti.bstr.byte1, ninto, eostr2);
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
 */
double eval1_10_cas(struct occstr ostr1, long long int xi, long long int xf,
		    struct occstr ostr2, int ne1, int ne2, double *moints1,
		    double *moints2, int ninto)
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
	nonzerobits(xi, ninto, &io);
	nonzerobits(xf, ninto, &fo);
#ifndef BIGCAS
	nonzerobits(ostr1.byte1, ninto, eostr1);
	nonzerobits(ostr2.byte1, ninto, eostr2);
#endif
	
	/* compute permuation index */
	pindx = pindex_single_rep(eostr1, io, fo, ne1);
	
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
 */
double eval2_11_cas(long long int axi, long long int axf, long long int bxi,
		    long long int bxf, double *moints2, int ninto,
                    long long int abyte, long long int bbyte)
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
	nonzerobits(axi, ninto, &aio);
	nonzerobits(axf, ninto, &afo);
	nonzerobits(bxi, ninto, &bio);
	nonzerobits(bxf, ninto, &bfo);
	
	/* compute permuation index */
        
#ifndef BIGCAS
	pindx = pindex_single_rep_cas(abyte, axi, axf, ninto);
        pindx = pindx * pindex_single_rep_cas(bbyte, bxi, bxf, ninto);
#endif
	/* compute matrix element */
	i1 = index2e(aio, bio, afo, bfo);
	val = pindx * (moints2[i1 - 1]);
	
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
double eval2_20_cas(long long int xi, long long int xf, double *moints2,
	int ninto, long long int str)
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
	nonzerobits(xi, ninto, init_orbs);
	nonzerobits(xf, ninto, finl_orbs);
	
	/* compute permuation index */
	pindx = pindex_double_rep_cas(str, init_orbs, finl_orbs, ninto);
	
	/* compute matrix element */
	i1 = index2e(init_orbs[0], init_orbs[1], finl_orbs[0], finl_orbs[1]);
	i2 = index2e(init_orbs[0], finl_orbs[0], init_orbs[1], finl_orbs[1]);
	val = pindx * (moints2[i1 - 1] - moints2[i2 - 1]);
	
	return val;
}

/* 
 * evaluate_dets_ncas: evaluate <i|H|j> for non-cas flagged determinants
 *
 * Input:
 *  ndiff   = number of differences (excitation level)
 *  deti    = determinant i
 *  detj    = determinant j
 *  numaxc  = number of alpha cas excitations
 *  numbxc  = number of beta  cas excitations
 *  numaxcv = number of alpha cas->virtual excitations
 *  numbxcv = number of beta  cas->virtual excitations
 *  numaxv  = number of alpha virtual excitaitons
 *  numbxv  = number of beta  virtual excitations
 *  axi     = alpha excitation initial orbitals
 *  axf     = alpha excitation final orbitals
 *  bxi     = beta  excitation initial orbitals
 *  bxf     = beta  excitation final orbitals
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 * Returns:
 *  value   = value of matrix element 
 */
double evaluate_dets_ncas(int ndiff, struct det deti, struct det detj,
			  int numaxc, int numbxc, int numaxcv, int numbxcv,
			  int numaxv, int numbxv, long long int axi,
			  long long int axf, long long int bxi, 
			  long long int bxf, int aelec, int belec, 
			  double *moints1, double *moints2, int ninto)
{
	double value = 0.0;

	/* check for inter-space interactions then number of differences */
	if (numaxcv + numbxcv == 0) {
		if (ndiff == 2) {
			if (numaxv == 1) {
				/* x,1;x,x */
				if (numbxv == 1) {
					/* 0,1;0,1 */
					value = eval2_ncas_c00cv00v11(
						deti.astr.virtx, detj.astr.virtx,
						detj.bstr.virtx, detj.bstr.virtx,
						deti.astr.nvrtx, deti.bstr.nvrtx,
						moints2);
				} else if (numaxc == 1) {
					/* 1,1;0,0 */
					value = eval2_ncas_c1cv0v1(
						deti.astr, axi, axf, 
						deti.astr.virtx,
						detj.astr.virtx, aelec, 
						deti.astr.nvrtx, moints2, ninto);
				} else if (numbxc == 1) {
					/* 0,1;1,0 */
					value = eval2_ncas_c1cv0v1(
						deti.bstr, bxi, bxf, 
						deti.astr.virtx,
						detj.astr.virtx, belec,
						deti.astr.nvrtx, moints2, ninto);
				}
			} else if (numbxv == 1) {
				/* x,x;x,1 */
				if (numbxc == 1) {
					/* 0,0;1,1 */
					value = eval2_ncas_c1cv0v1(
						deti.bstr, bxi, bxf, 
						deti.bstr.virtx, 
						detj.bstr.virtx, belec, 
						deti.bstr.nvrtx, moints2, ninto);
				} else if (numaxc == 1) {
					/* 1,0;0,1 */
					value = eval2_ncas_c1cv0v1(
						deti.astr, axi, axf, 
						deti.bstr.virtx,
						detj.bstr.virtx, aelec, 
						deti.bstr.nvrtx, moints2, ninto);
				}
			} else if (numaxv == 2) {
				/* 0,2;0,0 */
				value = eval2_ncas_c0cv0v2(
					deti.astr.virtx, detj.astr.virtx,
					moints2);
			} else if (numbxv == 2) {
				/* 0,0;0,2 */
				value = eval2_ncas_c0cv0v2(
					deti.bstr.virtx, detj.bstr.virtx,
					moints2);
			} else if (numaxc == 1) {
				/* 1,0;1,0 */
				value = eval2_11_cas(
					axi, axf, bxi, bxf, moints2, ninto,
                                        deti.astr.byte1, deti.bstr.byte1);
			} else if (numaxc == 2) {
				/* 2,0;0,0 */
				value = eval2_20_cas(
					axi, axf, moints2, ninto,
                                        deti.astr.byte1);
			} else {
				/* 0,0;2,0 */
				value = eval2_20_cas(
					bxi, bxf, moints2, ninto,
                                        deti.bstr.byte1);
			}
		} else if (ndiff == 1) {
			if (numaxv == 1) {
				value = eval1_ncas_c0cv0v1(
					deti.astr, detj.astr, aelec, deti.bstr,
					belec, moints1, moints2, ninto);
			} else if (numbxv == 1) {
				value = eval1_ncas_c0cv0v1(
					deti.bstr, detj.bstr, belec, deti.astr, 
					aelec, moints1, moints2, ninto);
			} else if (numaxc == 1) {
				value = eval1_ncas_c1cv0v0(
					deti.astr, axi, axf, deti.bstr, aelec,
					belec, moints1, moints2, ninto);
			} else {
				value = eval1_ncas_c1cv0v0(
					deti.bstr, bxi, bxf, deti.astr, belec,
					aelec, moints1, moints2, ninto);
			}
		} else {
			value = eval0_ncas(deti, aelec, belec, moints1, moints2, ninto);
		}
	} else if (numaxcv == 1 ) {
		/* x,1,x;x,x,x */
		if (ndiff == 2) {
			if (numaxv == 1) {
				/* 0,1,1;0,0,0 */
				value = eval2_ncas_c00cv10v10(
					axi, axf, deti.astr, detj.astr,
					moints2, ninto);
			} else if (numbxv == 1) {
				/* 0,1,0;0,0,1 */
				value = eval2_ncas_c00cv10v01(
					axi, axf, deti.astr.virtx, 
					detj.astr.virtx, deti.bstr.virtx, 
					detj.bstr.virtx, moints2, ninto,
                                        deti.astr.byte1);
			} else if (numaxc == 1) {
				/* 1,1,0;0,0,0 */
				value = eval2_ncas_c10cv10v00(
					deti.astr,
					axi, axf, deti.astr.virtx, 
					detj.astr.virtx, aelec, moints2, ninto);
			} else if (numbxc == 1) {
				/* 0,1,0;1,0,0 */
				value = eval2_ncas_c01cv10v00(
					deti.astr, deti.bstr, axi, axf, 
					deti.astr.virtx, detj.astr.virtx, 
					bxi, bxf, aelec, belec, moints2, ninto);
			} else if (numbxcv == 1) {
				/* 0,1,0;0,1,0 */
				value = eval2_ncas_c00cv11v00(
					axi, axf, deti.astr.virtx,
					detj.astr.virtx, bxi, bxf,
					deti.bstr.virtx, detj.bstr.virtx,
					moints2, ninto, deti.astr.byte1,
                                        deti.bstr.byte1);
			}
		} else if (ndiff == 1) {
			/* 0,1,0;0,0,0 */
			value = eval1_ncas_c0cv1v0(
				axi, axf, deti.astr, detj.astr, 
				deti.bstr, aelec, belec, 
				moints1, moints2, ninto);
		} 
	} else if (numbxcv == 1) {
		/* x,x,x;x,1,x */
		if (ndiff == 2) {
			if (numaxv == 1) {
				/* 0,0,1;0,1,0 */
				value = eval2_ncas_c00cv10v01(
					bxi, bxf, deti.bstr.virtx, 
					detj.bstr.virtx, deti.astr.virtx, 
					detj.astr.virtx, moints2, ninto,
                                        deti.bstr.byte1);
			} else if (numbxv == 1) {
				/* 0,0,0;0,1,1 */
				value = eval2_ncas_c00cv10v10(
					bxi, bxf, deti.bstr, detj.bstr,
					moints2, ninto);
			} else if (numaxc == 1) {
				/* 1,0,0;0,1,0 */
				value = eval2_ncas_c01cv10v00(
					deti.bstr, deti.astr, bxi, bxf, 
					deti.bstr.virtx, detj.bstr.virtx, 
					axi, axf, belec, aelec,moints2,
					ninto);
			} else if (numbxc == 1) {
				/* 0,0,0;1,1,0 */
				value = eval2_ncas_c10cv10v00(
					deti.bstr, bxi, bxf, deti.bstr.virtx, 
					detj.bstr.virtx, belec, moints2, ninto);
			}
		} else if (ndiff == 1) {
			/* 0,0,0;0,1,0 */
			value = eval1_ncas_c0cv1v0(
				bxi, bxf, deti.bstr, detj.bstr, 
				deti.astr, belec, aelec, 
				moints1, moints2, ninto);
		}	
	} else if (numaxcv == 2) {
		value = eval2_ncas_c0cv2v0(
			axi, axf, deti.astr.virtx, detj.astr.virtx, moints2,
                        ninto, deti.astr.byte1);
	} else if (numbxcv == 2) {
		value = eval2_ncas_c0cv2v0(
			bxi, bxf, deti.bstr.virtx, detj.bstr.virtx, moints2, 
                        ninto, deti.bstr.byte1);
	}
	
	return value;
}

/* 
 * eval0_ncas: evaluate matrix elements <i|H|i> with virtual occupations
 *
 * Input:
 *  deti    = determinant
 *  aelec   = alpha electrons
 *  belec   = beta  electrons
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 * Returns:
 *  value = <i|H|i>
 */
double eval0_ncas(struct det deti, int aelec, int belec, double *moints1,
		  double *moints2, int ninto)
{
	double val = 0.0;
	
	/* .. local scalars ..
	 * i1, i2 = integral indexes */
	int i1, i2;
	int i, j;

	/* .. local arrays ..
	 * eostr1 = electron orbital occupation string 1
	 * eostr2 = electron orbital occupation string 2 */
	int eostr1[aelec], eostr2[belec];
	
	/* form eostr1 and eostr2 */
	make_orbital_strings_virt(deti.astr, eostr1, aelec, ninto);
	make_orbital_strings_virt(deti.bstr, eostr2, belec, ninto);
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
 * eval1_ncas_c0cv1v0: evaluate single cas->virt replacement in matrix 
 *                     elements with non-cas-flagged determinants
 */
double eval1_ncas_c0cv1v0(long long int xi, long long int xf, 
			  struct occstr str1i, struct occstr str1j,
			  struct occstr str2i, int ne1, int ne2,
			  double *moints1, double *moints2, int ninto)
{
	double val = 0.0;
	int ifo[4] = {0};    /* initial, final orbitals */
	int i1, pindx; /* integral, permuational index */
	int eostr1[ne1], eostr2[ne2]; /* electron occupation string */
	
	/* get replacement information */
	cas_to_virt_replacements(1, 0, 0, xi, xf, str1i.virtx, str1j.virtx,
				 ifo, ninto);
        /* make orbital strings */
	make_orbital_strings_virt(str1i, eostr1, ne1, ninto);
	make_orbital_strings_virt(str2i, eostr2, ne2, ninto);
	/* compute permutational index and 1-e contribution */
	pindx = pindex_single_rep(eostr1, ifo[0], ifo[2], ne1);
	i1 = index1e(ifo[0], ifo[2]);
	val = pindx * moints1[i1 - 1];
	val = val + single_rep_2e_contribution(
		eostr1, ifo[0], ifo[2], pindx, eostr2, ne1, ne2, moints2);
	return val;
}

/*
 * eval1_ncas_c0cv0v1: evaluate single virtual replacement matrix elements
 *                     between non-cas-flagged determinants.
 */
double eval1_ncas_c0cv0v1(struct occstr ostr1i, struct occstr ostr1j,
			  int ne1, struct occstr ostr2i, int ne2,
			  double *moints1, double *moints2, int ninto)
{
	double val;
	int  pindx; /* permuational index */
	int ifo[2]; /* initial, final orbital */
	int i1, i2; /* integral indexes */
	int eostr1[ne1]; /* electron occupation string */
	int eostr2[ne2]; 
	int i;
	
	/* locate initial and final orbitals and construct 
	 * eostr1 and eostr2 */
	virtdiffs_single_rep(ostr1i.virtx, ostr1j.virtx, ifo);
	make_orbital_strings_virt(ostr1i, eostr1, ne1, ninto);
	make_orbital_strings_virt(ostr2i, eostr2, ne2, ninto);

	/* compute permuational index and 1-e contriubtion */
	pindx = pindex_single_rep(eostr1, ifo[0], ifo[1], ne1);
	i1 = index1e(ifo[0], ifo[1]);
	val = pindx * moints1[i1 - 1];
	
	val = val + single_rep_2e_contribution(
		eostr1, ifo[0], ifo[1], pindx, eostr2, ne1, ne2, moints2);

	return val;
}

/* 
 * eval1_ncas_c1cv0v0: evaluate single cas replacement matrix elements 
 *                     between non-cas-flagged determinants
 */
double eval1_ncas_c1cv0v0(struct occstr ostr1, long long int xi, 
			  long long int xf, struct occstr ostr2, int ne1, 
			  int ne2, double *moints1, double *moints2,
			  int ninto)
{
	double val; /* <i|H|j> */
	
	int pindx; /* permuational index */
	int io, fo; /* initial, final orbital */
	int i1; /* integral indexes */
	
	int eostr1[ne1]; /* electron occupation string */
	int eostr2[ne2]; 
	
	int i;
	
	/* locate the nonzero bits in xi and xf and form occ_str1, occ_str2 */
	nonzerobits(xi, ninto, &io);
	nonzerobits(xf, ninto, &fo);
	make_orbital_strings_virt(ostr1, eostr1, ne1, ninto);
	make_orbital_strings_virt(ostr2, eostr2, ne2, ninto);
	
	/* compute permuational index and 1-e contribution */
	pindx = pindex_single_rep(eostr1, io, fo, ne1);
	i1 = index1e(io, fo);
	val = pindx * moints1[i1 - 1];
	
	val = val + single_rep_2e_contribution(
		eostr1, io, fo, pindx, eostr2, ne1, ne2, moints2);
	
	return val;
}

/*
 * eval2_ncas_c0cv0v2: evaluate double virtual replacements for non-cas-flagged
 *                     determinants
 */
double eval2_ncas_c0cv0v2(int *vxi, int *vxj, double *moints2)
{
	double val;
	int i1, i2; /* integral indexes */
	
	i1 = index2e(vxi[0],vxj[0],vxi[1],vxj[1]);
	i2 = index2e(vxi[0],vxj[1],vxi[1],vxj[0]);
	val = (moints2[i1 - 1] - moints2[i2 - 1]);
	return val;
}

/* 
 * eval2_ncas_c0cv2v0: evaluate cas->virtual replacements for non-cas-flagged 
 *                     determinants
 */
double eval2_ncas_c0cv2v0(long long int xi, long long int xf, int *vxi, 
			  int *vxf, double *moints2, int ninto,
                          long long int str)
{
	double val;
	int i1,i2;         /* integral indexes */
        int pindx = 1;     /* permutational index */
	int ifo[4] = {0};  /* initial, final orbitals */
        long long int t;   /* pseudo-excitation byte */
        /* use pseudo-excitation byte to find number of bytes between
         * cas orbital in excitation and virtuals. */
#ifndef BIGCAS
	cas_to_virt_replacements(2, 0, 0, xi, xf, vxi, vxf, ifo, ninto);  
#endif
        for (i1 = 1; i1 >= 0; i1--) {
                t = ((long long int) 1) << (ifo[i1] - 1);
                pindx = pindx * pindex_single_rep_cas2virt(str, t, ninto);
        }
        pindx = pindx * (-1);
	i1 = index2e(ifo[0], ifo[2], ifo[1], ifo[3]);
	i2 = index2e(ifo[0], ifo[3], ifo[1], ifo[2]);
	val = pindx * (moints2[i1 - 1] - moints2[i2 - 1]);
	return val;
}

/* 
 * eval2_ncas_c1cv0v1: evaluate cas + virtual replacements for non-cas-flagged
 * determinants.
 */
double eval2_ncas_c1cv0v1(struct occstr ostr, long long int xi, long long int xf, 
			  int *vxi, int *vxj, int ne, int nvx, double *moints2,
			  int ninto)
{
	double val;
	int i1,i2, pindx; /* integral index, permuational index */
	int ifov[2], io, fo; /* initial, final orbitals */
	int estr[ne]; /* electron orbital index string */
	
	virtdiffs_single_rep(vxi, vxj, ifov);
#ifndef BIGCAS
	nonzerobits(xi, ninto, &io);
	nonzerobits(xf, ninto, &fo);
#endif
	make_orbital_strings_virt(ostr, estr, ne, ninto);
	pindx = pindex_single_rep(estr, io, fo, ne);
	pindx = pindx * pindex_single_rep(vxi, ifov[0], ifov[1], nvx);
	i1 = index2e(ifov[0], ifov[1], io, fo);
	i2 = index2e(ifov[0], fo, io, ifov[1]);
	val = pindx * (moints2[i1 - 1] - moints2[i2 - 1]);
	return val;
}
/*
 * eval2_ncas_c00cv00v11: evaluate single virtual replacements in both alpha
 *                        and beta strings of non-cas-flagged determinants
 */
double eval2_ncas_c00cv00v11(int *avirti, int *avirtj, int *bvirti, int *bvirtj,
			     int nvrtxa, int nvrtxb, double *moints2)
{
	double val;
	int i1, pindx; /* integral index, permutational index */
	int ifoa[2], ifob[2]; /* initial, final orbitals */
	
	virtdiffs_single_rep(avirti, avirtj, ifoa);
	virtdiffs_single_rep(bvirti, bvirtj, ifob);
	
	pindx = pindex_single_rep(avirti, ifoa[0], ifoa[1], nvrtxa);
	pindx = pindx * pindex_single_rep(bvirti, ifob[0], ifob[1], nvrtxb);
	i1 = index2e(ifoa[0],ifoa[1],ifob[0],ifob[1]);
	val = pindx * moints2[i1 - 1];
	return val;
}

/* 
 * eval2_ncas_c00cv10v01: evaluate double replacements for non-cas-flagged
 * determinants with one cas->virt replacement in one string and one
 * virt->virt replacement in the other string.
 */
double eval2_ncas_c00cv10v01(long long int xi, long long int xf,
			     int *vx1i, int *vx1j,
			     int *vx2i, int *vx2j, double *moints2,
			     int ninto, long long int stri1)
{
	double val = 0.0;
	int i1, pindx; /* integral index, permuational index */
	int ifo1[4] = { 0 }, ifo2[2];    /* initial, final orbitals */
	
	cas_to_virt_replacements(1,0,0, xi, xf, vx1i, vx1j, ifo1, ninto);
	virtdiffs_single_rep(vx2i, vx2j, ifo2);
        if (xi != 0x00) {
                pindx = pindex_single_rep_cas2virt(stri1, xi, ninto);
                pindx = pindx * pindex_single_rep_virt(ifo1[2], vx1j);
        } else {
                pindx = pindex_single_rep_cas2virt(stri1, xf, ninto);
                pindx = pindx * pindex_single_rep_virt(ifo1[2], vx1i);
        }
	i1 = index2e(ifo1[0], ifo1[2], ifo2[0], ifo2[1]);
	val = pindx * moints2[i1 - 1];
	return val;	
}

/* 
 * eval2_ncas_c00cv10v10: evaluate double replacements for non-cas-flagged 
 * determinants with one cas->virt replacement and one virt->virt replacement
 * in the same string.
 */
double eval2_ncas_c00cv10v10(long long int xi, long long int xf,
			     struct occstr stri, struct occstr strj,
			     double *moints2, int ninto)
{
	double val;
	int i1,i2, pindx; /* integral indexes, permuational index */
	int ifo[4] = {0};    /* initial, final orbital array       */
	
	cas_to_virt_replacements(1,0,1,xi,xf,stri.virtx,strj.virtx, ifo, ninto);
        if (xi != 0x00) {
#ifndef BIGCAS
                pindx = pindex_single_rep_cas2virt(stri.byte1, xi, ninto);
#endif
                pindx = pindx * pindex_single_rep_virt(ifo[2], strj.virtx);
        } else {
#ifndef BIGCAS
                pindx = pindex_single_rep_cas2virt(strj.byte1, xf, ninto);
#endif
                pindx = pindx * pindex_single_rep_virt(ifo[2], stri.virtx);
        }

	i1 = index2e(ifo[0], ifo[2], ifo[1], ifo[3]);
	i2 = index2e(ifo[0], ifo[3], ifo[1], ifo[2]);
	val = pindx * (moints2[i1 - 1] - moints2[i2 - 1]);
	return val;
}

/* 
 * eval2_ncas_c00cv11v00: evaluate double replacements for non-cas-flagged
 * determinants with one cas->virt replacement in each string. 
 */
double eval2_ncas_c00cv11v00(long long int xi1, long long int xf1,
			     int *vxi1, int *vxj1, long long int xi2,
			     long long int xf2, int *vxi2, int *vxj2,
			     double *moints2, int ninto,
                             long long int stri1, long long int stri2)
{
	double val;
	int i1, pindx; /* integral index, permuational index */
	int ifo1[4] = {0}, ifo2[4]={0};   /* initial, final orbital array */
	cas_to_virt_replacements(1,0,0, xi1, xf1, vxi1, vxj1, ifo1, ninto);
	cas_to_virt_replacements(1,0,0, xi2, xf2, vxi2, vxj2, ifo2, ninto);
	if (xi1 != 0x00) {
                //stri1 = stri1 | xi1;
                pindx = pindex_single_rep_cas2virt(stri1, xi1, ninto);
        } else {
                //stri1 = stri1 | xf1;
                pindx = pindex_single_rep_cas2virt(stri1, xf1, ninto);
        }
        if (xi2 != 0x00) {
                //stri2 = stri2 | xi2;
                pindx = pindx * pindex_single_rep_cas2virt(stri2, xi2, ninto);
        } else {
                //stri2 = stri2 | xf2;
                pindx = pindx * pindex_single_rep_cas2virt(stri2, xf2, ninto);
        }
        i1 = index2e(ifo1[0], ifo1[2], ifo2[0], ifo2[2]);
	val = pindx * moints2[i1 - 1];
	return val;
}
/* 
 * eval2_ncas_c01cv10v00: evaluate double replcaements for non-cas-flagged
 * determinants with one cas->virt replacement in one string and one
 * cas->cas replacement in the other string.
 */
double eval2_ncas_c01cv10v00(struct occstr str1, struct occstr str2,
			     long long int xi1, long long int xf1,
			     int *vx1i, int *vx1j, long long int xi2,
			     long long int xf2, int ne1, int ne2,
			     double *moints2, int ninto)
{
	double val;
	int i1, pindx; /* integral index, permuational index */
	int ifo1[4] = {0}, ifo2[2]; /* initial, final orbitals */

	cas_to_virt_replacements(1,0,0, xi1, xf1, vx1i, vx1j, ifo1, ninto);
	nonzerobits(xi2, ninto, &(ifo2[0]));
	nonzerobits(xf2, ninto, &(ifo2[1]));
#ifndef BIGCAS
	pindx = pindex_single_rep_cas(str2.byte1, xi2, xf2, ninto);
        if (xi1 != 0x00) {
                pindx = pindx * pindex_single_rep_cas2virt(
                                str1.byte1, xi1, ninto);
        } else {
                pindx = pindx * pindex_single_rep_cas2virt(
                                str1.byte1, xf1, ninto);
        }
#endif
        i1 = index2e(ifo1[0], ifo1[2], ifo2[0], ifo2[1]);
	val = pindx * moints2[i1 - 1];
	return val;
}

/*
 * eval2_ncas_c10cv10v00: evaluate double replacements for non-cas-flagged
 * determinants with one cas->virt replacement and one cas->cas replacement
 * int the same string.
 */
double eval2_ncas_c10cv10v00(struct occstr str, long long int xi, 
			     long long int xf, int *vxi,
			     int *vxj, int ne, double *moints2,
			     int ninto)
{
    double val;
    int i1, i2, pindx; /* integral index, permuational index */
    int ifo[4] = {0};    /* initial, final orbital array */
    long long int b1, b2, b3; /* cas orbital bytes */
    int estr[ne];  /* electron orbital index string */
    
    cas_to_virt_replacements(1, 1, 0, xi, xf, vxi, vxj, ifo, ninto);
    b1 = ((long long int) 1) << (ifo[1] - 1);
    b2 = ((long long int) 1) << (ifo[0] - 1);
    b3 = ((long long int) 1) << (ifo[2] - 1);
    make_orbital_strings_virt(str, estr, ne, ninto);
    pindx = pindex_single_rep_cas(str.byte1, b1, b3, ninto);
    i1 = index2e(ifo[0], ifo[2], ifo[1], ifo[3]);
    i2 = index2e(ifo[0], ifo[3], ifo[1], ifo[2]);
    val = pindx * (moints2[i2 - 1] - moints2[i1 - 1]);
    return val;
}

/* 
 * make_orbital_strings_virt: make orbital strings with virtual occupations.
 */
void make_orbital_strings_virt(struct occstr ostr1i, int *eostr1, int nelec1,
	int ninto)
{
	int i;

	init_int_array_0(eostr1, nelec1);

#ifndef BIGCAS
	nonzerobits(ostr1i.byte1, ninto, eostr1);
#endif
	if (eostr1[nelec1 - 2] == 0) {
		eostr1[nelec1 - 2] = ostr1i.virtx[0];
		eostr1[nelec1 - 1] = ostr1i.virtx[1];
	} else if (eostr1[nelec1 - 1] == 0) {
		eostr1[nelec1 - 1] = ostr1i.virtx[0];
	}
	return;
}

/*
 * pindex_double_rep_cas: compute permuational index for 2 replacements in
 * the cas orbitals.
 */
int pindex_double_rep_cas(long long int str, int *io, int *fo, int ninto)
{
        int pindx = 1; /* permutational index */
        long long int xi, xf; /* initial, final orbitals */
        int i;
        for (i = 0; i < 2; i++ ) {
                xi = ((long long int ) 1) << (io[0] - 1);
                xf = ((long long int ) 1) << (fo[0] - 1);
                pindx = pindx * pindex_single_rep_cas(str, xi, xf, ninto);
        }
        return pindx;
}

/*
 * pindex_single_rep_cas2virt: compute permutational index for single 
 * excitaiton from cas byte -> vitual orbitals.
 */
int pindex_single_rep_cas2virt(long long int stri, long long int xi,
        int ninto)
{
        int pindx; /* permuational index */
        int i;
        long long int xf; /* (psuedo final cas orbital */

        /* Turn on the (ninto + 1) bit, and treat this like a single cas
         * excitation.
         */
        xf= ((long long int) 1) << ninto;
        pindx = pindex_single_rep_cas(stri, xi, xf, ninto);
        return pindx;
}
        
/*
 * pindex_single_rep_cas: compute permuational index for single excitation
 * within CAS. This is done by counting occupations between the orbital in
 * xi and the orbital in xf.
 */
int pindex_single_rep_cas(long long int stri, long long int xi,
			  long long int xf, int ninto)
{
    int pindx = 1; /* permutational index */
    int i;
    long long int t; /* 2^64 */
    long long int x; /* xi ^ xf */
    /* Right shift $x and $stri until first non-zero bit of $x is 
     * right aligned. Now, left shift $x and $stri until first nonzero 
     * value is left aligned.
     */
    t = ((long long int) 1) << 63;
    x = xi ^ xf;
    for (i = 0; i < 64; i++) {
	if (x & 0x01) {
	    stri = stri >> 1;
	    x = x >> 1;
	    break;
	} else {
	    stri = stri >> 1;
	    x = x >> 1;
	}
    }
    for (i = 63; i >= 0; i--) {
	if (x & t) {
	    stri = stri << 1;
	    x = x << 1;
	    break;
	} else {
	    stri = stri << 1;
	    x = x << 1;
	}
    }
    /* count nonzero bits, there will be less than $ninto nonzero bits */
    for (i = 0; i <= ninto; i++) {
	if (stri & t) {
	    //printf(" \n PERMUTATION!\n");
	    pindx = pindx * (-1);
	    stri = stri << 1;
	} else {
	    stri = stri << 1;
	}
    }
    return pindx;
}

/* 
 * pindex_single_rep: compute permuational index for single excitation
 */
int pindex_single_rep(int *str, int io, int fo, int lstr)
{
	int pindx = 1;
	int pos, i;
	
	pos = find_pos_in_array_lnsrch(io, str, lstr);
	if (fo > str[pos + 1]) {
		i = pos + 1;
		while (fo > str[i] && i < lstr) {
			pindx = pindx * (-1);
			i++;
		}
	} else if (fo < str[pos - 1]) {
		i = pos - 1;
		while (fo < str[i] && i > 0) {
			pindx = pindx * (-1);
			i--;
		}
	}
	return pindx;
}
	
/*
 * pindex_single_rep_virt: excitation orbital, virtual orbitals.
 */
int pindex_single_rep_virt(int xorb, int *vorbs)
{
        int pindx = 1;
        /* Check the second virtual occupation. If this occupation is
         * zero, then we drop out; else, we check where xorb is.
         */
        if (vorbs[1] == 0) {
                return pindx;
        } else {
                if (xorb == vorbs[0]) {
                        return pindx;
                } else {
                        pindx = (-1) * pindx;
                        return pindx;
                }
        }
}
/*
 * single_rep_2e_contribution: compute contribution of 2e integrals to
 * single replacements.
 */
double single_rep_2e_contribution(int *eostr1, int io, int fo,
				  int pindx, int *eostr2, int ne1,
				  int ne2, double *moints2)
{
	double value = 0.0;
	int          i1,i2; /* integral indexes */
	int              i;
	
	for (i = 0; i < ne1; i++) {
		if (eostr1[i] != io || eostr1[i] != fo) {
			i1 = index2e(eostr1[i], eostr1[i], 
				     io, fo);
			i2 = index2e(eostr1[i], io, eostr1[i],
				     fo);
			value = value + pindx * (moints2[i1 - 1] -
						      moints2[i2 - 1]);
		}
	}
	for (i = 0; i < ne2; i++) {
		i1 = index2e(eostr2[i], eostr2[i], io, fo);
		value = value + pindx * moints2[i1 - 1];
	}
	return value;
}

/* 
 *virtdiffs_single_cas_to_virt: find location of cas->virt replacement
 */
void virtdiffs_single_cas_to_virt(int *vxi, int *vxj, int *repo)
{
	/* possible differences:
	 *  0 0  0 a  0 x  
	 *  0 x  a x  a a  */
	if (vxi[1] + vxj[1] == 0) {
		*repo = vxj[0];
	} else {
		if (vxi[0] == vxj[0] ) {
			*repo = vxj[1];
		} else {
			*repo = vxj[0];
		}
	}
	return;
}

/* 
 * virtdiffs_single_rep: find location of virtual orbital replacement
 *
 * Input:
 *  vxi = deti virtual orbitals
 *  vxj = detj virtual orbitals
 * Output:
 *  ifo = differences 
 */
void virtdiffs_single_rep(int *vxi, int *vxj, int *ifo)
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
