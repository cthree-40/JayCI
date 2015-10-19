/* davidson.c */
/*
 * Davidson algorithm and associated subfunctions.
 *
 * Written by:
 * Christopher L Malbon
 * Yarkony Group, Dept. of Chemistry, The Johns Hopkins University
 * (C) 2015 -
 *
 * Subfunctions contains:
 *  davidson_diagonalization_routine: davidson algorithm
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "action_util.h"
#include "arrayutil.h"
#include "binarystr.h"
#include "davidson.h"

/*
 * davidson_diagonalization_routine: performs davidson algorithm
 *
 * The algorithm procedes as follows:
 * 1) Perform Hv=c on initial space.
 * 2) diagonalize v.Hv
 * 3) compute residual
 * 4) generate new vector w
 * 5) perform Hw=c
 * 6) go to (2) with new space including w
 */
void davidson_diagonalization_routine(struct det *dlist, int ndets, 
				      double *moints1, double *moints2, 
				      int m1len, int m2len, int maxiter,
				      int aelec, int belec, int krymin,
				      int krymax, int nroots, double restol,
				      double *initv, double *civec, 
				      double *cienergy)
{
	double *hscr;  /* vhv scratch array */
	double *hvscr; /* hv scratch array */
	int niter = 1;
	
	hscr = (double *) malloc(krymax * krymax * sizeof(double));
	hvscr= (double *) malloc(krymax * ndets * sizeof(double));
	
	perform_hv_inititalspace(dlist, ndets, moints1, moints2,
				 aelec, belec, krymin, initv, hvscr);
	
		
	return;
}

/* 
 * perform_hv: perform hv=c
 */
void perform_hv(struct det *dlist, int ndets, double *moints1, double *moints2,
		int aelec, int belec, double *v, double *c)
{
	int i, j, k, l;
	double valij;
	
	init_dbl_array_0(c, ndets);

	for (i = 0; i < ndets; i++) {
		for (j = i; j < ndets; j++) {
			valij = hmatels(dlist[i], dlist[j], moints1, moints2,
				       aelec, belec);
			c[i] = c[i] + valij * v[j];
			c[j] = c[j] + valij * v[i];
		}
	}
	return;
}

/* 
 * perform_hv_initialspace: perform hv on inital vector space for davidson
 * algorithm
 */
void perform_hv_initialspace(struct det *dlist, int ndets, double *moints1,
			     double *moints2,  int aelec,
			     int belec, int idim, double *initv, double *finlv)
{
	int i;
	
	for (i = 0; i < idim; i++) {
		compute_hv(dlist, ndets, moints1, moints2, aelec, belec,
			   initv[i*ndets], finlv[i*ndets]);
	}
	return;
}
	
