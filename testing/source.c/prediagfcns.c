// File: prediagfcns.c
/*********************************************************************
 * prediagrcns: prediagonalization functions
 * -------------------------------------------------------------------
 * drefblock : diagonalize a reference space of the hamiltonian
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "binary.h"
#include "binarystr.h"
#include "mathutil.h"
#include "arrayutil.h"
#include "action_util.h"

/* drefblock: diagonalizes an explicitly constructed reference space 
 *            of the hamiltonian.
 * -------------------------------------------------------------------
 * Reference block is: H(1:refdim,1:refdim).
 *
 * Input:
 *  detlist = determinant list
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 *  m1len   = number of 1-e integrals
 *  m2len   = number of 2-e integrals
 *  aelec   = alpha electrons
 *  belec   = beta  electrons
 *  refdim = size of reference block
 * Output:
 *  evecs = eigenvectors of diagonalized block */
int drefblock(struct det *detlist,  double *moints1, double *moints2,
	      int m1len, int m2len, int aelec, int belec, int refdim, 
	      double *evecs, double frzce)
{
	int i, j, k, l;
	double *hmat;
	double *evals;
	int err;

	err = 0;
	
	hmat = (double *) malloc(refdim * refdim * sizeof(double));
	evals = (double *) malloc(refdim * sizeof(double));
	init_dbl_array_0(hmat, (refdim * refdim));
	init_dbl_array_0(evals,(refdim));
	for (i = 0; i < 1; i++) {
		for (j = i; j < refdim; j++) {
			
			k = (i * refdim) + j;
			l = (j * refdim) + i;
			hmat[k] = hmatels(detlist[i], detlist[j], moints1, moints2,
					  aelec, belec);
			printf("%15.8f\n", hmat[k]);
			hmat[l] = hmat[k];
		}
	}
	exit(1);
	/* diagonalize reference block */
	err = diagmat_dsyevr(hmat, refdim, evecs, evals);
        
	fprintf(stdout,
		"Lowest eigenvalue in reference space: %15.8lf\n", 
		(evals[0] + frzce));
	
	free(hmat);
	free(evals);
	
	return err;
	
}
