/* davidson.c */
/*
 * Davidson algorithm and associated subfunctions.
 *
 * Written by:
 * Christopher L Malbon
 * Yarkony Group, Dept. of Chemistry, The Johns Hopkins University
 * (C) 2015 -
 *
 * Subfunctions:
 *  diagonalize_subspace_ham: diagaonalize the v.Hv subspace
 *  davidson_diagonalization_routine: davidson algorithm
 *  perform_hv_initialspace: perform Hv=c on initial space
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

/* 
 * davidson_diagaonalization_routine: davidson algorithm.
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
	double *hvscr; /* hv vectors scratch array */
	double *bscr;  /* v  vectors scratch array */
	double *hevec, *heval; /* vhv eigenvetors and eigenvalues */
	double *rvec, *nvec;  /* residual vector, new vector */
	double rnorm;  /* norm of residual */
	int niter = 1, cdim, croot; /* iteration count, current dimension, 
				     * current root*/
	int i;
	
	/* allocate arrays. Then initialize them to 0d0. We then
	 * place the initial vectors in the basis vector scratch
	 * array and begin our loop */
	
	hscr = (double *) malloc(krymax * krymax * sizeof(double));
	hevec= (double *) malloc(krymax * krymax * sizeof(double));
	heval= (double *) malloc(krymax * sizeof(double));
	hvscr= (double *) malloc(krymax * ndets * sizeof(double));
	bscr = (double *) malloc(krymax * ndets * sizeof(double));
	rvec = (double *) malloc(ndets * sizeof(double));
	nvec = (double *) malloc(ndets * sizeof(double));
	init_dbl_array_0(hscr, (krymax * krymax));
	init_dbl_array_0(hvscr,(krymax * ndets));
	init_dbl_array_0(bscr, (krymax * ndets));
	init_dbl_array_0(hevec,(krymax * krymax));
	init_dbl_array_0(heval, krymax);
	init_dbl_array_0(rvec, ndets);
	init_dbl_array_0(nvec, ndets);
	for (i = 0; i < krymin * ndets; i++) {
		bscr[i] = initv[i];
	}
	
	cdim = krymin;
	croot = 1;

	while (niter <= maxiter) {
		
		perform_hv_inititalspace(dlist, ndets, moints1, moints2,
					 aelec, belec, krymin, initv, hvscr);
		make_subspace_hamiltonian(bscr, hvscr, ndets, cdim, hscr);
		diagonalize_subspace_ham(hscr, cdim, krymax, hevec, heval);
		
		while (cdim <= krymax) {
			
			generate_residual_vec(bscr, hvscr, ndets, cdim,
					      hevec, heval, croot, rvec);
			rnorm = compute_vector_norm(rvec, ndets);
			if (rnorm < res_tol) {
				printf(" Root %d converged! ");
				printf(" Eigenvalue = %15.8f\n");
				croot++;
				cdim+=100;
				if (croot >= nroot) nroot+=(maxiter*2);
				continue;
			}
			generate_new_vector(rnorm, dgls, heval(croot), ndets,
					   nvec);
			orthogonalize_vector(bscr, cdim, ndets, nvec);
			append_new_vector(bscr, cdim, ndets, nvec);
			cdim++;
			compute_hv(dlist, ndets, moints1, moints2, aelec,
				   belec, bscr[(cdim - 1) * ndets],
				   hvscr[(cdim - 1) * ndets]);
			init_dbl_array_0(hscr, (krymax * krymax));
			make_subspace_hamiltonian(bscr, hvscr, ndets, cdim,
						  hscr);
			diagonalize_subspace_ham(hscr, cdim, krymax, hevec,
						 heval);
		}
	}
	return;
}

/*
 * make_subspace_hamiltonian: build v.Hv subspace hamiltonian.
 */
void make_subspace_hamiltonian(double *v, double *hv, int ndt, int nv,
        double *vhv)
{
    int i, j;
    double *p;
    p = vhv;
    /* Compute v(i).Hv(j). vHv is stored in column-major order. */
    for (i = 0; i < nv; i++) {
        for (j = 0; j < nv; j++) {
            *p = dot_product(v[j * ndt], hv[i * ndt], ndt);
            p++;
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

