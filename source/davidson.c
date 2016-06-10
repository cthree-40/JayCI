// File: davidson.c

/*
 * Subfunctions for preforming the davidson algorithm to compute the lowest
 * roots of the CI Hamiltonian.
 *
 * Subfunctions:
 *  diag_subspacehmat: diagonalize subspace matrix, returning eigenvectors
 *  dvdalg: main davidson algorithm driver.
 *  make_subspacehmat: build subspace hamiltonian matrix, v.Hv
 *  perform_hv_initspace: perform Hv=c for all v in current space.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "errorlib.h"
#include "ioutil.h"
#include "arrayutil.h"
#include "allocate_mem.h"
#include "binarystr.h"
#include "cimapping.h"
#include "action_util.h"
#include "mathutil.h"
#include "initguess_sbd.h"
#include "davidson.h"

/*
 * append_vector_to_space: append new orthonormal basis vector to space.
 */
void append_vector_to_space(double **bvecs, int spsize, int vlen,
			    double *newvec)
 {
	 int i;
	 for (i = 0; i < vlen; i++) {
		 bvecs[spsize][i] = newvec[i];
	 }
	 return;
 }

/*
 * diag_subspacehmat: diagonalize subspace matrix, returning eigenvectors
 */
int diag_subspacehmat(double **hmat, double **hevec, double *heval, int dim,
		      double *hmat_scr, double *hevec_scr)
{
	int error = 0; /* error flag */
	/* copy 2d arrays into corresponding 1d arrays */
	cparray_2d1d(hmat, dim, dim, hmat_scr);
	/* diagonalize matrix */
	error = diagmat_dsyevr(hmat_scr, dim, hevec_scr, heval);
	if (error != 0) {
		error_flag(error, "diag_subspacehmat");
		return error;
	}
	/* copy 1d arrays into corresponding 2d arrays */
	cparray_1d2d(hevec_scr, hevec, dim, dim);
	/* zero out scratch arrays */
	init_dbl_array_0(hevec_scr, (dim * dim));
	init_dbl_array_0(hmat_scr, (dim * dim));
	return error;
}

/*
 * dvdalg: perform davidson algorithm on CI space to find lowest nroots
 * roots.
 *
 * dlist = determinant list
 * ndets = number of determinants
 * moints1 = 1-e integrals
 * moints2 = 2-e integrals
 * aelec = alpha electrons
 * belec = beta electrons
 * hdgls = <i|H|i>
 * ninto = (ndocc + nactv) CAS DOCC + CAS ACTIVE orbitals
 * totfrze = nucrep + frozen core
 * maxiter = maximum davidson iterations to exectute
 * krymin = minimum size of krylov space
 * krymax = maximum size of krylov space
 * nroots = number of roots to compute
 * restol = convergence tolerance of residual vector
 * hmap = map of nonzero matrix elements <i|H|j>
 * civec = civectors; civec[nroots][ndets] double array
 * cival = eigenvalues; cival[nroots] double array
 * predr = prediagonalization routine
 * plvl = print level
 */
int dvdalg(struct det *dlist, int ndets, double *moints1, double *moints2,
	   int aelec, int belec, double *hdgls, int ninto, double totfrze,
	   int maxiter, int krymin, int krymax, int nroots, double restol,
	   struct rowmap *hmap, double **civec, double *cival, int predr,
	   int plvl)
{
	int error = 0; /* error flag */
	double **hscr = NULL; /* krylov space hamiltonian v.Hv array */
	double *hscr_1d = NULL; /* one-dimensional form of v.Hv for
				 * interface with F90 */
	double **cscr = NULL; /* Hv=c vectors scratch array */
	double **vscr = NULL; /* v vectors scratch array */
	double *vscr_1d = NULL; /* one-dimensional form of basis vectors */
	double *tscr_1d = NULL; /* two-dimensional scratch for truncated vasis
				 * vectors */
	double **hevec = NULL; /* krylov eigenvectors array */
	double *hevec_1d = NULL; /* one-dimensional form of eigenvectors */
	double *heval = NULL; /* krylov eigenvalues array */
	double *rvec = NULL; /* residual vector: Hv - (v.Hv)c = r */
	double *nvec = NULL; /* new vector */
	double rnorm = 0.0; /* norm of residual vector */
	int croot = 0; /* current root solving for {1..nroot} */
	int ckdim = 0; /* current dimension of krylov space */
	int citer = 0; /* current iteration */
	
	/* allocate scratch arrays. size is constant throughout routine. */
	error = allocate_mem_double(&hscr, krymax, krymax);
	error = error + allocate_mem_double(&cscr, ndets, krymax);
	error = error + allocate_mem_double(&vscr, ndets, krymax);
	error = error + allocate_mem_double(&hevec, krymax, krymax);
	if (error != 0) {
		error_flag(error, "dvdalg");
		return error;
	}
	init_dbl_2darray_0(hscr, krymax, krymax);
	init_dbl_2darray_0(cscr, ndets, krymax);
	init_dbl_2darray_0(vscr, ndets, krymax);
	init_dbl_2darray_0(hevec, krymax, krymax);
	
	heval = (double *) malloc(krymax * sizeof(double));
	init_dbl_array_0(heval, krymax);
	
	rvec = (double *) malloc(ndets * sizeof(double));
	init_dbl_array_0(rvec, ndets);
	
	nvec = (double *) malloc(ndets * sizeof(double));
	init_dbl_array_0(nvec, ndets);
	
	hscr_1d = (double *) malloc(krymax * krymax * sizeof(double));
	init_dbl_array_0(hscr_1d, (krymax * krymax));
	
	hevec_1d = (double *) malloc(krymax * krymax * sizeof(double));
	init_dbl_array_0(hevec_1d, (krymax * krymax));
	
	vscr_1d = (double *) malloc(krymax * ndets * sizeof(double));
	init_dbl_array_0(vscr_1d, (krymax * krymax));
	
	tscr_1d = (double *) malloc(krymax * ndets * sizeof(double));
	init_dbl_array_0(tscr_1d, (krymax * krymax));
	
	/* Generate initial guess vectors by diagonalizing an upper block
	 * of the Hamiltonian. */
	error = initguess_sbd(dlist, ndets, moints1, moints2, aelec,
			      belec, ninto, krymin, vscr);
	if (error != 0) {
		error_flag(error, "dvdalg");
		return error;
	}
	
	/* Main Loop */
	citer = 1;
	croot = 1;
	ckdim = krymin;
	while (citer < maxiter && croot <= nroots) {
		/* Perform Hv=c on all basis vectors v. Then build the subspace
		 * Hamiltonian v.Hv = v.c */
		perform_hv_initspace(dlist, ndets, moints1, moints2, aelec,
				     belec, ninto, hmap, ckdim, vscr, cscr);
		make_subspacehmat(vscr, cscr, ndets, ckdim, hscr);
		if (plvl > 0) {
			fprintf(stdout, "\n Subspace hamiltonian:\n");
			print_array_2d(hscr, ckdim, ckdim);
		}
		/* Diagonalize v.Hv matrix, returning eigenvalues and
		 * eigenvectors */
		error = diag_subspacehmat(hscr, hevec, heval, ckdim, hscr_1d,
					  hevec_1d);
		if (error != 0) {
			error_flag(error, "dvdalg");
			return error;
		}
		
		/* Sub-loop */
		while (ckdim < krymax) {
			/* Generate the residual vector */
			generate_residual(vscr, cscr, ndets, ckdim, hevec,
					  heval, croot, rvec);
			rnorm = compute_vector_norm(rvec, ndets);
			/* Print out iteration information */
			print_iteration_info(heval, ckdim, croot,
					     rnorm, totfrze);
			/* test convergence */
			if (rnorm < restol) {
				fprintf(stdout, "\n Root %d converged!\n",
					croot);
				fprintf(stdout, " Eigenvalue = %15.8lf\n",
					(heval[croot - 1] + totfrze));
				croot++;
				ckdim+=100;
				if (croot == nroots) {
					fprintf(stdout,
						"\n ** CI CONVERGED! ** \n");
					return error;
				}
				continue; /* leave subloop */
			}
			
			/* generate new vector to add to space */
			generate_newvector(rvec, hdgls, heval[croot - 1],
					   ndets, nvec);
			orthonormalize_vector(vscr, ckdim, ndets, nvec);
			append_vector_to_space(vscr, ckdim, ndets, nvec);
			ckdim++;
			
			/* compute Hv=c for new vector */
			compute_hv(dlist, ndets, moints1, moints2, aelec,
				   belec, vscr[ckdim - 1], cscr[ckdim - 1],
				   ninto, hmap);
			
			/* make subspace hamiltonian */
			make_subspacehmat(vscr, cscr, ndets, ckdim, hscr);
			if (plvl > 0) {
				fprintf(stdout, "\n Subspace Hamiltonian:\n");
				print_array_2d(hscr, ckdim, ckdim);
			}
			
			/* diagonalize subspace hamiltonian */
			error = diag_subspacehmat(hscr, hevec, heval, ckdim,
						  hscr_1d, hevec_1d);
			if (error != 0) {
				error_flag(error, "dvdalg");
				return error;
			}
			citer++;
		}
		
		/* truncate space */
		error = truncate_kspace(vscr, ndets, krymin, krymax, croot,
					hevec, vscr_1d, hevec_1d, tscr_1d);
		if (error != 0) {
			error_flag(error, "dvdalg");
			return error;
		}
		ckdim = krymin;
		init_dbl_2darray_0(cscr, ndets, krymax);
	}
	fprintf(stdout, " Davidson algorithm finished.\n");
	
	return error;
}

/*
 * generate_newvector: build correction vector.
 */
void generate_newvector(double *rvec, double *dgls, double eval, int ndets,
			double *nvec)
{
	int i;
	for (i = 0; i < ndets; i++) {
		nvec[i] = (-1.0) * rvec[i] / (dgls[i] - eval);
	}
	return;
}

/*
 * generate_residual: generate residual vector for davidson algorithm.
 */
void generate_residual(double **vvecs, double **cvecs, int ndets, int nvec,
		       double **hevec, double *heval, int croot, double *rvec)
{
	int i, j;
	int root_id;
	root_id = croot - 1;
	/* zero out vector */
	init_dbl_array_0(rvec, ndets);
	for (i = 0; i < ndets; i++) {
		for (j = 0; j < nvec; j++) {
			rvec[i] += hevec[root_id][j] *
				(cvecs[j][i] - heval[root_id] *	vvecs[j][i]);
		}
	}
	return;
}

/*
 * make_subspacehmat: build krylov space hamiltonian v.Hv.
 */
void make_subspacehmat(double **v, double **hv, int ndets, int ckdim,
		       double **hmat)
{
	int i, j;
	/* Build v.Hv by taking dot product of v_i and (Hv)_j */
	for (i = 0; i < ckdim; i++) {
		for (j = 0; j < ckdim; j++) {
			hmat[i][j] = dot_product(v[i], hv[j], ndets);
		}
	}
	return;
}

/*
 * orthonormalize_vector: orthogonalize and normalize new vector to space
 * of basis vectors.
 */
void orthonormalize_vector(double **bvecs, int nbv, int vlen, double *nvec)
{
	
	double nvnorm; /* norm of orthogonalized vector */
	double ovrlp; /* overlap of new vector with basis vectors */
	int i, j;
	/* Loop over space removing the overlap of the new vector with each basis
	 * vector. */
	for (i = 0; i < nbv; i++) {
		ovrlp = dot_product(bvecs[i], nvec, vlen);
		for (j = 0; j < vlen; j++) {
			nvec[j] = nvec[j] - ovrlp * bvecs[i][j];
		}
	}
	/* Normalize this vector */
	nvnorm = compute_vector_norm(nvec, vlen);
	for (i = 0; i < vlen; i++) {
		nvec[i] = nvec[i] / nvnorm;
	}
	/* check orthogonalization */
	for (i = 0; i < nbv; i++) {
		ovrlp = dot_product(bvecs[i], nvec, vlen);
		if (ovrlp > 0.000001) {
			printf("*** Warning: ovrlp = %15.8lf ***\n", ovrlp);
		}
	}
	/* check norm */
	nvnorm = compute_vector_norm(nvec, vlen);
	if (pow(nvnorm - 1.0, 2) > 0.0000001) {
		printf("*** Warning: nvnorm = %15.8lf ***\n", nvnorm);
	}
	return;
}

/*
 * perform_hv_initspace: perform hv on inital vectors.
 */
void perform_hv_initspace(struct det *dlist, int ndets, double *moints1,
			  double *moints2, int aelec, int belec, int ninto,
			  struct rowmap *hmap, int nvec, double **vecs,
			  double **hvecs)
{
	int i, j, k;
	
#ifdef DEBUGGING
	double *v1, *c1;
	FILE *fptr;
	char  fname[30];
	
	v1 = (double *) malloc(ndets * sizeof(double));
	c1 = (double *) malloc(ndets * sizeof(double));
	for (j = 0; j < 1; j++) {
		for (i = 0; i < ndets; i++) {
			v1[i] = 0.0;
			c1[i] = 0.0;
		}
		v1[j] = 1.0;
		sprintf(fname, "file.%d", j);
		fptr = fopen(fname, "w");
		compute_hv(dlist, ndets, moints1, moints2, aelec, belec, v1, c1,
			   ninto, hmap);
		printf("Writing %d elements to file %d.\n", ndets, j);
		for (i = 0; i < ndets; i++) {
			k = fprintf(fptr, "%15.8lf\n", c1[i]);
		}
		fclose(fptr);
	}
	exit(1);
	
#endif
	/* Loop over vectors */
	for (i = 0; i < nvec; i++) {
		compute_hv(dlist, ndets, moints1, moints2, aelec, belec, vecs[i],
			   hvecs[i], ninto, hmap);
	}
	return;
}

/*
 * print_iteration_info: print davidson iteration information.
 */
void print_iteration_info(double *heval, int ckdim, int croot, double rnorm,
			  double totfrze)
{
	int i;
	fprintf(stdout, "\n");
	for (i = 0; i < 70; i++) {
		fprintf(stdout, "-");
	}
	fprintf(stdout, "\n  Eigenvalues: ");
	for (i = 0; i < ckdim; i++) {
		fprintf(stdout, " %15.8lf", (heval[i] + totfrze));
	}
	fprintf(stdout, "\n  ||r||: %15.8lf\n", rnorm);
	return;
}

/*
 * truncate_kspace: truncate the krylov space from krymax to krymin.
 */
int truncate_kspace(double **vscr, int vlen, int krymin, int krymax,
		    int croot, double **hevec, double *vscr1d, double *hevec1d,
		    double *scr)
{
	int error = 0; /* error flag */
	int i, j, k;
	/* copy arrays into proper 1d scratch arrays */
	cparray_2d1d(vscr, vlen, krymax, vscr1d);
	cparray_2d1d(hevec, krymax, krymax, hevec1d);
	/* perform matrix multiplication */
	error = matmul_dgemm(vscr1d, vlen, krymax, hevec1d, krymax, krymax, scr,
			     vlen, krymax);
	init_dbl_2darray_0(vscr, vlen, krymax);
	/* copy 1d array into 2d array */
	cparray_1d2d(scr, vscr, vlen, krymin);
	init_dbl_array_0(scr, (vlen * krymax));
	return error;
}
