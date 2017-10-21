/* File: mathutil.c */
/*
 * mathutil: math functions, some interfacing with LAPACK
 * -------------------------------------------------------------------
 * compute_vector_norm: compute norm of a vector
 * diagmat_dsyevr: diagonalizes a square matrix using dsyevr
 * dot_product: computes dot product of two vectors.
 * matmul_dgemm: perform A_ij B_jk = C_ik with dgemm
 * orthonormalize_vector: orthonormalize vector to space
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 */       

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arrayutil.h"
#include "mathutil.h"

/* 
 * fortran subroutines
 */
extern double _dlamch_fcn_();

/*
 * compute_vector_norm: compute euclidean norm of vector.
 */
double compute_vector_norm(double *v, int d)
{
    double nrm, dp = 0.0;
    dp = dot_product(v, v, d);
    nrm = sqrt(dp);
    return nrm;
}

/* 
 * diagmat_dsyevr: diagonalizes a square matrix using the dsyevr
 * LAPACK subroutine
 */
int diagmat_dsyevr(double *mat, int dim, double *evecs, double *evals)
{
    unsigned char jobz[1]  = "v"; /* return eigenvectors and eigenvalues */
    unsigned char range[1] = "a"; /* compute all eigenvectors */
    unsigned char uplo[1]  = "u"; /* upper triangle stored */
    long long int n;              /* order of matrix */
    long long int lda;            /* leading dimension of a */
    double vl, vu;                /* lower and upper bounds of eigenvalue 
                                   *  interval */
    double il, iu;                /* indices of smallest and largest 
                                   *  eigenvalues to find */
    double abstol;                /* absolute error tolerance */
    long long int m;              /* total number of eigenvalues found */
    long long int ldz;            /* leading dimension of eigenvector array */
    long long int *isuppz;        /* support array */
    double *work;                 /* workspace array */
    long long int lwork;          /* dimension of workspace array */
    long long int *iwork;         /* integer workspace array */
    long long int liwork;         /* dimension of integer workspace array */
    long long int *info;          /* information */

    int err = 0;

    /* set size */
    n = dim;
    lda = dim;
    vl = 1;
    vu = dim;
    il = 1;
    iu = dim;
    ldz= dim;
    abstol = dlamch_fcn_();
    lwork = 30 * dim + 10;
    liwork = 15 * dim;
    info = 0;
    isuppz = (long long int *) malloc(2 * dim * sizeof(long long int));
    work = (double *) malloc(lwork * sizeof(double));
    iwork = (long long int *) malloc(liwork * sizeof(double));

    dsyevr_(&jobz, &range, &uplo, &n, mat, &lda, &vl, &vu, &il, &iu,
	    &abstol, &m, evals, evecs, &ldz, isuppz, work, &lwork,
	    iwork, &liwork, &info);

    free(isuppz);
    free(work);
    free(iwork);

    if (info != 0) {
	fprintf(stderr,"*** ERROR in DSYEVR_! info = %lld ***\n", *info);
	err = -1;
    }
    return err;
}
    
/*
 * dot_product: compute the dot product of two vectors
 */
double dot_product(double *u, double *v, int d)
{
    int i;
    double dp = 0.0;
    for (i = 0; i < d; i++) {
        dp += u[i] * v[i];
    }
    return dp;
}

/*
 * gauss_fcn: gaussian function. exp(-1*alpha*x^2)
 */
double gauss_fcn(double alpha, double x)
{
        double result = 0.0;
        double negone = -1.0;
        double coeff = 0.0;
        coeff = negone * alpha * pow(x, 2);
        result = exp(coeff);
        return result;
}


/*
 * matmul_dgemm: matrix multiply two arrays with fortran routine DGEMM
 */
int matmul_dgemm(double *mata, int row_a, int col_a, double *matb, int row_b,
	   int col_b, double *matc, int row_c, int col_c)
{
    int error = 0; /* error flag */
    /* DGEMM variables */
    unsigned char transa[1] = "n"; /* Do not use transpose of matrix A */
    unsigned char transb[1] = "n"; /* Do not use transpose of matrix B */
    long long int m, n, k;
    long long int lda, ldb, ldc;
    double alpha, beta;
    /* set dgemm variables */
    m = (long long int) row_a;
    n = (long long int) col_c;
    k = (long long int) row_b;
    alpha = 1.0;
    beta = 0.0;
    lda = row_a;
    ldb = row_b;
    ldc = row_c;
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, mata, &lda, matb, &ldb, &beta,
	   matc, &ldc);
    return error;
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
