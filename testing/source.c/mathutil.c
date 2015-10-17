// File: mathutil.c
/*********************************************************************
 * mathutil: math functions interfacing with LAPACK
 * -------------------------------------------------------------------
 * diagmat_dsyevr: diagonalizes a square matrix using dsyevr
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathutil.h"

/** fortran subroutines **/
extern double _dlamch_fcn_();
/*************************/

/* diagmat_dsyevr: diagonalizes a square matrix using the dsyevr
 *                 LAPACK subroutine
 * -------------------------------------------------------------------
 * Input:
 *  mat  = matrix to diagonalize
 *  dim  = dimension of mat
 * Output:
 *  evec = eigenvectors of mat
 *  eval = eigenvalues of mat */
int diagmat_dsyevr(double *mat, int dim, double *evecs, double *evals)
{
    unsigned char jobz[1]  = "v"; /* return eigenvectors and eigenvalues */
    unsigned char range[1] = "a"; /* compute all eigenvectors */
    unsigned char uplo[1]  = "u"; /* upper triangle stored */
    long long int n; /* order of matrix */
    long long int lda; /* leading dimension of a */
    double vl, vu; /* lower and upper bounds of eigenvalue interval */
    double il, iu; /* indices of smallest and largest eigenvalues to find */
    double abstol; /* absolute error tolerance */
    long long int m; /* total number of eigenvalues found */
    long long int ldz; /* leading dimension of eigenvector array */
    long long int *isuppz; /* support array */
    double *work; /* workspace array */
    long long int lwork; /* dimension of workspace array */
    long long int *iwork; /* integer workspace array */
    long long int liwork; /* dimension of integer workspace array */
    long long int *info; /* information */

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
	fprintf(stderr,"*** ERROR in DSYEVR_! info = %d ***\n", info);
	err = -1;
    }
    return err;
}
    
