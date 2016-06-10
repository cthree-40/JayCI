// File: initguess_sbd.c

/*
 * Generate initial guess for davidson algorithm by diagonalizing a
 * subblock of the CI Hamiltonian.
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "arrayutil.h"
#include "iminmax.h"
#include "binarystr.h"
#include "cimapping.h"
#include "action_util.h"
#include "mathutil.h"
#include "det2string.h"
#include "initguess_sbd.h"

/*
 * initguess_sbd: diagonalize a sublock of CI Hamiltonian and return
 * krymin eigenvectors.
 */
int initguess_sbd(struct det *dlist, int ndets, double *moints1,
		  double *moints2, int aelec, int belec, int ninto,
		  int krymin, double **vscr)
{
    int error = 0; /* error flag */
    double *hmat = NULL; /* Hamiltonian subblock */
    double *hevec = NULL; /* eigenvectors of subblock */
    double *heval = NULL; /* eigenvalues of subblock */
    int sbsize = 0; /* subblock size */
    int i = 0, j = 0;
    int hp = 0, hp2 = 0;

    /* Get reference space size */
    sbsize = int_min(1000, ndets);

    /* Allocate arrays */
    hmat = (double *) malloc(sbsize * sbsize * sizeof(double));
    hevec= (double *) malloc(sbsize * sbsize * sizeof(double));
    heval= (double *) malloc(sbsize * sizeof(double));
    init_dbl_array_0(hmat, (sbsize * sbsize));
    init_dbl_array_0(hevec,(sbsize * sbsize));
    init_dbl_array_0(heval, sbsize);
    
    /* build subblock of Hamiltonian */
    hp = 0;
    for (i = 0; i < sbsize; i++) {
	for (j = i; j < sbsize; j++) {
	    hp = (i * sbsize) + j;
	    hmat[hp] = 0.0 + hmatels(dlist[i], dlist[j], moints1, moints2,
				     aelec, belec, ninto);
	    hp2 = (j * sbsize) + i;
	    hmat[hp2] = 0.0 + hmat[hp];
	}

    }

    /* diagonalize matrix */
    error = diagmat_dsyevr(hmat, sbsize, hevec, heval);
    if (error != 0) {
	error_flag(error, "initguess_sbd");
	return error;
    }
    fprintf(stdout, " Lowest eigenvalue in subblock = %15.8lf\n", heval[0]);
    
    /* set values for initial guess */
    for (i = 0; i < krymin; i++) {
	for (j = 0; j < sbsize; j++) {
	    vscr[i][j] = 0.0 + hevec[(i * sbsize) + j];
	}
    }

    free(hevec);
    free(hmat);
    free(heval);
    return error;
}
