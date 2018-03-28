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
		  int krymin, int rdim, double **vscr, double totfrze)
{
    int error = 0; /* error flag */
    double *hmat = NULL; /* Hamiltonian subblock */
    double *hevec = NULL; /* eigenvectors of subblock */
    double *heval = NULL; /* eigenvalues of subblock */
    int sbsize = 0; /* subblock size */
    int i = 0, j = 0;
    int hp = 0, hp2 = 0;

    /* Get reference space size */
    sbsize = int_min(rdim, ndets);

    /* Allocate arrays */
    hmat = (double *) malloc(sbsize * sbsize * sizeof(double));
    hevec= (double *) malloc(sbsize * sbsize * sizeof(double));
    heval= (double *) malloc(sbsize * sizeof(double));
    init_dbl_array_0(hmat, (sbsize * sbsize));
    init_dbl_array_0(hevec,(sbsize * sbsize));
    init_dbl_array_0(heval, sbsize);


    /*double tmp = 0;
    FILE *fptr;
    char fname[80];
    for (j = 0; j < 1000; j++) {
            sprintf(fname, "column%d.txt", (j + 1));
            fptr = fopen(fname, "w");
            for (i = 0; i < sbsize; i++) {
                    tmp = hmatels(dlist[i], dlist[j], moints1, moints2,
                                  aelec, belec, ninto);
                    fprintf(fptr, "%15.8lf", tmp);
                    fprintf(fptr, "\n");
            }
            fclose(fptr);
    }
    tmp = hmatels(dlist[2 - 1], dlist[228 - 1], moints1, moints2,
                  aelec, belec, ninto);
    printf("<2|H|228> = %15.8lf\n", tmp);
    //tmp = hmatels(dlist[254 - 1], dlist[26 - 1], moints1, moints2,
    //              aelec, belec, ninto);
    //printf("<254|H|26> = %15.8lf\n", tmp);
    
    //return (300);
    */

    
    /* build subblock of Hamiltonian */
    hp = 0;
    for (i = 0; i < sbsize; i++) {
            hp = (i * sbsize) + i;
            hmat[hp] = 0.0 + hmatels(dlist[i], dlist[i], moints1, moints2,
                                     aelec, belec, ninto);
	for (j = i + 1; j < sbsize; j++) {
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
    fprintf(stdout, " Lowest eigenvalue in subblock = %15.8lf\n",
            (heval[0] + totfrze));
    //return 1000;
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
