// File: initguess_sbd.h

/* Diagonalize subblock of CI hamiltonian for initial guess. */

#ifndef initguess_sbd_h
#define initguess_sbd_h

/*
 * initguess_sbd: diagonalize a sublock of CI Hamiltonian and return
 * krymin eigenvectors.
 */
int initguess_sbd(
    struct det *dlist, /* determinant lists */
    int ndets,         /* number of determinants */
    double *moints1,   /* 1-e integrals */
    double *moints2,   /* 2-e integrals */
    int aelec,         /* alpha electrons */
    int belec,         /* beta electrons */
    int ninto,         /* internal orbitals */
    int krymin,        /* minimum dimension of krylov space; also,
                        * output vector number */
    int rdim,          /* Reference space size */
    double **vscr,     /* output vectors */
    double totfrze     /* Total frozen core energy */
        );

#endif
