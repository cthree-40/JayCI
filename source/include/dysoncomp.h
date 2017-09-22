// File: dysoncomp.h

#ifndef dysoncomp_h
#define dysoncomp_h

/*
 * compute_dyson_orbital: compute the dyson orbital between one N+1 electron
 * wavefunction and on N electron wavefunction.
 */
void compute_dyson_orbital(
        struct det *dlist0, /* N+1 electron WVFCN determinants */
        double *civec0,     /* N+1 electron WVFCN civector */
        int ndets0,         /* Number of N+1 electron WVFCN determinants */
        struct det *dlist1, /* N electron   WVFCN determinants */
        double *civec1,     /* N electron   WVFCN civector */
        int ndets1,         /* Number of N   electron WVFCN determinants */
        int orbs,           /* Number of molecular orbitals */
        double *dysonorb    /* Dyson orbital */
        );

/*
 * compute_dyson_orbital_aeq: compute the dyson orbital between one N+1 electron
 * wavefunction and on N electron wavefunction when both wavefunctions have an
 * equal number of alpha electrons.
 */
void compute_dyson_orbital_aeq(struct det *dlist0, double *civec0, int ndets0,
                               struct det *dlist1, double *civec1, int ndets1,
                               int orbs, double *dysonorb);

/*
 * compute_dyson_orbital_beq: compute the dyson orbital between one N+1 electron
 * wavefunction and on N electron wavefunction when both wavefunctions have an
 * equal number of beta electrons.
 */
void compute_dyson_orbital_beq(struct det *dlist0, double *civec0, int ndets0,
                               struct det *dlist1, double *civec1, int ndets1,
                               int orbs, double *dysonorb);

/*
 * comparedets_dyson: compare two determinants for computation of dyson
 * orbital.
 */
int comparedets_dyson(
        struct det det0, /* N+1 electron WVFCN determinant */
        struct det det1  /* N   electron WVFCN determinant */
        );

/*
 * comparedets_dyson_beq: compare two determinants for computation of dyson
 * orbital with equal beta electrons.
 */
int comparedets_dyson_beq(struct det det0, struct det det1);

/*
 * comparedets_dyson_aeq: compare two determinants for computation of dyson
 * orbital with equal alpha electrons.
 */
int comparedets_dyson_aeq(struct det det0, struct det det1);

#endif
