// File: execute_ci_calculation.h

/* Perform CI calculation. */

#ifndef execute_ci_calculation_h
#define execute_ci_calculation_h

/*
 * execute_ci_calculation: perform CI calculation using davidson algorithm.
 */
int execute_ci_calculation(
        int aelec, /* alpha electrons */
        int belec, /* beta electrons */
        int orbs,  /* orbitals */
        int nastr, /* number of alpha strings */
        int nbstr, /* number of beta strings */
        int ndets, /* number of determinants */
        int ndocc, /* number of CAS docc orbitals */
        int nactv, /* number of CAS active orbitals */
        double *moints1, /* 1-e integrals */
        double *moints2, /* 2-e integrals */
        double nucrep, /* nuclear repulsion energy */
        double frzcore, /* frozen core energy */
        int plvl, /* print level */
        int pwvf  /* print wavefunction */
        );

/*
 * execute_ci_calculation_inpdlist: perform CI calculation using davidson
 * algorithm. This routine is 'given' the determinant list on input.
 */
int execute_ci_calculation_inpdlist(
        int aelec,
        int belec,
        int orbs,
        int nastr,
        int nbstr,
        int ndets,
        int ndocc,
        int nactv,
        double *moints1,
        double *moints2,
        double nucrep,
        double frzcore,
        struct det *detlist,
        int plvl,
        int pwvf);


#endif
