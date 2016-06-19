// File: genbindet.h

/*
 * Generate binary determinant list.
 */

#ifndef genbindet_h
#define genbindet_h

/* genbinarydetlist: generate list of determinants in binary format
 * -------------------------------------------------------------------
 * Input:
 *  aelec = alpha electrons
 *  belec = beta  electrons
 *  orbs  = orbitals
 *  ndocc = docc orbitals
 *  nactv = active orbitals
 * Output:
 *  dlist = determinant list
 * Returns:
 *  error = 0: no error, 1: error
 */
int genbinarydetlist(
    /* Determinant list */
    struct det *dlist,
    /* Alpha electrons */
    int aelec,
    /* Beta electrons */
    int belec,
    /* Number of orbitals */
    int orbs,
    /* Number of CAS DOCC orbitals */
    int ndocc,
    /* Number of CAS active orbitals */
    int nactv,
    /* Number of determinants */
        int ndets
    );
    
#endif

