//File: write_wavefunction.h
#ifndef write_wavefunction_h
#define write_wavefunction_h

/*
 * Write out a wavefunction: Expansion, vectors.
 */
void write_wavefunction(
        struct det *detlist, /* Expansion information */
        int ndets,           /* Number of determinants */
        int nroots,          /* Number of roots */
        double **civec,      /* CI Vectors */
        double *cival,       /* CI Eigenvalues */
        int orbs,            /* Number of molecular orbitals */
        int ninto,           /* Number of internal orbitals */
        int elec             /* Total electrons */
        );

#endif
