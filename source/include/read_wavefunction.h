// File: read_wavefunction.h

#ifndef read_wavefunction_h
#define read_wavefunction_h

/*
 * read_wavefunction: read a wavefunction.
 */
int read_wavefunction(
        struct det *detlist,
        int ndets,
        int nroots,
        double **civec,
        double *cival,
        int orbs,
        int elec,
        int ninto,
        char *file_name
        );

#endif
