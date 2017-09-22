// File: execute_dycicalc.h

#ifndef execute_dycicalc_h
#define execute_dycicalc_h

int execute_dycicalc(
        char *wvfcn_file0, /* Wavefunction file for N+1 electron system  */
        char *wvfcn_file1, /* Wavefunction file for N   electron system  */
        int nstates0,      /* Number of states to read in for N+1 system */
        int nstates1,      /* Number of states to read in for N   system */
        int nelecs0,       /* N+1 electrons */
        int nelecs1,       /* N   electrons */
        int orbitals,      /* Molecular orbitals */
        int ndets0,        /* N+1 electron wavefunction determinants */
        int ndets1         /* N   electron wavefunction determinants */
        );

#endif
