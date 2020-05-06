// File: pdyson_global.h
/*
 * Global variables for pdyson execution.
 */
#ifndef pdyson_global_h
#define pdyson_global_h

int mpi_num_procs;
int mpi_proc_rank;

extern const int mpi_root; /* GA Root process */

/*
 * Wavefunction information.
 * 0 refers to the N+1 electron wavefunction, and 1 refers to the
 * N electron wavefunction.
 */
int nstates0; /* Number of states */
int nelecs0;  /* Number of electrons */
int naelec0;  /* Number of alpha electrons */
int nbelec0;  /* Number of beta  electrons */
int nfrzc0;   /* Number of frozen core orbitals */
int ndocc0;   /* Number of DOCC orbitals */
int nactv0;   /* Number of CAS  orbitals */
int nfrzv0;   /* Number of frozen virtual orbitals */
int xlvl0;    /* CI excitation level */
int norbs0;   /* Number of orbitals */
int ciorbs0;  /* NUmber of ci orbitals */
int ciaelec0; /* Number of CI alpha electrons */
int cibelec0; /* Number of CI beta  electrons */

int nstates1; /* Number of states */
int nelecs1;  /* Number of electrons */
int naelec1;  /* Number of alpha electrons */
int nbelec1;  /* Number of beta  electrons */
int nfrzc1;   /* Number of frozen core orbitals */
int ndocc1;   /* Number of DOCC orbitals */
int nactv1;   /* Number of CAS  orbitals */
int nfrzv1;   /* Number of frozen virtual orbitals */
int xlvl1;    /* CI excitation level */
int norbs1;   /* Number of orbitals */
int ciorbs1;  /* NUmber of ci orbitals */
int ciaelec1; /* Number of CI alpha electrons */
int cibelec1; /* Number of CI beta  electrons */

#endif
