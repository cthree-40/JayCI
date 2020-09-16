// File: pdyson_global.h
/*
 * Global variables for pdyson execution.
 */
#ifndef pdyson_global_h
#define pdyson_global_h

/*
 * Wavefunction information.
 * 0 refers to the N+1 electron wavefunction, and 1 refers to the
 * N electron wavefunction.
 */
extern int nstates0;  /* Number of states */
extern int nelecs0;   /* Number of electrons */
extern int naelec0;   /* Number of alpha electrons */
extern int nbelec0;   /* Number of beta  electrons */
extern int nfrzc0;    /* Number of frozen core orbitals */
extern int ndocc0;    /* Number of DOCC orbitals */
extern int nactv0;    /* Number of CAS  orbitals */
extern int nfrzv0;    /* Number of frozen virtual orbitals */
extern int ninto0;    /* Number of internal orbitals */
extern int xlvl0;     /* CI excitation level */
extern int norbs0;    /* Number of orbitals */
extern int ciorbs0;   /* NUmber of ci orbitals */
extern int ciaelec0;  /* Number of CI alpha electrons */
extern int cibelec0;  /* Number of CI beta  electrons */
extern int pstr0_len; /* Number of alpha strings */
extern int qstr0_len; /* Number of beta  strings */

extern int nstates1;  /* Number of states */
extern int nelecs1;   /* Number of electrons */
extern int naelec1;   /* Number of alpha electrons */
extern int nbelec1;   /* Number of beta  electrons */
extern int nfrzc1;    /* Number of frozen core orbitals */
extern int ndocc1;    /* Number of DOCC orbitals */
extern int nactv1;    /* Number of CAS  orbitals */
extern int nfrzv1;    /* Number of frozen virtual orbitals */
extern int ninto1;    /* Number of internal orbitals */
extern int xlvl1;     /* CI excitation level */
extern int norbs1;    /* Number of orbitals */
extern int ciorbs1;   /* NUmber of ci orbitals */
extern int ciaelec1;  /* Number of CI alpha electrons */
extern int cibelec1;  /* Number of CI beta  electrons */
extern int pstr1_len; /* Number of alpha strings */
extern int qstr1_len; /* Number of beta  strings */


#endif
