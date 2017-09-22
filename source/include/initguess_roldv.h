// File: initguess_roldv.h

#ifndef initguess_roldv_h
#define initguess_roldv_h

/*
 * initguess_roldv: read old vectors for initial CI guess. This assumes
 * there is a file, civector.dat, present in the current working direcotry.
 * If this file is not found, this routine returns an error.
 */
int initguess_roldv(
        double **vscr, /* guess vectors */
        int krymin,    /* number of guess vectors */
        int orbs,      /* orbitals */
        int aelec,     /* alpha electrons */
        int belec,     /* beta electrons */
        int ndets      /* number of determinants */
        );

#endif
