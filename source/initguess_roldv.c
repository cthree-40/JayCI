// File: initguess_roldv.c

/*
 * Read old vectors for initial guess.
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "ioutil.h"
#include "arrayutil.h"
#include "mathutil.h"
#include "initguess_roldv.h"

/*
 * initguess_roldv: read old vectors for initial CI guess. This assumes
 * there is a file, civector.dat, present in the current working direcotry.
 * If this file is not found, this routine returns an error.
 */
int initguess_roldv(double **vscr, int krymin, int orbs, int aelec,
                    int belec, int ndets)
{
        int error = 0; /* error flag */
        int nroots = 0; /* roots to read */
        double *cival = NULL; /* Eigenvalues (throw away) */
        
        cival = (double *) malloc(krymin * sizeof(double));
        init_dbl_2darray_0(vscr, ndets, krymin);
        nroots = krymin;
        error = read_civectors(aelec, belec, orbs, vscr, ndets, nroots, cival);
        if (error != 0) {
                if (error == 1) {
                        return error;
                } else if (error == 2) {
                        return error;
                } else {
                        error_message("unknown error 1", "initguess_roldv");
                        error = 100;
                        return error;
                }
        }
        if (error != 0) {
                error_message("unknown error 2", "initguess_roldv");
                return error;
        }

        return error;
}
                        
                
