// File: read_wavefunction.c
/*
 * Read a wavefunction file: Expansion, civectors, civalues.
 */
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "binarystr.h"
#include "read_wavefunction.h"

/*
 * read_wavefunction: read a wavefunction.
 */
int read_wavefunction(struct det *detlist, int ndets, int nroots,
                      double **civec, double *cival, int ci_orbs, int ci_elec,
                      int ninto, char *file_name, int *orbitals, int *electrons,
                      int *nfrzc, int *ndocc, int *nactv, int *nfrzv, int *xlvl)
{
        int error = 0;
        int fl_ndets, fl_nroots, fl_orbs, fl_elec; /* File values */
        int fl_ninto;
        FILE *fptr = NULL;
        char emsg[100];
        int i, j;
        
        /* Open file */
        fptr = fopen(file_name, "r");
        if (fptr == NULL) {
                error = 1;
                sprintf(emsg, "Could not open file: %s", file_name);
                error_flag(error, "read_wavefunction");
                error_message(emsg, "read_wavefunction");
                return error;
        }
        /* Read first lines */
        fscanf(fptr, "%d %d\n", electrons, orbitals);
        fscanf(fptr, "%d %d %d %d\n", nfrzc, ndocc, nactv, nfrzv);
        fscanf(fptr, "%d\n", xlvl);
        fscanf(fptr, "%d %d %d %d\n\n", &fl_elec, &fl_orbs, &fl_ninto,
               &fl_ndets);
        if (fl_elec != ci_elec || fl_orbs != ci_orbs || fl_ndets != ndets) {
                if (fl_elec != ci_elec) {
                        printf("fl_elec != ci_elec\n");
                }
                if (fl_orbs != ci_orbs) {
                        printf("fl_orbs != ci_orbs\n");
                }
                if (fl_ndets != ndets) {
                        printf("fl_ndets != ndets\n");
                }
                fprintf(stderr, "fl_elec,   elec = %9d, %9d\n", fl_elec,ci_elec);
                fprintf(stderr, "fl_orbs,   orbs = %9d, %9d\n", fl_orbs,ci_orbs);
                fprintf(stderr, "fl_ndets, ndets = %9d, %9d\n", fl_ndets, ndets);
                error = -1;
                error_flag(error, "read_wavefunction");
                error_message("Incorrect parameters.", "read_wavefunction");
                return error;
        }
        /* Read determinants. */
        for (i = 0; i < ndets; i++) {
                fscanf(fptr, "%lld %lld %d %d %d  %lld %lld %d %d %d  %d\n",
                        &detlist[i].astr.byte1, &detlist[i].astr.byte2,
                        &detlist[i].astr.virtx[0], &detlist[i].astr.virtx[1],
                        &detlist[i].astr.nvrtx,
                        &detlist[i].bstr.byte1, &detlist[i].bstr.byte2,
                        &detlist[i].bstr.virtx[0], &detlist[i].bstr.virtx[1],
                        &detlist[i].bstr.nvrtx,
                        &detlist[i].cas);
        }
        fscanf(fptr,"\n\n");
        fscanf(fptr,"%d\n\n", &fl_nroots);
        if (fl_nroots != nroots) {
                error = -2;
                error_flag(error, "read_wavefunction");
                error_message("Incorrect number of roots", "read_wavefunction");
                return error;
        }
        for (i = 0; i < nroots; i++) {
                fscanf(fptr, "%lf\n\n", &(cival[i]));
                for (j = 0; j < ndets; j++) {
                        fscanf(fptr, "%lf\n", &(civec[i][j]));
                }
                fscanf(fptr, "\n");
        }
        fclose(fptr);
        return error;
}
