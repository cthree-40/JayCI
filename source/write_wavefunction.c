// File: write_wavefunction.c
/*
 * Write out a wavefunction: Expansion, vectors.
 */
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "ioutil.h"
#include "binarystr.h"
#include "write_wavefunction.h"

void write_wavefunction(struct det *detlist, int ndets, int nroots,
                        double **civec, double *cival, int orbs, int ninto,
                        int elec)
{
        int error = 0;
        FILE *fptr = NULL;
        int electrons, orbitals;
        int nfrzncore, ndocc, nactive, nfrznvirt;
        int xlevel;
        int printlvl, printwvf;
        int i, j;

        /* Read &general input file to get wavefunction information */
        readgeninput(&electrons, &orbitals, &nfrzncore, &ndocc, &nactive,
                     &xlevel, &nfrznvirt, &printlvl, &printwvf, &error);
        if (error != 0) {
                error_message("Could not input file.", "write_wavefunction");
                return;
        }
        
        /* Open file and print wavefunction information */
        fptr = fopen("ciwvfcn.dat", "w");
        if (fptr == NULL){
                error_message("Could not open ciwvfcn.dat!",
                              "write_wavefunction");
                return;
        }
        fprintf(fptr, "%d %d\n", electrons, orbitals);
        fprintf(fptr, "%d %d %d %d\n", nfrzncore, ndocc, nactive, nfrznvirt);
        fprintf(fptr, "%d\n", xlevel);
        fprintf(fptr, "%d %d %d %d\n\n", elec, orbs, ninto, ndets);
        /* Loop over determinants, printing information. Alpha string is
         * printed first. */
        for (i = 0; i < ndets; i++) {
                fprintf(fptr, "%lld %lld %d %d %d  %lld %lld %d %d %d  %d\n",
                        detlist[i].astr.byte1, detlist[i].astr.byte2,
                        detlist[i].astr.virtx[0], detlist[i].astr.virtx[1],
                        detlist[i].astr.nvrtx,
                        detlist[i].bstr.byte1, detlist[i].bstr.byte2,
                        detlist[i].bstr.virtx[0], detlist[i].bstr.virtx[1],
                        detlist[i].bstr.nvrtx,
                        detlist[i].cas);
        }
        /* Print vectors */
        fprintf(fptr, "\n\n");
        fprintf(fptr, "%d\n\n", nroots);
        for (i = 0; i < nroots; i++) {
                fprintf(fptr, "%15.8lf\n\n", cival[i]);
                for (j = 0; j < ndets; j++) {
                        fprintf(fptr, "%15.12lf\n", civec[i][j]);
                }
                fprintf(fptr, "\n");
        }
        fclose(fptr);
        return;
}
        
                
