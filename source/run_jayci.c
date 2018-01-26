// File: run_jayci.c

/*
 * Main driver for jayci.x
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errorlib.h"
#include "ioutil.h"
#include "arrayutil.h"
#include "abecalc.h"
#include "moindex.h"
#include "binarystr.h"
#include "citruncate.h"
#include "execute_ci_calculation.h"
#include "run_jayci.h"

/*
 * run_jayci: Executes a determinant-based CI calculation.
 */
int run_jayci()
{
        int error = 0; /* Error flag */

        int m1len = 0; /* number of 1-e integrals */
        int m2len = 0; /* number of 2-e integrals */
        int itype = 1; /* integral type. Always 1 */
        double nucrep_energy = 0.0;  /* nuclear repulsion energy */
        double frzcore_energy = 0.0; /* frozen core orbital energy */
        double *moints1 = NULL; /* 1-e integrals */
        double *moints2 = NULL; /* 2-e integrals */
        char moflname[FLNMSIZE] = {""}; /* SIFS integral filename */

        int electrons = 0; /* total number of electrons */
        int orbitals = 0;  /* total number of orbitals */
        int nfrzc = 0;     /* number of frozen core orbitals */
        int nfrzv = 0;     /* number of frozen virtual orbitals */
        int nactv = 0;     /* number of CAS active orbitals */
        int ndocc = 0;     /* number of CAS doubly-occupied orbitals */
        int xlvl = 0;      /* CI expansion excitation level. must be <= 2 */
        int printlvl = 0;  /* print level */
        int printwvf = 0;  /* print wavefunction info */

        int aelec = 0;    /* Alpha electrons */
        int belec = 0;    /* Beta  eletrons */
        int ci_aelec = 0; /* CI alpha electrons */
        int ci_belec = 0; /* CI beta  electrons */
        int ci_orbs = 0;  /* active CI orbitals */

        int astr_len;  /* Number of alpha strings */
        int bstr_len;  /* Number of beta  strings */
        int dstr_len;  /* Number of determinants */
        struct det *detlist = NULL; /* Determinant list */
        
        /* Read &general namelist input. This will be system information, and
           expansion information. Check that the CAS orbital space is not
           larger than 64, as we are restricted to 64-bit terms. Generate the
           CI expansion wavefunction, and read in the molecular orbitals.
           Execute Davidson algorithm. */
	readgeninput(&electrons, &orbitals, &nfrzc, &ndocc, &nactv, &xlvl,
		     &nfrzv, &printlvl, &printwvf, &error);
	if (error != 0) {
		error_flag(error, "run_jayci");
		return error;
	}
	if ((ndocc + nactv) > 64) {
		error = ndocc + nactv;
		error_flag(error, "run_jayci");
		fprintf(stderr, "CI-active orbitals = %d.", error);
		fprintf(stderr,
			"CI-active orbital number must be less than 64.\n");
		return error;
	}
	abecalc(electrons, &aelec, &belec);
        detlist = citrunc_rtnlist(aelec, belec, orbitals, nfrzc, ndocc, nactv,
                                  nfrzv, xlvl, &astr_len, &bstr_len, &dstr_len,
                                  &ci_orbs, &ci_aelec, &ci_belec);

	m1len = index1e(orbitals, orbitals);
	m2len = index2e(orbitals, orbitals, orbitals, orbitals);
	strncpy(moflname, "moints", FLNMSIZE);
	moints1 = (double *) malloc(m1len * sizeof(double));
	moints2 = (double *) malloc(m2len * sizeof(double));
	init_dbl_array_0(moints1, m1len);
	init_dbl_array_0(moints2, m2len);
	readmointegrals(moints1, moints2, itype, ci_orbs, moflname, m1len, m2len,
			&nucrep_energy, &frzcore_energy);
	fprintf(stdout, "Nuclear repulsion energy = %15.8lf\n", nucrep_energy);
	if (nfrzc > 0) {
		fprintf(stdout,
			"Frozen core energy = %15.8lf\n", frzcore_energy);
	}
	
        /* Execute CI calculation */
	printf("Executing CI calculation.\n");
	error = execute_ci_calculation_inpdlist(ci_aelec, ci_belec, ci_orbs,
                                                astr_len, bstr_len, dstr_len,
                                                ndocc, nactv, moints1, moints2,
                                                nucrep_energy, frzcore_energy,
                                                detlist, printlvl, printwvf);
	if (error != 0) {
		error_flag(error, "run_jayci");
		return error;
	}
	
	return error;

        
	return error;
}
    
    
