// File: run_jayci.c

/*
 * Main driver for jayci.x
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "ioutil.h"
#include "moindex.h"
#include "run_jayci.h"

/*
 * run_jayci: Executes a determinant-based CI calculation.
 */
int run_jayci(long long int memory)
{
    /* .. local scalars .. */
    int error = 0; /* error handling flag */
    int m1len; /* number of 1-e integrals */
    int m2len; /* number of 2-e integrals */
    int itype = 1; /* integral type. Always 1 */
    double nucrep_energy; /* nuclear repulsion energy */
    double frzcore_energy; /* frozen core orbital energy */

    /* .. &general scalars .. */
    int electrons; /* total number of electrons */
    int orbitals; /* total number of orbitals */
    int nfrzc; /* number of frozen core orbitals */
    int nfrzv; /* number of frozen virtual orbitals */
    int nactv; /* number of CAS active orbitals */
    int ndocc; /* number of CAS doubly-occupied orbitals */
    int xlvl; /* excitation level. must be <= 2 */
    int printlvl; /* print level */

    /* .. input.jayci scalars .. */
    int ci_aelec; /* ci alpha electrons */
    int ci_belec; /* ci beta electrons */
    int ci_orbs; /* active ci orbitals */
    int nastr; /* number of alpha strings */
    int nbstr; /* number of beta strings */
    int ndets; /* number of determinants */
    
    /* Read &general namelist input. */
    readgeninput(&electrons, &orbitals, &nfrzc, &ndocc, &nactv, &xlvl,
		 &nfrzv, &printlvl, &error);
    if (error != 0) {
	error_flag(error, "run_jayci");
	return err;
    }
    /* Check size of CAS orbital space. */
    if ((ndocc + nactv) > 64) {
	error = ndocc + nactv;
	error_flag(error, "run_jayci");
	fprintf(stderr, " CI-active orbitals = %d.", error);
	fprintf(stderr, " CI-active orbital number must be less than 64.\n");
	return error;
    }

    /* Read input.jayci input file. This contains expansion information. */
    err = readinputjayci(&ci_aelec, &ci_belec, &ci_orbs, &nastr, &nbstr,
			 &ndets);
    if (error != 0) {
	error_flag(error, "run_jayci");
	return error;
    }

    /* Read molecular orbitals */
    m1len = index1e(orbitals, orbitals);
    m2len = index2e(orbitals, orbitals, orbitals, orbitals);
    strncpy(moflname, "moints", FLNMSIZE);
    moints1 = (double *) malloc(m1len * sizeof(double));
    moints2 = (double *) malloc(m2len * sizeof(double));
    readmointegrals(moints1, moints2, itype, ci_orbs, moflname, m1len, m2len,
		    &nucrep_energy, &frzcore_energy);
    fprintf(stdout, " Nuclear repulsion energy = %15.8lf\n", nucrep_energy);
    if (nfrzc > 0) {
	fprintf(stdout, " Frozen core energy = %15.8lf\n", frzcore_energy);
    }

    /* Execute CI calculation */
    error = execute_ci_calculation(ci_aelec, ci_belec, ci_orbs, nastr, nbstr,
				   ndets, ndocc, nactv, moints1, moints2,
				   nucrep_energy, frzcore_energy, printlvl);
    if (error != 0) {
	error_flag(error, "run_jayci");
	return error;
    }

    return error;
}
    
    
