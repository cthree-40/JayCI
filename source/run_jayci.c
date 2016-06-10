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
#include "moindex.h"
#include "execute_ci_calculation.h"
#include "run_jayci.h"

/*
 * run_jayci: Executes a determinant-based CI calculation.
 */
int run_jayci()
{
	/* .. local scalars .. */
	int error = 0; /* error handling flag */
	int m1len = 0; /* number of 1-e integrals */
	int m2len = 0; /* number of 2-e integrals */
	int itype = 1; /* integral type. Always 1 */
	double nucrep_energy = 0.0; /* nuclear repulsion energy */
	double frzcore_energy = 0.0; /* frozen core orbital energy */
	
	/* .. local arrays .. */
	double *moints1 = NULL; /* 1-e integrals */
	double *moints2 = NULL; /* 2-e integrals */
	char moflname[FLNMSIZE] = {""}; /* SIFS integral filename */
	
	/* .. &general scalars .. */
	int electrons = 0; /* total number of electrons */
	int orbitals = 0; /* total number of orbitals */
	int nfrzc = 0; /* number of frozen core orbitals */
	int nfrzv = 0; /* number of frozen virtual orbitals */
	int nactv = 0; /* number of CAS active orbitals */
	int ndocc = 0; /* number of CAS doubly-occupied orbitals */
	int xlvl = 0; /* excitation level. must be <= 2 */
	int printlvl = 0; /* print level */
	
	/* .. input.jayci scalars .. */
	int ci_aelec = 0; /* ci alpha electrons */
	int ci_belec = 0; /* ci beta electrons */
	int ci_orbs = 0; /* active ci orbitals */
	int nastr = 0; /* number of alpha strings */
	int nbstr = 0; /* number of beta strings */
	int ndets = 0; /* number of determinants */
	
	/* Read &general namelist input. */
	readgeninput(&electrons, &orbitals, &nfrzc, &ndocc, &nactv, &xlvl,
		     &nfrzv, &printlvl, &error);
	if (error != 0) {
		error_flag(error, "run_jayci");
		return error;
	}
	/* Check size of CAS orbital space. */
	if ((ndocc + nactv) > 64) {
		error = ndocc + nactv;
		error_flag(error, "run_jayci");
		fprintf(stderr, "CI-active orbitals = %d.", error);
		fprintf(stderr,
			"CI-active orbital number must be less than 64.\n");
		return error;
	}
	
	/* Read input.jayci input file. This contains expansion information. */
	error = readinputjayci(&ci_aelec, &ci_belec, &ci_orbs, &nastr, &nbstr,
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
	error = execute_ci_calculation(ci_aelec, ci_belec, ci_orbs, nastr, nbstr,
				       ndets, ndocc, nactv, moints1, moints2,
				       nucrep_energy, frzcore_energy, printlvl);
	if (error != 0) {
		error_flag(error, "run_jayci");
		return error;
	}
	
	return error;
}
    
    
