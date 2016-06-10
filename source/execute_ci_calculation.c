// File: execute_ci_calculation.c

/*
 * Perform CI calculation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "errorlib.h"
#include "allocate_mem.h"
#include "binarystr.h"
#include "ioutil.h"
#include "genbindet.h"
#include "cimapping.h"
#include "action_util.h"
#include "davidson.h"
#include "det2string.h"
#include "execute_ci_calculation.h"

/*
 * execute_ci_calculation: perform CI calculation using davidson algorithm.
 */
int execute_ci_calculation(int aelec, int belec, int orbs, int nastr, int nbstr,
			   int ndets, int ndocc, int nactv, double *moints1,
			   double *moints2, double nucrep, double frzcore,
			   int plvl)
{
	int error = 0; /* error flag */
	double totfrze = 0.0; /* nucrep + frzcore */
	int ninto = 0; /* internal orbitals (docc + active) */
	
	int maxiter = 0; /* maximum number of davidson algorithm iterations */
	int krymin = 0; /* minimum size of the krylov space */
	int krymax = 0; /* maximum size of the krylov space */
	int nroots = 0; /* number of roots to solve for */
	int prediag_routine = 0; /* prediagonalization routine */
	int refdim = 0; /* prediagonalization reference space */
	double restol = 0.0; /* residual norm converegence tolerance */
	
	struct det *detlist = NULL; /* determinant list */
	double *hdgls = NULL; /* diagonal elements of hamiltonian */
	double **civec = NULL; /* final CI eigenvectors */
	double *cival = NULL; /* final ci eigenvalues */
	struct rowmap *hmap = NULL; /* valid <i|H|j> combinations */
	
	clock_t curr_time, prev_time, start_time;
	
	start_time = clock();
	fprintf(stdout, " Start time: %10.5lf\n", (double) start_time);
	
	/* Compute total contribution of nuclei and frozen core. Compute the
	   number of internal orbitals. */
	totfrze = nucrep + frzcore;
	ninto = ndocc + nactv;
	
	/* Read &dgalinfo namelist. */
	readdaiinput(&maxiter, &krymin, &krymax, &nroots, &prediag_routine,
		     &refdim, &restol, &error);
	if (error != 0) {
		error_flag(error, "execute_ci_calculation");
		return error;
	}
	
	/* Generate binary determinant list. */
	prev_time = clock();
	if (plvl > 0) fprintf(stdout, " Generating binary determinant list.\n");
	detlist = (struct det *) malloc(ndets * sizeof(struct det));
	init_detlist(detlist, ndets);
	error = genbinarydetlist(detlist, aelec, belec, orbs, ndocc, nactv,
				 ndets);
	if (error != 0) {
		error_flag(error, "execute_ci_calculation");
		return error;
	}
	curr_time = clock();
	if (plvl > 1) fprintf(stdout, " Total time reading list: %10.5f sec\n",
			      (double) (curr_time - prev_time)/CLOCKS_PER_SEC);
	
	/* Generate CI map of valid determinants. */
	prev_time = clock();
	if (plvl > 0) fprintf(stdout, " Generating cimap.\n");
	hmap = (struct rowmap *) malloc(ndets * sizeof(struct rowmap));
	error = generate_cimap(detlist, ndets, ninto, hmap);
	if (error != 0) {
		error_flag(error, "execute_ci_calculation");
		return error;
	}
	curr_time = clock();
	if (plvl > 1) fprintf(stdout,
			      " Total time generating ci map: %10.5f sec\n",
			      (double) (curr_time - prev_time)/CLOCKS_PER_SEC);
	
	/* Compute diagonal elements <i|H|i>. */
	if (plvl > 0) fprintf(stdout, " Computing diagonal elements.\n");
	hdgls = (double *) malloc(ndets * sizeof(double));
	compute_hdgls(detlist, ndets, moints1, moints2, aelec, belec, hdgls,
		      ninto);
	if (plvl > 0) fprintf(stdout, " Hartree Fock Energy = %15.8lf\n",
			      (hdgls[0] + totfrze));
	
	/* Allocate civec(m,n) and cival(m).*/
	error = allocate_mem_double(&civec, ndets, nroots);
	if (error != 0) {
		error_flag(error, "execute_ci_calculation");
		return error;
	}
	cival = (double *) malloc(nroots * sizeof(double));
	/* Davidson algorithm */
	error = dvdalg(detlist, ndets, moints1, moints2, aelec, belec, hdgls,
		       ninto, totfrze, maxiter, krymin, krymax, nroots, restol,
		       hmap, civec, cival, prediag_routine, plvl);
	return error;
}
