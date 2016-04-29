// File: execute_ci_calculation.c

/*
 * Perform CI calculation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "errorlib.h"
#include "binarystr.h"
#include "ioutil.h"
#include "genbindet.h"
#include "cimapping.h"
#include "execute_ci_calculation.c"

/*
 * execute_ci_calculation: perform CI calculation using davidson algorithm.
 */
int execute_ci_calculation(int aelec, int belec, int orbs, int nastr, int nbstr,
			   int ndets, int ndocc, int nactv, double *moints1,
			   double *moints2, double nucrep, double frzcore,
			   int plvl)
{
    int error = 0; /* error flag */
    double totfrze; /* nucrep + frzcore */

    int maxiter; /* maximum number of davidson algorithm iterations */
    int krymin; /* minimum size of the krylov space */
    int krymax; /* maximum size of the krylov space */
    int nroots; /* number of roots to solve for */
    int prediag_routine; /* prediagonalization routine */
    int refdim; /* prediagonalization reference space */
    double restol; /* residual norm converegence tolerance */

    struct det *detlist; /* determinant list */
    double *hdgls; /* diagonal elements of hamiltonian */
    double *civec; /* final CI eigenvectors */
    double *cival; /* final ci eigenvalues */
    struct rowmap *hmap; /* valid <i|H|j> combinations */
    
    clock_t curr_time, prev_time, start_time;

    start_time = clock();
    fprintf(stdout, " Start time: %10.5f\n", start_time);

    /* Read &dgalinfo namelist. */
    readdaiinput(&maxiter, &krymin, &krymax, &nroots, &prediag_routine, &refdim,
		 &restol, &error);
    if (error != 0) {
	error_flag(error, "execute_ci_calculation");
	return err;
    }

    /* Generate binary determinant list. */
    prev_time = clock();
    if (plvl > 0) fprintf(stdout, " Generating binary determinant list.\n");
    detlist = (struct det *) malloc(ndets * sizeof(struct det));
    error = genbinarydetlist(detlist, aelec, belec, orbs, ndocc, nactv, ndets);
    if (error != 0) {
	error_flag(error, "execute_ci_calculation");
	return err;
    }
    curr_time = clock();
    if (plvl > 1) fprintf(stdout, " Total time reading list: %10.5f sec\n",
			  (double) (curr_time - prev_time)/CLOCKS_PER_SEC);

    /* Generate CI map of valid determinants. */
    hmap = (struct rowmap *) malloc(ndets * sizeof(struct rowmap));
    error = generate_cimap(detlist, ndets, hmap);
    if (error != 0) {
	error_flag(error, "execute_ci_calculation");
	return err;
    }

    return err;
}
