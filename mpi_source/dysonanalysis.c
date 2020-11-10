// File: dysonanalysis.c
/* Analyze dyson orbitals */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dysonio.h"

double *allocate_mem_double(double ***a, int n, int m);
void    analyze_dyson_orbitals(int ndyo, int norbs, double **dyorb);
double  compute_vector_norm(double *v, int l);
void    get_maxcoeff_info(double *dyorb, int norbs, double *maxval, int *indx);

/* Main routine */
int main ()
{
#define MAXSTATES 10
    
    int error = 0;   /* Error flag */
    
    int num_dyo = 0; /* Number of dyson orbitals */
    int ndyst0  = 0; /* Number of anion dyson orbital states */
    int ndyst1  = 0; /* Number of neutral dyson orbital states */
    int norbs   = 0; /* Number of molecular orbitals */
    int *dyst0  = NULL; /* Anion dyson orbital states */
    int *dyst1  = NULL; /* Neutral dyson orbital states */

    double **dyson = NULL; /* Dyson orbitals */
    double *dyson1d= NULL; /* 1-d memory block for dyson orbitals */

    double dnorm = 0.0;
    
    int i;

    /* Get dyson orbital information. */
    dyst0 = malloc(sizeof(int) * MAXSTATES);
    dyst1 = malloc(sizeof(int) * MAXSTATES);
    error = read_dysonorbital_info("dysonorb.dat", &num_dyo, &ndyst0, &ndyst1,
				   &norbs, dyst0, dyst1);
    print_dysonorbital_info(num_dyo, ndyst0, ndyst1, norbs, dyst0, dyst1);

    /* Allocate dyson orbital memory array and read in dyson orbitals */
    dyson1d = allocate_mem_double(&dyson, norbs, num_dyo);
    error = read_dysonorbitals_from_file("dysonorb.dat", num_dyo, norbs, dyson);

    /* Analyze dyson orbitals */
    analyze_dyson_orbitals(num_dyo, norbs, dyson);

    return 0;
    
}

/* Subfunctions */
/*
 * 2-d memory allocator
 */
double *allocate_mem_double(double ***array, int n, int m)
{
    double *ptr = NULL; /* Data */
    int i = 0;
    if ((ptr = (double *) malloc(n * m * sizeof(double))) == NULL) {
	return ptr;
    }
    *array = (double **) malloc(m * sizeof(double *));
    for (i = 0; i < m; i++) {
	(*array)[i] = ptr + i * n;
    }
    return ptr;
}

/*
 * analyze_dyson_orbitals: perform analysis of dyson orbitals.
 */
void analyze_dyson_orbitals(int ndyo, int norbs, double **dyson)
{
    double *dnorms = NULL; /* Dyson orbital norms */
    double *maxcval= NULL; /* Max coefficient value */
    int    *maxcidx= NULL; /* Max coefficient index */
    double *percoef= NULL; /* % of coefficient */
    int i;

    /* Allocate arrays */
    dnorms = malloc(sizeof(double) * ndyo);
    maxcval= malloc(sizeof(double) * ndyo);
    maxcidx= malloc(sizeof(int) * ndyo);
    percoef= malloc(sizeof(double) * ndyo);

    /* Compute norms */
    for (i = 0; i < ndyo; i++) {
	dnorms[i] = compute_vector_norm(dyson[i], norbs);
    }
    printf(" Norms: ");
    for (i = 0; i < ndyo; i++) {
	printf(" %8.5lf", dnorms[i]);
    }
    printf("\n");

    /* Get max coefficient values and indices */
    for (i = 0; i < ndyo; i++) {
	get_maxcoeff_info(dyson[i], norbs, &maxcval[i], &maxcidx[i]);
	percoef[i] = maxcval[i] / dnorms[i];
    }
    printf(" Max coef: ");
    for (i = 0; i < ndyo; i++) {
	printf(" %8.5lf (%5.3lf)", maxcval[i], percoef[i]);
    }
    printf("\n");

    
    return;
}

/*
 * compute_vector_norm: compute the norm of a vector
 */
double compute_vector_norm(double *v, int l)
{
    double nrm = 0.0;
    for (int i = 0; i < l; i++) {
	nrm += v[i] * v[i];
    }
    nrm = sqrt(nrm);
    return nrm;
}

/*
 * get_maxcoeff_info: get maximum coefficient information from dyson orbital
 */
void get_maxcoeff_info(double *dyorb, int norbs, double *maxval, int *indx)
{
    int    curr_idx = 0;
    double curr_val = 0.0;
    int i;
    curr_idx = 0;
    curr_val = dyorb[0];
    for (i = 1; i < norbs; i++) {
	if (curr_val < dyorb[i]) {
	    curr_val = dyorb[i];
	    curr_idx = i;
	}
    }
    *maxval = curr_val;
    *indx = curr_idx;
    return;
}
