// File: cimapping.c

/*
 * Generate linked-list for evaluation of CI Hamiltonian.
 */

#include <stdio.h>
#include <stdlib.h>
#include "binarystr.h"
#include "cimapping.h"

/* -- MPI options -- */
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif
/*------------------ */
  
/*
 * convertrow2map: convert a row of bits into rowmap structure
 */
void convertrow2map(int *row, int n2v, struct rowmap *v2v)
{
	int i, j;
	int currstat; /* current status: 0 or 1 */
	struct section *csec; /* current section */
	struct section *tmp;
	/* allocate first section */
	v2v->sec = (struct section *) malloc(sizeof(struct section));
	v2v->nsec = 1; /* at least one section */
	csec = v2v->sec;
	v2v->sec->next = NULL;
	
	/* loop over row */
	csec->first = 0;  /* diagonal is always on */
	currstat = 1;    /* and our current status is "ON" */
	for (i = 1; i < n2v; i++) {
		if (currstat == row[i]) continue;
		/* new section found */
		if (row[i-1] == 1) {
			csec->last = i - 1;
			csec->next =(struct section *)
				malloc(sizeof(struct section));
			csec->next->prev = csec;
			csec = csec->next;
			csec->next = NULL;
			currstat = 0;
		} else {
			csec->first = i;
			v2v->nsec++;
			currstat = 1;
		}
	}
	if (currstat == 1){
		csec->last = n2v - 1;
		csec->next = NULL;
	} else {
		tmp = csec;
		csec = csec->prev;
		csec->next = NULL;
		free(tmp);
	}
	return;
}

/*
 * generate_cimap: generate map for evaluating nonzero matrix elements of
 * the CI Hamiltonian.
 */
int generate_cimap(struct det *dlist, int ndets, int nactv, struct rowmap *hmap)
{
	int i, j = 0;
	int v;
	int ddiff;
	int naxc, nbxc, naxv, nbxv, naxcv, nbxcv;
	long long int axi, axf, bxi, bxf;
	int row[ndets];

	/* OMP Section */
#pragma omp parallel \
	shared(ndets, dlist, nactv, hmap) \
	private(i, j, v, naxc, nbxc, naxv, nbxv, naxcv, nbxcv, \
		axi, bxi, axf, bxf, row, ddiff)
	{
#pragma omp for schedule(dynamic)
	/* Loop over each row in H. Flag nonzero matrix elements, creating
	 * an array of elements 1 or 0. */
	for (i = 0; i < ndets; i++) {
		v = 0;
		for (j = i; j < ndets; j++) {
			/* Check if determinants are CAS-flagged. */
			if (dlist[i].cas + dlist[j].cas < 2) {
				ddiff = comparedets_ncas(
					dlist[i], dlist[j],
					&naxc, &nbxc, &naxv, &nbxv, &naxcv,
					&nbxcv, &axi, &axf, &bxi, &bxf, nactv);
			} else {
				ddiff = comparedets_cas(
					dlist[i], dlist[j],
					&naxc, &nbxc, &axi, &axf, &bxi, &bxf,
					nactv);
			}
			if (ddiff > 2) {
				row[v] = 0;
			} else {
				row[v] = 1;
			}
			v++;
		}
		/* Convert row of 1,0,1,1,... into rowmap structure */
		convertrow2map(row, v, &(hmap[i]));
	}
	} /* End of OMP Section */
	return 0;
}
