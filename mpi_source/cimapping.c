// File: cimapping.c
/*
 * Hamiltonian mapping utilities.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "binarystr.h"
#include "cimapping.h"

/*
 * convertrow2map: convert row of bits into rowmap structure.
 */
void convertrow2map(int *row, int n2v, struct rowmap *map)
{
	int i = 0;
	int currstat; /* current status: 0 or 1 */
	struct sectionll *head; /* head of linked list */
	struct sectionll *csec; /* current section */
	struct sectionll *tmp; /* temporary section */
	
        /* allocate first section */
	head = (struct sectionll *) malloc(sizeof(struct sectionll));
	csec = head;
	csec->next = NULL;
	map->nsec = 1; /* at least one section */

	/* loop over row */
	csec->first = 0; /* diagonal is always on */
	currstat = 1; /* and our current status is "ON" */
	for (i = 1; i < n2v; i++) {
		if (currstat == row[i]) continue;
		/* new section found */
		if (row[i-1] == 1) {
			csec->last = i - 1;
			csec->next = (struct sectionll *)
				malloc(sizeof(struct sectionll));
			csec->next->prev = csec;
			csec = csec->next;
			csec->next = NULL;
			currstat = 0;
		} else {
			csec->first = i;
			map->nsec++;
			currstat = 1;
		}
	}
	if (currstat == 1) {
		csec->last = n2v - 1;
		csec->next = NULL;
	} else {
		tmp = csec;
		csec = csec->prev;
		csec->next = NULL;
		free(tmp);
	}
	/* Convert sections to arrays. */
	map->sec = (struct section *) malloc(map->nsec * sizeof(struct section));
	csec = head;
	for (i = 0; i < map->nsec; i++) {
		map->sec[i].first = csec->first;
		map->sec[i].last = csec->last;
		csec = csec->next;
	}
	/* Deallocate linked list */
	csec = head;
	while (csec->next != NULL) {
		tmp = csec;
		csec = csec->next;
		free(tmp);
	}
	free(csec);
	return;
}

/*
 * pci_generate_cimap: generate rowmap for <a|H|i>..<b|H|i>, where a = lwrbnd,
 * and b = upprbnd. Note: hmap is a linked list.
 */
int pci_generate_cimap(struct det *dlist, int ndets, int nactv, int lwrbnd,
		       int upprbnd, struct rowmap *hmap, int nrows)
{
	int error = 0;
	int i, j = 0;
	int v;
	int ddiff;
	int naxc, nbxc, naxv, nbxv, naxcv, nbxcv;
	long long int axi, axf, bxi, bxf;
	int row[ndets];
	struct rowmap **hmap_index = NULL; /* pointers to hmap rows. This will
					   * allow for threaded loop over rows.*/
	struct rowmap *curr_hmap = NULL; /* pointer to hmap row. */

        /* Allocate hmap linked list and hmap_index arrays */
	hmap_index = (struct rowmap **) malloc(nrows * sizeof(struct rowmap *));
	curr_hmap = hmap;
	for (i = 0; i < nrows; i++) {
		/* Set index. Then allocate next node */
		hmap_index[i] = curr_hmap;
		curr_hmap->next = (struct rowmap *) malloc(sizeof(struct rowmap));
		curr_hmap->next->next = NULL;
		curr_hmap = curr_hmap->next;
	}

	/* -- OMP Section -- */
	/* Loop over hmap_index entries. These are pointers to hmap structures. */
	/* Loop over each row in H, flagging nonzero matrix elements. */
#pragma omp parallel \
	shared(ndets, dlist, nactv, hmap_index, upprbnd, lwrbnd) \
	private(i, j, v, naxc, nbxc, naxv, nbxv, naxcv, nbxcv,   \
		axi, bxi, axf, bxf, row, ddiff)
	{
#pragma omp for schedule(guided)
		
	for (i = lwrbnd; i < upprbnd; i++) {
		v = 0;
		hmap_index[(i - lwrbnd)]->row = i;
		hmap_index[(i - lwrbnd)]->nelm = 0;
		for (j = i; j < ndets; j++) {
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
				/* increment total number of nonzero elements. */
				hmap_index[(i - lwrbnd)]->nelm++;
			}
			v++;
		}
		/* Convert row of 1,0,1,1,... into rowmap structure */
		convertrow2map(row, v, hmap_index[(i - lwrbnd)]);
	}
	} /* -- End OMP Section -- */
	return error;
}

/*
 * get_total_nelm: get total number of 'on' matrix elements for an hmap.
 */
int get_total_nelm(struct rowmap *hmap, int nrows, int *total)
{
	int error = 0;
	struct rowmap *curr_hmap;
        /* Set total to zero. Sum over all rows */
	*total = 0;
	curr_hmap = hmap;
	while(curr_hmap->next != NULL) {
		*total = *total + curr_hmap->nelm;
		curr_hmap = curr_hmap->next;
	}
	return error;
}

/*
 * get_upptri_row: get row of upper triangle matrix.
 */
int get_upptri_row(int element, int matdim)
{
	int result = 0;
	int i = 0;
	for (i = 0; i < matdim; i++) {
		if (upptri_rowpack_index(i, i, matdim) > element) break;
	}
	result = i - 1;
	return result;
}

/*
 * get_upptri_size: get size of upper triangle of square matrix.
 */
int get_upptri_size(int matdim)
{
	int result = 0;
	/* Upper triangle elements are indexed as follows:
	 *  M(i,j) = j(i-1)/2 + i */
	result = ((matdim * matdim) + matdim) / 2;
	return result;
}

/*
 * upptri_rowpack_index: return index of uppertriangle element packed by rows.
 */
int upptri_rowpack_index(int i, int j, int matdim)
{
	int result = 0;
	/* M(i,j) = j + (2*M - i)*(i - 1)/2 */
	result = (2 * matdim - i) * (i - 1);
	result = result / 2;
	result = result + j;
	return result;
}
