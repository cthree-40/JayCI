/* File: genbindet.c */

/*
 * Generate binary determinant list.
 */

#include <stdio.h>
#include <stdlib.h>
#include "straddress.h"
#include "binarystr.h"
#include "genbindet.h"

/* genbinarydetlist: generate list of determinants in binary format
 * -------------------------------------------------------------------
 * Input:
 *  aelec = alpha electrons
 *  belec = beta  electrons
 *  orbs  = orbitals
 *  ndocc = docc orbitals
 *  nactv = active orbitals
 * Output:
 *  dlist = determinant list
 * Returns:
 *  error = 0: no error, 1: error
 */
int genbinarydetlist(struct det *dlist, int aelec, int belec, int orbs,
		     int ndocc, int nactv, int ndets)
{
    int i;
    int ai, bi;
    int *aistr, *bistr, *strscr;
    int error;
    FILE *dlistfl;

    error = 0;

    aistr = (int *) malloc(aelec * sizeof(int));
    bistr = (int *) malloc(belec * sizeof(int));
    strscr= (int *) malloc(aelec * sizeof(int));

    dlistfl = fopen("det.list", "r");
    if (dlistfl == NULL) {
	fprintf(stderr, "*** ERROR: Cannot open det.list! ***\n");
	error = 1;
	return error;
    }
    i = 0;
    printf("Reading file.\n");
    while (fscanf(dlistfl, " %d %d\n", &ai, &bi) != EOF) {
	str_adr2str(ai, strscr, aelec, orbs, aistr);
	str_adr2str(bi, strscr, belec, orbs, bistr);
	dlist[i].astr = str2occstr(aistr, aelec, ndocc, nactv, orbs);
	dlist[i].bstr = str2occstr(bistr, belec, ndocc, nactv, orbs);
	/* is this determinant a CAS determinant? */
	if (dlist[i].astr.virtx[0] == 0 && dlist[i].bstr.virtx[0] == 0) {
	    dlist[i].cas = 1;
	} else {
	    dlist[i].cas = 0;
	    dlist[i].nvo = dlist[i].astr.nvrtx + dlist[i].bstr.nvrtx;
	}
	i++;
    }

    fclose(dlistfl);
    free(aistr);
    free(bistr);
    free(strscr);

    fprintf(stdout, "Read in %d determinants.\n", i);
    if (i != ndets) {
	fprintf(stderr, "*** ERROR: Incorrect number of determinants! ***\n");
	error = 1;
	return error;
    }
    return error;
}
