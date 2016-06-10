// File: det2string.c

/* Convert struct det to orbital strings. Results are printed. */

#include <stdio.h>
#include <stdlib.h>
#include "binarystr.h"
#include "bitutil.h"
#include "det2string.h"

void det2string(struct det d)
{
	int nzb[64];
	int i;
	
	printf("ALPHA STRING = ");
	nonzerobits(d.astr.byte1, 64, nzb);
	for (i = 0; i < 4; i++) {
		if (nzb[i] > 0) printf(" %d ", nzb[i]);
	}
	for (i = 0; i < d.astr.nvrtx; i++) {
		printf(" %d ", d.astr.virtx[i]);
	}
	printf("\n");
	printf("BETA STRING = ");
	nonzerobits(d.bstr.byte1, 64, nzb);
	for (i = 0; i < 4; i++) {
		if (nzb[i] > 0) printf(" %d ", nzb[i]);
	}
	for (i = 0; i < d.bstr.nvrtx; i++) {
		printf(" %d ", d.bstr.virtx[i]);
	}
	printf("\n\n");
	return;
}
