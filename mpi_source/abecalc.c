// File: abecalc.c
/*********************************************************************
 * abecalc
 * -------
 * Calculate number of alpha and beta electrons.
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include "abecalc.h"

void abecalc(int tot_elec, int *aelec, int *belec)
/* Input:
 *  tot_elec = total number of electrons
 * Output:
 *  aelec = alpha electrons 
 *  belec = beta  electrons
 */
{
	/* test if even number of electrons */
	if (tot_elec % 2 == 0) {
		/* even */
		*aelec = tot_elec / 2;
		*belec = *aelec;
	} else {
		/* odd */
		*aelec = (tot_elec + 1) / 2;
		*belec = (tot_elec - 1) / 2;
	}
	
	return;
}
