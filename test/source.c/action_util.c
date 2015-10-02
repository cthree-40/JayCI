// File: action_util.c
/*********************************************************************
 * action_util
 * -----------
 * Utilities for performing Hv=c
 *
 * Subfunctions:
 *  hmatels: compute matrix element <i|H|j>
 *
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include "bitutil.h"
#include "binarystr.h"
#include "action_util.h"

/* hmatels: compute matrix element <i|H|j>
 * -------------------------------------------------------------------
 * Input:
 *  deti = determinant i
 *  detj = determinant j
 *  moints1 = 1-e integrals
 *  moints2 = 2-e integrals
 *  m1len   = length of 1-e integrals
 *  m2len   = length of 2-e integrals
 * Returns:
 *  val = <i|H|j> */
double hmatels(struct det deti, struct det detj, double *moints1,
	       double *moints2, int m1len, int m2len)
{
    double val;
    
    /* .. local scalars ..
     * detdiff = differences between determinants
     * numax = number of alpha excitations
     * numbx = number of beta excitations
     * axi = alpha excitations initial orbitals
     * axf = alpha excitations final orbitals
     * bxi = beta excitaitons initial orbitals
     * bxf = beta excitations final orbitals */
    
    int detdiff; 
    int numax, numbx;
    long long int axi, axf, bxi, bxf;
    
    val = 0.0;
    /* test if determinants are CAS-flagged */
    if (deti.cas + detj.cas == 0) {
	/* non-CAS determinants  */
	detdiff = comparedets_virt(deti, detj, &numax, &numbx,
				   &axi, &axf, &bxi, &bxf);
	if (detdiff > 2) return val;
    } else {
	/* CAS-flagged 
	 * Note: return from comparedets_cas must be divided by 2 */
	detdiff = comparedets_cas(deti, detj, &numax, &numbx, 
				  &axi, &axf, &bxi, &bxf);
	detdiff = detdiff / 2;
	if (detdiff > 2) return val;
    }
    
    val = evaluate_dets(detdiff, deti, detj, numax, numbx,
			axi, axf, bxi, bxf);
    return val;
}
/* evaluate_dets: evaluate matrix element <i|H|j>.
 * -------------------------------------------------------------------
 * Input:
 *  ndiff = number of differences (excitation level)
 *  deti  = determinant i
 *  detj  = determinant j
 *  numax = number of alpha excitations
 *  numbx = number of beta  excitations
 *  axi   = alpha excitation initial orbitals
 *  axf   = alpha excitation final orbitals
 *  bxi   = beta  excitation initial orbitals
 *  bxf   = beta  excitation final orbitals
 * Returns:
 *  value = value of matrix element */
double evaluate_dets(int ndiff, struct det deti, struct det detj, int numax, 
		     int numbx, long long int axi, long long int axf, 
		     long long int bxi, long long int bxf)
{
    double value;
    

    value = 1.0;
    
    return value;
}
