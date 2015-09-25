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

    val = 1.0;

    return val;
}
