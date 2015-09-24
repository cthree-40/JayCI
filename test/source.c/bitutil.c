// File: bitutil.c
/*********************************************************************
 * bitutil: utilities for byte to byte bitwise comparisons
 * -------------------------------------------------------------------
 *
 * ndiffbytes: compute number of differences between two bytes
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include "bitutil.h"

/* ndiffbytes: compute number of differences between two bytes
 * -------------------------------------------------------------------
 * Input:
 *  byte1 = 64 bit byte
 *  byte2 = 64 bit byte
 * Output:
 *  diffs = number of differences between two bytes */
int ndiffbytes(long long int byte1, long long int byte2)
{
    int i, ndiffs;
    long long int xorbit;

    xorbit = byte1 ^ byte2;
    ndiffs = 0;
    for (i = 0; i < 64; i++) {
	if ((xorbit >> i) & 0x01)
	    ndiffs = ndiffs + 1;
    }
    
    return ndiffs;
}
