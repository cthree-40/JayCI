// File: bitutil.c
/*********************************************************************
 * bitutil: utilities for byte to byte bitwise comparisons
 * -------------------------------------------------------------------
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include "bitutil.h"

/* ndiffsbyte: compute number of differences between two bytes
 * -------------------------------------------------------------------
 * Input:
 *  byte1 = 64 bit byte
 *  byte2 = 64 bit byte
 * Output:
 *  diffs = number of differences between two bytes */
void ndiffsbyte(long long int byte1, long long int byte2, int diffs)
{
    long long int xorbit;

    xorbit = byte1 ^ byte2;
    
    
}
