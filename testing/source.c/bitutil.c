// File: bitutil.c
/*
 * bitutil: utilities for byte to byte bitwise comparisons
 * -------------------------------------------------------------------
 *
 * ndiffbytes: compute number of differences between two bytes
 * nsamebytes: compute number of similarities between two bytes
 *
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 */

#include <stdio.h>
#include <math.h>
#include "bitutil.h"

/* 
 * ndiffbytes: compute number of differences between two bytes
 * 
 * Input:
 *  byte1 = 64 bit byte
 *  byte2 = 64 bit byte
 * Output:
 *  ndiff = number of differences between two bytes 
 *  diffs  = location of each difference as a byte 
 */
int ndiffbytes(long long int byte1, long long int byte2, 
	       int bl, long long int *diffs)
{
	int i, ndiff;
	long long int xorbit;

	xorbit = byte1 ^ byte2;
	ndiff = 0;
	*diffs = 0;
	for (i = 0; i < bl; i++) {
		if (xorbit & 0x01) {
			*diffs = *diffs + (1 << i);
			ndiff = ndiff + 1;
		}
		xorbit = xorbit >> 1;
	}
	
	return ndiff;
}

/* 
 * nsamebytes: compute number of similarities between two bytes
 * 
 * Input:
 *  byte1 = 64 bit byte
 *  byte2 = 64 bit byte
 * Output:
 *  nsame = number of similar bits
 *  sames = location of each similar bit 
 */
int nsamebytes(long long int byte1, long long int byte2, 
	       int bl, long long int *sames)
{
	int i, nsame = 0;
	long long int andbit;
	
	andbit = byte1 & byte2;
	//nsame = 0;
	*sames = 0;
        if (andbit & 0x01) {
                *sames = *sames + (1 << 0);
                nsame = nsame + 1;
        }
        andbit = andbit >> 1;
        for (i = 1; i < bl; i++) {
                if (andbit & 0x01) {
                        *sames = *sames + (1 << i);
                        nsame = nsame + 1;
                }
                andbit = andbit >> 1;
        }
	return nsame;
}

/* 
 * nonzerobits: find nonzero bits of 64 bit byte
 */
void nonzerobits(long long int byt, int bl, int *nzb)
{
	int i;
	
	/* check for nonzero bits, iterate if one is found */
	for (i = 0; i < bl; i++) {
		if (byt & 0x01) {
			*nzb = i+1;
			nzb++;
		}
		byt = byt >> 1;
	}
	return;
}
