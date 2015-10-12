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
int ndiffbytes(long long int byte1, long long int byte2, long long int *diffs)
{
	int i, ndiff;
	long long int xorbit;

	xorbit = byte1 ^ byte2;
	ndiff = 0;
	*diffs = 0;
	for (i = 0; i < 64; i++) {
		if ((xorbit >> i) & 0x01) {
			*diffs = *diffs + pow(2,i);
			ndiff = ndiff + 1;
		}
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
int nsamebytes(long long int byte1, long long int byte2, long long int *sames)
{
	int i, nsame;
	long long int andbit;
	
	andbit = byte1 & byte2;
	nsame = 0;
	*sames = 0;
	for (i = 0; i < 64; i++) {
		if ((andbit >> i) & 0x01) {
			*sames = *sames + pow(2,i);
			nsame = nsame + 1;
		}
	}
	
	return nsame;
}

/* 
 * nonzerobits: find nonzero bits of 64 bit byte
 * 
 * Input:
 *  byt = 64 bit byte
 * Output:
 *  nzb = non zero bits array 
 * Notes: 
 *  This does NOT check the bounds of nzb 
 */
void nonzerobits(long long int byt, int *nzb)
{
	int i;
	
	/* check for nonzero bits, iterate if one is found */
	for (i = 0; i < 64; i++) {
		if ((byt >> i) & 0x01) {
			*nzb = i+1;
			nzb++;
		}
	}
	return;
}
