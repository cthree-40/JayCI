// File: bitutil.h
#ifndef bitutil_h
#define bitutil_h

int ndiffbytes(long long int byte1,
	       long long int byte2,
	       int bl,                /* bit length */
	       long long int *diffs);

int nsamebytes(long long int byte1,
	       long long int byte2,
	       int bl,                 /* bit length */
	       long long int *sames);

/*
 * nonzerobits: find nonzero bits of 64 bit byte
 */
void nonzerobits(
	long long int byt,  /* 64 bit byte */
	int bl,             /* bit length */
	int *nzb);          /* numer of nonzero bits */

#endif
