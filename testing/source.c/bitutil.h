// File: bitutil.h
#ifndef bitutil_h
#define bitutil_h

/*
 * ndiffbytes: compute number of differences between two bytes
 */
int ndiffbytes(
        long long int           byte1,
        long long int           byte2,
        int           length_of_bytes,
        long long int       *xor_byte
        );

/*
 * nsamebytes: compute number of similarities between two bytes
 */
int nsamebytes(
        long long int           byte1,
        long long int           byte2,
        int           length_of_bytes,                 /* bit length */
        long long int       *and_byte
        );

/*
 * nonzerobits: find nonzero bits of 64 bit byte
 */
void nonzerobits(
	long long int byt,  /* 64 bit byte */
	int bl,             /* bit length */
	int *nzb);          /* numer of nonzero bits */

#endif
