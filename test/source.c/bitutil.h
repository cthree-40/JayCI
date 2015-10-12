// File: bitutil.h
#ifndef bitutil_h
#define bitutil_h

int ndiffbytes(long long int byte1,
	       long long int byte2,
	       long long int *diffs);

int nsamebytes(long long int byte1,
	       long long int byte2,
	       long long int *sames);

void nonzerobits(long long int byt,
		 int *nzb);
#endif
