// File: prediagfcns.h
#ifndef prediagfcns_h
#define prediagfcns_h

int drefblock(
	struct det *detlist,
	double *moints1,
	double *moints2,
	int m1len,
	int m2len,
	int aelec,
	int belec,
	int refdim,           /* dimension of reference space */
	double *evecs,        /* output eigenvectors */
	double frzce,         /* frozen core energy */ 
	int ninto           /* internal orbitals  */
	);
	
#endif
