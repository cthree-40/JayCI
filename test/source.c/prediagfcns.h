// File: prediagfcns.h
#ifndef prediagfcns_h
#define prediagfcns_h

int drefblock(struct det *detlist,
	      double *moints1,
	      double *moints2,
	      int m1len,
	      int m2len,
	      int aelec,
	      int belec,
	      int refdim,
	      double *evecs);

#endif
