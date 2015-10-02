// File: action_util.h
#ifndef action_util_h
#define action_util_h

double hmatels(struct det deti,
	       struct det detj,
	       double *moints1,
	       double *moints2,
	       int m1len,
	       int m2len);

double evaluate_dets(int ndiff,
		     struct det deti,
		     struct det detj,
		     int numax,
		     int numbx,
		     long long int axi,
		     long long int axf,
		     long long int bxi,
		     long long int bxf);
#endif
