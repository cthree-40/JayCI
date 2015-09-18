// FILE: ioutil.h
#ifndef ioutil_h_
#define ioutil_h_


void readinput1();
void readmointegrals(double *moints1, double *moints2,
		     int type,
		     int orbitals, unsigned char *moflname,
		     int m1len,
		     int m2len,
		     double *nuc_rep, double *fcenergy);

#endif
