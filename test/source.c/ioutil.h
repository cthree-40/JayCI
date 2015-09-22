// FILE: ioutil.h
#ifndef ioutil_h_
#define ioutil_h_

#define FLNMSIZE           30
#define MAX_LINE_SIZE     300
#define MAX_NAMELIST_SIZE  30

int checkinputfiles();

int geninput(int dlen,
	     int alen,
	     int blen,
	     int aelec,
	     int belec,
	     int orbs,
	     int nfrzc,
	     int nfrzv);

void readgeninput(int *elec,
                  int *orbs,
		  int *nfrozen,
		  int *ndocc,
		  int *nactive,
		  int *xlevel,
		  int *nfrzvirt,
		  int *printlvl,
		  int *err);

void readmointegrals(double *moints1,
		     double *moints2,
		     int type,
		     int orbitals,
		     unsigned char *moflname,
		     int m1len,
		     int m2len,
		     double *nuc_rep,
		     double *fcenergy);

#endif
