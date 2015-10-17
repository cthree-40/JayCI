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

void readdaiinput(int *maxiter,
		  int *krymin,
		  int *krymax,
		  int *nroots,
		  int *prediagr,
		  int *refdim,
		  double *restol,
		  int *err);

void readgeninput(int *elec,
                  int *orbs,
		  int *nfrozen,
		  int *ndocc,
		  int *nactive,
		  int *xlevel,
		  int *nfrzvirt,
		  int *printlvl,
		  int *err);

int readinputjayci(int *ci_aelec,
		   int *ci_belec,
		   int *ci_orbs,
		   int *nastr,
		   int *nbstr,
		   int *ndets);

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
