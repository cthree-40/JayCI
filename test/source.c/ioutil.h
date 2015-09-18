// FILE: ioutil.h
#ifndef ioutil_h_
#define ioutil_h_


void readgeninput(int *elec,
                  int *orbs,
		  int *nfrozen,
		  int *ndocc,
		  int *nactive,
		  int *xlevel,
		  int *nfrzvirt,
		  int *printlvl,
		  int *err
 );
void readmointegrals(double *moints1,
		     double *moints2,
		     int type,
		     int orbitals,
		     unsigned char *moflname,
		     int m1len,
		     int m2len,
		     double *nuc_rep,
		     double *fcenergy);

/* file name size */
#define FLNMSIZE 30

/* max line size */
#define MAX_LINE_SIZE 300

/* max namelist size */
#define MAX_NAMELIST_SIZE 30

#endif
