// FILE: ioutil.h
#ifndef ioutil_h_
#define ioutil_h_

#define FLNMSIZE           30
#define MAX_LINE_SIZE     300
#define MAX_NAMELIST_SIZE  30

/* check_for_file
 * --------------
 * Check if file exists and can be opened for desired operation.
 */
int check_for_file(
        char *filename,     /* name of file */
        char *fileoperation /* mode */
        );

int checkinputfiles();

int find_str_count_in_file(
        char *string,
        FILE *fptr
        );

FILE *find_str_line(
        char *string,
        FILE *fptr
        );

int geninput(int dlen,
	     int alen,
	     int blen,
	     int aelec,
	     int belec,
	     int orbs,
	     int nfrzc,
	     int nfrzv);
/*
 * print_array_2d: print 2d array
 */
void print_array_2d(
    double **array,
    int rows,
    int cols
    );

/*
 * print_civectors: print out ci vectors
 */
int print_civectors(
        int aelec, /* alpha electrons */
        int belec, /* beta  electrons */
        int orbs, /* number of molecular orbitals */
	double **civec,
	int ndets,
	int roots,
	double *cival
	);

/*
 * read_civectors: read in ci information and vectors
 */
int read_civectors(
        int aelec,
        int belec,
        int orbs,
        double **civec,
        int ndets,
        int roots,
        double *cival
        );

void readdaiinput(int *maxiter,
		  int *krymin,
		  int *krymax,
		  int *nroots,
		  int *prediagr,
		  int *refdim,
		  double *restol,
		  int *err);

void read_dysonorb_input(
        char *wvfcn_file0,
        char *wvfcn_file1,
        int *nstates0,
        int *nstates1,
        int *nelecs0,
        int *nelecs1,
        int *norbs0,
        int *norbs1,
        int *ndets0,
        int *ndets1,
        int *ninto0,
        int *ninto1,
        int *err);

void readgeninput(int *elec,
                  int *orbs,
		  int *nfrozen,
		  int *ndocc,
		  int *nactive,
		  int *xlevel,
		  int *nfrzvirt,
		  int *printlvl,
                  int *printwvf,
                  int *err);

int readinputjayci(int *ci_aelec,
		   int *ci_belec,
		   int *ci_orbs,
		   int *nastr,
		   int *nbstr,
		   int *ndets);

/* readmocoeffs: subroutine to read molecular coefficient file.
 * -------------------------------------------------------------------
 * Calls fortran subroutine readmocoef()
 *
 *
 */
void readmocoeffs(double *c, int clen);

void readmointegrals(double *moints1,
		     double *moints2,
		     int type,
		     int orbitals,
		     char *moflname,
		     int m1len,
		     int m2len,
		     double *nuc_rep,
		     double *fcenergy);

/*
 * substring: gets substring from string and returns pointer to said
 * substring.
 */
char *substring(
        char *string,
        int position,
        int length);

#endif
