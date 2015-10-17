// File: binarystr.h
#ifndef binarystr_h
#define binarystr_h

/* occstr: determinant alpha/beta string occupation information */
struct occstr {
	long long int byte1; /* occupation for DOCC+CAS orbitals 1 -> 64  */
	long long int byte2; /* occupation for DOCC+CAS orbitals 65-> 128 */
	int virtx[2];        /* virtual orbital occupations               */
	int nvrtx;           /* number of virtual orbitals                */
};

/* det: determinant composed of alpha and beta occupation strings */
struct det {
	struct occstr astr; /* alpha string                       */
	struct occstr bstr; /* beta  string                       */
	int cas;            /* 1 = no occupations in virtuals, 
		             * 0 = virtual occupations            */
};

/* comparedets_cas: compare two determinants, returning differences */
int comparedets_cas(struct det deti,
		    struct det detj,
		    int *numax,
		    int *numbx,
		    long long int *axi,
		    long long int *axf,
		    long long int *bxi,
		    long long int *bxf);

/* comparedets_virt: compare two determinants with virtual orbital excitations */
int comparedets_ncas(struct det deti,
		     struct det detj,
		     int *numaxc,
		     int *numbxc,
		     int *numaxv,
		     int *numbxv,
		     int *numaxcv,
		     int *numbxcv,
		     long long int *axi,
		     long long int *axf,
		     long long int *bxi,
		     long long int *bxf);

/* str2occstr: convert orbital index string -> occstr type */
struct occstr str2occstr(int *istr, /* orbital index string */
			 int  elec, /* number of electrons  */
			 int ndocc, /* number of doubly-occupied orbitals */
			 int nactv); /* number of active orbitals */


#endif
