// File: binarystr.h
#ifndef binarystr_h
#define binarystr_h

/* occstr: determinant alpha/beta string occupation information */
struct occstr {
     long long int byte1; /* occupation for DOCC+CAS orbitals 1 -> 64 */
     long long int byte2; /* occupation for DOCC+CAS orbitals 65-> 128*/
     int virtx[2]; /* virtual orbital excitations */
};

/* str2occstr: convert orbital index string -> occstr type */
struct occstr str2occstr(int *istr, /* orbital index string */
			 int  elec, /* number of electrons  */
			 int ndocc, /* number of doubly-occupied orbitals */
			 int nactv); /* number of active orbitals */


#endif
