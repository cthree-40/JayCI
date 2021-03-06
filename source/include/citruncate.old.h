// FILE: citruncate.h
#ifndef citruncate_h
#define citruncate_h

int citrunc(int aelec,
	    int belec,
	    int orbs,
	    int nfrzc,
	    int ndocc,
	    int nactv,
	    int nfrzv,
	    int xlvl,
	    int *astr_len,
	    int *bstr_len,
	    int *dtrm_len);

int str_enfactv(int *str,
		int elec,
		int ndocc,
	        int nactv);

int str_enfdocc(int *str,
		int elec,
		int ndocc,
		int nactv);

#endif

