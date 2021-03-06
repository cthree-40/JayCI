// FILE: citruncate.c
/*********************************************************************
 * citruncate: Subfunctions for truncating ci expansion.
 * -----------------------------------------------------
 *  citrunc: performs truncation
 *  str_enfdocc: enforces DOCC orbital restrictions strings
 *  str_enfactv: enforces ACTIVE orbital restrictions on strings
 *  
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ioutil.h"
#include "binary.h"
#include "citruncate.h"
#include "combinatorial.h"
#include "straddress.h"

struct validstr {
	int index;
	int xdocc;
	int xactv;
	struct validstr *next;
};

struct validdet {
        int pindx;
        int qindx;
        struct validdet *next;
};


/* citrunc: truncate ci expansion
 * -------------------------------------------------------------------
 * Input:
 *  aelec = alpha electrons
 *  belec = beta  electrons
 *  orbs  = orbitals
 *  nfrzc = number of frozen core orbitals
 *  ndocc = number of doubly-occupied orbitals
 *  nactv = number of active orbitals
 *  nfrzv = number of frozen virtual orbitals
 *  xlvl  = excitation level
 * Output:
 *  astr_len = number of alpha strings
 *  bstr_len = number of beta  strings
 *  dtrm_len = number of determinants */
int citrunc(int aelec, int belec, int orbs, int nfrzc, int ndocc,
	    int nactv, int nfrzv, int xlvl, int *astr_len,
	    int *bstr_len, int *dtrm_len)
{
	/* output scalar 
	 * err = error handling */
	int err;
	
	/* local scalars 
	 * ci_orbs  = orbitals - (frozen core) - (frozen virtual)
	 * ci_aelec = aelec - (frozen core)
	 * ci_belec = belec - (frozen core)
	 * axdocc   = alpha string docc   orbital excitations
	 * axactv   = alpha string active orbital excitations
	 * bxdocc   = beta  string docc   orbital excitations
	 * bxactv   = beta  string active orbital excitations
	 * astrings = binomial_coef(ci_orbs, ci_aelec)
	 * bstrings = binomial_coef(ci_orbs, ci_belec)
	 * tot_dets = astrings * bstrings 
	 * aindx    = alpha string index
	 * bindx    = beta  string index */
	int  ci_orbs, ci_aelec, ci_belec;
	int   axdocc,   axactv,   bxdocc,  bxactv;
	int astrings, bstrings, tot_dets;
	int aindx, bindx;
	
	/* local arrays 
	 * astr1      = alpha string 1
	 * astr2      = alpha string 2
	 * bstr1      = beta  string 1
	 * bstr2      = beta  string 2 
	 * bstr_hf    = beta  string start (Hartree-Fock)
	 * astr_hf    = alpha string start (Hartree-Fock)
	 * detlstfl   = determinant list file name
	 * strlstfl   = string list file name  */
	int *astr1, *astr2, *bstr1, *bstr2, *bstr_hf, *astr_hf;
	char detlstfl[FLNMSIZE], strlstfl[FLNMSIZE];
	
	/* local linked-list
	 * qindxlist = valid beta strings
         * pindxlist = valid alpha strings
         * detlist   = valid determinants */
	struct validstr *qindxlist;
	struct validstr *qptr, *tmp;
        struct validstr *pindxlist;
        struct validstr *pptr;
        struct validdet *detlist;
        struct validdet *dptr;
	
	/* local files
	 * detfileptr = determinant list file pointer
	 * strfileptr = string list file pointer */
	FILE *detfileptr, *strfileptr;
	
	int i, j;
	
	/* initialize error flag */
	err = 0;
	/* set determinant and string list file names */
	strncpy(detlstfl, "det.list", FLNMSIZE);
	strncpy(strlstfl, "str.list", FLNMSIZE);
	
	/* compute electron numbers in truncated space */
	ci_aelec = aelec - nfrzc;
	ci_belec = belec - nfrzc;
	ci_orbs  = orbs - nfrzc - nfrzv;
	
	fprintf(stdout, " CI Expansion: %d %d %d\n", ci_aelec, ci_belec, ci_orbs);
	
	/* check input */
	/* ndocc values */
	if (ndocc > ci_belec) {
		fprintf(stderr,
			"*** ERROR: Unoccupied DOCC orbitals! ***\n");
		err = -100;
		return err;
	}
	if (ci_orbs == 0 ) {
		fprintf(stderr,
			"WARNING: No virtual orbitals!\n");
		err = 10;
		return err;
	}
	
	/* compute total, untruncated - excepting frozen core truncation -
	 * number of alpha/beta/determinant strings */
	astrings = binomial_coef(ci_orbs, ci_aelec);
	bstrings = binomial_coef(ci_orbs, ci_belec);
	tot_dets = astrings * bstrings;
	fprintf(stdout, "Untruncated expansion size = %10d\n", tot_dets);
	
	/* allocate arrays */
	astr1 = (int *) malloc(ci_aelec * sizeof(int));
	astr2 = (int *) malloc(ci_aelec * sizeof(int));
	bstr1 = (int *) malloc(ci_belec * sizeof(int));
	bstr2 = (int *) malloc(ci_belec * sizeof(int));
	astr_hf = (int *) malloc(ci_aelec * sizeof(int));
	bstr_hf = (int *) malloc(ci_belec * sizeof(int));
	
	/* form HF strings */
	for (i = 0; i < ci_aelec; i++) {
		astr_hf[i] = i + 1;
		astr1[i] = i + 1;
	}
	for (i = 0; i < ci_belec; i++) {
		bstr_hf[i] = i + 1;
		bstr1[i] = i + 1;
	}
	
	/* open determinant list file */
	detfileptr = fopen("det.list", "w");
	if (detfileptr == NULL) {
		fprintf(stderr,
			"*** ERROR: det.list could not be opened! ***\n");
		exit(1);
	}
	/* open string list file */
	strfileptr = fopen("str.list", "w");
	if (strfileptr == NULL) {
		fprintf(stderr,
			"*** ERROR: str.list could not be opened! ***\n");
		exit(1);
	}
	
	/* first strings are included */
	aindx = str_adrfind(astr1, ci_aelec, ci_orbs);
	bindx = str_adrfind(bstr1, ci_belec, ci_orbs);
	fprintf(detfileptr, " %14d %14d\n", aindx, bindx);
	fprintf(strfileptr, "BETA\n");
	fprintf(strfileptr, " %14d\n", bindx);
	
	*astr_len = 1;
	*bstr_len = 1;
	*dtrm_len = 1;
	
	/* allocate first node of valid q string linked list */
	qindxlist = malloc(sizeof(struct validstr));
	qindxlist->next = NULL;
	qptr = qindxlist;
	
	qptr->index = bindx;
	qptr->xdocc = 0;
	qptr->xactv = 0;
	
	qptr->next = malloc(sizeof(struct validstr));
	qptr->next->next = NULL;
	qptr = qptr->next;
        
	/* allocate first node of valid p string linked list */
	pindxlist = malloc(sizeof(struct validstr));
	pindxlist->next = NULL;
	pptr = pindxlist;
	
	pptr->index = aindx;
	pptr->xdocc = 0;
	pptr->xactv = 0;
	
	pptr->next = malloc(sizeof(struct validstr));
	pptr->next->next = NULL;
	pptr = pptr->next;

        /* allocate first node of valid determinant list */
        detlist = malloc(sizeof(struct validdet));
        detlist->next = NULL;
        dptr = detlist;
        dptr->pindx = 1;
        dptr->qindx = 1;
        dptr->next = malloc(sizeof(struct validdet));
        dptr->next->next = NULL;
        dptr = dptr->next;
        
	/* first loop: |p,q> for p=1; q=1,..,bstrings */
	for (i = 2; i <= bstrings; i++) {
		/* generate beta electron string */
		str_strfind1(bstr1, ci_belec, ci_orbs, bstr2);
		
		/* test docc restrictions */
		bxdocc = str_enfdocc(bstr2, ci_belec, ndocc, nactv);
		if (bxdocc > xlvl) {
			for (j = 0; j < ci_belec; j++) {
				bstr1[j] = bstr2[j];
			}
			continue;
		}
		
		/* test active orbital restrictions */
		bxactv = str_enfactv(bstr2, ci_belec, ndocc, nactv);
		if (bxactv > xlvl) {
			for (j = 0; j < ci_belec; j++) {
				bstr1[j] = bstr2[j];
			}
			continue;
		}
		
		/* if here, string is valid */
		*bstr_len = *bstr_len + 1;
		//fprintf(strfileptr, " %14d\n", i);
		
		*dtrm_len = *dtrm_len + 1;
                dptr->pindx = 1;
                dptr->qindx = i;
                //fprintf(detfileptr, " %14d %14d\n", 1, i);
		dptr->next = malloc(sizeof(struct validdet));
                dptr->next->next = NULL;
                dptr = dptr->next;
                
		qptr->index = i;
		qptr->xdocc = bxdocc;
		qptr->xactv = bxactv;
		
		qptr->next = malloc(sizeof(struct validstr));
		qptr->next->next = NULL;
		qptr = qptr->next;
		
		/* increment string */
		for (j = 0; j < ci_belec; j++) {
			bstr1[j] = bstr2[j];
		}
		
	}
	
	/* second loop: |p,q> for p=1,..,astrings; q=1,..,bstrings */
	fprintf(strfileptr, "ALPHA\n");
	fprintf(strfileptr, " %14d\n", aindx);
	for (i = 2; i <= astrings; i++) {
		
		/* generate alpha electron string */
		str_strfind1(astr1, ci_aelec, ci_orbs, astr2);
		
		/* test docc restrictions */
		axdocc = str_enfdocc(astr2, ci_aelec, ndocc, nactv);
		if (axdocc > xlvl) {
			for (j = 0; j < ci_aelec; j++) {
				astr1[j] = astr2[j];
			}
			continue;
		}
		
		/* test active orbital restrictions */
		axactv = str_enfactv(astr2, ci_aelec, ndocc, nactv);
		if (axactv > xlvl) {
			for (j = 0; j < ci_aelec; j++) {
				astr1[j] = astr2[j];
			}
			continue;
		}
		
		/* if here, |p, 1> is valid (p = i) */
		//fprintf(strfileptr, " %14d\n", i);
                pptr->index = i;
                pptr->xdocc = axdocc;
                pptr->xactv = axactv;
                pptr->next = malloc(sizeof(struct validstr));
                pptr->next->next = NULL;
                pptr = pptr->next;
                
                *astr_len = *astr_len + 1;
		
		for (j = 0; j < ci_belec; j++) {
			bstr1[j] = bstr_hf[j];
		}
		bindx = 1;
		
		/* now loop over valid beta strings */
		qptr = qindxlist;
		while (qptr->next != NULL) {
			
			/* generate orbital string */
			str_strfind2(bstr1, bindx, ci_belec, ci_orbs, qptr->index,
				     bstr2);
			bindx = qptr->index;
			for (j = 0; j < ci_belec; j++) {
				bstr1[j] = bstr2[j];
			}
			
			/* test docc restrictions */
			if ((axdocc + qptr->xdocc) > xlvl) {
				for (j = 0; j < ci_belec; j++) {
					bstr1[j] = bstr2[j];
				}
				qptr = qptr->next;
				continue;
			}
			
			/* test active orbital restrictions */
			if ((axactv + qptr->xactv) > xlvl) {
				for (j = 0; j < ci_belec; j++) {
					bstr1[j] = bstr2[j];
				}
				qptr = qptr->next;
				continue;
			}
			
			/* test together */
			if (((axactv + qptr->xactv) + (axdocc + qptr->xdocc))
			    > xlvl) {
				for (j = 0; j < ci_belec; j++) {
					bstr1[j] = bstr2[j];
				}
				qptr = qptr->next;
				continue;
			}
			
			/* if here, write determinant */
			//fprintf(detfileptr, " %14d %14d\n", i, qptr->index);
                        dptr->pindx = i;
                        dptr->qindx = qptr->index;
                        dptr->next = malloc(sizeof(struct validdet));
                        dptr->next->next = NULL;
                        dptr = dptr->next;
                        
                        *dtrm_len = *dtrm_len + 1;
			
			for (j = 0; j < ci_belec; j++) {
				bstr1[j] = bstr2[j];
			}
			qptr = qptr->next;
		}
		
		for (j = 0; j < ci_aelec; j++) {
			astr1[j] = astr2[j];
		}
	}
	
	fprintf(stdout, "Deallocating electron arrays.\n");
	free(astr1);
	free(astr2);
	free(bstr1);
	free(bstr2);
	free(bstr_hf);
	free(astr_hf);
	
	/* free linked list */
	qptr = qindxlist;
	fprintf(stdout, "Deallocating linked list.\n");
	while (qptr->next != NULL) {
		tmp = qptr;
		qptr = qptr->next;
		
		free(tmp);
	}
	
	fclose(detfileptr);
	fclose(strfileptr);
	return err;
}
/* str_enfactv: enforce ACTIVE orbital restrictions
 * -------------------------------------------------------------------
 * Input:
 *  str   = orbital index string
 *  elec  = electrons
 *  ndocc = number of docc orbitals
 *  nactv = number of active orbitals */
int str_enfactv(int *str, int elec, int ndocc, int nactv)
{
	int i;
	int x;
	
	x = 0;
	
	if (nactv == 0) return x;
	
	for (i = 0; i < elec; i++ ) {
		if (str[i] > (ndocc + nactv))
			x = x + 1;
	}
	
	return x;
	
}
/* str_enfdocc: enforce DOCC orbital restrictions
 * -------------------------------------------------------------------
 * Input:
 *  str   = orbital index string
 *  elec  = electrons
 *  ndocc = number of docc orbitals
 *  nactv = number of active orbitals */
int str_enfdocc(int *str, int elec, int ndocc, int nactv)
{
	int i;
	int x;
	
	x = 0;
	
	if (ndocc == 0) return x;
	
	for (i = 0; i < ndocc; i++) {
		if (str[i] > ndocc)
			x = x + 1;
	}
	
	return x;
}


