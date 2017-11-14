// File: citruncate_2.c
/*
 * Truncate the CI expansion.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ioutil.h"
#include "iminmax.h"
#include "binary.h"
#include "errorlib.h"
#include "arrayutil.h"
#include "allocate_mem.h"
#include "combinatorial.h"
#include "straddress.h"
#include "citruncate_2.h"

/*
 * allocate_strings_array: allocate the strings array.
 */
struct eostring *allocate_strings_array(int nstr, int nelecs)
{
        struct eostring *ptr = NULL;
        int i;
        ptr = (struct eostring *) malloc(sizeof(struct eostring) * nstr);
        for (i = 0; i < nstr; i++) {
                ptr[i].string = (int *) malloc(sizeof(int) * nelecs);
                ptr[i].index = 0;
        }
        return ptr;
}

/*
 * check_ci_input_values: check the validity of truncation values before
 * execution.
 */
int check_ci_input_values(int cia, int cib, int cio, int nd)
{
        int err = 0;
        if (nd > cib || nd > cia) {
                error_message("Unoccupied DOCC orbitals!",
                              "check_ci_input_values");
                err++;
                return err;
        }
        if (cio == 0) {
                error_message("No virtual orbitals!",
                              "check_ci_input_values");
                err++;
                return err;
        }
        return err;
}

/*
 * citrunc2: truncate the ci expansion.
 */
int citrunc(int aelec, int belec, int orbs, int nfrzc, int ndocc,
	    int nactv, int nfrzv, int xlvl, int *astr_len,
	    int *bstr_len, int *dtrm_len)
{
        int error = 0; /* Error flag */

        int ci_orbs;   /* CI orbitals (frozen core/virt removed) */
        int ci_aelec;  /* CI alpha electrons (frozen core removed) */
        int ci_belec;  /* CI beta  electrons (frozen core removed) */

        char detlstfl[FLNMSIZE]; /* Determinant list file */
        char strlstfl[FLNMSIZE]; /* String listS file */
        FILE *dtflptr = NULL;    /* Determinant list file pointer */
        FILE *stflptr = NULL;    /* STring listS file pointer */

        struct eostring *pstrings; /* Alpha electron strings */
        struct eostring *qstrings; /* Beta  electron strings */
        
        /* Compute electron and orbital numbers in truncated CI space. Then
         * compute the number of valid alpha/beta strings. Allocate valid string
         * arrays */
        ci_aelec = aelec - nfrzc;
        ci_belec = belec - nfrzc;
        ci_orbs  = orbs - nfrzc - nfrzv;
        error = check_ci_input_values(ci_aelec, ci_belec, ci_orbs, ndocc);
        if (error != 0) return error;
        fprintf(stdout, " CI Expansion: %d electrons in %d orbitals\n",
                (ci_aelec + ci_belec), ci_orbs);
        *astr_len = compute_stringnum(ci_orbs, ci_aelec, ndocc, nactv, xlvl);
        *bstr_len = compute_stringnum(ci_orbs, ci_belec, ndocc, nactv, xlvl);
        printf(" Total number of alpha strings = %d\n", *astr_len);
        printf(" Total number of beta  strings = %d\n", *bstr_len);
        pstrings = allocate_strings_array(*astr_len, ci_aelec);
        qstrings = allocate_strings_array(*bstr_len, ci_belec);
        if (pstrings == NULL || qstrings == NULL) {
                error_message("Failed allocating electron string arrays",
                              "citrunc");
                error = 100;
                return error;
        }

        /* Generate string lists */
        generate_string_list(pstrings, *astr_len, ci_orbs, ci_aelec, ndocc,
                             nactv, xlvl);
        generate_string_list(qstrings, *bstr_len, ci_orbs, ci_aelec, ndocc,
                             nactv, xlvl);
        error = 100;
        return error;
}

/*
 * compute_eostrings: compute electron occupation strings for all string
 * combinations of DOCC, ACTV, and VIRT electrons.
 */
void compute_eostrings(struct eostring *strlist, int *pos, int ci_orbs,
                       int ndocc, int nactv, int docc_elec, int actv_elec,
                       int virt_elec, int *docc, int *actv, int *virt,
                       int *doccscr, int *actvscr, int *virtscr)
{
        int i, j, k;
        int virt_orbs;
        int *ptr;
        int vstr, astr, dstr; /* VIRT, ACTV, DOCC string number */
        virt_orbs = ci_orbs - ndocc - nactv;
        vstr = binomial_coef(virt_orbs, virt_elec);
        astr = binomial_coef(nactv, actv_elec);
        dstr = binomial_coef(ndocc, docc_elec);
        *pos = 0;
        /* Loop over DOCC strings */
        for (i = 1; i <= dstr; i++) {
                str_adr2str(i, doccscr, docc_elec, ndocc,
                            &(strlist[*pos].string[0]));
                /* Loop over ACTV strings */
                for (j = 1; j <= astr; j++) {
                        str_adr2str(j, actvscr, actv_elec, nactv,
                                    &(strlist[*pos].string[docc_elec]));
                        if (vstr == 0) (*pos)++;
                        /* Loop over VIRT strings */
                        for (k = 1; k <= vstr; k++) {
                                str_adr2str(k, virtscr, virt_elec, virt_orbs,
                                            &(strlist[*pos].string[(docc_elec + actv_elec)]));
                                (*pos)++;
                        }
                }
        }
        
        return;
}

/*
 * compute_stringnum: compute number of valid strings in expansion.
 */
int compute_stringnum(int orbs, int elecs, int ndocc, int nactv, int xlvl)
{
        int num = 0; /* Number of strings */
        int min_docc_elecs = 0;
        int max_docc_elecs = 0;
        int max_actv_elecs = 0;
        int min_actv_elecs = 0;
        int count = 0;
        int i, j, k;

        min_docc_elecs = ndocc - int_min(ndocc, xlvl);
        max_docc_elecs = ndocc;
        min_actv_elecs = (elecs - ndocc) - int_min((elecs - ndocc), xlvl);
        max_actv_elecs = int_min(((elecs - ndocc) + xlvl), nactv);
        printf("  Min DOCC elecs: %d, Min ACTV elecs: %d, Max ACTV elecs: %d\n",
               min_docc_elecs, min_actv_elecs, max_actv_elecs);
        /* Loop over holes in DOCC space */
        for (i = min_docc_elecs; i <= max_docc_elecs; i++) {
                /* Loop over ACTV space */
                for (j = min_actv_elecs; j <= max_actv_elecs; j++) {
                        /* Loop over virtual space */
                        for (k = 0; k <= xlvl; k++) {
                                num = num + string_number(i, j, k, ndocc, nactv,
                                                          orbs, elecs);
                        }
                }
        }
        
        return num;
}

/*
 * generate_string_list: generate full *valid* alph/beta string lists.
 */
void generate_string_list(struct eostring *strlist, int nstr, int orbs,
                          int elecs, int ndocc, int nactv, int xlvl)
{
        int min_docc_elecs = 0;
        int max_docc_elecs = 0;
        int max_actv_elecs = 0;
        int min_actv_elecs = 0;
        int *docc = NULL, *doccscr = NULL;
        int *actv = NULL, *actvscr = NULL;
        int *virt = NULL, *virtscr = NULL;
        int count = 0;
        int i, j, k;
        min_docc_elecs = ndocc - int_min(ndocc, xlvl);
        max_docc_elecs = ndocc;
        min_actv_elecs = (elecs - ndocc) - int_min((elecs - ndocc), xlvl);
        max_actv_elecs = int_min(((elecs - ndocc) + xlvl), nactv);
        docc = (int *) malloc(sizeof(int) * ndocc);
        doccscr = (int *) malloc(sizeof(int) * ndocc);
        actv = (int *) malloc(sizeof(int) * nactv);
        actvscr = (int *) malloc(sizeof(int) * nactv);
        virt = (int *) malloc(sizeof(int) * xlvl);
        virtscr = (int *) malloc(sizeof(int) * xlvl);
        /* Loop over DOCC occupations */
        for (i = min_docc_elecs; i <= max_docc_elecs; i++) {
                /* Loop over ACTV occupations */
                for (j = min_actv_elecs; j <= max_actv_elecs; j++) {
                        /* Loop over VIRT occupations */
                        for (k = 0; k <= xlvl; k++) {
                                if ((i + j + k ) == elecs) continue;
                                // COMPUTE STRINGS
                                compute_eostrings(strlist, &count, orbs, ndocc,
                                                  nactv, i, j, k, docc, actv,
                                                  virt, doccscr,actvscr,virtscr);
                        }
                }
        }
        /* Set index */
        for (i = 0; i < nstr; i++) {
                strlist[i].index = str_adrfind(strlist[i].string, elecs, orbs);
        }
        return;
}
        
/*
 * generate_actv_strings: generate ACTV orbital strings array.
 */
int generate_actv_strings(int **qactv_array, int nstr, int nactv, int elecs,
                          int ndocc, int xlvl)
{
        int error = 0;
        int maxelecs;
        int minelecs;
        int *scr = NULL;
        int cnt = 0;
        int i, j;
        minelecs = (elecs - ndocc) - int_min((elecs - ndocc), xlvl);
        maxelecs = int_min(((elecs - ndocc) + xlvl), nactv);
        scr = (int *) malloc(sizeof(int) * maxelecs);
        for (i = minelecs; i <= maxelecs; i++) {
                for (j = 1; j <= binomial_coef(nactv, i); j++) {
                        str_adr2str(j, scr, i, nactv, qactv_array[cnt]);
                        cnt++;
                }
        }
        free(scr);
        return error;
}

/*
 * generate_docc_strings: generate DOCC orbital strings array.
 */
int generate_docc_strings(int **qdocc_array, int nstr, int orbs, int elec,
                          int xlvl)
{
        int error = 0; /* Error flag */
        int maxholes = 0;
        int *scr = NULL;
        int cnt = 0;
        int i, j;
        maxholes = int_min(orbs, xlvl);
        scr = (int *) malloc(sizeof(int) * elec);
        for (i = 0; i <= maxholes; i++) {
                for (j = 1; j <= binomial_coef(orbs, (elec - i)); j++) {
                        str_adr2str(j, scr, (elec - i), orbs, qdocc_array[cnt]);
                        cnt++;
                }
        }
        free(scr);
        return error;
}

/*
 * generate_virt_strings: generate VIRT orbital strings array.
 */
int generate_virt_strings(int **virt_array, int nstr, int orbs, int elec,
                          int ndocc, int nactv, int xlvl)
{
        int error = 0;
        int max_elecs = 0;
        int virt_orbs = 0;
        int *scr = NULL;
        int cnt = 0;
        int i, j;
        max_elecs = xlvl;
        virt_orbs = orbs - ndocc - nactv;
        scr = (int *) malloc(sizeof(int) * xlvl);
        for (i = 0; i <= xlvl; i++) {
                for (j = 1; j <= binomial_coef(virt_orbs, i); j++) {
                        str_adr2str(j, scr, i, virt_orbs, virt_array[cnt]);
                        cnt++;
                }
        }
        free(scr);
        return error;
}

/*
 * max_actv_space_strings: return the maximum number of ACTV space strings.
 */
int max_actv_space_strings(int orbs, int elecs, int ndocc, int xlvl)
{
        int val = 0;
        int minelecs;
        int maxelecs;
        int i;
        minelecs = (elecs - ndocc) - int_min((elecs - ndocc), xlvl);
        maxelecs = int_min(((elecs - ndocc) + xlvl), orbs);
        for (i = minelecs; i <= maxelecs; i++) {
                val = val + binomial_coef(orbs, i);
        }
        return val;
}
        

/*
 * max_docc_space_strings: return the maximum number of DOCC space strings.
 */
int max_docc_space_strings(int orbs, int elecs, int xlvl)
{
        int val = 0;
        int maxholes;
        int i;
        maxholes = int_min(orbs, xlvl);
        for (i = 0; i <= maxholes; i++) {
                val = val + binomial_coef(orbs, (elecs - i));
        }
        return val;
}

/*
 * string_number: generate number of possible *valid* strings given the number
 * of DOCC electrons, ACTV electrons, and VIRT electrons.
 */
int string_number(int nde, int nae, int nve, int ndocc, int nactv, int orbs,
                  int elecs)
{
        int total = 0;
        if ((nde + nae + nve) != elecs) return total;
        total = binomial_coef(ndocc, nde) * binomial_coef(nactv, nae) *
                binomial_coef((orbs - (ndocc + nactv)), nve);
        return total;
}
