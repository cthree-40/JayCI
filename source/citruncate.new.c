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
 * allocate_xmap: allocate 2d excitation map.
 */
struct xstrmap **allocate_xmap(int xlvl)
{
        struct xstrmap **ptr = NULL;
        int i;
        ptr = (struct xstrmap **) malloc(sizeof(struct xstrmap) * (xlvl + 1));
        for (i = 0; i <= xlvl; i++) {
                ptr[i] = (struct xstrmap *) malloc(sizeof(struct xstrmap) * (xlvl + 1));
        }
        return ptr;
}

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
                init_int_array_0(ptr[i].string, nelecs);
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

        int single_ref = 0; /* Single reference flag */
        
        int ci_orbs;   /* CI orbitals (frozen core/virt removed) */
        int ci_aelec;  /* CI alpha electrons (frozen core removed) */
        int ci_belec;  /* CI beta  electrons (frozen core removed) */

        char detlstfl[FLNMSIZE]; /* Determinant list file */
        char strlstfl[FLNMSIZE]; /* String listS file */
        FILE *dtflptr = NULL;    /* Determinant list file pointer */
        FILE *stflptr = NULL;    /* STring listS file pointer */
        
        struct eostring *pstrings;  /* Alpha electron strings */
        struct eostring *qstrings;  /* Beta  electron strings */

        struct xstrmap **pxmap; /* Alpha string excitations map */
        struct xstrmap **qxmap; /* Beta  string excitations map */
        
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

        if (nactv == 0) {
                single_ref = 1;
                printf("Single reference CI. Setting single_ref = %d\n",
                       single_ref);
        }
        
        /* Generate string lists */
        pxmap = allocate_xmap(xlvl);
        qxmap = allocate_xmap(xlvl);
        if (single_ref == 1) {
                generate_string_list_sr(pstrings, *astr_len, ci_orbs, ci_aelec,
                                        ndocc, nactv, xlvl, pxmap);
        } else {
                generate_string_list(pstrings, *astr_len, ci_orbs, ci_aelec,
                                     ndocc, nactv, xlvl, pxmap);
        }
        printf("Finished generating alpha strings.\n");
        for (int i = 0; i <= 2; i++) {
                for (int j = 0; j <= 2; j++) {
                        printf("%d %d: start=%d finish=%d\n", i, j,
                               pxmap[i][j].start, pxmap[i][j].finish);
                }
        }
        if (single_ref == 1) {
                generate_string_list_sr(qstrings, *bstr_len, ci_orbs, ci_belec,
                                        ndocc, nactv, xlvl, qxmap);
        } else {
                generate_string_list(qstrings, *bstr_len, ci_orbs, ci_belec,
                                     ndocc, nactv, xlvl, qxmap);
        }
        for (int i = 0; i <= 2; i++) {
                for (int j = 0; j <= 2; j++) {
                        printf("%d %d: start=%d finish=%d\n", i, j,
                               qxmap[i][j].start, qxmap[i][j].finish);
                }
        }
        printf("Finished generating beta  strings.\n");

        /* Generate determinant list */
        printf("Generating determinant list...\n");
        *dtrm_len = generate_determinant_list(pstrings, *astr_len,
                                              qstrings, *bstr_len, ci_orbs,
                                              ci_aelec, ci_belec, ndocc, nactv,
                                              xlvl, pxmap, qxmap);
        printf("determinants = %d\n", *dtrm_len);
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
        int i, j, k, l;
        int virt_orbs;
        int *ptr;
        int vstr, astr, dstr; /* VIRT, ACTV, DOCC string number */
        virt_orbs = ci_orbs - ndocc - nactv;
        vstr = binomial_coef2(virt_orbs, virt_elec);
        if (virt_elec == 0) vstr = 0;
        astr = binomial_coef2(nactv, actv_elec);
        if (actv_elec == 0) astr = 0;
        dstr = binomial_coef2(ndocc, docc_elec);
        if (docc_elec == 0) dstr = 0;

        /* Loop over DOCC strings */
        printf("dstr = %d\n", dstr);
        printf("astr = %d\n", astr);
        printf("vstr = %d\n", vstr);
        for (i = 1; i <= dstr; i++) {
                strlist[*pos].doccx = ndocc - docc_elec;
                str_adr2str(i, doccscr, docc_elec, ndocc, docc);
                for (j = 0; j < docc_elec; j++) {
                        strlist[*pos].string[j] = docc[j];
                }
                /* Loop over ACTV strings */
                for (j = 1; j <= astr; j++) {
                        strlist[*pos].actvx = virt_elec;
                        str_adr2str(j, actvscr, actv_elec, nactv, actv);
                        for (k = 0; k < docc_elec; k++) {
                                strlist[*pos].string[k] = docc[k];
                        }
                        for (k = docc_elec; k < docc_elec + actv_elec; k++) {
                                strlist[*pos].string[k] = actv[k - docc_elec] + ndocc;
                        }
                        if (vstr == 0) (*pos)++;
                        /* Loop over VIRT strings */
                        for (k = 1; k <= vstr; k++) {
                                str_adr2str(k, virtscr, virt_elec, virt_orbs, virt);
                                for (l = 0; l < docc_elec; l++) {
                                        strlist[*pos].string[l] = docc[l];
                                }
                                for (l = docc_elec; l < docc_elec + actv_elec; l++) {
                                        strlist[*pos].string[l] = actv[l - docc_elec] + ndocc;
                                }
                                for (l = actv_elec; l < docc_elec + actv_elec + virt_elec; l++) {
                                        strlist[*pos].string[l] =
                                                virt[l - actv_elec] + (ndocc + nactv);
                                }
                                (*pos)++;
                        }
                }
        }

        return;
}

/*
 * compute_eostrings_sr: compute electron occupation strings for all string
 * combinations of DOCC, and VIRT electrons.
 */
void compute_eostrings_sr(struct eostring *strlist, int *pos, int ci_orbs,
                          int ndocc, int docc_elec, int virt_elec, int *docc,
                          int *virt, int *doccscr, int *virtscr)
{
        int i, j, k, l;
        int virt_orbs;
        int *ptr;
        int vstr, dstr; /* VIRT, ACTV, DOCC string number */
        virt_orbs = ci_orbs - ndocc;
        vstr = binomial_coef2(virt_orbs, virt_elec);
        if (virt_elec == 0) vstr = 0;
        dstr = binomial_coef2(ndocc, docc_elec);
        if (docc_elec == 0) dstr = 0;

        /* Loop over DOCC strings */
        printf("dstr = %d\n", dstr);
        printf("vstr = %d\n", vstr);
        for (i = 1; i <= dstr; i++) {
                strlist[*pos].doccx = ndocc - docc_elec;
                str_adr2str(i, doccscr, docc_elec, ndocc, docc);
                for (j = 0; j < docc_elec; j++) {
                        strlist[*pos].string[j] = docc[j];
                }
                if (vstr == 0) (*pos)++;
                /* Loop over VIRT strings */
                for (k = 1; k <= vstr; k++) {
                        str_adr2str(k, virtscr, virt_elec, virt_orbs, virt);
                        for (l = 0; l < docc_elec; l++) {
                                strlist[*pos].string[l] = docc[l];
                        }
                        for (l = docc_elec; l < docc_elec + virt_elec; l++) {
                                strlist[*pos].string[l] =
                                        virt[l - docc_elec] + (ndocc);
                        }
                        (*pos)++;
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
 * generate_determinant_list: generate full *valid* determinant list.
 */
int generate_determinant_list(struct eostring *pstrings, int plen,
                              struct eostring *qstrings, int qlen,
                              int orbs, int aelec, int belec,
                              int ndocc, int nactv, int xlvl,
                              struct xstrmap **pxmap,
                              struct xstrmap **qxmap)
{
        int dlen = 0; /* Number of valid determinants */
        int i, j, k, l;

        /* 0,0 matches with all determinants */
        for (i = pxmap[0][0].start; i < pxmap[0][0].finish; i++) {
                dlen = dlen + qlen;
        }
        printf("%d %d\n", pxmap[0][0].start, pxmap[0][0].finish);
        printf("%d |1,n> determinants\n", dlen);
        
        printf("Single excitations...\n");
        /* 1,0 matches with 0,0; 1,0; and 0,1 */
        for (i = pxmap[1][0].start; i < pxmap[1][0].finish; i++) {
                for (j = qxmap[0][0].start; j < qxmap[0][0].finish; j++) {
                        dlen++;
                }
                if (xlvl == 1) continue; /* only 1 excitation */
                for (j = qxmap[1][0].start; j < qxmap[1][0].finish; j++) {
                        dlen++;
                }
                for (j = qxmap[0][1].start; j < qxmap[0][1].finish; j++) {
                        dlen++;
                }
        }
        /* 0,1 matches with 0,0; 1,0; and 0,1 */
        for (i = pxmap[0][1].start; i < pxmap[0][1].finish; i++) {
                for (j = qxmap[0][0].start; j < qxmap[0][0].finish; j++) {
                        dlen++;
                }
                if (xlvl == 1) continue; /* only 1 excitation */
                for (j = qxmap[1][0].start; j < qxmap[1][0].finish; j++) {
                        dlen++;
                }
                for (j = qxmap[0][1].start; j < qxmap[0][1].finish; j++) {
                        dlen++;
                }
        }
        if (xlvl == 1) return dlen;
        
        printf("Double excitations...\n");
        /* 2,0 matches with 0,0 only */
        for (i = pxmap[2][0].start; i < pxmap[2][0].finish; i++) {
                for (j = qxmap[0][0].start; j < qxmap[0][0].finish; j++) {
                        dlen++;
                }
        }
        for (i = pxmap[0][2].start; i < pxmap[0][2].finish; i++) {
                for (j = qxmap[0][0].start; j < qxmap[0][0].finish; j++) {
                        dlen++;
                }
        }
        /* 1,1 matches with 0,0 only */
        for (i = pxmap[1][1].start; i < pxmap[1][1].finish; i++) {
                for (j = qxmap[0][0].start; j < qxmap[0][0].finish; j++) {
                        dlen++;
                }
        }
        for (i = pxmap[1][1].start; i < pxmap[1][1].finish; i++) {
                for (j = qxmap[0][0].start; j < qxmap[0][0].finish; j++) {
                        dlen++;
                }
        }
                
        return dlen;
}

/*
 * generate_string_list: generate full *valid* alph/beta string lists.
 */
void generate_string_list(struct eostring *strlist, int nstr, int orbs,
                          int elecs, int ndocc, int nactv, int xlvl,
                          struct xstrmap **xmap)
{
        int min_docc_elecs = 0;
        int max_docc_elecs = 0;
        int max_actv_elecs = 0;
        int min_actv_elecs = 0;
        int *docc = NULL, *doccscr = NULL;
        int *actv = NULL, *actvscr = NULL;
        int *virt = NULL, *virtscr = NULL;
        int count = 0;
        int xmi, xmj;
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
        for (i = max_docc_elecs; i >= min_docc_elecs; i--) {
                xmi = max_docc_elecs - i;
                /* Loop over VIRT occupations */
                for (k = 0; k <= xlvl; k++) {
                        xmj = k;
                        xmap[xmi][xmj].start = count;
                        /* Loop over ACTV occupations */
                        for (j = max_actv_elecs; j >= min_actv_elecs; j--) {
                                printf("docc=%d actv=%d virt=%d\n", i, j, k);
                                if ((i + j + k ) != elecs) continue;

                                // COMPUTE STRINGS
                                compute_eostrings(strlist, &count, orbs, ndocc,
                                                  nactv, i, j, k, docc, actv,
                                                  virt, doccscr,actvscr,virtscr);
                                printf("count = %d\n", count);
                        }
                        xmap[xmi][xmj].finish = count;
                }
        }
        /* Set index */
        for (i = 0; i < nstr; i++) {
                strlist[i].index = str_adrfind(strlist[i].string, elecs, orbs);
        }
        return;
}

/*
 * generate_string_list_sr: generate full *valid* alph/beta string lists for
 * single reference CI.
 */
void generate_string_list_sr(struct eostring *strlist, int nstr, int orbs,
                             int elecs, int ndocc, int nactv, int xlvl,
                             struct xstrmap **xmap)
{
        int min_docc_elecs = 0;
        int max_docc_elecs = 0;
        int *docc = NULL, *doccscr = NULL;
        int *virt = NULL, *virtscr = NULL;
        int count = 0;
        int xmi, xmj;
        int i, j, k;

        /* Should be no active orbitals */
        if (nactv != 0) {
                error_message("NACTV != 0! Single reference calc.",
                              "generate_string_list_sr");
                exit(1);
        }
        
        min_docc_elecs = ndocc - int_min(ndocc, xlvl);
        max_docc_elecs = ndocc;
        docc = (int *) malloc(sizeof(int) * ndocc);
        doccscr = (int *) malloc(sizeof(int) * ndocc);
        virt = (int *) malloc(sizeof(int) * xlvl);
        virtscr = (int *) malloc(sizeof(int) * xlvl);
        /* Loop over DOCC occupations */
        for (i = max_docc_elecs; i >= min_docc_elecs; i--) {
                xmi = max_docc_elecs - i;
                /* Loop over VIRT occupations */
                for (k = 0; k <= xlvl; k++) {
                        xmj = k;
                        xmap[xmi][xmj].start = count;
                        if ((i + k) != elecs) continue;
                        compute_eostrings_sr(strlist, &count, orbs, ndocc, i, k,
                                             docc, virt, doccscr, virtscr);
                }
                xmap[xmi][xmj].finish = count;
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
                for (j = 1; j <= binomial_coef2(nactv, i); j++) {
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
                for (j = 1; j <= binomial_coef2(orbs, (elec - i)); j++) {
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
                for (j = 1; j <= binomial_coef2(virt_orbs, i); j++) {
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
                val = val + binomial_coef2(orbs, i);
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
                val = val + binomial_coef2(orbs, (elecs - i));
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
        total = binomial_coef2(ndocc, nde) * binomial_coef2(nactv, nae) *
                binomial_coef2((orbs - (ndocc + nactv)), nve);
        return total;
}
