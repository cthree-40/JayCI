// File: citruncate.c
/*
 * Truncate a Full CI expansion.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pjayci_global.h"
#include "mpi_utilities.h"
#include "ioutil.h"
#include "iminmax.h"
#include "binary.h"
#include "bitutil.h"
#include "errorlib.h"
#include "arrayutil.h"
#include "allocate_mem.h"
#include "combinatorial.h"
#include "straddress.h"
#include "binarystr.h"
#include "action_util.h"
#include "citruncate.h"
#include <mpi.h>
#include <ga.h>
#include <macdecls.h>

/*--------------------------------------------------------------------------*/
/* -- MAIN DRIVER --                                                        */
/*--------------------------------------------------------------------------*/

/*
 * citrunc: truncate the CI expansion.
 * Input:
 *  aelec = alpha electrons
 *  belec = beta  electrons
 *  orbs  = molecular orbitals
 *  nfrzc = number of frozen core orbitals
 *  ndocc = number of doubly-occupied reference orbitals
 *  nactv = number of active reference orbitals
 *  nfrzv = number of frozen virtual orbitals
 *  xlvl  = excitation level (DOCC->ACTV & (DOCC + ACTV)->VIRT)
 *  astrings = alpha electron strings
 *  astr_len = alpha string number
 *  bstrings = beta  electron strings
 *  bstr_len = beta  string number
 *  peosp    = alpha electron space array
 *  pegrps   = alpha electron groupings
 *  qeosp    = beta electron space array
 *  qegrps   = beta electron groupings
 *  astr_len = alpha string number
 *  bstr_len = beta  string number
 *  dtrm_len = determinant number
 *  pq_spaces = pq pairings
 *  num_pq   = number of valid pq pairings
 */
int citrunc(int aelec, int belec, int orbs, int nfrzc, int ndocc,
	    int nactv, int nfrzv, int xlvl, struct occstr *astrings,
            int astr_len, struct occstr *bstrings, int bstr_len,
            struct eospace *peosp, int pegrps,
            struct eospace *qeosp, int qegrps, int *dtrm_len,
            int **pq_spaces, int *num_pq)
{
        int error = 0; /* Error flag */
        int ci_orbs;   /* CI orbitals (frozen core/virt removed) */
        int ci_aelec;  /* CI alpha electrons (frozen core removed) */
        int ci_belec;  /* CI beta  electrons (frozen core removed) */

        struct eostring *pstrings;  /* Alpha electron strings */
        struct eostring *qstrings;  /* Beta  electron strings */

        /* Compute electron and orbital number in CI expansion space.
         * Compute alpha and beta string numbers */
        ci_aelec = aelec - nfrzc;
        ci_belec = belec - nfrzc;
        ci_orbs  = orbs - nfrzc - nfrzv;

        /* Allocate the electron string arrays. */
        pstrings = allocate_strings_array(astr_len, ci_aelec);
        qstrings = allocate_strings_array(bstr_len, ci_belec);
        if (pstrings == NULL || qstrings == NULL ||
            astrings == NULL || bstrings == NULL) {
                error_message(mpi_proc_rank,
                              "Failed allocating electron string arrays",
                              "citrunc");
                error = 100;
                return error;
        }

        /* Generate alpha/beta string lists */
        generate_string_list(pstrings, astr_len, ci_orbs, ci_aelec, ndocc,
                             nactv, xlvl, peosp, pegrps);
        generate_binstring_list(pstrings, astr_len, ci_aelec, ndocc, nactv,
                                astrings);
        generate_string_list(qstrings, bstr_len, ci_orbs, ci_belec, ndocc,
                             nactv, xlvl, qeosp, qegrps);
        generate_binstring_list(qstrings, bstr_len, ci_belec, ndocc, nactv,
                                bstrings);

        /* Compute number of determinants */
        *dtrm_len = compute_detnum(peosp, pegrps, qeosp, qegrps, ndocc,
                                   nactv, xlvl, pq_spaces, num_pq);
        
        deallocate_eostrings_array(pstrings, astr_len);
        deallocate_eostrings_array(qstrings, bstr_len);
        return error;
}

/* -------------------------------------------------------------------------- */
/* -- SUBROUTINES --                                                          */
/* -------------------------------------------------------------------------- */

/*
 * allocate_eospace_array: allocate the electron number space array.
 */
struct eospace *allocate_eospace_array(int nelec, int norbs, int ndocc,
                                       int nactv, int xlvl,  int *ngrp)
{
        struct eospace *ptr = NULL;
        int min_docc_elecs = 0;
        int max_docc_elecs = 0;
        int min_actv_elecs = 0;
        int max_actv_elecs = 0;
        int vorbs = 0;
        int cnt = 0;                 /* Counter */
        int i, j, k;

        vorbs = norbs - ndocc - nactv;
        min_docc_elecs = ndocc - int_min(ndocc, xlvl);
        max_docc_elecs = ndocc;
        min_actv_elecs = (nelec - ndocc) - int_min((nelec - ndocc), xlvl);
        //max_actv_elecs = int_min(((nelec - ndocc) + xlvl), nactv);
        max_actv_elecs = int_min(int_min(ndocc, xlvl) + (nelec - ndocc), nactv);
        /* Loop over DOCC, ACTV, and VIRT occupations. */
        *ngrp = 0;
        for (i = max_docc_elecs; i >= min_docc_elecs; i--) {
                for (j = max_actv_elecs; j >= min_actv_elecs; j--) {
                        for (k = 0; k <= int_min(vorbs, xlvl); k++) {
                                if ((i + j + k) != nelec) continue;
                                (*ngrp)++;
                        }
                }
        }

        /* Allocate groupings list */
        ptr = (struct eospace *) malloc(sizeof(struct eospace) * (*ngrp));

        /* Loop over DOCC, ACTV, and VIRT occupations. */
        cnt = 0;
        for (i = max_docc_elecs; i >= min_docc_elecs; i--) {
                for (j = max_actv_elecs; j >= min_actv_elecs; j--) {
                        for (k = 0; k <= int_min(vorbs, xlvl); k++) {
                                if ((i + j + k) != nelec) continue;
                                ptr[cnt].docc = i;
                                ptr[cnt].actv = j;
                                ptr[cnt].virt = k;
                                cnt++;
                        }
                }
        }
        return ptr;
}

/*
 * allocate_occstr_arrays: allocate binary electron orbital occupation
 * string arrays.
 */
struct occstr *allocate_occstr_arrays(int nstr)
{
        struct occstr *ptr = NULL;
        ptr = (struct occstr *) malloc(sizeof(struct occstr) * nstr);
        for (int i = 0; i < nstr; i++) {
                ptr[i].byte1 = 0x0;
                ptr[i].virtx[0] = 0;
                ptr[i].virtx[1] = 0;
                ptr[i].nvrtx = 0;
        }
        return ptr;
}                

/*
 * allocate_strings_array: allocate the strings array.
 */
struct eostring *allocate_strings_array(int nstr, int nelecs)
{
        struct eostring *ptr = NULL;
        int i, j;
        ptr = (struct eostring *) malloc(sizeof(struct eostring) * nstr);
        for (i = 0; i < nstr; i++) {
                ptr[i].string = (int *) malloc(sizeof(int) * nelecs);
                //init_int_array_0(ptr[i].string, nelecs);
                for (j = 0; j < nelecs; j++) {
                        ptr[i].string[j] = 0;
                }
                ptr[i].index = 0;
                ptr[i].doccx = 0;
                ptr[i].actvx = 0;
        }
        return ptr;
}

/*
 * allocate_xmap: allocate 2d excitation map.
 */
struct xstrmap **allocate_xmap(int xlvl)
{
        struct xstrmap **ptr = NULL;
        int i;
        ptr = (struct xstrmap **) malloc(sizeof(struct xstrmap) * (xlvl + 1));
        for (i = 0; i <= xlvl; i++) {
                ptr[i] = (struct xstrmap *)
                        malloc(sizeof(struct xstrmap) * (xlvl + 1));
        }
        return ptr;
}

/*
 * construct_xlist: construct list of X excitations for electron strings
 */
int **construct_xlist(struct occstr *strlist, int nstr, int intorb, int x, int *max)
{
        int **ptr = NULL;
        int *data = NULL;
        int cnt;
        int numx;
        int numxv, numxcv, numxc;
        int samei, samej;
        long long int diffs, axi, axj;
        int i, j;
        *max = 0;
        for (i = 0; i < nstr; i++) {
                cnt = 0;
                for (j = 0; j < nstr; j++) {
                        numxv = compute_virt_diffs(strlist[i],strlist[j]);
                        numxcv= abs(strlist[i].nvrtx - strlist[j].nvrtx);
                        numxc = ndiffbytes(strlist[i].byte1,strlist[j].byte1,
                                           intorb, &diffs);
                        samei = nsamebytes(strlist[i].byte1,diffs,intorb,&axi);
                        samej = nsamebytes(strlist[j].byte1,diffs,intorb,&axj);
                        numxc = int_min(samei,samej);
                        numxv = numxv - numxcv;
                        numx  = numxc + numxcv + numxv;
                        if (numx != x) continue;
                        cnt++;
                }
                if (cnt > *max) *max = cnt;
        }
        *max = *max + 1; /* First slot is for how many strings are in row */
        data = allocate_mem_int_cont(&ptr, *max, nstr);
        for (i = 0; i < nstr; i++) {
                cnt = 1;
                for (j = 0; j < nstr; j++) {
                        numxv = compute_virt_diffs(strlist[i],strlist[j]);
                        numxcv= abs(strlist[i].nvrtx - strlist[j].nvrtx);
                        numxc = ndiffbytes(strlist[i].byte1,strlist[j].byte1,
                                           intorb, &diffs);
                        samei = nsamebytes(strlist[i].byte1,diffs,intorb,&axi);
                        samej = nsamebytes(strlist[j].byte1,diffs,intorb,&axj);
                        numxc = int_min(samei,samej);
                        numxv = numxv - numxcv;
                        numx  = numxc + numxcv + numxv;
                        if (numx != x) continue;
                        ptr[i][cnt] = j;
                        cnt++;
                }
                ptr[i][0] = cnt;
        }

        return ptr;
}

/*
 * compute_ci_elecs_and_orbitals: compute number of ci electrons and orbitals
 */
void compute_ci_elecs_and_orbitals(int aelec, int belec, int orbitals, int nfrzc,
                                   int nfrzv, int *ci_aelec, int *ci_belec,
                                   int *ci_orbitals)
{
        *ci_aelec = aelec - nfrzc;
        *ci_belec = belec - nfrzc;
        *ci_orbitals  = orbitals - nfrzc - nfrzv;
}        

/*
 * compute_detnum: compute the number of determinants in expansion.
 */
int compute_detnum(struct eospace *peosp, int pegrps, struct eospace *qeosp,
                   int qegrps, int ndocc, int nactv, int xlvl,
                   int **pq_spaces, int *num_pq)
{
        int dcnt = 0;    /* Determinant count. */
        int doccmin = 0; /* Minimum determinant DOCC occupations. */
        int xcnt = 0;
        int i, j;

        /* Set min DOCC occupation numbers of alpha + beta strings. The
         * excitation level is the maximum occupation for VIRT orbitals. */
        doccmin = 2 * (ndocc - int_min(ndocc, xlvl));

        /* Loop over p (alpha) string groups. */
        *num_pq = 0;
        for (i = 0; i < pegrps; i++) {
                /* Loop over q (beta) string groups. */
                for (j = 0; j < qegrps; j++) {
                        if ((peosp[i].docc + qeosp[j].docc) < doccmin) continue;
                        if ((peosp[i].virt + qeosp[j].virt) > xlvl) continue;
                        pq_spaces[*num_pq][0] = i;
                        pq_spaces[*num_pq][1] = j;
                        (*num_pq)++;
                        dcnt = dcnt + peosp[i].nstr * qeosp[j].nstr;
                }
        }

        /* Generate list of valid space pairings for each electron space
         * group. */
        for (i = 0; i < pegrps; i++) {
            xcnt = 0;
            for (j = 0; j < (*num_pq); j++) {
                if (pq_spaces[j][0] == i) {
                    peosp[i].pairs[xcnt] = pq_spaces[j][1];
                    xcnt++;
                }
            }
            peosp[i].npairs = xcnt;
        }
        for (i = 0; i < qegrps; i++) {
            xcnt = 0;
            for (j = 0; j < (*num_pq); j++) {
                if (pq_spaces[j][1] == i) {
                    qeosp[i].pairs[xcnt] = pq_spaces[j][0];
                    xcnt++;
                }
            }
            qeosp[i].npairs = xcnt;
        }

        return dcnt;
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
        int ci_elecs;
        int ptr = 0;
        int elecs[3]; /* Number of electrons in DOCC, ACTV, VIRT */
        int orbs[3];  /* Number of orbitals in DOCC, ACTV, VIRT */
        int nstr[3];  /* Number of strings in DOCC, ACTV, VIRT */
        int pegs[3];  /* Orbital index "pegs" for each space. */
        int nspcs = 0;  /* DOCC or ACTV or VIRT */
        int vstr, astr, dstr; /* VIRT, ACTV, DOCC string number */

        virt_orbs = ci_orbs - ndocc - nactv;
        vstr = binomial_coef2(virt_orbs, virt_elec);
        if (virt_elec == 0) vstr = 0;
        astr = binomial_coef2(nactv, actv_elec);
        if (actv_elec == 0) astr = 0;
        dstr = binomial_coef2(ndocc, docc_elec);
        if (docc_elec == 0) dstr = 0;

        setup_eostrings_compute(elecs, orbs, nstr, pegs, &nspcs, docc_elec,
                                actv_elec, virt_elec, dstr, vstr, astr, ndocc,
                                nactv, virt_orbs);

        ptr = *pos; /* Set counter to start of list */
        ci_elecs = elecs[0] + elecs[1] + elecs[2];

        /* Loop over first spaces's strings */
        for (i = 1; i <= nstr[0]; i++) {
                strlist[ptr].doccx = ndocc - docc_elec;
                strlist[ptr].actvx = virt_elec;

                str_adr2str(i, doccscr, elecs[0], orbs[0], docc);
                for (j = 0; j < elecs[0]; j++) {
                        strlist[ptr].string[j] = docc[j] + pegs[0];
                }
                if (nspcs == 1) {
                        /* Next string */
                        strlist[ptr].index = str_adrfind(strlist[ptr].string,
                                                         ci_elecs, ci_orbs);
                        ptr++;
                        continue;
                }
                /* Loop over second space's strings */
                for (j = 1; j <= nstr[1]; j++) {
                        str_adr2str(j, actvscr, elecs[1], orbs[1], actv);
                        for (k = 0; k < elecs[0]; k++) {
                                strlist[ptr].string[k] = docc[k] + pegs[0];
                        }
                        for (k = elecs[0]; k < elecs[0] + elecs[1]; k++) {
                                strlist[ptr].string[k] =
                                        actv[k - elecs[0]] + pegs[1];
                        }
                        if (nspcs == 2) {
                                /* Next string */
                                strlist[ptr].index = str_adrfind(
                                        strlist[ptr].string,
                                        ci_elecs, ci_orbs);
                                ptr++;
                                continue;
                        }
                        /* Loop over third space's strings */
                        for (k = 1; k <= nstr[2]; k++) {
                                str_adr2str(k, virtscr, elecs[2], orbs[2], virt);
                                for (l = 0; l < elecs[0]; l++) {
                                        strlist[ptr].string[l] = docc[l] + pegs[0];
                                }
                                for (l = elecs[0];
                                     l < elecs[0] + elecs[1]; l++) {
                                        strlist[ptr].string[l] =
                                                actv[l - elecs[0]] + pegs[1];
                                }
                                for (l = elecs[0] + elecs[1];
                                     l < elecs[0] + elecs[1] + elecs[2]; l++) {
                                        strlist[ptr].string[l] =
                                                virt[l - elecs[0] - elecs[1]] + pegs[2];
                                }
                                /* Next string */
                                strlist[ptr].index = str_adrfind(
                                        strlist[ptr].string,
                                        ci_elecs, ci_orbs);

                                ptr++;
                        }
                }
        }
        *pos = ptr;
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
        int vorbs = 0;
        int i, j, k;

        vorbs = orbs - ndocc - nactv;
        min_docc_elecs = ndocc - int_min(ndocc, xlvl);
        max_docc_elecs = ndocc;
        min_actv_elecs = (elecs - ndocc) - int_min((elecs - ndocc), xlvl);
        //max_actv_elecs = int_min(((elecs - ndocc) + xlvl), nactv);
        max_actv_elecs = int_min(int_min(ndocc, xlvl) + (elecs - ndocc), nactv);
        /* Loop over holes in DOCC space */
        for (i = min_docc_elecs; i <= max_docc_elecs; i++) {
                /* Loop over ACTV space */
                for (j = min_actv_elecs; j <= max_actv_elecs; j++) {
                        /* Loop over virtual space */
                        for (k = 0; k <= int_min(vorbs, xlvl); k++) {
                                num = num + string_number(i, j, k, ndocc, nactv,
                                                          orbs, elecs);
                        }
                }
        }

        return num;
}

/*
 * deallocate_eostrings_array: deallocate a struct eostring *array.
 */
void deallocate_eostrings_array(struct eostring *array, int nstr)
{
        struct eostring *ptr = NULL;
        ptr = array;
        for (int i = 0; i < nstr; i++) {
                free(ptr[i].string);
        }
        free(ptr);
}

/*
 * generate_binstring_list: convert the integer string lists to struct occstr
 */
void generate_binstring_list(struct eostring *str, int nstr, int elec,
                             int ndocc, int nactv, struct occstr *binstr)
{
        for (int i = 0; i < nstr; i++) {
                binstr[i] = str2occstr(str[i].string, elec, ndocc, nactv);
                for (int j = 0; j < elec; j++) {
                    binstr[i].istr[j] = str[i].string[j];
                }
        }
        return;
}

/*
 * generate_determinant_list: generate determinant list.
 */
void generate_determinant_list(struct eostring *pstrlist, int npstr, int aelec,
                               struct eostring *qstrlist, int nqstr, int belec,
                               struct eospace *peosp, int pegrps,
                               struct eospace *qeosp, int qegrps,
                               int ndocc, int nactv, int xlvl, int *dtrm_len)
{
        int dcnt = 0; /* Determinant count */
        int doccmin = 0; /* Minimum occupation for DOCC orbitals */
        int virtmax = 0; /* Maximum occupation ofr VIRT orbitals */
        FILE *fptr = NULL;
        int i, j;

        /* Set max/min occupation numbers of alpha + beta strings */
        doccmin = 2 * (ndocc - int_min(ndocc, xlvl));
        virtmax = xlvl;

        /* Open determinant file */
        fptr = fopen("det.list","w");
        
        /* Loop over p string groups */
        for (i = 0; i < pegrps; i++) {
                /* Loop over q string groups */
                for (j = 0; j < qegrps; j++) {
                        if ((peosp[i].docc + qeosp[j].docc) < doccmin) continue;
                        if ((peosp[i].virt + qeosp[j].virt) > virtmax) continue;

                        write_determinant_strpairs(fptr, peosp[i].start,
                                                   peosp[i].nstr, qeosp[j].start,
                                                   qeosp[j].nstr, &dcnt,
                                                   pstrlist, qstrlist);
                }
        }
        *dtrm_len = dcnt;
        return;
}

/*
 * generate_determinant_list_rtnlist: generate determinant list, returning
 * the list.
 */
void generate_determinant_list_rtnlist(struct eostring *pstrlist, int npstr,
                                       int aelec, struct eostring *qstrlist,
                                       int nqstr, int belec,
                                       struct eospace *peosp, int pegrps,
                                       struct eospace *qeosp, int qegrps,
                                       int ndocc, int nactv, int xlvl,
                                       int dtrm_len, struct det *dtlist)
{
        int dcnt = 0;     /* Determinant count */
        int doccmin = 0;  /* Minimum occupation for DOCC orbitals */
        int virtmax = 0;  /* Maximum occupation ofr VIRT orbitals */
        int cflag = 0;    /* CAS flag */
        int i, j;

        /* Set max/min occupation numbers of alpha + beta strings */
        doccmin = 2 * (ndocc - int_min(ndocc, xlvl));
        virtmax = xlvl;

        /* Loop over p string groups */
        for (i = 0; i < pegrps; i++) {
                /* Loop over q string groups */
                for (j = 0; j < qegrps; j++) {
                        if ((peosp[i].docc + qeosp[j].docc) < doccmin) continue;
                        if ((peosp[i].virt + qeosp[j].virt) > virtmax) continue;
                        if ((peosp[i].virt + qeosp[j].virt) == 0) {
                                cflag = 1;
                        } else {
                                cflag = 0;
                        }

                        write_determinant_strpairs_dtlist(peosp[i].start,
                                                          peosp[i].nstr,
                                                          qeosp[j].start,
                                                          qeosp[j].nstr,
                                                          pstrlist, qstrlist,
                                                          aelec, belec, ndocc,
                                                          nactv, dtlist, &dcnt,
                                                          cflag);
                }
        }
        return;
}

/*
 * generate_single_excitations: generate the single excitations in
 * EOSPACE for an input string.
 */
int generate_single_excitations(struct occstr str, struct eospace eosp,
                                int nelec, int ndocc, int nactv,
                                int intorb, int vorbs,
                                struct xstr *singlex,
                                int *elecs, int *orbsx)
{
    int n1x = 0;
    int nvo = 0;  /* Number of virtual orbitals */
    /* str eospace information */
    int str_docc = 0;
    int str_actv = 0;
    int str_virt = 0;
    
    get_string_eospace_info(str, ndocc, nactv, &str_docc, &str_actv, &str_virt);
    if (abs(str_docc - eosp.docc) > 1) return n1x;
    if (abs(str_actv - eosp.actv) > 1) return n1x;
    if (abs(str_virt - eosp.virt) > 1) return n1x;

    /* Get available orbitals */
    get_available_orbital_list(str, intorb, vorbs, orbsx, &nvo);
    
    /* Internal excitations */
    generate_doccx(1, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                   vorbs, nelec, orbsx, singlex, &n1x, nvo);
    generate_actvx(1, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                   vorbs, nelec, orbsx, singlex, &n1x, nvo);
    generate_virtx(1, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                   vorbs, nelec, orbsx, singlex, &n1x, nvo);
    /* External excitations */
    generate_docc2virtx(1, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                        vorbs, nelec, orbsx, singlex, &n1x, nvo);
    generate_docc2actvx(1, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                        vorbs, nelec, orbsx, singlex, &n1x, nvo);
    generate_actv2virtx(1, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                        vorbs, nelec, orbsx, singlex, &n1x, nvo);

    return n1x;
}

/*
 * generate_double_excitations: generate the double excitations in
 * EOSPACE for an input string.
 */
int generate_double_excitations(struct occstr str, struct eospace eosp,
                                int nelec, int ndocc, int nactv,
                                int intorb, int vorbs,
                                struct xstr *doublex,
                                int *elecs, int *orbsx)
{
    int n2x = 0;
    int nvo = 0;
    
    /* str eospace information */
    int str_docc = 0;
    int str_actv = 0;
    int str_virt = 0;

    get_string_eospace_info(str, ndocc, nactv, &str_docc, &str_actv, &str_virt);
    if (abs(str_docc - eosp.docc) > 2) return n2x;
    if (abs(str_actv - eosp.actv) > 2) return n2x;
    if (abs(str_virt - eosp.virt) > 2) return n2x;

    /* Get available orbitals */
    get_available_orbital_list(str, intorb, vorbs, orbsx, &nvo);

    /* Internal excitations */
    generate_doccx(2, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                   vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_actvx(2, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                   vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_virtx(2, str, str_docc, str_actv, str_virt, eosp, ndocc, nactv,
                   vorbs, nelec, orbsx, doublex, &n2x, nvo);
    /* External excitations */
    generate_docc2virtx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                        nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_docc2actvx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                        nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_actv2virtx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                        nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    /* External + External excitations */
    generate_docc2actvvirtx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                            nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_actv2doccvirtx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                            nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_virt2doccactvx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                            nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    /* Internal + External excitations */
    generate_actv2virtx_actv1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_actv2virtx_docc1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_actv2virtx_virt1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_docc2actvx_actv1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_docc2actvx_docc1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_docc2actvx_virt1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_docc2virtx_docc1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_docc2virtx_actv1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_docc2virtx_virt1(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                              nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    /* Internal + Internal excitations */
    generate_doccx_actvx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                         nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_doccx_virtx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                         nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    generate_actvx_virtx(2, str, str_docc, str_actv, str_virt, eosp, ndocc,
                         nactv, vorbs, nelec, orbsx, doublex, &n2x, nvo);
    return n2x;
}

/*
 * generate_actvx: generate excitations within the ACTV space given a
 * string.
 */
void generate_actvx(int nrep, struct occstr str, int str_docc, int str_actv,
                    int str_virt, struct eospace eosp, int ndocc, int nactv,
                    int nvirt, int elec, int *scr, struct xstr *xlist,
                    int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    //long long int ibyte = 0x00; /* Internal orbitals */
    //long long int xbyte = 0x00; /* Excitation orbitals */
    int intorb;
    int i, j, k, l;
    intorb = ndocc + nactv;

    /* Check if excitations within ACTV are possible */
    if (str_actv < nrep) return;
    if (str_actv != eosp.actv) return;
    if (str_actv == nactv) return;
    if (str_actv == 0) return;
    if (str_docc != eosp.docc) return;
    if (str_virt != eosp.virt) return;
    
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    
    /* Loop over occupied orbitals */
    newstr = str;
    switch (nrep) {
    case 1: /* Single replacements. Skip DOCC. */
        for (i = (str_actv + str_docc - 1); i >= str_docc; i--) {
            for (j = (ndocc - str_docc);
                 j < (intorb - str_actv - str_docc); j++) {
                newstr = str;
                newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv, nvirt,
                                                    elec, elecs);
                xlist[*numx].io[0] = str.istr[i];
                xlist[*numx].io[1] = 0;
                xlist[*numx].fo[0] = scr[j];
                xlist[*numx].fo[1] = 0;
                
                xlist[*numx].permx = pindex_single_rep(str.istr, str.istr[i],
                                                       scr[j], elec);
                
                (*numx)++;
            }
        }
        break;
    case 2: /* Double replacements */
        if ((nactv - str_actv) < 2) return; /* Can't have 2x */
        for (i = (str_actv + str_docc - 1); i >= str_docc + 1; i--) {
            for (j = (ndocc - str_docc);
                 j < (intorb - str_actv - str_docc) - 1; j++) {
                for (k = (i - 1); k >= str_docc; k--) {
                    for (l = (j + 1); l < (intorb - str_actv - str_docc); l++) {
                        newstr = str;
                        newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                        newstr.byte1 = newstr.byte1 - pow(2, (str.istr[k] - 1));
                        newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                        newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                        xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                            nvirt, elec, elecs);
                        xlist[*numx].io[0] = str.istr[i];
                        xlist[*numx].io[1] = str.istr[k];
                        xlist[*numx].fo[0] = scr[j];
                        xlist[*numx].fo[1] = scr[l];
                        
                        xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                          
                          
                          ;
                        (*numx)++;
                    }
                }
            }
        }
        break;
    default:
        printf(" case != {1, 2}.\n");
        break;
    }
    return;
}

/*
 * generate_actv2virtx: generate replacements between ACTV and VIRT spaces.
 */
void generate_actv2virtx(int nrep, struct occstr str, int str_docc, int str_actv,
                         int str_virt, struct eospace eosp, int ndocc, int nactv,
                         int nvirt, int elec, int *scr, struct xstr *xlist,
                         int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    int xtype = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int tmp;
    //int nvo = 0;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (abs(str_actv - eosp.actv) != nrep) return;
    if (abs(str_virt - eosp.virt) != nrep) return;
    if (str_docc != eosp.docc) return;
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation ACTV -> VIRTX or ACTV <- VIRTX */
    xtype = str_actv - eosp.actv;
    if (xtype > 0) {
        /* ACTV(str) -> VIRT(eosp) */
        newstr = str;
        switch (nrep) {
        case 1: /* Single replacements */
            /* Loop over occupied ACTV orbitals. Removing them */
            for (i = (str_docc + str_actv - 1); i >= str_docc; i--) {
                /* Add new orbital in virtual orbital space */
                for (j = intorb; j < (nvo + intorb); j++) {
                    newstr = str;
                    newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                    newstr.virtx[0] = str.virtx[0];
                    newstr.virtx[1] = str.virtx[1];
                    newstr.virtx[str.nvrtx] = scr[j];
                    newstr.nvrtx = str.nvrtx + 1;
                    if (newstr.virtx[0] > newstr.virtx[1] && newstr.virtx[1] != 0) {
                        tmp = newstr.virtx[0];
                        newstr.virtx[0] = newstr.virtx[1];
                        newstr.virtx[1] = tmp;
                    }
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                        nvirt, elec, elecs);
                    xlist[*numx].io[0] = str.istr[i];
                    xlist[*numx].io[1] = 0;
                    xlist[*numx].fo[0] = scr[j];
                    xlist[*numx].fo[1] = 0;
                    
                    xlist[*numx].permx = pindex_single_rep(str.istr, str.istr[i],
                                                       scr[j], elec);;
                                    
                    (*numx)++;
                }
            }
            break;
        case 2:
            /* Loop over occupied ACTV orbitals. Removing them. */
            for (i = (str_docc + str_actv - 1); i >= str_docc + 1; i--) {
                for (j = (i - 1); j >= str_docc; j--) {
                    /* Add new orbitals */
                    for (k = intorb; k < (nvo + intorb - 1); k++) {
                        for (l = k + 1; l < (nvo + intorb); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            /* If replacing double excitation into virt space,
                             * then str_virt == 0. */
                            newstr.virtx[0] = scr[k];
                            newstr.virtx[1] = scr[l];
                            newstr.nvrtx = 2;
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                                nvirt, elec, elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];
                            
                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* ACTV(str) <- VIRT(eosp) */
        newstr = str;
        switch(nrep) {
        case 1: /* Single replacements */
            /* Loop over unoccupied ACTV orbitals. Turning them "on" */
            for (i = (intorb - str_docc - str_actv - 1);
                 i >= (ndocc - str_docc);
                 i--){
                /* Remove one occupied virtual orbital */
                for (j = 0; j < str_virt; j++) {
                    newstr = str;
                    newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                    newstr.virtx[0] = str.virtx[0];
                    newstr.virtx[1] = str.virtx[1];
                    newstr.virtx[j] = 0;
                    newstr.nvrtx = str.nvrtx - 1;
                    /* Make sure virtual orbitals are ordered correctly after
                     * removal of virtx[j] */
                    if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                        tmp = newstr.virtx[0];
                        newstr.virtx[0] = newstr.virtx[1];
                        newstr.virtx[1] = tmp;
                    }
                    if (newstr.virtx[0] > newstr.virtx[1] && newstr.virtx[1] != 0) {
                        tmp = newstr.virtx[0];
                        newstr.virtx[0] = newstr.virtx[1];
                        newstr.virtx[1] = tmp;
                    }
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                        nvirt, elec, elecs);
                    xlist[*numx].io[0] = str.virtx[j];
                    xlist[*numx].io[1] = 0;
                    xlist[*numx].fo[0] = scr[i];
                    xlist[*numx].fo[1] = 0;
                    
                    xlist[*numx].permx = pindex_single_rep(str.istr, str.virtx[j],
                                                       scr[i], elec);;
                    (*numx)++;
                }
            }
            break;
        case 2: /* Double replacements */
            /* Loop over unoccupied ACTV orbitals. Turning them "on" */
            for (i = (intorb - str_docc - str_actv - 1); i >= (ndocc - str_docc + 1); i--){
                for (j = (i - 1); j >= (ndocc - str_docc); j--) {
                    newstr = str;
                    /* Remove occupied virtual orbitals. If double replacement,
                     * both occupations are removed. */
                    newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                    newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                    newstr.virtx[0] = 0;
                    newstr.virtx[1] = 0;
                    newstr.nvrtx = 0;
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                        nvirt, elec, elecs);
                    xlist[*numx].io[0] = str.virtx[0];
                    xlist[*numx].io[1] = str.virtx[1];
                    xlist[*numx].fo[0] = scr[i];
                    xlist[*numx].fo[1] = scr[j];
                    
                    xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                    (*numx)++;
                }
            }
            break;
        default:
            printf("WARNING: Unknown case. nrep != 1, 2\n");
            break;
        }
    }
    
    return;
}

/*
 * generate_actv2virtx_actv1: generate replacements between ACTV and VIRT spaces.
 * This is a special case where one replacement occurs within ACTV, in addition
 * to the ACTV <-> VIRT replacement.
 */
void generate_actv2virtx_actv1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    int xtype = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int tmp;
//    int nvo = 0;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (str_actv == eosp.actv) return;
    if (abs(str_virt - eosp.virt) != 1) return;
    if (str_docc != eosp.docc) return;
        
    /* Generate internal orbital byte */
//    for (i = 0; i < intorb; i++) {
//        ibyte = ibyte + pow(2, i);
//    }
    /* Get possible excitations */
//    xbyte = ibyte ^ str.byte1;
//    nonzerobits(xbyte, intorb, scr);
//    for (i = intorb + 1; i <= (intorb + nvirt); i++) {
//        if (i != str.virtx[0] && i != str.virtx[1]) {
//            scr[intorb + nvo] = i;
//            nvo++;
//        }
//    }

    /* Is the excitation ACTV -> VIRTX or ACTV <- VIRTX */
    xtype = str_actv - eosp.actv;
    if (xtype > 0) {
        /* ACTV(str) -> VIRT(eosp) + Internal ACTV replacement*/
        switch (nrep) {
        case 2:
            /* Loop over occupied ACTV orbitals. Removing them. */
            for (i = (str_docc + str_actv - 1); i >= str_docc + 1; i--) {
                for (j = (i - 1); j >= str_docc; j--) {
                    /* Add new orbitals */
                    /* ACTV replacement */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_docc - str_actv);
                         k++) {
                        for (l = intorb; l < (nvo + intorb); l++) {
                            newstr = str;
                            /* Remove */
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            /* Add */
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[str_virt] = scr[l];
                            newstr.nvrtx = str_virt + 1;
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                          nactv, nvirt, elec,
                                                          elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];
                            
                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* ACTV(str) <- VIRT(eosp) */
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over unoccupied ACTV orbitals. Turning them "on" */
            for (i = (intorb - str_docc - str_actv - 1);
                 i >= (ndocc - str_docc + 1);
                 i--){
                for (j = (i - 1); j >= (ndocc - str_docc); j--) {
                    /* Loop over occupied ACTV orbitals. Turning bits off */
                    for (k = str_docc; k < str_docc + str_actv; k++) {
                        /* loop over virtual orbitals, removing them */
                        for (l = 0; l < str_virt; l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[l] = 0;
                            newstr.nvrtx = str.nvrtx - 1;
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }

                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[k];
                            xlist[*numx].io[1] = str.virtx[l];
                            xlist[*numx].fo[0] = scr[i];
                            xlist[*numx].fo[1] = scr[j];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            printf("WARNING: Unknown case. nrep != 1, 2\n");
            break;
        }
    }
    
    return;
}

/*
 * generate_actv2virtx_virt1: generate replacements between ACTV and VIRT spaces.
 * This is a special case where one replacement occurs within VIRT, in addition
 * to the ACTV <-> VIRT replacement.
 */
void generate_actv2virtx_virt1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    int xtype = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, k, l;
    //int nvo = 0;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (abs(str_actv - eosp.actv) != 1) return;
    if (abs(str_virt - eosp.virt) != 1) return;
    if (str_virt == 0 || eosp.virt == 0) return;
    if (str_virt != 2 && eosp.virt != 2) return;
    if (str_docc != eosp.docc) return;

    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation ACTV -> VIRTX + VIRTX or ACTV <- VIRTX + VIRTX */
    xtype = str_actv - eosp.actv;
    if (xtype > 0) {
        /* ACTV(str) -> VIRT(eosp) + Internal VIRT replacement*/
        newstr = str;
        switch (nrep) {
        case 2:
            /* Loop over occupied ACTV orbitals. Removing 1. */
            for (i = (str_docc + str_actv - 1); i >= str_docc; i--) {
                /* VIRT replacement */
                for (k = intorb; k < (nvo + intorb - 1); k++) {
                    for (l = (k + 1); l < (nvo + intorb); l++) {
                        newstr = str;
                        /* Remove */
                        newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                        /* Replace virtual orbital and add excitation */
                        newstr.virtx[0] = scr[k];
                        newstr.virtx[1] = scr[l];
                        newstr.nvrtx = 2;
                        /* Automatically ordred. */
                        xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                            nactv, nvirt, elec,
                                                            elecs);
                        xlist[*numx].io[0] = str.istr[i];
                        xlist[*numx].io[1] = str.virtx[0];
                        xlist[*numx].fo[0] = scr[k];
                        xlist[*numx].fo[1] = scr[l];

                        xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                        (*numx)++;
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* ACTV(str) <- VIRT(eosp)  + VIRT*/
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over unoccupied ACTV orbitals. Turning them "on" */
            for (i = (intorb - str_docc - str_actv - 1);
                 i >= (ndocc - str_docc);
                 i--){
                /* Loop over unoccupied Virtuals for replacement */
                for (l = intorb; l < (intorb + nvo); l++) {
                    newstr = str;
                    newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                    newstr.virtx[0] = scr[l];
                    newstr.virtx[1] = 0;
                    newstr.nvrtx = 1;
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                        nactv, nvirt, elec,
                                                        elecs);
                    xlist[*numx].io[0] = str.virtx[0];
                    xlist[*numx].io[1] = str.virtx[1];
                    xlist[*numx].fo[0] = scr[i];
                    xlist[*numx].fo[1] = scr[l];

                    xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                    (*numx)++;
                }
            }
            break;
        default:
            printf("WARNING: Unknown case. nrep != 1, 2\n");
            break;
        }
    }
    
    return;
}

/*
 * generate_actv2virtx_docc1: generate replacements between DOCC, ACTV and VIRT
 * spaces. This is a special case where one replacement occurs within DOCC, in
 * addition to the ACTV <-> VIRT replacement.
 */
void generate_actv2virtx_docc1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    int xtype = 0;
    int tmp;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (nrep < 2) return;
    if (str_docc != eosp.docc) return;
    if (str_actv == eosp.actv) return;
    if (str_virt == eosp.virt) return;
    if (abs(str_virt - eosp.virt) != 1) return;
    
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    /* VIRT */
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation DOCC, ACTV->VIRT or DOCC, ACTV<-VIRT */
    xtype = str_actv - eosp.actv;
    if (xtype > 0) {
        /* DOCC (str) --> (ACTV(eosp), VIRT(eosp)) */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied DOCC orbitals. Removing one */
            for (i = (str_docc - 1); i >= 0; i--) {
                /* Loop over occupied ACTV orbitals. Removing one. */
                for (j = (str_docc + str_actv - 1); j >= str_docc; j--) {
                    /* Add new occupation to DOCC space */
                    for (k = 0;
                         k < (ndocc - str_docc);
                         k++) {
                        /* Add new orbital to VIRT space */
                        for (l = intorb; l < (nvo + intorb); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            /* If adding occupation to VIRT, there must be
                             * str.nvirt = 0 or 1. */
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[str.nvrtx] = scr[l];
                            newstr.nvrtx = str.nvrtx + 1;
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str) (ACTV(eosp) <-- VIRT(eosp)) */
        /* Need to add occupations to ACTV and remove them from VIRT */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied DOCC orbitals. Turning them "off" */
            for (i = (str_docc - 1); i >= 0; i--) {
                /* Loop over unoccupied DOCC orbitals. Turning them "on" */
                for (j = 0; j < (ndocc - str_docc); j++) {
                    /* Loop over unoccupied ACTV orbitals, turning them "on" */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_docc - str_actv);
                         k++) {
                        /* Loop over occupied VIRT orbitals, turning them "off"*/
                        for (l = 0; l < str_virt; l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.nvrtx = str_virt - 1;
                            newstr.virtx[l] = 0;
                            /* Make sure virtual orbitals are ordered correctly after
                             * removal of virtx[j] */
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.virtx[l];
                            xlist[*numx].fo[0] = scr[j];
                            xlist[*numx].fo[1] = scr[k];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    }
    return;

}


/*
 * generate_actv2doccvirtx: generate replacements between ACTV and (DOCC,VIRT)
 * spaces.
 */
void generate_actv2doccvirtx(int nrep, struct occstr str, int str_docc,
                             int str_actv, int str_virt, struct eospace eosp,
                             int ndocc, int nactv, int nvirt, int elec, int *scr,
                             struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    int xtype = 0;
    //long long ibyte = 0x00;
    //long long xbyte = 0x00;
    int tmp;
    //int nvo = 0;
    int i, j, k, l;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (nrep < 2) return;
    if (str_actv == eosp.actv) return;
    if (str_docc == eosp.docc) return;
    if (str_virt == eosp.virt) return;
    if (abs(str_actv - eosp.actv) != nrep) return;

    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    /* VIRT */
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation ACTV -> (DOCC, VIRT) or ACTV <- (DOCC, VIRT) */
    xtype = str_actv - eosp.actv;
    if (xtype > 0) {
        /* ACTV -> (DOCC, VIRT) */
        newstr = str;
        switch (nrep) {
        case 2: /* Double replacements */
            /* Loop over ACTV orbitals, removing them */
            for (i = (str_docc + str_actv - 1); i >= str_docc + 1; i--) {
                for (j = (i - 1); j >= str_docc; j--) {
                    /* Add new orbital to DOCC */
                    for (k = 0; k < (ndocc - str_docc); k++) {
                        /* Add new orbital to VIRT */
                        for (l = intorb; l < (nvo + intorb); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[str.nvrtx] = scr[l];
                            newstr.nvrtx = str.nvrtx + 1;
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* ACTV <- (DOCC, VIRT) */
        newstr = str;
        switch (nrep) {
        case 2: /* Double replacements */
            /* Loop over unoccupied ACTV orbitals, adding them */
            for (i = (intorb - str_docc - str_actv - 1);
                 i >= (ndocc - str_docc + 1);
                 i--) {
                for (j = (i - 1); j >= (ndocc - str_docc); j--) {
                    /* Loop over occupied DOCC orbitals, removing them */
                    for (k = 0; k < str_docc; k++) {
                        /* Loop over occupied VIRT orbitals, removing them */
                        for (l = 0; l < str_virt; l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.nvrtx = str.nvrtx - 1;
                            newstr.virtx[l] = 0;
                            /* Make sure virtual orbitals are ordered correctly after
                             * removal of virtx[j] */
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[k];
                            xlist[*numx].io[1] = str.virtx[l];
                            xlist[*numx].fo[0] = scr[i];
                            xlist[*numx].fo[1] = scr[j];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    }
    return;
}

/*
 * generate_docc2actvx: generate replacements between DOCC and ACTV spaces.
 */
void generate_docc2actvx(int nrep, struct occstr str, int str_docc, int str_actv,
                         int str_virt, struct eospace eosp, int ndocc, int nactv,
                         int nvirt, int elec, int *scr, struct xstr *xlist,
                         int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    int xtype = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;
    /* Check if excitations are possible */
    if (abs(str_docc - eosp.docc) != nrep) return;
    if (abs(str_actv - eosp.actv) != nrep) return;
    if (str_virt != eosp.virt) return;
    
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);

    /* Is the excitation DOCC -> ACTV or DOCC <- ACTV */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC(str) -> ACTV(eosp) */
        /* Need to remove occupations from DOCC, and add to ACTV */
        newstr = str;
        switch(nrep) {
        case 1: /* Single replacements */
            /* Loop over occupied docc orbitals. Removing them */
            for (i = (str_docc - 1); i >= 0; i--) {
                /* Add new orbital in actv space */
                for (j = (ndocc - str_docc);
                     j < (intorb - str_actv - str_docc);
                     j++) {
                    newstr = str;
                    newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                    newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                        nvirt, elec, elecs);
                    xlist[*numx].io[0] = str.istr[i];
                    xlist[*numx].io[1] = 0;
                    xlist[*numx].fo[0] = scr[j];
                    xlist[*numx].fo[1] = 0;

                    xlist[*numx].permx = pindex_single_rep(str.istr, str.istr[i],
                                                           scr[j], elec);;
                    (*numx)++;
                }
            }
            break;
        case 2: /* Double replacements */
            /* Loop over occupied docc orbitals. Removing them. */
            for (i = (str_docc - 1); i >= 1; i--) {
                for (j = (i - 1); j >= 0; j--) {
                    /* Add new orbitals to actv space */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_actv - str_docc - 1);
                         k++){
                        for (l = (k + 1);
                             l < (intorb - str_actv - str_docc);
                             l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str) <- ACTV(eosp) */
        /* Need to add occupations to DOCC, and remove from ACTV */
        newstr = str;
        switch(nrep) {
        case 1: /* Single replacements */
            /* Loop over occupied actv orbitals. Removing them */
            for (i = (str_actv + str_docc - 1); i >= str_docc; i--) {
                /* Add new orbital in docc space */
                for (j = 0; j < (ndocc - str_docc); j++) {
                    newstr = str;
                    newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                    newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                        nvirt, elec, elecs);
                    xlist[*numx].io[0] = str.istr[i];
                    xlist[*numx].io[1] = 0;
                    xlist[*numx].fo[0] = scr[j];
                    xlist[*numx].fo[1] = 0;

                    xlist[*numx].permx = pindex_single_rep(str.istr, str.istr[i],
                                                           scr[j], elec);;
                    (*numx)++;
                }
            }
            break;
        case 2: /* Double replacements */
            /* Loop over occupied actv orbitals. Removing them */
            for (i = (str_actv + str_docc - 1); i >= str_docc + 1; i--) {
                for (j = (i - 1); j >= str_docc; j--) {
                    /* Loop over unoccupied DOCC orbitals. Turning them on. */
                    for (k = 0; k < (ndocc - str_docc - 1); k++) {
                        for (l = (k + 1); l < (ndocc - str_docc); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }   
    }
    return;
}

/*
 * generate_docc2actvx_actv1: generate excitations from DOCC -> ACTV with
 * a replacement within the ACTV.
 */
void generate_docc2actvx_actv1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    int xtype = 0;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (abs(str_docc - eosp.docc) != 1) return;
    if (str_virt != eosp.virt) return;
    if (str_actv == eosp.actv) return;
    
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation DOCC -> ACTV or ACTV -> DOCC */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC(str) -> ACTV(eosp) */
        /* Need to remove 1 occupation from DOCC, and add to ACTV */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied docc orbitals. Removing them. */
            for (i = (str_docc - 1); i >= 0; i--) {
                /* Loop over occupied ACTV orbitals. Turning them "off" */
                for (j = str_docc + str_actv - 1; j >= str_docc; j--) {
                    /* Add new orbital to actv space */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_actv - str_docc - 1);
                         k++){
                        /* Add new ACTV orbital to ACTV space */
                        for (l = k+1; l < (intorb - str_actv - str_docc);
                             l++) {
                            newstr = str;
                            newstr.nvrtx = str.nvrtx;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str) <- ACTV(eosp) */
        /* Need to add occupation to DOCC, and remove from ACTV */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied ACTV orbitals. Removing 2 */
            for (i = (str_actv + str_docc - 1); i >= str_docc + 1; i--) {
                for (j = (i - 1); j >= str_docc; j--) {
                    /* Loop over unoccupied ACTV orbitals. Turning 1 on. */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_actv - str_docc);
                         k++) {
                        /* Loop over unoccupied DOCC orbitals Turn 1 on */
                        for (l = 0; l < (ndocc - str_docc); l++) {
                            newstr = str;
                            newstr.nvrtx = str.nvrtx;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }   
    }
    
    return;
}


/*
 * generate_docc2actvx_docc1: generate excitations from DOCC -> ACTV with
 * a replacement within the DOCC.
 */
void generate_docc2actvx_docc1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    int xtype = 0;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (abs(str_docc - eosp.docc) != 1) return;
    if (str_virt != eosp.virt) return;
    if (str_actv == eosp.actv) return;
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //     if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation DOCC -> ACTV or ACTV <- DOCC */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC(str) -> ACTV(eosp) */
        /* Need to remove 1 occupation from DOCC, and add to ACTV */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied docc orbitals. Removing them. */
            for (i = (str_docc - 1); i >= 1; i--) {
                for (j = (i - 1); j >= 0; j--) {
                    /* Add new orbital to actv space */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_actv - str_docc);
                         k++){
                        /* Add new DOCC orbital to DOCC space */
                        for (l = 0;
                             l < (ndocc - str_docc);
                             l++) {
                            newstr = str;
                            newstr.nvrtx = str.nvrtx;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str) <- ACTV(eosp) */
        /* Need to add occupation to DOCC, and remove from ACTV */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied actv orbitals. Removing 1 */
            for (i = (str_actv + str_docc - 1); i >= str_docc; i--) {
                /* Loop over occupied DOCC orbital. Removing 1 */
                for (j = (str_docc - 1); j >= 0; j--) {
                    /* Loop over unoccupied DOCC orbitals. Turning them on. */
                    for (k = 0; k < (ndocc - str_docc - 1); k++) {
                        for (l = (k + 1); l < (ndocc - str_docc); l++) {
                            newstr = str;
                            newstr.nvrtx = str.nvrtx;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }   
    }
    
    return;
}

/*
 * generate_docc2actvx_virt1: generate excitations from DOCC -> ACTV with
 * a replacement within the VIRT.
 */
void generate_docc2actvx_virt1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    int xtype = 0;
    int tmp;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (abs(str_docc - eosp.docc) != 1) return;
    if (str_virt != eosp.virt) return;
    if (str_actv == eosp.actv) return;

    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation DOCC -> ACTV or ACTV <- DOCC */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC(str) -> ACTV(eosp) */
        /* Need to remove 1 occupation from DOCC, and add to ACTV */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied docc orbitals. Removing one. */
            for (i = (str_docc - 1); i >= 0; i--) {
                /* Loop over virtual orbitals. Swapping one */
                for (j = 0; j < str_virt; j++) {
                    /* Add new orbital to actv space */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_actv - str_docc);
                         k++){
                        /* Add new virtual orbitual space */
                        for (l = intorb; l < (intorb + nvo); l++) {
                            newstr = str;
                            newstr.nvrtx = str.nvrtx;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[j] = scr[l];
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.virtx[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str) <- ACTV(eosp) */
        /* Need to add occupation to DOCC, and remove from ACTV */
        /* And perform VIRT replacement. */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied ACTV orbitals. Removing 1 */
            for (i = (str_actv + str_docc - 1); i >= str_docc; i--) {
                /* Loop over occupied VIRT orbital. Swapping orbital */
                for (j = (str_virt - 1); j >= 0; j--) {
                    /* Loop over unoccupied DOCC orbitals. Turning them on. */
                    for (k = 0; k < (ndocc - str_docc); k++) {
                        /* Add new virtual orbital */
                        for (l = intorb; l < (intorb + nvo); l++) {
                            newstr = str;
                            newstr.nvrtx = str.nvrtx;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[j] = scr[l];
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.virtx[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }   
    }
    
    return;
    
    return;
}

/*
 * generate_docc2virtx_docc1: generate excitations from DOCC -> VIRT with a
 * replacement within the DOCC space.
 */
void generate_docc2virtx_docc1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    int xtype = 0;
    int tmp;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (nrep < 2) return;
    if (str_docc == eosp.docc) return;
    if (str_actv != eosp.actv) return;
    if (str_virt == eosp.virt) return;
    if (abs(str_virt - eosp.virt) != 1) return;

    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    /* VIRT */
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation DOCC, DOCC->VIRT or DOCC, DOCC<-VIRT */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC (str) --> VIRT(eosp) */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied DOCC orbitals. Removing two */
            for (i = (str_docc - 1); i >= 1; i--) {
                for (j = (i - 1); j >= 0; j--) {
                    /* Add new occupation to DOCC space */
                    for (k = 0; k < (ndocc - str_docc); k++) {
                        /* Add new orbital to VIRT space */
                        for (l = intorb; l < (nvo + intorb); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            /* If adding occupation to VIRT, there must be
                             * str.nvirt = 0 or 1. */
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[str.nvrtx] = scr[l];
                            newstr.nvrtx = str.nvrtx + 1;
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str,eosp) (DOCC(eosp) <-- VIRT(str)) */
        /* Need to add occupations to DOCC and remove them from VIRT */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied DOCC orbitals. Turning one "off" */
            for (i = (str_docc - 1); i >= 0; i--) {
                /* Loop over unoccupied DOCC orbitals. Turning two "on" */
                for (j = 0; j < (ndocc - str_docc - 1); j++) {
                    for (k = (j + 1); k < (ndocc - str_docc); k++) {
                        /* Loop over occupied VIRT orbitals, turning them "off"*/
                        for (l = 0; l < str_virt; l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.nvrtx = str_virt - 1;
                            newstr.virtx[l] = 0;
                            /* Make sure virtual orbitals are ordered correctly after
                             * removal of virtx[j] */
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.virtx[l];
                            xlist[*numx].fo[0] = scr[j];
                            xlist[*numx].fo[1] = scr[k];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    }
    return;

}

/*
 * generate_docc2virtx_actv1: generate excitations from DOCC -> VIRT with
 * a replacement within the ACTV.
 */
void generate_docc2virtx_actv1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    int xtype = 0;
    int tmp;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;
    
    /* Check if excitations are possible */
    if (nrep < 2) return;
    if (str_docc == eosp.docc) return;
    if (str_actv != eosp.actv) return;
    if (str_virt == eosp.virt) return;
    if (abs(str_virt - eosp.virt) != 1) return;

    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    /* VIRT */
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation ACTV, DOCC->VIRT or ACTV, DOCC<-VIRT */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC (str) --> VIRT(eosp) */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied DOCC orbitals. Removing one */
            for (i = (str_docc - 1); i >= 0; i--) {
                /* Loop over occupied ACTV orbitals. Removing one. */
                for (j = (str_docc + str_actv - 1); j >= str_docc; j--) {
                    /* Add new occupation to ACTV space */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_docc - str_actv); k++) {
                        /* Add new orbital to VIRT space */
                        for (l = intorb; l < (nvo + intorb); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            /* If adding occupation to VIRT, there must be
                             * str.nvirt = 0 or 1. */
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[str.nvrtx] = scr[l];
                            newstr.nvrtx = str.nvrtx + 1;
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[(*numx)].index = occstr2address(newstr, eosp, ndocc,
                                                                  nactv, nvirt, elec,
                                                                  elecs);
                            xlist[(*numx)].io[0] = str.istr[i];
                            xlist[(*numx)].io[1] = str.istr[j];
                            xlist[(*numx)].fo[0] = scr[k];
                            xlist[(*numx)].fo[1] = scr[l];

                            xlist[(*numx)].permx = pindex_double_rep_str(newstr.istr,
                                                                         xlist[(*numx)].io[0],
                                                                         xlist[(*numx)].fo[0],
                                                                         xlist[(*numx)].io[1],
                                                                         xlist[(*numx)].fo[1],
                                                                         elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* ACTV, DOCC(eosp) <-- VIRT(str) */
        /* Need to add occupation to DOCC and remove one from VIRT */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over ACTV occupied  orbitals. Turning one "off" */
            for (i = (str_docc + str_actv - 1); i >= str_docc; i--) {
                /* Loop over unoccupied ACTV orbitals. Turning one "on" */
                for (j = (ndocc - str_docc);
                     j < (intorb - str_actv - str_docc); j++) {
                    /* Loop over unoccupied DOCC orbitals. Turning one "on" */
                    for (k = 0; k < (ndocc - str_docc); k++) {
                        /* Loop over occupied VIRT orbitals, turning them "off"*/
                        for (l = 0; l < str_virt; l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.nvrtx = str_virt - 1;
                            newstr.virtx[l] = 0;
                            /* Make sure virtual orbitals are ordered correctly after
                             * removal of virtx[j] */
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.virtx[l];
                            xlist[*numx].fo[0] = scr[j];
                            xlist[*numx].fo[1] = scr[k];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    }
    return;

}

/*
 * generate_docc2virtx_virt1: generate excitations from DOCC -> VIRT with
 * a replacement within the virt.
 */
void generate_docc2virtx_virt1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    int xtype = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, k, l;
    //int tmp;
    //int nvo = 0;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (abs(str_docc - eosp.docc) != 1) return;
    if (abs(str_virt - eosp.virt) != 1) return;
    if (str_virt == 0 || eosp.virt == 0) return;
    if (str_virt != 2 && eosp.virt != 2) return;
    if (str_actv != eosp.actv) return;
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation DOCC -> VIRTX + VIRTX or DOCC <- VIRTX + VIRTX */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC(str) -> VIRT(eosp) + Internal VIRT replacement*/
        newstr = str;
        switch (nrep) {
        case 2:
            /* Loop over occupied DOCC orbitals. Removing 1. */
            for (i = str_docc - 1; i >= 0; i--) {
                /* VIRT replacements */
                for (k = intorb; k < (nvo + intorb - 1); k++) {
                    for (l = (k + 1); l < (nvo + intorb); l++) {
                        newstr = str; 
                        /* Remove */
                        newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                        /* Replace virtual orbital and add excitation */
                        newstr.virtx[0] = scr[k];
                        newstr.virtx[1] = scr[l];
                        newstr.nvrtx = 2;
                        /* Automatically ordred. */
                        xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                            nactv, nvirt, elec,
                                                            elecs);
                        xlist[*numx].io[0] = str.istr[i];
                        xlist[*numx].io[1] = str.virtx[0];
                        xlist[*numx].fo[0] = scr[k];
                        xlist[*numx].fo[1] = scr[l];

                        xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                        (*numx)++;
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str) <- VIRT(eosp)  + VIRT*/
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over unoccupied DOCC orbitals. Turning them "on" */
            for (i = (ndocc - str_docc - 1); i >= 0; i--){
                /* Loop over unoccupied Virtuals for replacement */
                for (l = intorb; l < (intorb + nvo); l++) {
                    newstr = str;
                    newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                    newstr.virtx[0] = scr[l];
                    newstr.virtx[1] = 0;
                    newstr.nvrtx = 1;
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                        nactv, nvirt, elec,
                                                        elecs);
                    xlist[*numx].io[0] = str.virtx[0];
                    xlist[*numx].io[1] = str.virtx[1];
                    xlist[*numx].fo[0] = scr[i];
                    xlist[*numx].fo[1] = scr[l];

                    xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                    (*numx)++;
                }
            }
            break;
        default:
            printf("WARNING: Unknown case. nrep != 1, 2\n");
            break;
        }
    }

    return;
    
}

/*
 * generate_doccx_actvx: generate replacements  DOCC + ACTV.
 */
void generate_doccx_actvx(int nrep, struct occstr str, int str_docc,
                          int str_actv, int str_virt, struct eospace eosp,
                          int ndocc, int nactv, int nvirt, int elec,
                          int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    //int xtype = 0;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;
    
    /* Check if excitations are possible */
    if (str_docc != eosp.docc) return;
    if (str_virt != eosp.virt) return;
    if (str_actv != eosp.actv) return;
    if (nrep != 2) return;
    
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    //* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);

    newstr = str;
    switch(nrep) {
    case 2: /* Double replacements */
            /* Loop over occupied docc orbitals. Removing one. */
        for (i = (str_docc - 1); i >= 0; i--) {
            /* Loop over occupied actv orbitals. removing one */
            for (j = (str_docc + str_actv - 1); j >= str_docc; j--) {
                /* Add new orbital to actv space */
                for (k = (ndocc - str_docc);
                     k < (intorb - str_actv - str_docc);
                     k++){
                    /* Add new DOCC orbital to DOCC space */
                    for (l = 0;
                         l < (ndocc - str_docc);
                         l++) {
                        newstr = str;
                        newstr.nvrtx = str.nvrtx;
                        newstr.virtx[0] = str.virtx[0];
                        newstr.virtx[1] = str.virtx[1];
                        newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                        newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                        newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                        newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                        xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                            nactv, nvirt, elec,
                                                            elecs);
                        xlist[*numx].io[0] = str.istr[i];
                        xlist[*numx].io[1] = str.istr[j];
                        xlist[*numx].fo[0] = scr[k];
                        xlist[*numx].fo[1] = scr[l];

                        xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                        (*numx)++;
                    }
                }
            }
        }
        break;
    default:
        break;
    }
    return;
}

/*
 * generate_doccx_virtx: generate replacements  DOCC + VIRT.
 */
void generate_doccx_virtx(int nrep, struct occstr str, int str_docc,
                          int str_actv, int str_virt, struct eospace eosp,
                          int ndocc, int nactv, int nvirt, int elec,
                          int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    //int xtype = 0;
    int tmp;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;
    
    /* Check if excitations are possible */
    if (str_docc != eosp.docc) return;
    if (str_virt != eosp.virt) return;
    if (str_actv != eosp.actv) return;
    if (nrep != 2) return;
    
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    newstr = str;
    switch(nrep) {
    case 2: /* Double replacements */
            /* Loop over occupied docc orbitals. Removing one. */
        for (i = (str_docc - 1); i >= 0; i--) {
            /* Loop over unoccupied docc orbitals. Adding one */
            for (j = 0; j < (ndocc - str_docc); j++) {
                /* Add swap occupation in virt */
                for (k = 0; k < str_virt; k++){
                    /* Add new DOCC orbital to DOCC space */
                    for (l = intorb; l < (intorb + nvo); l++) {
                        newstr = str;
                        newstr.nvrtx = str.nvrtx;
                        newstr.virtx[0] = str.virtx[0];
                        newstr.virtx[1] = str.virtx[1];
                        newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                        newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                        newstr.virtx[k] = scr[l]; // Swapped
                        if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                            tmp = newstr.virtx[0];
                            newstr.virtx[0] = newstr.virtx[1];
                            newstr.virtx[1] = tmp;
                        }
                        if (newstr.virtx[0] > newstr.virtx[1] &&
                            newstr.virtx[1] != 0) {
                            tmp = newstr.virtx[0];
                            newstr.virtx[0] = newstr.virtx[1];
                            newstr.virtx[1] = tmp;
                        }
                        xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                            nactv, nvirt, elec,
                                                            elecs);
                        xlist[*numx].io[0] = str.istr[i];
                        xlist[*numx].io[1] = str.virtx[k];
                        xlist[*numx].fo[0] = scr[j];
                        xlist[*numx].fo[1] = scr[l];

                        xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                        (*numx)++;
                    }
                }
            }
        }
        break;
    default:
        break;
    }
    return;
}

/*
 * generate_actvx_virtx: generate replacements  ACTV + VIRT.
 */
void generate_actvx_virtx(int nrep, struct occstr str, int str_docc,
                          int str_actv, int str_virt, struct eospace eosp,
                          int ndocc, int nactv, int nvirt, int elec,
                          int *scr, struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    int xtype = 0;
    int tmp;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (str_docc != eosp.docc) return;
    if (str_virt != eosp.virt) return;
    if (str_actv != eosp.actv) return;
    if (eosp.virt == 0) return;
    if (nrep != 2) return;
    
    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    newstr = str;
    switch(nrep) {
    case 2: /* Double replacements */
            /* Loop over occupied actv orbitals. Removing one. */
        for (i = (str_docc + str_actv - 1); i >= str_docc; i--) {
            /* Loop over unoccupied actv orbitals. Adding one */
            for (j = (ndocc - str_docc);
                 j < (intorb - str_docc - str_actv);
                 j++) {
                /* Add swap occupation in virt */
                for (k = 0; k < str_virt; k++){
                    for (l = intorb; l < (intorb + nvo); l++) {
                        newstr = str;
                        newstr.nvrtx = str.nvrtx;
                        newstr.virtx[0] = str.virtx[0];
                        newstr.virtx[1] = str.virtx[1];
                        newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                        newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                        newstr.virtx[k] = scr[l]; // Swapped
                        if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                            tmp = newstr.virtx[0];
                            newstr.virtx[0] = newstr.virtx[1];
                            newstr.virtx[1] = tmp;
                        }
                        if (newstr.virtx[0] > newstr.virtx[1] &&
                            newstr.virtx[1] != 0) {
                            tmp = newstr.virtx[0];
                            newstr.virtx[0] = newstr.virtx[1];
                            newstr.virtx[1] = tmp;
                        }
                        xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                            nactv, nvirt, elec,
                                                            elecs);
                        xlist[*numx].io[0] = str.istr[i];
                        xlist[*numx].io[1] = str.virtx[k];
                        xlist[*numx].fo[0] = scr[j];
                        xlist[*numx].fo[1] = scr[l];

                        xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                        (*numx)++;
                    }
                }
            }
        }
        break;
    default:
        break;
    }
    return;
}

/*
 * generate_docc2actvvirtx: generate excitations from DOCC -> (ACTV, VIRT)
 * for a given string.
 */
void generate_docc2actvvirtx(int nrep, struct occstr str, int str_docc,
                             int str_actv, int str_virt, struct eospace eosp,
                             int ndocc, int nactv, int nvirt, int elec, int *scr,
                             struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr; /* New string. */
    int elecs[20];
    int xtype = 0;
    int tmp;
    //int nvo = 0;
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (nrep < 2) return;
    if (str_docc == eosp.docc) return;
    if (str_actv == eosp.actv) return;
    if (str_virt == eosp.virt) return;
    if (abs(str_docc - eosp.docc) != nrep) return;


    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    /* VIRT */
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is the excitation DOCC -> (ACTV, VIRT) or DOCC <- (ACTV, VIRT) */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC (str) --> (ACTV(eosp), VIRT(eosp)) */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied DOCC orbitals. Removing them */
            for (i = (str_docc - 1); i >= 1; i--) {
                for (j = (i - 1); j >= 0; j--) {
                    /* Add new orbital to ACTV space */
                    for (k = (ndocc - str_docc);
                         k < (intorb - str_actv - str_docc);
                         k++) {
                        /* Add new orbital to VIRT space */
                        for (l = intorb; l < (nvo + intorb); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[k] - 1));
                            /* If adding occupation to VIRT, there must be
                             * str.nvirt = 0 or 1. */
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.virtx[str.nvrtx] = scr[l];
                            newstr.nvrtx = str.nvrtx + 1;
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str) <-- (ACTV(eosp), VIRT(eosp)) */
        /* Need to add occupations to DOCC and remove them from ACTV and VIRT */
        newstr = str;
        switch(nrep) {
        case 2: /* Double replacements */
            /* Loop over unoccupied DOCC orbitals. Turning them "on" */
            for (i = (ndocc - str_docc - 1); i >= 1; i--) {
                for (j = (i - 1); j >= 0; j--) {
                    /* Loop over occupied ACTV orbitals, turning them "off" */
                    for (k = str_docc; k < (str_docc + str_actv); k++) {
                        /* Loop over occupied VIRT orbitals, removing them */
                        for (l = 0; l < str_virt; l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                            newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, ((str.istr[k] - 1)));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.nvrtx = str_virt - 1;
                            newstr.virtx[l] = 0;
                            /* Make sure virtual orbitals are ordered correctly after
                             * removal of virtx[j] */
                            if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            if (newstr.virtx[0] > newstr.virtx[1] &&
                                newstr.virtx[1] != 0) {
                                tmp = newstr.virtx[0];
                                newstr.virtx[0] = newstr.virtx[1];
                                newstr.virtx[1] = tmp;
                            }
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[k];
                            xlist[*numx].io[1] = str.virtx[l];
                            xlist[*numx].fo[0] = scr[i];
                            xlist[*numx].fo[1] = scr[j];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    }
    return;
}
    
/*
 * generate_docc2virtx: generate excitations from DOCC -> VIRT spaces for
 * a given string.
 */
void generate_docc2virtx(int nrep, struct occstr str, int str_docc, int str_actv,
                         int str_virt, struct eospace eosp, int ndocc, int nactv,
                         int nvirt, int elec, int *scr, struct xstr *xlist,
                         int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    //long long int ibyte = 0x00;
    //long long int xbyte = 0x00;
    int i, j, k, l;
    int tmp;
    //int nvo = 0;
    int xtype = 0;
    int intorb = ndocc + nactv;
    
    /* Check if excitations are possible */
    if (abs(str_docc - eosp.docc) != nrep) return;
    if (str_actv != eosp.actv) return;

    /* Generate DOCC orbital byte */
    //for (i = 0; i < ndocc; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, ndocc, scr);
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}
    /* Is excitation DOCC -> VIRT or DOCC <- VIRT? */
    xtype = str_docc - eosp.docc;
    if (xtype > 0) {
        /* DOCC(str) -> VIRT(eosp) */
        /* Need to remove occupations from DOCC, and add to VIRT */
        newstr = str;
        switch(nrep) {
        case 1: /* Single replacements */
            /* Loop over occupied docc orbitals. Removing them */
            for (i = (str_docc - 1); i >= 0; i--) {
                /* Add new orbital */
                for (j = intorb; j < (nvo + intorb); j++) {
                    newstr = str;
                    newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                    newstr.virtx[0] = str.virtx[0];
                    newstr.virtx[1] = str.virtx[1];
                    newstr.virtx[str.nvrtx] = scr[j];
                    newstr.nvrtx = str.nvrtx + 1;
                    if (newstr.virtx[0] > newstr.virtx[1] && newstr.virtx[1] != 0) {
                        tmp = newstr.virtx[0];
                        newstr.virtx[0] = newstr.virtx[1];
                        newstr.virtx[1] = tmp;
                    }
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                        nvirt, elec, elecs);
                    xlist[*numx].io[0] = str.istr[i];
                    xlist[*numx].io[1] = 0;
                    xlist[*numx].fo[0] = scr[j];
                    xlist[*numx].fo[1] = 0;

                    xlist[*numx].permx = pindex_single_rep(str.istr, str.istr[i],
                                                           scr[j], elec);;
                    (*numx)++;
                }
            }
            break;
        case 2: /* Double replacements */
            /* Loop over occupied docc orbitals, removing them. */
            for (i = (str_docc - 1); i >= 1; i--) {
                for (j = (i - 1); j >= 0; j--) {
                    /* Add new orbitals to VIRT space. If 2 excitations, both
                     * virtx[] elements are replaced. */
                    for (k = intorb; k < (nvo + intorb - 1); k++) {
                        for (l = (k + 1); l < (nvo + intorb); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.virtx[0] = str.virtx[0];
                            newstr.virtx[1] = str.virtx[1];
                            newstr.nvrtx = 2;
                            newstr.virtx[0] = scr[k];
                            newstr.virtx[1] = scr[l];
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* DOCC(str) <- VIRT(eosp) */
        /* Need to add occupations to DOCC, and remove from VIRT */
        newstr = str;
        switch(nrep) {
        case 1: /* Single replacements */
            /* Loop over unoccupied docc orbitals. Adding one. */
            for (i = (ndocc - str_docc - 1); i >= 0; i--) {
                /* Remove one occupied virtual orbital */
                for (j = str_virt - 1; j >= 0; j--) {
                    newstr = str;
                    newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                    newstr.virtx[0] = str.virtx[0];
                    newstr.virtx[1] = str.virtx[1];
                    newstr.virtx[j] = 0;
                    newstr.nvrtx = str.nvrtx - 1;
                    /* Make sure virtual orbitals are ordered correctly after
                     * removal of virtx[j] */
                    if (newstr.virtx[0] == 0 && newstr.virtx[1] != 0) {
                        tmp = newstr.virtx[0];
                        newstr.virtx[0] = newstr.virtx[1];
                        newstr.virtx[1] = tmp;
                    }
                    if (newstr.virtx[0] > newstr.virtx[1] && newstr.virtx[1] != 0) {
                        tmp = newstr.virtx[0];
                        newstr.virtx[0] = newstr.virtx[1];
                        newstr.virtx[1] = tmp;
                    }
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                        nvirt, elec, elecs);
                    xlist[*numx].io[0] = str.virtx[j];
                    xlist[*numx].io[1] = 0;
                    xlist[*numx].fo[0] = scr[i];
                    xlist[*numx].fo[1] = 0;

                    xlist[*numx].permx = pindex_single_rep(str.istr, str.virtx[j],
                                                           scr[i], elec);;
                    (*numx)++;
                }
            }
            break;
        case 2: /* Double replacements */
            /* Loop over unoccupied docc orbitals, adding two */
            for (i = (ndocc - str_docc - 1); i >= 1; i--) {
                for (j = (i - 1); j >= 0; j--) {
                    /* Remove virtual orbitals. This is a double replacement
                     * so both are removed. */
                    newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                    newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                    newstr.nvrtx = 0;
                    newstr.virtx[0] = 0;
                    newstr.virtx[1] = 0;
                    xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                        nvirt, elec, elecs);
                    xlist[*numx].io[0] = str.virtx[0];
                    xlist[*numx].io[1] = str.virtx[1];
                    xlist[*numx].fo[0] = scr[i];
                    xlist[*numx].fo[1] = scr[j];

                    xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                    (*numx)++;
                }
            }
            break;
        default:
            break;
        }
    }
    return;
}

/*
 * generate_doccx: generate excitations within the DOCC space given a
 * string.
 */
void generate_doccx(int nrep, struct occstr str, int str_docc, int str_actv,
                    int str_virt, struct eospace eosp, int ndocc, int nactv,
                    int nvirt, int elec, int *scr, struct xstr *xlist,
                    int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    //long long int ibyte = 0x00; /* Internal orbitals */
    //long long int xbyte = 0x00; /* Excitation orbitals */
    //int intorb;
    int i, j, k ,l;
    //intorb = ndocc + nactv;
    
    /* Check if excitations are possible */
    if (str_docc != eosp.docc) return;
    if (str_docc == ndocc) return;
    if (str_actv != eosp.actv) return;
    if (str_virt != eosp.virt) return;

    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);

    /* Loop over occupied orbitals */
    newstr = str;
    switch (nrep) {
    case 1: /* Single replacements */
        for (i = (str_docc - 1); i >= 0; i--) {
            for (j = 0; j < (ndocc - str_docc); j++) {
                newstr = str;
                newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv, nvirt,
                                                    elec, elecs);
                xlist[*numx].io[0] = str.istr[i];
                xlist[*numx].io[1] = 0;
                xlist[*numx].fo[0] = scr[j];
                xlist[*numx].fo[1] = 0;

                xlist[*numx].permx = pindex_single_rep(str.istr, str.istr[i],
                                                       scr[j], elec);;
                (*numx)++;
            }
        }
        break;
    case 2: /* Double replacements */
        if ((ndocc - str_docc) != 2) return; /* Can't have 2x */
        for (i = (str_docc - 1); i >= 1; i--) {
            for (j = 0; j < (ndocc - str_docc - 1); j++) {
                for (k = (i - 1); k >= 0; k--) {
                    for (l = (j + 1); l < (ndocc - str_docc); l++) {
                        newstr = str;
                        newstr.byte1 = str.byte1 - pow(2, (str.istr[k] - 1));
                        newstr.byte1 = newstr.byte1 - pow(2, (str.istr[i] - 1));
                        newstr.byte1 = newstr.byte1 + pow(2, (scr[l] - 1));
                        newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                        xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                            nactv, nvirt, elec, elecs);
                        xlist[*numx].io[0] = str.istr[k];
                        xlist[*numx].io[1] = str.istr[i];
                        xlist[*numx].fo[0] = scr[l];
                        xlist[*numx].fo[1] = scr[j];

                        xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                        (*numx)++;
                    }
                }
            }
        }
        break;
    default:
        printf(" case != {1, 2}.\n");
        break;
    }
    return;
}

/*
 * generate_virtx: generate excitations within the VIRT space given a
 * string.
 */
void generate_virtx(int nrep, struct occstr str, int str_docc, int str_actv,
                    int str_virt, struct eospace eosp, int ndocc, int nactv,
                    int nvirt, int elec, int *scr, struct xstr *xlist,
                    int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    int i, j, k, l;
    int tmp;
    //int nvo = 0;
    int intorb = ndocc + nactv;

    /* Check if excitations are possible */
    if (str_virt != eosp.virt) return;
    if (str_virt == nvirt) return;
    if (str_virt == 0) return;
    if (str_actv != eosp.actv) return;
    if (str_docc != eosp.docc) return;
    
    /* Generate possible orbitals for replacement */
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[nvo] = i;
    //        nvo++;
    //    }
    //}
    /* Loop over occupied orbitals */
    newstr = str;
    switch (nrep) {
    case 1: /* Single replacements */
        for (i = str_virt - 1; i >= 0; i--) {
            for (j = intorb; j < (nvo + intorb); j++) {
                newstr = str;
                newstr.virtx[0] = str.virtx[0];
                newstr.virtx[1] = str.virtx[1];
                newstr.nvrtx = str.nvrtx;
                newstr.virtx[i] = scr[j];
                /* Ensure orbitals are ordered properly */
                if (newstr.virtx[0] > newstr.virtx[1] && newstr.virtx[1] != 0) {
                    tmp = newstr.virtx[0];
                    newstr.virtx[0] = newstr.virtx[1];
                    newstr.virtx[1] = tmp;
                }
                xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv, nvirt,
                                                    elec, elecs);
                xlist[*numx].io[0] = str.virtx[i];
                xlist[*numx].io[1] = 0;
                xlist[*numx].fo[0] = scr[j];
                xlist[*numx].fo[1] = 0;

                xlist[*numx].permx = pindex_single_rep(str.istr, str.virtx[i],
                                                       scr[j], elec);;
                (*numx)++;
            }
        }
        break;
    case 2: /* Double replacements */
        if ((nvirt - str_virt) < 2) return; /* Can't have 2x */
        for (i = str_virt - 1; i >= 1; i--) {
            for (j = intorb; j < (intorb + nvo - 1); j++) {
                for (k = (i - 1); k >= 0; k--) {
                    for (l = (j + 1); l < (intorb + nvo); l++) {
                        newstr = str;
                        newstr.virtx[0] = str.virtx[0];
                        newstr.virtx[1] = str.virtx[1];
                        newstr.virtx[i] = scr[j];
                        newstr.virtx[k] = scr[l];
                        newstr.nvrtx = str.nvrtx;
                        /* Ensure orbitals are ordered properly */
                        if (newstr.virtx[0] > newstr.virtx[1] &&
                            newstr.virtx[1] != 0) {
                            tmp = newstr.virtx[0];
                            newstr.virtx[0] = newstr.virtx[1];
                            newstr.virtx[1] = tmp;
                        }
                        xlist[*numx].index = occstr2address(newstr, eosp, ndocc, nactv,
                                                            nvirt, elec, elecs);
                        xlist[*numx].io[0] = str.virtx[0];
                        xlist[*numx].io[1] = str.virtx[1];
                        xlist[*numx].fo[0] = scr[j];
                        xlist[*numx].fo[1] = scr[l];
                        
                        xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                        (*numx)++;
                    }
                }
            }
        }
        break;
    default:
        printf(" nrep != 1 or 2\n");
        break;
    }
    return;
}

/*
 * generate_virt2doccactvx: generate replacements between VIRT and (DOCC,ACTV)
 * spaces.
 */
void generate_virt2doccactvx(int nrep, struct occstr str, int str_docc,
                             int str_actv, int str_virt, struct eospace eosp,
                             int ndocc, int nactv, int nvirt, int elec, int *scr,
                             struct xstr *xlist, int *numx, int nvo)
{
    struct occstr newstr;
    int elecs[20];
    int xtype = 0;
    long long ibyte = 0x00;
    long long xbyte = 0x00;
    int tmp;
    //int nvo = 0;
    int i, j, k, l;
    int intorb = ndocc + nactv;
    
    /* Check if excitations are possible */
    if (nrep < 2) return;
    if (str_virt == eosp.virt) return;
    if (str_docc == eosp.docc) return;
    if (str_actv == eosp.actv) return;
    if (abs(str_virt - eosp.virt) != nrep) return;


    /* Generate internal orbital byte */
    //for (i = 0; i < intorb; i++) {
    //    ibyte = ibyte + pow(2, i);
    //}
    /* Get possible excitations */
    //xbyte = ibyte ^ str.byte1;
    //nonzerobits(xbyte, intorb, scr);
    /* VIRT */
    //for (i = intorb + 1; i <= (intorb + nvirt); i++) {
    //    if (i != str.virtx[0] && i != str.virtx[1]) {
    //        scr[intorb + nvo] = i;
    //        nvo++;
    //    }
    //}

    /* Is this excitation VIRT -> (DOCC, ACTV) or VIRT <- (DOCC, ACTV)? */
    xtype = str_virt - eosp.virt;
    if (xtype > 0) {
        /* VIRT -> (DOCC, ACTV) */
        newstr = str;
        switch (nrep) {
        case 2: /* Double replacements */
            /* Loop over unoccupied DOCC orbitals, turning them "on" */
            for (i = 0; i < (ndocc - str_docc); i++) {
                /* Loop over unoccupied ACTV orbitals, turning them "on" */
                for (j = (ndocc - str_docc);
                     j < (intorb - str_actv - str_docc);
                     j++) {
                    newstr = str;
                    newstr.byte1 = str.byte1 + pow(2, (scr[i] - 1));
                    newstr.byte1 = newstr.byte1 + pow(2, (scr[j] - 1));
                    newstr.virtx[0] = 0;
                    newstr.virtx[1] = 0;
                    newstr.nvrtx = 0;
                    xlist[*numx].index =  occstr2address(newstr, eosp, ndocc,
                                                         nactv, nvirt, elec,
                                                         elecs);
                    xlist[*numx].io[0] = str.virtx[0];
                    xlist[*numx].io[1] = str.virtx[1];
                    xlist[*numx].fo[0] = scr[i];
                    xlist[*numx].fo[1] = scr[j];

                    xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                    (*numx)++;
                }
            }
            break;
        default:
            break;
        }
    } else {
        /* VIRT <- (DOCC, ACTV) */
        newstr = str;
        switch (nrep) {
        case 2: /* Double replacements */
            /* Loop over occupied DOCC orbitals, turning them "off" */
            for (i = 0; i < str_docc; i++) {
                /* Loop over occupied ACTV orbitals, turning them "off" */
                for (j = str_docc; j < str_docc + str_actv; j++) {
                    /* Loop over available VIRT orbitals, adding them */
                    for (k = intorb; k < (nvo + intorb - 1); k++) {
                        for (l = (k + 1); l < (nvo + intorb); l++) {
                            newstr = str;
                            newstr.byte1 = str.byte1 - pow(2, (str.istr[i] - 1));
                            newstr.byte1 = newstr.byte1 - pow(2, (str.istr[j] - 1));
                            newstr.nvrtx = str_virt + 2;
                            newstr.virtx[0] = scr[k];
                            newstr.virtx[1] = scr[l];
                            /* No need to check ordering */
                            xlist[*numx].index = occstr2address(newstr, eosp, ndocc,
                                                                nactv, nvirt, elec,
                                                                elecs);
                            xlist[*numx].io[0] = str.istr[i];
                            xlist[*numx].io[1] = str.istr[j];
                            xlist[*numx].fo[0] = scr[k];
                            xlist[*numx].fo[1] = scr[l];

                            xlist[*numx].permx = pindex_double_rep_str(newstr.istr,
                                                                   xlist[*numx].io[0],
                                                                   xlist[*numx].fo[0],
                                                                   xlist[*numx].io[1],
                                                                   xlist[*numx].fo[1],
                                                                   elec);
                            (*numx)++;
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    }
    return;
}


/*
 * generate_string_list: generate full *valid* alpha/beta string lists.
 */
void generate_string_list(struct eostring *strlist, int nstr, int orbs,
                          int elecs, int ndocc, int nactv, int xlvl,
                          struct eospace *eosp, int egrps)
{

        int count = 0; /* String counter */
        int *docc = NULL, *doccscr = NULL;
        int *actv = NULL, *actvscr = NULL;
        int *virt = NULL, *virtscr = NULL;
        int max_space_size = 20;
        int i;
        //docc = (int *) malloc(sizeof(int) * max_space_size);
        //doccscr = (int *) malloc(sizeof(int) * max_space_size);
        docc =  malloc(sizeof(int) * max_space_size);
        doccscr =  malloc(sizeof(int) * max_space_size);
        init_int_array_0(docc, max_space_size);
        init_int_array_0(doccscr, max_space_size);
        //actv = (int *) malloc(sizeof(int) * max_space_size);
        //actvscr = (int *) malloc(sizeof(int) * max_space_size);
        actv =  malloc(sizeof(int) * max_space_size);
        actvscr =  malloc(sizeof(int) * max_space_size);
        init_int_array_0(actv, max_space_size);
        init_int_array_0(actvscr, max_space_size);
        virt =  malloc(sizeof(int) * max_space_size);
        virtscr = malloc(sizeof(int) * max_space_size);
        init_int_array_0(virt, max_space_size);
        init_int_array_0(virtscr, max_space_size);

        /* Loop over electron groupings */
        for (i = 0; i < egrps; i++) {
                eosp[i].nstr  = count;
                eosp[i].start = count;
                compute_eostrings(strlist, &count, orbs, ndocc, nactv,
                                  eosp[i].docc, eosp[i].actv, eosp[i].virt,
                                  docc, actv, virt, doccscr, actvscr, virtscr);
                eosp[i].nstr = count - eosp[i].nstr;
        }
        free(actv);
        free(actvscr);
        free(virt);
        free(virtscr);
        free(docc);
        free(doccscr);
        return;
}

/*
 * get_available_orbital_list: get list of available orbitals for input str
 * Input:
 *  str = input occupation string
 *  intorb = internal orbitals
 *  nvirt = number of virtual orbitals
 * Output:
 *  orbsx = orbitals available for excitations
 *  nvo   = number of available orbitals.
 */
void get_available_orbital_list(struct occstr str, int intorb, int nvirt,
                                int *orbsx, int *nvo)
{
    /* Replacement byte and internal orbital byte */
    long long int xbyte = 0x00, ibyte = 0x00;
    int i = 0;
    for (i = 0; i < intorb; i++) {
        ibyte = ibyte + pow(2, i);
    }
    xbyte = ibyte ^ str.byte1;
    nonzerobits(xbyte, intorb, orbsx);
    *nvo = 0;
    for (i = intorb + 1; i <= (intorb + nvirt); i++) {
        if (i != str.virtx[0] && i != str.virtx[1]) {
            orbsx[intorb + (*nvo)] = i;
            (*nvo)++;
        }
    }
    return;
}

/*
 * get_eospace_detrange: get starting and final determinant indexes
 * for an eospace pairing.
 */
void get_eospace_detrange(int **pq, int npq, int pair, struct eospace *peosp,
                          int pegrps, struct eospace *qeosp, int qegrps,
                          int *start, int *final)
{
    int i, j;
    *start = 0;
    *final = 0;
    for (i = 0; i < pair; i++) {
        *start = *start + peosp[pq[i][0]].nstr * qeosp[pq[i][1]].nstr;
    }
    *final = *start + peosp[pq[i][0]].nstr * qeosp[pq[i][1]].nstr - 1;
    return;
}
/*
 * get_string_eospace: get eospace of input string.
 */
int get_string_eospace(struct occstr str, int ndocc, int nactv,
                       struct eospace *esp, int egrps)
{
    int spindx = 0; /* Index of eospace for input string */
    /* str eospace info */
    int str_docc = 0;
    int str_actv = 0;
    int str_virt = 0;
    int i;
    get_string_eospace_info(str, ndocc, nactv, &str_docc, &str_actv, &str_virt);
    for (i = 0; i < egrps; i++) {
        if (esp[i].docc == str_docc &&
            esp[i].actv == str_actv &&
            esp[i].virt == str_virt) {
            spindx = i;
            break;
        }
    }
    return spindx;
}
    
/*
 * occstr2address: compute the string index of given an occupation string.
 */
int occstr2address(struct occstr str, struct eospace eosp, int ndocc, int nactv,
                   int nvirt, int nelec, int *elecs)
{
    int i, j;
    int dstr = 0, astr = 0, vstr = 0;
    int daddr = 0, aaddr = 0, vaddr = 0;
    int addr = 0;
    
    /* Create electron string necessary for address look up. */
    nonzerobits(str.byte1, (ndocc + nactv), elecs);
    for (i = 0; i < str.nvrtx; i++) {
        //elecs[nelec - str.nvrtx + i] = str.virtx[i];
        elecs[nelec - str.nvrtx + i] = str.virtx[i] - ndocc - nactv;
    }
    /* Give orbitals in string "relative" values for space */
    /* DOCC space can be skipped. */
    for (i = eosp.docc; i < eosp.docc + eosp.actv; i++) {
        elecs[i] = elecs[i] - ndocc;
    }
    //for (i = eosp.docc + eosp.actv; i < nelec; i++) {
    //    elecs[i] = elecs[i] - ndocc - nactv;
    //}
    
#ifdef DEBUGGING
    for (i = 0; i < nelec; i++) {
        if (elecs[i] <= 0) {
            printf("Error! ");
            for (j = 0; j < nelec; j++) {
                printf(" %d", elecs[j]);
            }
            printf("\n");
            print_occstring(str, nelec, ndocc, nactv);
            return -10;
        }
    }
#endif
    
    if ((eosp.docc * eosp.actv * eosp.virt) != 0) {
        /* Electrons present in all spaces */
//            dstr = binomial_coef2(ndocc, eosp.docc);
        astr = binomial_coef3(nactv, eosp.actv);
        vstr = binomial_coef3(nvirt, eosp.virt);
        
        daddr = str_adrfind_fast(&(elecs[0]),eosp.docc,ndocc);
        aaddr = str_adrfind_fast(&(elecs[eosp.docc]),eosp.actv,nactv);
        vaddr = str_adrfind_fast(&(elecs[eosp.docc + eosp.actv]),eosp.virt,nvirt);
        
        addr = (daddr - 1) * astr * vstr + (aaddr - 1) * vstr + vaddr;
        
    } else if (eosp.actv != 0 && eosp.virt != 0) {
        /* (0, 1, 1) */
        astr = binomial_coef3(nactv, eosp.actv);
        vstr = binomial_coef3(nvirt, eosp.virt);
        
        aaddr = str_adrfind_fast(&(elecs[eosp.docc]),eosp.actv,nactv);
        vaddr = str_adrfind_fast(&(elecs[eosp.docc + eosp.actv]),eosp.virt,nvirt);
        
        addr = (aaddr - 1) * vstr + vaddr;
        
    } else if (eosp.docc != 0 && eosp.virt != 0) {
        /* (1, 0, 1) */
//            dstr = binomial_coef2(ndocc, eosp.docc);
        vstr = binomial_coef3(nvirt, eosp.virt);
        
        daddr = str_adrfind_fast(&(elecs[0]),eosp.docc,ndocc);
        vaddr = str_adrfind_fast(&(elecs[eosp.docc + eosp.actv]),eosp.virt,nvirt);
        
        addr = (daddr - 1) * vstr + vaddr;
        
    } else if (eosp.docc != 0 && eosp.actv != 0) {
        /* (1, 1, 0) */
//            dstr = binomial_coef2(ndocc, eosp.docc);
        astr = binomial_coef3(nactv, eosp.actv);
        
        daddr = str_adrfind_fast(&(elecs[0]),eosp.docc,ndocc);
        aaddr = str_adrfind_fast(&(elecs[eosp.docc]),eosp.actv,nactv);
        
        addr = (daddr - 1) * astr + aaddr;
        
    } else if (eosp.docc != 0) {
        /* (1, 0, 0) */
        dstr = binomial_coef3(ndocc, eosp.docc);
        
        addr = str_adrfind_fast(&(elecs[0]),eosp.docc,ndocc);
        
    } else if (eosp.actv != 0) {
        /* (0, 1, 0) */
        astr = binomial_coef3(nactv, eosp.actv);
        
        addr = str_adrfind_fast(&(elecs[eosp.docc]),eosp.actv,nactv);
        
    } else if (eosp.virt != 0) {
        /* (0, 0, 1) */
        vstr = binomial_coef3(nvirt, eosp.virt);
        
        addr = str_adrfind_fast(&(elecs[eosp.docc + eosp.actv]),eosp.virt,nvirt);
        
    }
    
    addr = addr + eosp.start - 1;
    if (addr > (eosp.start + eosp.nstr)) {
        printf("Error! addr = %d, eosp.start = %d, eosp.nstr = %d\n",
               addr, eosp.start, eosp.nstr);
        printf("eosp.docc = %d, eosp.actv = %d, eosp.virt = %d\n",
               eosp.docc, eosp.actv, eosp.virt);
        printf("Relative occupations: ");
        for (i = 0; i < (eosp.docc + eosp.actv + eosp.virt); i++) {
            printf(" %d", elecs[i]);
        }
        printf("\n");
    }
    return addr;
}

/*
 * remove_grt_xstr: remove replacements with indices greater than
 * val. (Modification of remove_leq_int in arrayutil.c)
 */
void remove_grt_xstr(int val, struct xstr *list, int *nx, struct xstr *scr)
{
    int num = 0;
    int i;
    if (*nx == 0) return;
    for (i = 0; i < *nx; i++) {
        /* If index is less than or equal to val add to scr */
        if (list[i].index <= val) {
            scr[num] = list[i];
            num++;
        }
    }
    /* Copy over elements */
    for (i = 0; i < num; i++) {
        list[i] = scr[i];
    }
    /* Reset *nx to refelect new list */
    *nx = num;
    return;
}

/*
 * remove_leq_xstr: remove replacements with indices less than, or equal to
 * val. (Modification of remove_leq_int in arrayutil.c)
 */
void remove_leq_xstr(int val, struct xstr *list, int *nx, struct xstr *scr)
{
    int num = 0;
    int i;
    if (*nx == 0) return;
    for (i = 0; i < *nx; i++) {
        /* If index is greater than val add to scr */
        if (list[i].index > val) {
            scr[num] = list[i];
            num++;
        }
    }
    /* Copy over elements */
    for (i = 0; i < num; i++) {
        list[i] = scr[i];
    }
    /* Reset *nx to refelect new list */
    *nx = num;
    return;
}

/*
 * setup_eostrings_compute: set up a electron, orbital, strings arrays for
 * compuation of the electron occupation strings.
 */
void setup_eostrings_compute(int *elecs, int *orbs, int *nstr, int *pegs,
                             int *nspcs, int docc_elec, int actv_elec,
                             int virt_elec, int dstr, int vstr, int astr,
                             int ndocc, int nactv, int vorbs)
{
        int typ = 0;

        init_int_array_0(elecs, 3);
        init_int_array_0(orbs, 3);
        init_int_array_0(nstr, 3);
        init_int_array_0(pegs, 3);
        
        if ((vstr * astr * dstr) == 0) {
                if (dstr == 0) typ = 1;
                if (astr == 0) typ = 2;
                if (vstr == 0) typ = 3;
                if (dstr == 0 && astr == 0) typ = 4;
                if (dstr == 0 && vstr == 0) typ = 5;
                if (astr == 0 && vstr == 0) typ = 6;
                if (dstr == 0 && astr == 0 && vstr == 0) return;
        }
        switch (typ) {
        case 0: // DOCC,ACTV,VIRT
                *nspcs = 3;
                elecs[0] = docc_elec;
                elecs[1] = actv_elec;
                elecs[2] = virt_elec;
                orbs[0] = ndocc;
                orbs[1] = nactv;
                orbs[2] = vorbs;
                nstr[0] = dstr;
                nstr[1] = astr;
                nstr[2] = vstr;
                pegs[0] = 0;
                pegs[1] = ndocc;
                pegs[2] = ndocc + nactv;
                break;
        case 1: // ACTV, VIRT
                *nspcs = 2;
                elecs[0] = actv_elec;
                elecs[1] = virt_elec;
                orbs[0] = nactv;
                orbs[1] = vorbs;
                nstr[0] = astr;
                nstr[1] = vstr;
                pegs[0] = ndocc;
                pegs[1] = ndocc + nactv;
                break;
        case 2: // DOCC, VIRT
                *nspcs = 2;
                elecs[0] = docc_elec;
                elecs[1] = virt_elec;
                orbs[0] = ndocc;
                orbs[1] = vorbs;
                nstr[0] = dstr;
                nstr[1] = vstr;
                pegs[0] = 0;
                pegs[1] = ndocc + nactv;
                break;
        case 3: // DOCC, ACTV
                *nspcs = 2;
                elecs[0] = docc_elec;
                elecs[1] = actv_elec;
                orbs[0] = ndocc;
                orbs[1] = nactv;
                nstr[0] = dstr;
                nstr[1] = astr;
                pegs[0] = 0;
                pegs[1] = ndocc;
                break;
        case 4: // VIRT
                *nspcs = 1;
                elecs[0] = virt_elec;
                orbs[0] = vorbs;
                nstr[0] = vstr;
                pegs[0] = ndocc + nactv;
                break;
        case 5: // ACTV
                *nspcs = 1;
                elecs[0] = actv_elec;
                orbs[0] = nactv;
                nstr[0] = astr;
                pegs[0] = ndocc;
                break;
        case 6: // DOCC
                *nspcs = 1;
                elecs[0] = docc_elec;
                orbs[0] = ndocc;
                nstr[0] = dstr;
                pegs[0] = 0;
                break;
        default:
                error_message(mpi_proc_rank,
                              "case != {0..6}", "compute_eostrings");
                break;
        }
        
        return;
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

/*
 * write_determinant_strpairs: write determinant alpha/beta string pairs for
 * valid determinants to file.
 */
void write_determinant_strpairs(FILE *fptr, int pstart, int pnstr, int qstart,
                                int qnstr, int *cnt,
                                struct eostring *pstr, struct eostring *qstr)
{
        FILE *ptr; /* New file pointer */
        int i, j;
        ptr = fptr;
        for (i = pstart; i < (pstart + pnstr); i++) {
                for (j = qstart; j < (qstart + qnstr); j++) {
                        fprintf(ptr, " %14d %14d\n",
                                pstr[i].index, qstr[j].index);
                        (*cnt)++;
                }
        }
        fptr = ptr;
        return;
}

/*
 * write_determinant_strpairs_dtlist: write determinant alpha/beta string
 * pairs for valid determinants to dtlist.
 */
void write_determinant_strpairs_dtlist(int pstart, int pnstr, int qstart,
                                       int qnstr, struct eostring *pstr,
                                       struct eostring *qstr, int aelec,
                                       int belec, int ndocc, int nactv,
                                       struct det *dlist, int *dptr,
                                       int cflag)
{
        
        int i, j;
        
        if (cflag == 1) {
                for (i = pstart; i < (pstart + pnstr); i++) {
                        for (j = qstart; j < (qstart + qnstr); j++) {
                                dlist[*dptr].astr = str2occstr(pstr[i].string,
                                                               aelec,
                                                               ndocc, nactv);
                                dlist[*dptr].bstr = str2occstr(qstr[j].string,
                                                               belec,
                                                               ndocc, nactv);
                                dlist[*dptr].cas = 1;
                                (*dptr)++;
                        }
                }
        } else {
                for (i = pstart; i < (pstart + pnstr); i++) {
                        for (j = qstart; j < (qstart + qnstr); j++) {
                                dlist[*dptr].astr = str2occstr(pstr[i].string,
                                                               aelec,
                                                               ndocc, nactv);
                                dlist[*dptr].bstr = str2occstr(qstr[j].string,
                                                               belec,
                                                               ndocc, nactv);
                                dlist[*dptr].cas = 0;
                                (*dptr)++;
                        }
                }
        }
        return;
}
