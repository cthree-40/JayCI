// File: citruncate.c
/*
 * Truncate a Full CI expansion.
 */
#include <stdio.h>
#include <stdlib.h>
#include "pjayci_global.h"
#include "ioutil.h"
#include "iminmax.h"
#include "binary.h"
#include "errorlib.h"
#include "arrayutil.h"
#include "allocate_mem.h"
#include "combinatorial.h"
#include "straddress.h"
#include "binarystr.h"
#include "citruncate.h"
#include <mpi.h>

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
        int i;
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
        int typ = 0; /* 0: normal, 1: No DOCC, 2: No CAS, 3: No Virt,
                      * 4: No DOCC+CAS, 5: No DOCC+VIRT, 6: No CAS+VIRT */
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
        int count = 0;
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
        docc = (int *) malloc(sizeof(int) * max_space_size);
        doccscr = (int *) malloc(sizeof(int) * max_space_size);
        init_int_array_0(docc, max_space_size);
        init_int_array_0(doccscr, max_space_size);
        actv = (int *) malloc(sizeof(int) * max_space_size);
        actvscr = (int *) malloc(sizeof(int) * max_space_size);
        init_int_array_0(actv, max_space_size);
        init_int_array_0(actvscr, max_space_size);
        virt = (int *) malloc(sizeof(int) * max_space_size);
        virtscr = (int *) malloc(sizeof(int) * max_space_size);
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
        int array[4];
        
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
