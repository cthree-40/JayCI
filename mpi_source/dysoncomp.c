// File: dysoncomp.c
/*
 * Functions to compute dyson orbital.
 */
#include <stdio.h>
#include <stdlib.h>
#include "mpi_utilities.h"
#include "errorlib.h"
#include "allocate_mem.h"
#include "ioutil.h"
#include "binary.h"
#include "bitutil.h"
#include "binarystr.h"
#include "citruncate.h"
#include "action_util.h"
#include "iminmax.h"
#include "dysoncomp.h"

#include <ga.h>
#include <macdecls.h>

#include <mpi.h>

/* -- OpenMP options -- */
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif
/* -------------------- */

/*
 * build_ppo_triples: build (p0, p1, o) triples for dyson evaluation.
 */
int build_ppo_triples(int **ppo, struct occstr *str0, int nstr0,
		      struct occstr *str1, int nstr1, int ninto)
{
    int nppo = 0;
    int oindex = 0;
    int i, j;
    for (i = 0; i < nstr0; i++) {
	for (j = 0; j < nstr1; j++) {
	    
	    oindex = comparestrings_dyson(str0[i], str1[j], ninto);
	    if (oindex == 0) continue;
	    ppo[nppo][0] = i;
	    ppo[nppo][1] = j;
	    ppo[nppo][2] = oindex;
	    nppo++;
	}
    }
    return nppo;
}

/*
 * comparestrings_dyson: compare the strings of N and N+1 determinants
 * to get the index of orbital contribution.
 */
int comparestrings_dyson(struct occstr str0, struct occstr str1, int ninto)
{
        int oindex = 0;      /* Orbital index */
        int numxv  = 0;      /* Number of virtual orbital differences */
        int numxc  = 0;      /* Number of cas orbital differences */
        int ifo[2] = {0, 0}; /* Initial, final virtual orbitals */
        long long int diffsb = 0; /* difference 64-bit byte. */

        /* Test virtual orbital blocks for differences. There can be no
         * more than one difference. */
        numxv = compute_virt_diffs(str0, str1);
        if (numxv > 1) return oindex;

        /* If there is one virtual orbital difference, there can be no
         * differences in CAS byte. */
        if (numxv == 1) {
                if ((str0.byte1 ^ str1.byte1) != 0x00) return oindex;
                virtdiffs_single_rep(str0.virtx, str1.virtx, ifo);
                oindex = ifo[0];
                return oindex;
        }

        /* If here, there are no differences in the virtual orbital blocks.
         * The difference must be in the active space bytes. Again, there can
         * only be one difference. */
        numxc = ndiffbytes(str0.byte1, str1.byte1, ninto, &diffsb);
        if (numxc > 1) return oindex;

        /* There must be one difference */
        if (numxc == 1) {
                nonzerobits(diffsb, ninto, &oindex);
        }
                
        return oindex;
}

/*
 * compute_dyson_orbital: compute the dyson orbital between electronic
 * states of N+1 and N electron wavefuntions by comparing alpha/beta strings.
 * Input:
 *
 */
void compute_dyson_orbital(int v0_hndl, int v1_hndl, int w0_hndl, int w1_hndl,
			   struct occstr *pstr0, struct eospace *peosp0, int npe0,
			   struct occstr *qstr0, struct eospace *qeosp0, int nqe0,
			   struct occstr *pstr1, struct eospace *peosp1, int npe1,
			   struct occstr *qstr1, struct eospace *qeosp1, int nqe1,
			   int norbs, int ndocc, int nactv, int ndyst0,
			   int *dysnst0, int ndyst1, int *dysnst1, int ndets0,
			   int ndets1, int sp1, int sp2)
{
#define MAXBUFFER 10000
    /* Local V and W buffers */
    double **v0 = NULL, *v0data = NULL;
    int v0cols = 0, v0rows = 0;
    double **v1 = NULL, *v1data = NULL;
    int v1cols = 0, v1rows = 0;
    int **w0 = NULL, *w0data = NULL;
    int **w1 = NULL, *w1data = NULL;
    /* GA dimensions for buffers */
    int v0_lo[2] = {0, 0}, v0_hi[2] = {0, 0}, v0_ld[1] = {0};
    int v1_lo[2] = {0, 0}, v1_hi[2] = {0, 0}, v1_ld[1] = {0};
    int w0_lo[2] = {0, 0}, w0_hi[2] = {0, 0}, w0_ld[1] = {0};
    int w1_lo[2] = {0, 0}, w1_hi[2] = {0, 0}, w1_ld[1] = {0};

    int buflen;
    int i, imax, j, jmax;
    int nvirt;

    nvirt = norbs - ndocc - nactv;
    
    /* Compute dimensions of V0. Allocate array. Get buffer. */
    NGA_Distribution(v0_hndl, mpi_proc_rank, v0_lo, v0_hi);
    v0rows = v0_hi[1] - v0_lo[1] + 1;
    v0cols = v0_hi[0] - v0_lo[0] + 1;
    v0_ld[0]= v0rows;
    v0data = allocate_mem_double_cont(&v0, v0rows, v0cols);
    NGA_Get(v0_hndl, v0_lo, v0_hi, v0data, v0_ld);
    /* Get W0 wavefunction information */
    w0_lo[0] = v0_lo[1];
    w0_lo[1] = 0;
    w0_hi[0] = v0_hi[1];
    w0_hi[1] = 2;
    w0data = allocate_mem_int_cont(&w0, 3, v0rows);
    NGA_Get(w0_hndl, w0_lo, w0_hi, w0data, w0_ld);

    buflen = MAXBUFFER; /* Set buffer length */

    /* Get state number for v1 by getting number of columns
     * on this process. (Array is distributed by row.) */
    NGA_Distribution(v1_hndl, mpi_proc_rank, v1_lo, v1_hi);
    v1cols = v1_hi[0] - v1_lo[0] + 1;
    v1rows = buflen;
    v1_ld[0] = v1rows;
    v1_lo[0] = 0;
    v1_hi[0] = 0;
    v1data = allocate_mem_double_cont(&v1, v1rows, v1cols);
    w1_lo[1] = 0;
    w1_hi[1] = 2;
    w1data = allocate_mem_int_cont(&w1, 3, v1rows);
    
    GA_Sync();
    
    for (i = 0; i < ndets1; i += buflen) {
	imax = int_min((i + buflen - 1), (ndets1 - 1));
	v1_lo[1] = i;
	v1_hi[1] = imax;
	v1rows = v1_hi[1] - v1_lo[1] + 1;
	NGA_Get(v1_hndl, v1_lo, v1_hi, v1data, v1_ld);
	w1_lo[0] = i;
	w1_hi[0] = imax;
	NGA_Get(w1_hndl, w1_lo, w1_hi, w1data, w1_ld);

	compute_det_contributions2(w0, v0, v0rows, v0cols, w1, v1, v1rows,
				   v1cols, pstr0, peosp0, qstr0, qeosp0,
				   pstr1, peosp1, qstr1, qeosp1, ndocc,
				   nactv, norbs, sp1, sp2, nvirt);
    }

    GA_Sync();
    deallocate_mem_cont_int(&w0, w0data);
    deallocate_mem_cont_int(&w1, w1data);
    deallocate_mem_cont(&v0, v0data);
    deallocate_mem_cont(&v1, v1data);
    return;
}

/*
 * compute_det_contributions2: compute determinant contributions to dyson
 * orbitals between two buffers v0 and v1.
 */
void compute_det_contributions2(int **w0, double **v0, int v0rows, int v0cols,
				int **w1, double **v1, int v1rows, int v1cols,
				struct occstr *str0, struct eospace *eosp0, int ne0,
				struct occstr *str1, struct eospace *eosp1, int ne1,
				int ndocc, int nactv, int norbs, int nelec,
				int sp1, int sp2, int nvirt)
{
    int elecx[20];
    struct occstr newstr;
    int naddr, neosp1;
    int neosp0;
    int nstrd0, nstra0, nstrv0;
    int i, j;
    
    for (i = 0; i < v0rows; i++) {

	newstr = str0[w0[i][sp1]];
	neosp0 = get_string_eospace(&newstr, ndocc, nactv, eosp0, ne0);
	nstrd0 = eosp0[neosp].docc;
	nstra0 = eosp0[neosp].actv;
	nstrv0 = eosp0[neosp].virt;
	/* Remove internal orbitals */
	for (j = 0; j < (nstrd0 + nstra0); j++) {
	    newstr.byte1 = newstr.byte1 - pow(2, (newstr.istr[j] - 1));
	    neosp1 = get_string_eospace(&newstr, ndocc, nactv, eosp1, ne1);
	    naddr = occstr2address(&newstr, eosp1[neosp1], ndocc, nactv, 
				   nvirt, nelec, elecx);
	    
	    newstr.byte1 = newstr.byte1 + pow(2, (newstr.istr[j] - 1));
	}
	/* Remove external orbitals */
	for (j = 0; j < nstrv0; j++) {
	    newstr.virtx[j] = 0;
	    if (nstrv == 2 && j == 0) {
		newstr.virtx[0] = newstr.virtx[1];
		newstr.virtx[1] = 0;
	    }
	    

	}
    }
    return;
}

/*
 * compute_det_contributions: compute determinant contributions to dyson
 * orbitals.
 * Input:
 *  w0      = p, q, CAS
 *  v0      = CI vectors
 *  v0_rows = rows of CI vectors
 *  v0_cols = columns of CI vectors
 *  w1      = p, q, CAS
 *  v1      = CI vectors
 *  v1_rows = rows of CI vectors
 *  v1_cols = columns of CI vectors
 *  sp      = spin index (0 = alpha, 1 = beta)
 *  spx     = opposite spin index (1 = beta, 0 = alpha)
 *  str0    = p/q strings
 *  str1    = p/q strings
 *  dyst0   = dyson orbital N+1 electron states
 *  ndyst0  = number of dyson orbital N+1 electron states
 *  dyst1   = dyson orbital N electron states
 *  ndyst1  = number of dyson orbital N electron states
 *  ninto   = number of internal orbitals
 * Output:
 *  dyorb   = dyson orbitals
 */
void compute_det_contributions(int **w0, double **v0, int v0_rows, int v0_cols,
                               int **w1, double **v1, int v1_rows, int v1_cols,
                               int sp, int spx,
                               struct occstr *str0, struct occstr *str1,
                               int *dyst0, int ndyst0, int *dyst1, int ndyst1,
                               double **dyorb, int ninto)
{
    int i = 0, j = 0;
    int k = 0, l = 0;
    int dyind = 0;      /* dyson orbital index */
    int ndiff = 0;      /* number of differences */
    int orbindx = 0;    /* orbital index of contribution */
    
    /* Loop over determinants. N electrons must be in the same slots for
     * both determinants. The N+1 electron, in a unique slot, is the MO
     * index to which the matrix element contributes to the dyson orbital. */
    
#pragma omp parallel					      \
    default(none)					      \
    shared(w0,v0,v0_rows,v0_cols,			      \
	   w1,v1,v1_rows,v1_cols,			      \
	   sp, spx, str0, str1, dyst0, ndyst0, dyst1, ndyst1, \
	   dyorb,ninto)					      \
    private(i, j, orbindx, ndiff, dyind, k, l)
    {
#pragma omp for schedule(runtime)
	for (i = 0; i < v0_rows; i++) {
	    for (j = 0; j < v1_rows; j++) {
		/* If the other spin strings are different, skip. */
		if (w0[i][spx] != w1[j][spx]) continue;
		orbindx = comparestrings_dyson(str0[w0[i][sp]],
					       str1[w1[j][sp]], ninto);
		if (orbindx == 0) continue;
		
		/* Compute contributions to orbital (orbindx) for all
		 * CI vectors involved in dyson orbitals */
		dyind = 0;
		for (k = 0; k < ndyst0; k++) {
		    for (l = 0; l < ndyst1; l++) {
			dyorb[dyind][(orbindx - 1)] =
			    dyorb[dyind][(orbindx - 1)] +
			    v0[dyst0[k]][i]*
			    v1[dyst1[l]][j];
			dyind++;
		    }
		}
	    }
        }
    }
    return;
}
