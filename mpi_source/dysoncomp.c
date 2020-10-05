// File: dysoncomp.c
/*
 * Functions to compute dyson orbital.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi_utilities.h"
#include "errorlib.h"
#include "allocate_mem.h"
#include "arrayutil.h"
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
 * compute_dyson_orbital_b: compute the dyson orbital between electronic
 * states of N+1 and N electron wavefuntions by comparing beta strings.
 * Input:
 *
 */
void compute_dyson_orbital_b(int v0_hndl, int v1_hndl, int w0_hndl, int w1_hndl,
                             struct occstr *pstr0, struct eospace *peosp0, int npe0,
                             struct occstr *qstr0, struct eospace *qeosp0, int nqe0,
                             struct occstr *pstr1, struct eospace *peosp1, int npe1,
                             struct occstr *qstr1, struct eospace *qeosp1, int nqe1,
                             int **pq1, int npq1,
                             int norbs, int ndocc, int nactv, int ndyst0,
                             int *dysnst0, int ndyst1, int *dysnst1, int ndets0,
                             int ndets1, int **strcont, int nbelec0,
                             double **dyorb)
{
#define MAXBUFFER 10000
    /* (N+1) and (N) Vector arrays */
    double **v0 = NULL,  *v01d = NULL; 
    double **v1 = NULL,  *v11d = NULL; 
    int v0_lo[2], v0_hi[2], v0_ld[1];
    int v1_lo[2], v1_hi[2], v1_ld[1];
    int v0_rows, v0_cols;
    int v1_rows, v1_cols;
    /* (N+1) and (N) wavefunction arrays */
    int **w0 = NULL, *w01d = NULL;
    int **w1 = NULL, *w11d = NULL;
    int w0_lo[2], w0_hi[2], w0_ld[1];
    int w1_lo[2], w1_hi[2], w1_ld[1];
    /* Local variables */
    int p0, q0; /* w0 values */
    int p1, q1; /* w1 values */
    int *j0indx = NULL;  /* Determinant index of |p0,q0> */
    int *j1indx = NULL;  /* Determinant index of |p1,q1> */
    int *o1indx = NULL;  /* Orbital difference of <p0,q0|p1,q1> */
    int nj1 = 0;
    int **vindx = NULL, *vindx1d = NULL;
    int i, j, k;
    int cnt;
    
    /* Get local distribution of V0. Allocate array and read in values. */
    NGA_Distribution(v0_hndl, mpi_proc_rank, v0_lo, v0_hi);
    v0_rows = v0_hi[1] - v0_lo[1] + 1;
    v0_cols = v0_hi[0] - v0_lo[0] + 1;
    v0_ld[0]= v0_rows;
    v01d = allocate_mem_double_cont(&v0, v0_rows, v0_cols);
    NGA_Get(v0_hndl, v0_lo, v0_hi, v01d, v0_ld);
    /* Get corresponding wavefunction information */
    w01d = allocate_mem_int_cont(&w0, 3, v0_rows);
    w0_lo[0] = v0_lo[1];
    w0_lo[1] = 0;
    w0_hi[0] = v0_hi[1];
    w0_hi[1] = 2;
    w0_ld[0] = 3;
    NGA_Get(w0_hndl, w0_lo, w0_hi, w01d, w0_ld); 

    /* Allocate J index array and j-index array to grab V1 values from
     * global array. */
    nj1 = v0_rows * nbelec0;
    j0indx = malloc(sizeof(int) * nj1);
    j1indx = malloc(sizeof(int) * nj1);
    o1indx = malloc(sizeof(int) * nj1);
    vindx1d = allocate_mem_int_cont(&vindx, 2, (ndyst1 * nj1));

    /* Allocate local buffer of V1 */
    v11d = allocate_mem_double_cont(&v1, nj1, ndyst1);
    
    /* Loop over w0 = (p,q), finding q'/p' in w1 */
    cnt = 0;
    for (i = 0; i < v0_rows; i++) {
        p0 = w0[i][0];
        q0 = w0[i][1];
        /* Loop over q' */
        for (j = 0; j < nbelec0; j++) {
            p1 = p0;
            q1 = strcont[q0][j]; /* New string */
            j1indx[cnt] = string_info_to_determinant(p1, q1, peosp1, npe1,
                                                     qeosp1, nqe1, pq1, npq1);
            o1indx[cnt] = qstr0[q0].istr[j] - 1;
            j0indx[cnt] = i;

            cnt++;
        }
    }

    /* Get indexes for all columns of J that we need */
    set_ga_det_indexes_spec(j1indx, cnt, ndyst1, dysnst1, vindx);

    /* Get V1 elements */
    NGA_Gather(v1_hndl, v11d, vindx, (cnt * ndyst1));

    /* Loop over contributions */
    for (i = 0; i < cnt; i++) {

        for (j = 0; j < ndyst0; j++) {
            for (k = 0; k < ndyst1; k++) {
                dyorb[j * ndyst1 + k][o1indx[i]] += v0[j][j0indx[i]] *
                    v1[k][i];
            }
        }
    }
    
               
    /* Deallocate arrays */
    deallocate_mem_cont(&v0, v01d);
    deallocate_mem_cont(&v1, v11d);
    deallocate_mem_cont_int(&w0, w01d);
    deallocate_mem_cont_int(&vindx, vindx1d);
    free(j0indx);
    free(j1indx);
    free(o1indx);
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

/*
 * generate_strcontlist: generate contribution list for each string.
 * Input:
 *  str    = string list
 *  nstr   = number of alpha/beta strings
 *  eosp0  = electron orbital space list (N+1 electrons)
 *  ne0    = number of electron orbital spaces (N+1 electrons)
 *  ndocc  = number of docc orbitals
 *  nactv  = number of actv orbitals
 *  nvirt  = number of virt orbitals
 *  nelec1 = number of N-electron alpha/beta electrons
 *  eosp1  = electron orbital space list (N-electron)
 *  ne1    = number of electron orbital spaces (N-electron)
 * Output:
 *  strcont = string contribution list
 */
void generate_strcontlist(struct occstr *str, int nstr, struct eospace *eosp0,
			  int ne0, int ndocc, int nactv, int nvirt,
			  int **strcont, int nelec1, struct eospace *eosp1,
			  int ne1)
{
    int elecx[20];           /* Scratch electron array */
    struct occstr newstr;    /* New string */
    int naddr = 0;           /* New string address */
    int neospx = 0;          /* New string electron orbital space index */
    int nelec0 = 0;          /* N + 1 electrons */
    int cnt = 0;             
    int i, j;
    
    nelec0 = nelec1 + 1;

    /* Loop over strings */
    for (i = 0; i < nstr; i++) {
	cnt = 0; // Counter for new strings
        /* Loop over internal orbital electrons, removing them */
	for (j = 0; j < (nelec0 - str[i].nvrtx); j++) {
	    newstr.byte1 = str[i].byte1 - pow(2, (str[i].istr[j] - 1));
	    newstr.nvrtx = str[i].nvrtx;
	    newstr.virtx[0] = str[i].virtx[0];
	    newstr.virtx[1] = str[i].virtx[1];

#ifdef DEBUGGING
	    //printf(" ORIG: ");
	    //print_occstring(&(str[i]), nelec0, ndocc, nactv);
	    //printf(" NEW:  ");
	    //print_occstring(&(newstr), nelec1, ndocc, nactv);
#endif
	    neospx = get_string_eospace(&newstr, ndocc, nactv, eosp1, ne1);
	    strcont[i][cnt] = occstr2address(&newstr, eosp1[neospx], ndocc,
					     nactv, nvirt, nelec1, elecx);
	    cnt++;
	}
	/* Loop over external orbital electrons, removing them */
	for (j = 0; j < str[i].nvrtx; j++) {
            newstr.byte1 = str[i].byte1;
            newstr.nvrtx = str[i].nvrtx - 1;
            newstr.virtx[0] = str[i].virtx[0];
            newstr.virtx[1] = str[i].virtx[1];
            newstr.virtx[j] = 0;
            if (j == 0 && str[i].nvrtx == 2) {
                newstr.virtx[0] = newstr.virtx[1];
                newstr.virtx[1] = 0;
            }
#ifdef DEBUGGING
	    //printf(" ORIG: ");
	    //print_occstring(&(str[i]), nelec0, ndocc, nactv);
	    //printf(" NEW:  ");
	    //print_occstring(&(newstr), nelec1, ndocc, nactv);
#endif
	    neospx = get_string_eospace(&newstr, ndocc, nactv, eosp1, ne1);
	    strcont[i][cnt] = occstr2address(&newstr, eosp1[neospx], ndocc,
					     nactv, nvirt, nelec1, elecx);
	    cnt++;
            
	}
        //printf("Count: %d\n", cnt);
    }
    return;
}

/*
 * string_info_to_determinant: compute the determinant index given
 * the p and q string information.
 * Input:
 *  pval    = p string
 *  qval    = q string
 *  peosp   = alpha electron orbital spaces
 *  pegrps  = number of alpha electron orbital spaces
 *  qeosp   = beta  electron orbital spaces
 *  qegrps  = number of beta  electron orbital spaces
 *  pq      = (p,q)-space pairings
 *  npq     = number of (p,q)-space pairings
 * Output:
 *  detindx = determinant index in expansion.
 */
int string_info_to_determinant(int pval, int qval, struct eospace *peosp,
                               int pegrps, struct eospace *qeosp, int qegrps,
                               int **pq, int npq)
{
        int i, j, k;
        int jstart, jmax, kstart, kmax;
        int iloc, jloc;
        int cntr = 0; /* Counter */
        /* Find space combo for p,q */
        for (i = 0; i < pegrps - 1; i++) {
                if (pval >= peosp[i].start && pval < peosp[i + 1].start) break;
        }
        iloc = i;
        for (j = 0; j < qegrps - 1; j++) {
                if (qval >= qeosp[j].start && qval < qeosp[j + 1].start) break;
        }
        jloc = j;
        for (i = 0; i < npq; i++) {
                if (pq[i][0] == iloc && pq[i][1] == jloc) break;
        }
        iloc = i;
        /* p,q pairing is found: i */
        for (i = 0; i < iloc; i++) {
                cntr = cntr + peosp[(pq[i][0])].nstr * qeosp[(pq[i][1])].nstr;
        }
        /* Search for index */
        jstart = peosp[(pq[iloc][0])].start;
        kstart = qeosp[(pq[iloc][1])].start;
        jmax   = jstart + peosp[(pq[iloc][0])].nstr;
        kmax   = kstart + qeosp[(pq[iloc][1])].nstr;
        for (j = jstart; j < jmax; j++) {
                for (k = kstart; k < kmax; k++) {
                        if (j == pval && k == qval) {
                                k = kmax;
                                j = jmax;
                                break;
                        }
                        cntr++;
                }
        }
        return cntr;
}

/*
 * set_ga_det_indexes_spec: set the array of indices to gather from global array.
 * Input:
 *  jindx   = list of row numbers in vector V
 *  num     = number of row numbers
 *  cols    = number of columns of V
 *  colindx = column indices
 * Output:
 *  vindx   = indices of global array V to gather
 */
void set_ga_det_indexes_spec(int *jindx, int num, int cols, int *colindx,
                             int **vindx)
{
    int cntr = 0; /* Index counter */
    /* Arrays in global arrays are stored like: [column][row] */
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < num; j++) {
            vindx[cntr][0] = colindx[i];
            vindx[cntr][1] = jindx[j];
            cntr++;
        }
    }
    return;
}
