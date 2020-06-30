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
 *  w0_hndl = N+1 wavefunction handle
 *  dlen0   = number of N+1 determinants
 *  w1_hndl = N   wvaefunction handle
 *  dlen1   = number of N determinants
 *  v0_hndl = N+1 CI vectors handle
 *  v1_hndl = N   CI vectors handle
 *  str0    = N+1 alpha/beta strings
 *  str1    = N   alpha/beta strings
 *  spindx1 = alpha: 0,  beta: 1
 *  spindx2 = beta:  1, alpha: 0
 *  ninto0  = N+1 internal orbitals
 *  ninto1  = N   internal orbitals
 *  norbs   = total number of orbitals
 *  dysnst0 = dyson orbital N+1 states
 *  ndyst0  = number of N+1 states for dyson orbitals
 *  dysnst1 = dyson orbital N   states
 *  ndyst1  = number of N   states for dyson orbitals
 *  ndyorbs = number of dyson orbitals
 *  dyorb   = dyson orbital
 */
void compute_dyson_orbital(int w0_hndl, int dlen0, int w1_hndl, int dlen1,
                           int v0_hndl, int v1_hndl,
                           struct occstr *str0, struct occstr *str1,
                           int spindx1, int spindx2, int ninto0, int ninto1,
                           int norbs,  int *dysnst0, int ndyst0, int *dysnst1,
                           int ndyst1, int ndyorbs, double **dyorb)
{
#define BUFSIZE 1000
        int buflen = 0;         /* Buffer length */
        
        double **v0buff = NULL; /* N+1 electron ci buffer */
        double *v0data  = NULL;
        double **v1buff = NULL; /* N   electron ci buffer */
        double *v1data  = NULL;
        int **w0buff    = NULL; /* N+1 electron wf buffer */
        int  *w0data    = NULL;
        int **w1buff    = NULL; /* N   electron wf buffer */
        int  *w1data    = NULL;

        int v0_lo[2] = {0, 0};
        int v0_hi[2] = {0, 0};
        int v0_ld[1] = {0};
        int v0_rows = 0;
        int v0_cols = 0;        /* Rows and columns of v0 buffer */
        int v1_lo[2] = {0, 0};
        int v1_hi[2] = {0, 0};
        int v1_ld[1] = {0};
        int v1_rows = 0;
        int v1_cols = 0;
        
        int w0_lo[2] = {0, 0};
        int w0_hi[2] = {0, 0};
        int w0_ld[1] = {0};
        int w1_lo[2] = {0, 0};
        int w1_hi[2] = {0, 0};
        int w1_ld[1] = {0};
        int w1len = 0;
        
        int i0 = 0, i1 = 0;
        int i1max = 0;
        
        int ninto = 0; /* Number of internal orbitals */

        ninto = int_max(ninto0, ninto1);
        if (mpi_proc_rank == mpi_root) {
                printf(" Internal orbitals: %d\n", ninto);
        }

        /* Determine which blocks of data are locally owned for the N+1
         * electron vector. This will be the "static" buffer for which
         * we compute the dyson orbital contributions. */
        NGA_Distribution(v0_hndl, mpi_proc_rank, v0_lo, v0_hi);
        v0_cols = v0_hi[0] - v0_lo[0] + 1;
        v0_rows = v0_hi[1] - v0_lo[1] + 1;
        v0_ld[0]= v0_rows;
        /* Allocate local array */
        v0data = allocate_mem_double_cont(&v0buff, v0_rows, v0_cols);
        NGA_Get(v0_hndl, v0_lo, v0_hi, v0data, v0_ld);
        /* Allocate corresponding W0 array and get W0 data */
        w0data = allocate_mem_int_cont(&w0buff, 3, v0_rows);
        w0_lo[0] = v0_lo[1];
        w0_lo[1] = 0;
        w0_hi[0] = v0_hi[1];
        w0_hi[1] = 2;
        w0_ld[0] = 3;
        NGA_Get(w0_hndl, w0_lo, w0_hi, w0data, w0_ld);
        
        buflen = BUFSIZE; /* Set buffer length */

        /* Get state number for v1 by getting number of columns
         * on this process. (Array is distributed by row.) */
        NGA_Distribution(v1_hndl, mpi_proc_rank, v1_lo, v1_hi);
        v1_cols = v1_hi[0] - v1_lo[0] + 1;
        /* Set values for v1_lo and v1_hi */
        v1_lo[0] = 0;
        v1_lo[1] = 0;
        v1_hi[0] = v1_cols - 1;
        v1_hi[1] = dlen1 - 1;
        v1_rows = buflen;
        v1_ld[0] = v1_rows;
        /* Allocate local array buffer (buflen x v1_cols) */
        v1data = allocate_mem_double_cont(&v1buff, v1_rows, v1_cols);

        /* Allocate W1 array */
        w1data = allocate_mem_int_cont(&w1buff, 3, v1_rows);
        w1_lo[1] = 0;
        w1_hi[1] = 2;
        w1_ld[0] = 3;

        GA_Sync();

        for (i1 = 0; i1 < dlen1; i1 += buflen) {

                i1max = int_min((i1 + buflen - 1), (dlen1 - 1));
                
                /* Get patch of V1 vectors and corresponding W1 info */
                v1_lo[1] = i1;
                v1_hi[1] = i1max;
                v1_rows = v1_hi[1] - v1_lo[1] + 1;
                NGA_Get(v1_hndl, v1_lo, v1_hi, v1data, v1_ld);
                w1_lo[0] = i1;
                w1_hi[0] = i1max;
                //w1len = i1max - i1 + 1;
                NGA_Get(w1_hndl, w1_lo, w1_hi, w1data, w1_ld);
                
                compute_det_contributions(w0buff,  v0buff, v0_rows,
                                          v0_cols, w1buff, v1buff,
                                          v1_rows, v1_cols, spindx1, spindx2,
                                          str0, str1, dysnst0, ndyst0,
                                          dysnst1, ndyst1, dyorb, ninto);
        }
        GA_Sync();
        printf("Finished computing dyson orbitals.\n");
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
        
        printf("Evaluating det contributions.\n");

        /* Loop over determinants. N electrons must be in the same slots for
         * both determinants. The N+1 electron, in a unique slot, is the MO
         * index to which the matrix element contributes to the dyson orbital. */
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
                                        printf("%d: dyst1[%d] = %d\n",
                                               mpi_proc_rank, l, dyst1[l]);
                                        dyorb[dyind][(orbindx - 1)] =
                                                v0[dyst0[k]-1][i]*
                                                   v1[dyst1[l]-1][j];
                                        dyind++;
                                }
                        }
                }
        }

        return;
};
