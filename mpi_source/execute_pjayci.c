// File: execute_pjayci.c
/*
 * Execute parallel, determinant-based CI algorithm.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pjayci_global.h"
#include "errorlib.h"
#include "arrayutil.h"
#include "allocate_mem.h"
#include "ioutil.h"
#include "combinatorial.h"
#include "moindex.h"
#include "abecalc.h"
#include "mpi_utilities.h"
#include "binarystr.h"
#include "citruncate.h"
#include "iminmax.h"
#include "action_util.h"
#include "execute_pjayci.h"
#include <ga.h>
#include <macdecls.h>
#include <mpi.h>

/*
 * execute_pjayci: execute parallel CI.
 */
int execute_pjayci ()
{
        int error = 0;     /* Error flag */

        int electrons = 0; /* total number of electrons */
        int orbitals = 0;  /* total number of orbitals */
        int nfrzc = 0;     /* number of frozen core orbitals */
        int nfrzv = 0;     /* number of frozen virtual orbitals */
        int nactv = 0;     /* number of CAS active orbitals */
        int ndocc = 0;     /* number of CAS doubly-occupied orbitals */
        int intorb = 0;    /* number of CI internal orbitals: ndocc + nactv */
        int xlvl = 0;      /* CI expansion excitation level. must be <= 2 */
        int printlvl = 0;  /* print level */
        int printwvf = 0;  /* print wavefunction info */

        int aelec = 0;    /* Alpha electrons */
        int belec = 0;    /* Beta  eletrons */
        int ci_aelec = 0; /* CI alpha electrons */
        int ci_belec = 0; /* CI beta  electrons */
        int ci_orbs = 0;  /* active CI orbitals */

        struct occstr  *pstrings;    /* Alpha electron strings */
        struct occstr  *qstrings;    /* Beta  electron strings */
        struct eospace *peospace;    /* Alpha electron string spaces */
        struct eospace *qeospace;    /* Beta  electron string spaces */
        int **pq_space_pairs = NULL; /* Valid (p,q) space pairs. */
        int num_pq = 0;              /* Number of valid (p,q) space pairs */
        int pegrps = 0;              /* Alpha electron occupation groups */
        int qegrps = 0;              /* Beta  electron occupation groups */
        int pstr_len = 0;            /* Number of alpha strings */
        int qstr_len = 0;            /* Number of beta  strings */
        int dtrm_len = 0;            /* Total number of determinants */

        int m1len = 0;                  /* Number of 1-e integrals */
        int m2len = 0;                  /* Number of 2-e integrals */
        double *moints1 = NULL;         /* 1-e integrals */
        double *moints2 = NULL;         /* 2-e integrals */
        double frzcore_e = 0.0;         /* Frozen core energy */
        double nucrep_e  = 0.0;         /* Nuclear repulsion energy */
        char moflname[FLNMSIZE] = {""}; /* SIFS integral filename */
        int itype = 1;                  /* Integral type. Always 1 */

        int hmat_chunk  = 0;              /* Size of Hv=c chunks per process */
        int row_lwrbnd = 0;               /* Lower bound of Hv=c rows */
        int row_uprbnd = 0;               /* Upper bound of Hv=c rows */
        int pq_start = 0, pq_final = 0;   /* p,q pairings start/final */
        int pstart = 0, pfinal = 0;       /* p string start/final */
        int qstart = 0, qfinal = 0;       /* q string start/final */
        
        double *hdgls_local = NULL;  /* <i|H|i> local elements */
        int hdgls_handle = 0;        /* <i|H|i> global array handle */
        int hdgls_dims[1];           /* <i|H|i> global array dimensions */
        int hdgls_chunk[1];          /* <i|H|i> global array blocking */
        int hdgls_lo[1];             /* <i|H|i> first index (i) */
        int hdgls_hi[1];             /* <i|H|i> final index (i) */
        
        /* Read in the &general namelist. Ensure that the expansion's
         * CAS space is not greater than 64 orbitals. */
        if (mpi_proc_rank == mpi_root) {
                printf("Reading &general input\n");
                readgeninput(&electrons, &orbitals, &nfrzc, &ndocc, &nactv,
                             &xlvl, &nfrzv, &printlvl, &printwvf, &error);
                if (error != 0) {
                        error_flag(mpi_proc_rank, error, "execute_pjayci");
                        return error;
                }
                if ((ndocc + nactv) > 64) {
                        error = ndocc + nactv;
                        error_flag(mpi_proc_rank, error, "execute_pjayci");
                        error_message(mpi_proc_rank,
                                      "CI-active orbital number must be < 64.\n",
                                      "execute_pjayci");
                        return error;
                }
        }
        /* Broadcast &general values. These are needed on all processes. */
        MPI_Bcast(&electrons, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&orbitals,  1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzc,     1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&ndocc,     1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nactv,     1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&xlvl,      1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzv,     1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&printlvl,  1, MPI_INT, mpi_root, MPI_COMM_WORLD);

        /* Get number of alpha/beta electrons */
        abecalc(electrons, &aelec, &belec);
        
        /* Generate pstrings and qstrings. */
        compute_ci_elecs_and_orbitals(aelec, belec, orbitals, nfrzc, nfrzv,
                                      &ci_aelec, &ci_belec, &ci_orbs);
        pstr_len = compute_stringnum(ci_orbs, ci_aelec, ndocc, nactv, xlvl);
        qstr_len = compute_stringnum(ci_orbs, ci_belec, ndocc, nactv, xlvl);
        pstrings = allocate_occstr_arrays(pstr_len);
        qstrings = allocate_occstr_arrays(qstr_len);
        peospace = allocate_eospace_array(ci_aelec, ci_orbs, ndocc, nactv, xlvl,
                                       &pegrps);
        qeospace = allocate_eospace_array(ci_belec, ci_orbs, ndocc, nactv, xlvl,
                                       &qegrps);
        num_pq = pegrps * qegrps;
        error = allocate_mem_int(&pq_space_pairs, 2, num_pq);
        error = citrunc(aelec, belec, orbitals, nfrzc, ndocc, nactv, nfrzv,
                        xlvl, pstrings, pstr_len, qstrings, qstr_len,
                        peospace, pegrps, qeospace, qegrps, &dtrm_len,
                        pq_space_pairs, &num_pq);
        intorb = ndocc + nactv;
        MPI_Barrier(MPI_COMM_WORLD);
        if (printlvl > 0 && mpi_proc_rank == mpi_root) {
                printf("Determinants = %15d\n", dtrm_len);
                fflush(stdout);
        }
        
        /* Read the molecular orbitals */
        m1len = index1e(orbitals, orbitals);
        m2len = index2e(orbitals, orbitals, orbitals, orbitals);
        strncpy(moflname, "moints", FLNMSIZE);
        moints1 = malloc(sizeof(double) * m1len);
        moints2 = malloc(sizeof(double) * m2len);
        init_dbl_array_0(moints1, m1len);
        init_dbl_array_0(moints2, m2len);
        if (mpi_proc_rank == mpi_root) {
                readmointegrals(moints1, moints2, itype, ci_orbs, moflname,
                                m1len, m2len, &nucrep_e, &frzcore_e);
        }
        MPI_Bcast(moints1, m1len, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(moints2, m2len, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nucrep_e, 1, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&frzcore_e,1, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);

        /* Compute chunk size and lower bounds for each chunk. */
        compute_cimat_chunks(dtrm_len, &hmat_chunk, &row_lwrbnd, &row_uprbnd);
        
        if (printlvl > 0 && mpi_proc_rank == mpi_root) {
                printf("Chunksize = %15d\n", hmat_chunk);
                fflush(stdout);
        }
        if (printlvl > 1) {
                printf("proc %2d: lower bound = %10d, upper bound = %10d\n",
                       mpi_proc_rank, row_lwrbnd, row_uprbnd);
                fflush(stdout);
        }
        
        /* Compute starting values and finishing values for row chunks of
         * Hmat. */
        determinant_string_info(row_lwrbnd, peospace, pegrps, qeospace, qegrps,
                                pq_space_pairs, num_pq, &pq_start, &pstart, &qstart);
        determinant_string_info(row_uprbnd, peospace, pegrps, qeospace, qegrps,
                                pq_space_pairs, num_pq, &pq_final, &pfinal, &qfinal);
        MPI_Barrier(MPI_COMM_WORLD);

        /* Compute diagonal matrix elements. Set global array variables, and
         * allocate global array and local buffer. */
        hdgls_dims[0] = dtrm_len;
        mpi_split_work_array_1d(dtrm_len, &hdgls_chunk[0], &hdgls_lo[0], &hdgls_hi[0]);

        hdgls_handle = NGA_Create(C_DBL, 1, hdgls_dims, "<i|H|i>", hdgls_chunk);
        if (!hdgls_handle) GA_Error("Create failed: <i|H|i>",0);
        if (printlvl > 0 && mpi_proc_rank == mpi_root)
                printf("Global <i|H|i> successfully created.\n");
        hdgls_local = malloc(sizeof(double) * (hdgls_chunk[0] + 10));
        if (hdgls_local == NULL) {
                error_message(mpi_proc_rank, "Could not allocate hdgls_local.",
                              "execute_pjayci");
                error = -100;
                return error;
        }
        
        compute_diagonal_matrix_elements(hdgls_local, hdgls_lo[0], hdgls_hi[0],       
                                         moints1, moints2, ci_aelec, ci_belec,
                                         peospace, pegrps, qeospace, qegrps,
                                         pq_space_pairs, num_pq, pstrings,
                                         qstrings, intorb);  
        if (printlvl > 0 && mpi_proc_rank == mpi_root) {
                printf("<1|H|1> = ");
                printf("%lf\n", (hdgls_local[0] + nucrep_e + frzcore_e));
        }


        
        GA_Destroy(hdgls_handle);
        
        free(moints1);
        free(moints2);
        return error;
}

/*
 * compute_cimat_chunks: compute chunksize of bounds of H for evaluation.
 * Only elements in the upper triangle are computed.
 * Input:
 *  dlen = number of determinants
 */
void compute_cimat_chunks (int dlen, int *chunk, int *lwrbnd, int *uppbnd)
{
        int tmp;
        *chunk = get_upptri_size(dlen);
        *chunk = *chunk / mpi_num_procs;
        tmp = mpi_proc_rank * (*chunk);
        *lwrbnd = get_upptri_element_rownumber(tmp, dlen);
        tmp = (mpi_proc_rank + 1) * (*chunk);
        *uppbnd = get_upptri_element_rownumber(tmp, dlen) - 1;
        /* Upper bound of last chunk should be last row. */
        if (mpi_proc_rank == (mpi_num_procs - 1)) {
                *uppbnd = dlen - 1;
        }
        return;
}

/*
 * compute_diagonal_matrix_elements: compute <i|H|i> elements.
 * Input:
 *  hdgls   = digonal elements <i|H|i> for (i = start,...,final)
 *  start   = starting determinant index
 *  final   = final determinant index
 *  mo1     = 1-e integrals
 *  mo2     = 2-e integrals
 *  aelec   = CI alpha electrons
 *  belec   = CI beta  electrons
 *  peosp   = alpha electron orbital spaces
 *  pegrps  = number of alpha electron orbital spaces
 *  qeosp   = beta  electron orbital spaces
 *  qegrps  = number of beta  electron orbital spaces
 *  pq      = (p,q)-space pairings
 *  npq     = number of (p,q)-space pairings
 *  pstrings= alpha strings
 *  qstrings= beta  strings
 *  intorb  = internal orbitals (DOCC + CAS)
 */
void compute_diagonal_matrix_elements(double *hdgls, int start, int final,
                                      double *mo1, double *mo2, int aelec,
                                      int belec, struct eospace *peosp,
                                      int pegrps,struct eospace *qeosp,
                                      int qegrps, int **pq, int npq,
                                      struct occstr *pstrings,
                                      struct occstr *qstrings,
                                      int intorb)
{
        int pq_start = 0; /* Starting pq pairing index */
        int pq_final = 0; /* Final pq pairing index */
        int p_start  = 0; /* Starting pstring index */
        int p_final  = 0; /* Final pstring index */
        int q_start  = 0; /* Starting qstring index */
        int q_final  = 0; /* Final qstring index */
        int i, j, k;
        int j_start, j_max;
        int k_start, k_max;
        struct det deti;
        int cnt = 0; /* diagonal counter */

        /* Get starting/ending values for determinant-generating loop */
        determinant_string_info(start, peosp, pegrps, qeosp, qegrps, pq,
                                npq, &pq_start, &p_start, &q_start);
        determinant_string_info(final, peosp, pegrps, qeosp, qegrps, pq,
                                npq, &pq_final, &p_final, &q_final);

        /* Each determinant is associated with a triple: (pq, p, q).
         * Because of the loop structure this triple is dependent on we
         * have to test for whether the iteration is the first/last
         * in an block, to determine what values are set. */

        /* Loop over p,q string pairings */
        cnt = 0;
        for (i = pq_start; i <= pq_final; i++) {
                if (i == pq_start) {
                        /* If first pq pair, then we skip ahead */
                        j_start = p_start;
                } else {
                        /* If not first pair, start at beginning p string. */
                        j_start = peosp[(pq[i][0])].start;
                }
                if (i == pq_final) {
                        /* If last pq pair, we end prematurely */
                        j_max = p_final + 1;
                } else {
                        /* If not last pair, we end at last p string */
                        j_max = peosp[(pq[i][0])].start +
                                peosp[(pq[i][0])].nstr;
                }

                /* Set cas flags for electron space */
                if ((peosp[(pq[i][0])].virt + qeosp[(pq[i][1])].virt) == 0) {
                        deti.cas = 1;
                } else {
                        deti.cas = 0;
                }

                /* Loop over alpha (p) strings */ 
                for (j = j_start; j < j_max; j++) {
                        /* Set alpha (p) strings */
                        deti.astr = pstrings[j];

                        /* First loop, we skip ahead */
                        if (i == pq_start &&
                            j_start == p_start && j == p_start) {
                                k_start = q_start;
                        } else {
                                k_start = qeosp[(pq[i][1])].start;
                        }
                        /* Last loop */
                        if (i == pq_final &&
                            j_max == (p_final + 1) && j == (p_final)) {
                                k_max = q_final;
                        } else {
                                k_max = qeosp[(pq[i][1])].start +
                                        qeosp[(pq[i][1])].nstr;
                        }
                        
                        /* Loop over beta (q) strings */
                        for (k = k_start; k < k_max; k++) {
                                        
                                deti.bstr = qstrings[k];
                                hdgls[cnt] = hmatels(deti, deti, mo1, mo2, aelec,
                                                     belec, intorb);
                                cnt++;
                        }
                }
        }
        return;
}
        
/*
 * determinant_string_info: compute string information given
 * a determinant.
 * Input:
 *  detindx = determinant index (position in list)
 *  peosp   = alpha electron orbital spaces
 *  pegrps  = number of alpha electron orbital spaces
 *  qeosp   = beta  electron orbital spaces
 *  qegrps  = number of beta  electron orbital spaces
 *  pq      = (p,q)-space pairings
 *  npq     = number of (p,q)-space pairings
 * Output:
 *  pqval = (p,q) pairing start
 *  pval  = p string of (p,q)-space pairing
 *  qval  = q string of (p,q)-space pairing
 */
void determinant_string_info(int detindx, struct eospace *peosp, int pegrps,
                             struct eospace *qeosp, int qegrps, int **pq,
                             int npq, int *pqval, int *pval, int *qval)
{
        int i, j, k;
        int jstart, jmax, kstart, kmax;
        int cntr = 0; /* counter */
        /* Loop over p,q string pairings */
        for (i = 0; i < npq; i++) {
                jstart = peosp[(pq[i][0])].start;
                kstart = qeosp[(pq[i][1])].start;
                jmax   = jstart + peosp[(pq[i][0])].nstr;
                kmax   = kstart + qeosp[(pq[i][1])].nstr;
                for (j = jstart; j < jmax; j++) {
                        for (k = kstart; k < kmax; k++) {
                                if (cntr == detindx) {
                                        *pqval = i;
                                        *pval = j;
                                        *qval = k;
                                        return;
                                }
                                cntr++;
                        }
                }
        }
        error_message(mpi_proc_rank, "Determinant not found",
                      "determinant_string_info");
        error_flag(mpi_proc_rank, detindx, "determinat_string_info");
        return;
}

/*
 * get_upptri_element_index: get index of an element (i,j) in
 * list of elements in upper triangle of H matrix (n x n).
 */
int get_upptri_element_index (int i, int j, int n)
{
        int result = 0;
        result = (2 * n - i) * (i - 1);
        result = result / 2;
        result = result + j;
        return result;
}

/*
 * get_upptri_element_rownumber: get row number of upper triangle
 * matrix list element.
 */
int get_upptri_element_rownumber (int element, int n)
{
        int row = 0;
        int i = 0;
        for (i = 0; i <= n; i++) {
                if (get_upptri_element_index(i, i, n) > element) break;
        }
        row = i - 1;
        return row;
}

/*
 * get_upptri_size: compute the size of H matrix upper triangle.
 */
int get_upptri_size (int n)
{
        int result = 0;
        result = (n * (n + 1)) / 2;
        return result;
}

