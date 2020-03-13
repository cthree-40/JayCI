// File: pdavidson.c
/*
 * Parallel davidson algorithm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pjayci_global.h"
#include "errorlib.h"
#include "iminmax.h"
#include "arrayutil.h"
#include "mathutil.h"
#include "allocate_mem.h"
#include "mpi_utilities.h"
#include "binarystr.h"
#include "citruncate.h"
#include "action_util.h"
#include "pdavidson.h"

#include <ga.h>
#include <macdecls.h>

#include <mpi.h>


/*
 * pdavidson: parallel implementation of davidson algorithm.
 * NOTE: Global arrays are transposed in C: [columns, rows]
 * dim: [0,1] = [columns,rows]
 */
int pdavidson(struct occstr *pstrings, struct eospace *peospace, int pegrps,
              struct occstr *qstrings, struct eospace *qeospace, int qegrps,
              int **pq_space_pairs, int num_pq, double *moints1, double *moints2,
              int aelec, int belec, int intorb, int ndets, double nucrep_e,
              double frzcore_e, int printlvl, int maxiter, int krymin,
              int krymax, int nroots, int prediagr, int refdim, int restol)
{
        int v_hndl = 0;           /* GLOBAL basis vectors, V */
        int v_dims[2]  = {0, 0};  /* GLOBAL basis vectors dimensions */
        int v_chunk[2] = {0, 0};  /* GLOBAL basis vectors chunk sizes */
        int c_hndl = 0;           /* GLOBAL Hv=c  vectors, C */
        int n_hndl = 0;           /* GLOBAL new vector, N */
        int n_dims[1]  = {0};     /* GLOBAL new vector dimensions */
        int n_chunk[1] = {0};     /* GLOBAL new vector chunk size */
        int r_hndl = 0;           /* GLOBAL residual vector, R */
        int d_hndl = 0;           /* GLOBAL <i|H|i> vector, D */

        double *d_local = NULL;   /* LOCAL <i|H|i> array. */

        double **vhv = NULL;      /* LOCAL v.Hv array */
        double *vhv_data = NULL;  /* LOCAL v.Hv memory block */
        double *vhv_scr  = NULL;  /* LOCAL v.Hv scratch array */
        
        double **hevec = NULL;    /* LOCAL v.Hv eigenvectors */
        double *hevec_data = NULL;/* LOCAL v.Hv eigenvectors memory block */
        double *heval = NULL;     /* LOCAL v.Hv eigenvalues */
        double *hevec_scr= NULL;  /* LOCAL v.Hv eigenvector scratch array */
        
        
        int citer = 0; /* current iteration */
        int croot = 0; /* current root */
        int ckdim = 0; /* current dimension of krylo space */
        
        int lo[2] = {0, 0}; 
        int hi[2] = {0, 0};
        int ld_ndets[1] = {0};

        double totcore_e = 0.0;
        
        int error = 0;
        
        totcore_e = nucrep_e + frzcore_e;
        ld_ndets[0] = ndets;

        /* Allocate GLOBAL arrays: V, Hv=c, N, R, and D */
        if (mpi_proc_rank == mpi_root) printf("Creating global arrays...\n");
        v_dims[0]  = krymax;
        v_dims[1]  = ndets;
        v_chunk[0] = krymax;
        v_chunk[1] = -1; // Distribute evenly
        v_hndl = NGA_Create(C_DBL, 2, v_dims, "Basis vectors", v_chunk);
        if (!v_hndl) GA_Error("Create failed: Basis vectors", 2);
        c_hndl = NGA_Duplicate(v_hndl, "Hv=c vectors");
        if (!c_hndl) GA_Error("Duplicate failed: Hv=c vectors", 2);
        n_dims[0]  = ndets;
        n_chunk[0] = -1; // Distribute evenly
        n_hndl = NGA_Create(C_DBL, 1, n_dims, "New vector", n_chunk);
        if (!n_hndl) GA_Error("Create failed: New vector", 1);
        r_hndl = NGA_Duplicate(n_hndl, "Residual vector");
        if (!r_hndl) GA_Error("Duplicate failed: Residual vector", 1);
        d_hndl = NGA_Duplicate(n_hndl, "Diagonal vector");
        if (!d_hndl) GA_Error("Duplicate failed: Diagonal vectors", 1);
        
        if (mpi_proc_rank == mpi_root) printf("Global arrays created.\n\n");

        /* Allocate local arrays: d, vhv, hevec, heval */
        d_local = malloc(((ndets / mpi_num_procs) + 10) * sizeof(double));
        vhv_data = allocate_mem_double_cont(&vhv, krymax, krymax);
        hevec_data = allocate_mem_double_cont(&hevec, krymax, krymax);
        heval = malloc(sizeof(double) * krymax);
        /* Allocate 1-d scratch arrays. */
        vhv_scr = malloc(krymax * krymax * sizeof(double));
        hevec_scr = malloc(krymax * krymax * sizeof(double));
        
        GA_Sync();
        
        /* Compute diagonal matrix elements <i|H|i> */
        NGA_Zero(d_hndl);
        NGA_Distribution(d_hndl, mpi_proc_rank, lo, hi);
        compute_diagonal_matrix_elements(d_local, lo[0], hi[0], moints1, moints2,
                                         aelec, belec, peospace, pegrps,
                                         qeospace, qegrps, pq_space_pairs,
                                         num_pq, pstrings, qstrings, intorb);
        NGA_Put(d_hndl, lo, hi, d_local, ld_ndets);
        if (printlvl > 0 && lo[0] == 0) {
                printf("<1|H|1> = ");
                printf("%lf\n", (d_local[0] + totcore_e));
        }

        if (mpi_proc_rank == mpi_root)
                printf("Beginning Davidson algorithm...\n");
        
        /* Build initial guess basis vectors. */
        build_init_guess_vectors(prediagr, v_hndl, refdim, krymin, ndets,
                                 pstrings, peospace, pegrps, qstrings,
                                 qeospace, qegrps, pq_space_pairs, num_pq,
                                 moints1, moints2, aelec, belec, intorb);
        if (mpi_proc_rank == mpi_root)
                printf(" Initial guess vectors set.\n");

        
        GA_Sync();
        print_vector_space(v_hndl, 6, ndets);

        /* .. MAIN LOOP .. */
        citer = 1; croot = 1; ckdim = krymin;
        while (citer < maxiter && croot <= nroots) {
                perform_hv_initspace(pstrings, peospace, pegrps, qstrings,
                                     qeospace, qegrps, pq_space_pairs, num_pq,
                                     moints1, moints2, aelec, belec, intorb,
                                     ndets, totcore_e, ckdim, krymax, v_hndl, d_hndl,
                                     c_hndl);
                print_vector_space(c_hndl, 6, ndets);
                make_subspacehmat_ga(v_hndl, c_hndl, ndets, ckdim, vhv);
                print_subspacehmat(vhv, ckdim);
                error = diag_subspacehmat(vhv, hevec, heval, ckdim, krymax,
                                          vhv_scr, hevec_scr);
                if (error != 0)  return error;
                print_subspace_eigeninfo(hevec, heval, ckdim, totcore_e);
                
                citer = maxiter + 1;
        }
        
        return error;
}

/*
 * add_new_vector: add a new vector to basis space.
 */
void add_new_vector(int v_hndl, int cdim, int len, int n_hndl)
{
        int nlo[1] = {0};
        int nhi[1] = {0};
        int nlen = 0;
        int vlo[2] = {0, 0};
        int vhi[2] = {0, 0};
        char trans = "N";

        /*nlo[0] = 0;
        nhi[0] = len - 1;
        vlo[0] = 0;
        vlo[1] = cdim;
        vhi[0] = len - 1;
        vhi[1] = cdim;*/
        nlo[0] = 0;
        nhi[0] = len - 1;
        vlo[0] = 0;
        vhi[0] = len - 1;
        vlo[1] = cdim;
        vhi[1] = cdim;
        int type, ndim, dims[2];
        NGA_Inquire(v_hndl, &type, &ndim, dims);
        printf("cdim = %d\n", cdim);
        printf("len  = %d\n", (len - 1));
        printf("type, ndim, dims = %d, %d, %d %d\n", type, ndim, dims[0], dims[1]);
        NGA_Copy_patch(trans, n_hndl, nlo, nhi, v_hndl, vlo, vhi);

        
        return;
}

/*
 * build_init_guess_vectors: build initial guess basis vectors for davidson
 * procedure.
 * Input:
 *  n     = input routine: 1 = subblock diagonalization, 2 = unit vectors
 *  v     = basis vectors handle (GLOBAL ARRAY)
 *  dim   = dimension of subspace
 *  kmin  = krylov space minimum dimension
 *  ndets = number of determinants
 */
void build_init_guess_vectors(int n, int v, int dim, int kmin, int ndets,
                              struct occstr *pstr, struct eospace *peosp,
                              int pegrps, struct occstr *qstr,
                              struct eospace *qeosp, int qegrps,
                              int **pqs, int num_pq, double *m1, double *m2,
                              int aelec, int belec, int intorb)
{
        double **refspace = NULL; /* Reference space */
        double *rdata = NULL;     /* Reference space data */
        double **initvecs = NULL; /* Initial vectors */
        double *ivdata = NULL;    /* Initial vectors data */
        int error = 0;
        int i = 0;
        int j = 0;
        int lo[2]={0,0}; int hi[2]={0,0};
        int ld[1] = {0};
        
        /* Set all vectors to 0.0 */
        NGA_Zero(v);

        /* Only root process will perform initial guess vector generation. */
        if (mpi_proc_rank == mpi_root) {
                rdata  = allocate_mem_double_cont(&refspace, dim, dim);
                ivdata = allocate_mem_double_cont(&initvecs, ndets, kmin);
                init_dbl_2darray_0(refspace, dim, dim);
                init_dbl_2darray_0(initvecs, ndets, kmin);
        }
        
        if (mpi_proc_rank == mpi_root && n == 1) {
                /* Diagonalize subspace of H */
                init_diag_H_subspace(pstr, peosp, pegrps, qstr, qeosp, qegrps,
                                     pqs, num_pq, m1, m2, aelec, belec, intorb,
                                     ndets, dim, refspace);
                
        } else if (mpi_proc_rank == mpi_root && n == 2) {
                /* Unit vectors for initial guess */
                
                for (i = 0; i < dim; i++) {
                        refspace[i][i] = 1.0;
                }
        }

        if (mpi_proc_rank == mpi_root) {
                /* Set initial guess vectors */
                for (i = 0; i < kmin; i++) {
                        for (j = 0; j < dim; j++) {
                                        initvecs[i][j] = refspace[i][j];
                        }
                }

                ld[0] = ndets; /* (rows in C/C++) N of M X N matrix */
                lo[0] = 0; lo[1] = 0;
                hi[1] = ndets-1; hi[0] = kmin-1;

                NGA_Put(v, lo, hi, ivdata, ld);
                
                deallocate_mem_cont(&refspace, rdata);
                deallocate_mem_cont(&initvecs, ivdata);
        }
        
        return;
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
        /* Upper bound of last chunk should be last row. Because we compute
         * diagonal elements separately, the last row computed is actually
         * the n-1 row, if n is the size of the matrix; so the index for the
         * upperbound is n-2.*/
        if (mpi_proc_rank == (mpi_num_procs - 1)) {
                *uppbnd = dlen - 2;
        }
        return;
}

/*
 * compute_cimat_chunks_2: compute chunksize of bounds of H for evaluation.
 * Only elements in the upper triangle are computed.
 * Input:
 *  dlen = number of determinants
 */
void compute_cimat_chunks_2 (int dlen, int *chunk, int *row_lo, int *row_hi,
                             int *col_lo, int *col_hi)
{
        int tmp;
        *chunk = get_upptri_size(dlen);
        *chunk = *chunk / mpi_num_procs;
        tmp = mpi_proc_rank * (*chunk);
        //*lwrbnd = get_upptri_element_rownumber(tmp, dlen);
        get_upptri_element_position(tmp, dlen, row_lo, col_lo);
        tmp = (mpi_proc_rank + 1) * (*chunk);
        //*uppbnd = get_upptri_element_rownumber(tmp, dlen) - 1;
        get_upptri_element_position(tmp, dlen, row_hi, col_hi);
        /* Upper bound of last chunk should be last row. Because we compute
         * diagonal elements separately, the last row computed is actually
         * the n-1 row, if n is the size of the matrix; so the index for the
         * upperbound is n-2.*/
        //if (mpi_proc_rank == (mpi_num_procs - 1)) {
        //        *uppbnd = dlen - 2;
        //}
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
         * Because of the loop structure this triple is dependent on, we
         * have to test for whether the iteration is the first/last
         * in an block to determine what values are set. */

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
                                k_max = q_final + 1;
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
 * compute_hv_newvector: compute Hv=c for newest vector in basis space.
 * Input:
 *  v_hndl = GA handle for basis vectors (ckdim is new vector)
 *  c_hndl = GA handle for Hv=c vectors
 *  ckdim  = current dimension of space
 */
void compute_hv_newvector(int v_hndl, int c_hndl, int ckdim, struct occstr *pstr,
                          struct eospace *peosp, int pegrps, struct occstr *qstr,
                          struct eospace *qeosp, int qegrps, int **pqs,
                          int num_pq, double *m1, double *m2, int aelec,
                          int belec, int intorb, int ndets, int kmax)
{
        int error = 0; /* Error flag */

        double **c_local = NULL; /* Local c array */
        double *cdata = NULL;
        int c_rows = 0;         /* Local c rows  */
        int c_cols = 0;         /* Local c columns (=1) */

        double **v_local = NULL; /* Local v array */
        double *vdata = NULL;
        int v_rows = 0;         /* Local v rows  */
        int v_cols = 0;         /* Local v columns (=1) */

        /*
         * The following convention is used:
         *    H(i,j)*V(j,k)=C(i,k)
         */
        int start_det_i = 0;    /* Starting determinant index of i */
        int final_det_i = 0;    /* Ending determinant index of i   */
        int start_det_j = 0;    /* Starting determinant index of j */
        int final_det_j = 0;    /* Ending determinant index of j   */
        int ndetsi = 0;         /* Number of i determinants */
        
        int pq_start_i = 0;     /* Starting pq-pair index of i */
        int pstart_i   = 0;     /* Starting alpha string of i  */
        int qstart_i   = 0;     /* Starting beta  string of i  */
        int pq_final_i = 0;     /* Ending pq-pair index of i   */
        int pfinal_i   = 0;     /* Ending alpha string of i    */
        int qfinal_i   = 0;     /* Ending beta  string of i    */
        int pq_start_j = 0;     /* Starting pq-pair index of j */
        int pstart_j   = 0;     /* Starting alpha string of j  */
        int qstart_j   = 0;     /* Starting beta  string of j  */
        int pq_final_j = 0;     /* Ending pq-pair index of j   */
        int pfinal_j   = 0;     /* Ending alpha string of j    */
        int qfinal_j   = 0;     /* Ending beta  string of j    */

        
        int c_lo[2] = {0, 0}; /* starting indices for memory block */
        int c_hi[2] = {0, 0}; /* ending indices of memory block */
        int v_lo[2] = {0, 0}; /* starting indices of memory block */
        int v_hi[2] = {0, 0}; /* ending indices of memory block */

        int v_ld[1] = {0}; /* Leading dimensions of V local buffer */
        int c_ld[1] = {0}; /* Leading dimensions of C local buffer */

        double alpha[1] = {1.0}; /* Scale factor for c_local into c_global. */
        
        /*
         * Determine which block of data is locally owned. And get the blocks
         * of V that are required to compute c.
         * Recall: C is split by rows, and not columns, across processes.
         *         We only require the last columns of both V and C.
         */
        NGA_Distribution(c_hndl, mpi_proc_rank, c_lo, c_hi);
        c_rows = c_hi[0] - c_lo[0] + 1;
        c_cols = c_hi[1] - c_lo[1] + 1;
        c_cols = int_min(c_cols, ckdim);
        c_hi[1] = c_cols - 1;
        c_lo[1] = c_cols - 1; // Last column is all that is required.
        c_ld[0] = 1; //c_hi[1] - c_lo[1];
        c_cols = c_ld[0];
        c_ld[0] = 0;
        cdata = allocate_mem_double_cont(&c_local, c_rows, c_cols);
        v_lo[0] = 0;
        v_lo[1] = c_lo[1];
        v_hi[0] = ndets - 1;
        v_hi[1] = c_hi[1];
        v_rows = v_hi[0] - v_lo[0] + 1;
        v_cols = v_hi[1] - v_lo[1] + 1;
        v_ld[0] = 0;
        vdata = allocate_mem_double_cont(&v_local, v_rows, v_cols);


        NGA_Get(v_hndl, v_lo, v_hi, vdata, v_ld);
        NGA_Get(c_hndl, c_lo, c_hi, cdata, c_ld);
        
        
        GA_Sync();

        /* Get values for determinnat index at beginning/ending of block */
        start_det_i = c_lo[0];
        final_det_i = c_hi[0];
        start_det_j = v_lo[0];
        final_det_j = v_hi[0];
        determinant_string_info(start_det_i, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_start_i, &pstart_i, &qstart_i);
        determinant_string_info(final_det_i, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_final_i, &pfinal_i, &qfinal_i);
        determinant_string_info(start_det_j, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_start_j, &pstart_j, &qstart_j);
        determinant_string_info(final_det_j, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_final_j, &pfinal_j, &qfinal_j);

        /* Evaluate block of H(i,j)*V(j,k)=C(i,k). This is done via the        
         * space indexes. Each block is over all V vectors, k.*/
        evaluate_hdblock_ij(pq_start_i, pstart_i, qstart_i,
                            pq_final_i, pfinal_i, qfinal_i,
                            pq_start_j, pstart_j, qstart_j,
                            pq_final_j, pfinal_j, qfinal_j,
                            v_rows, v_cols, v_local,
                            c_rows, c_cols, c_local,
                            start_det_i, final_det_i,
                            start_det_j, final_det_j,
                            ndets, peosp, pegrps, pstr, qeosp, qegrps, qstr,
                            pqs, num_pq, m1, m2, aelec, belec, intorb);

        NGA_Acc(c_hndl, c_lo, c_hi, cdata, c_ld, alpha);
        GA_Sync();

        deallocate_mem_cont(&v_local, vdata);
        deallocate_mem_cont(&c_local, cdata);
        
        return;
}

/*
 * compute_residual_norm: compute the norm of the residual vector.
 */
void compute_residual_norm (int r_hndl,  double *norm)
{
        double *r_local = NULL; /* r local array */
        int r_lo[1] = {0};
        int r_hi[1] = {0};
        int r_ld[1] = {0};
        int r_len;
        double local_norm;
        double global_norm;
        int i;

        /* Find what chunk of R is held locally. Allocate local
         * array of R. Get chunk of R from global array. */
        NGA_Distribution(r_hndl, mpi_proc_rank, r_lo, r_hi);
        r_len = r_hi[0] - r_lo[0] + 1;
        r_local = malloc(sizeof(double) * r_len);
        init_dbl_array_0(r_local, r_len);
        if (r_local == NULL) {
                error_message(mpi_proc_rank, "Error allocating memory",
                              "compute_residual_norm");
                return;
        }
        NGA_Get(r_hndl, r_lo, r_hi, r_local, r_ld);

        /* Sum squares of elements */
        local_norm = dot_product(r_local, r_local, r_len);
        /* Sum all partial sums and take square root */
        MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        *norm = sqrt(global_norm);

        free(r_local);
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
 * diag_subspacehmat: diagonalize the v.Hv subspace hamiltonian, returning
 * eigenvectors and eigenvalues.
 * Input:
 *  vhv = v.Hv
 *  hevec = Eigenvectors of v.Hv
 *  heval = Eigenvalues of v.Hv
 *  cdim  = current dimension of v.Hv
 *  maxdim = maximum dimension of v.Hv
 * Output:
 *  error = error flag
 */
int diag_subspacehmat(double **vhv, double **hevec, double *heval, int cdim,
                      int maxdim, double *vhv_scr, double *hevec_scr)
{
        int error = 0; /* Error flag */
        /* Copy 2d arrays into corresponding 1d arrays */
        cparray_2d1d(vhv, cdim, cdim, vhv_scr);
        error = diagmat_dsyevr(vhv_scr, cdim, hevec_scr, heval);
        if (error != 0) {
                error_flag(mpi_proc_rank, error, "diag_subspacehmat");
                return error;
        }
        cparray_1d2d(hevec_scr, hevec, cdim, cdim);
        init_dbl_array_0(vhv_scr, (maxdim * maxdim));
        init_dbl_array_0(hevec_scr, (maxdim * maxdim));
        return error;
}

/*
 * evaluate_hdblock_ij: evaluate a block of H:
 *  H(i,j)*V(j,k)=C(i,k)
 *  Input:
 *   vrows  = columns of C
 *   vcols  = rows of C
 *   crows  = columns of C
 *   ccols  = rows of C
 *   v      = V(j,k)
 *   c      = C(i,k)
 *   starti = starting index i
 *   finali = ending index i
 *   startj = starting index j
 *   finalj = ending index j
 *   ndets  = number of determinants
 *   peosp  = alpha electron orbital spaces
 *   pegrps = number of alpha string orbital spaces
 *   pstr   = alpha electron strings
 *   qeosp  = beta electron orbital spaces
 *   qegrps = number of beta string orbital spaces
 *   qstr   = beta electron strings
 *   pq     = (p,q)-space pairings
 *   npq    = number of (p, q)-space pairings
 *   mo1    = 1-e integrals
 *   mo2    = 2-e integrals
 *   aelec  = alpha electrons
 *   belec  = beta  electrons
 *   intorb = internal orbitals (DOCC + CAS)
 *  Output:
 *   C      = C(i,k)
 */
void evaluate_hdblock_ij(int pq_start_i, int pstart_i, int qstart_i,
                         int pq_final_i, int pfinal_i, int qfinal_i,
                         int pq_start_j, int pstart_j, int qstart_j,
                         int pq_final_j, int pfinal_j, int qfinal_j,
                         int vrows, int vcols, double **v,
                         int crows, int ccols, double **c,
                         int starti, int finali, int startj, int finalj,
                         int ndets, struct eospace *peosp, int pegrps,
                         struct occstr *pstr,
                         struct eospace *qeosp, int qegrps,
                         struct occstr *qstr,  int **pq,
                         int npq, double *mo1, double *mo2, int aelec,
                         int belec, int intorb)
{
        int cnti = 0, cntj = 0; /* i,j determinant counters */

        int **d_triplet = NULL; /* |i> = (p, q, CAS-flag) list*/
        int *d_trip_dat = NULL; /* d_triplet memory block */
        int **dtj = NULL;
        int *dtjdat=NULL;
        int ndetj;
        
        struct det deti;        /* Determinant i */
        struct det detj;        /* Determinant j */
        int ndeti;              /* Number of i determinants (chunk size) */

        double hijval = 0.0;    /* <i|H|j> value */
        
        int pstart = 0, pmax = 0;
        int qstart = 0, qmax = 0;
        int i = 0, j = 0, k = 0, l = 0;
        int ii= 0;
        
        printf(" start i = %d, final i = %d\n", starti, finali);
        printf(" start j = %d, final j = %d\n", startj, finalj);
        /*
         * Each determinant is associated with a triple: (p, q, flag).
         * First, we generate this list.
         */
        ndeti = finali - starti + 1;
        ndetj = finalj - startj + 1;
        d_trip_dat = allocate_mem_int_cont(&d_triplet, 3, ndeti);
        dtjdat = allocate_mem_int_cont(&dtj, 3, ndetj);
        generate_det_triples(ndeti, d_triplet, pq_start_i, pstart_i, qstart_i,
                             pq_final_i, pfinal_i, qfinal_i, pq, npq,
                             peosp, pegrps, qeosp, qegrps);
        generate_det_triples(ndetj, dtj, pq_start_j, pstart_j, qstart_j,
                             pq_final_j, pfinal_j, qfinal_j, pq, npq,
                             peosp, pegrps, qeosp, qegrps);

        /* Loop through list of triplets for determinants |i>. */
        for (i = 0; i < ndeti; i++) {
                deti.astr = pstr[d_triplet[i][0]];
                deti.bstr = qstr[d_triplet[i][1]];
                deti.cas = d_triplet[i][2];
                /* Loop over determinants |j> */
                for (j = 0; j < ndetj; j++) {
                        detj.astr = pstr[dtj[j][0]];
                        detj.bstr = qstr[dtj[j][1]];
                        detj.cas  = dtj[j][2];
                        
                        hijval = hmatels(deti, detj, mo1, mo2,
                                         aelec, belec, intorb);
                        /* H_ij*v_jl = c_il */
                        for (l = 0; l < vcols; l++) {
                                c[l][i] = c[l][i] +
                                        hijval * v[l][j];
                        }
                }
        }
        
        return;
                /* Loop over determinants |j> */
//                for (ii = 0; ii < npq; ii++) {
//                        pstart = peosp[(pq[ii][0])].start;
//                        pmax   = pstart + peosp[(pq[ii][0])].nstr;
//                        qstart = qeosp[(pq[ii][1])].start;
//                        qmax   = qstart + qeosp[(pq[ii][1])].nstr;
//                        detj.cas = 0;
//                        if ((peosp[(pq[ii][0])].virt +
//                             qeosp[(pq[ii][0])].virt) == 0) {
//                                detj.cas = 1;
//                        }
//                        cntj = 0;                
//                        /* Loop over alpha strings. */
//                        for (j = pstart; j < pmax; j++) {
//                                detj.astr = pstr[j];
//                                /* Loop over beta strings */
//                                for (k = qstart; k < qmax; k++) {
//                                        detj.bstr = qstr[k];
//                                        hijval = hmatels(deti, detj, mo1, mo2,
//                                                         aelec, belec, intorb);
//                                        /* H_ij*v_jk = c_ik */
//                                        for (l = 0; l < vcols; l++) {
//                                                c[l][starti + i] = c[l][starti + i] +
//                                                        hijval * v[l][cntj];
//                                        }
//                                        cntj++;
//                                }
//                        }
//                        
}

/*
 * generate_det_triples: generate list of triplets for each determinant:
 *  |i> = |(pq, p, q)>
 */
void generate_det_triples (int ndeti, int **d_triplet, int pq_start,
                           int p_start, int q_start, int pq_final,
                           int p_final, int q_final, int **pq, int npq,
                           struct eospace *peosp, int pegrps,
                           struct eospace *qeosp, int qegrps)
{
        int i, j, k;
        int j_start, j_max;
        int k_start, k_max;
        int cnt;

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

                /* Loop over alpha (p) strings */ 
                for (j = j_start; j < j_max; j++) {

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
                                k_max = q_final + 1;
                        } else {
                                k_max = qeosp[(pq[i][1])].start +
                                        qeosp[(pq[i][1])].nstr;
                        }
                        
                        /* Loop over beta (q) strings */
                        for (k = k_start; k < k_max; k++) {
                                d_triplet[cnt][0] = j;
                                d_triplet[cnt][1] = k;
                                /* Set cas flags for electron space */
                                if ((peosp[(pq[i][0])].virt +
                                     qeosp[(pq[i][1])].virt) == 0) {
                                        d_triplet[cnt][2] = 1;
                                } else {
                                        d_triplet[cnt][2] = 0;
                                }
                                cnt++;
                        }
                }
        }
        return;

}

/*
 * generate_newvector: generate new vector to be added to space.
 *  n = - r / (d - E)
 * Input:
 *  r_hndl = GA handle for residual vector
 *  d_hndl = GA handle for diagonal elements
 *  n_hndl = GA handle for new vector
 */
void generate_newvector (int r_hndl, int d_hndl, double eval, int ndets,
                         int n_hndl)
{
        int n_lo[1] = {0};  /* starting indices of memory block */
        int n_hi[1] = {0};  /* ending indices of memory block */
        int n_ld[1] = {0};  /* leading dimension of local N buffer */

        double *r_local = NULL;
        double *d_local = NULL;
        double *n_local = NULL;

        int nlen;
        int i;

        /* Zero-out new vector GA. Determine which blocks of N are
         * locally owned. Read-in this data. Use these dimensions for
         * R and D. */
        NGA_Zero(n_hndl);
        NGA_Distribution(n_hndl, mpi_proc_rank, n_lo, n_hi);
        nlen = n_hi[0] - n_lo[0] + 1;
        printf("nlen = %d\n", nlen);
        /* Allocate local arrays */
        n_local = malloc(sizeof(double) * nlen);
        r_local = malloc(sizeof(double) * nlen);
        d_local = malloc(sizeof(double) * nlen);
        /* Get data for each local array */
        NGA_Get(n_hndl, n_lo, n_hi, n_local, n_ld);
        NGA_Get(r_hndl, n_lo, n_hi, r_local, n_ld);
        NGA_Get(d_hndl, n_lo, n_hi, d_local, n_ld);

        /* Compute new vector */
        for (i = 0; i < nlen; i++) {
                n_local[i] = (-1.0) * r_local[i] / (d_local[i] - eval);
        }

        /* Put new vector into global array */
        NGA_Put(n_hndl, n_lo, n_hi, n_local, n_ld);
        
        free(n_local);
        free(r_local);
        free(d_local);
        
        return;
}
/*
 * generate_residual: generate residual vector.
 *  r_i = v_aj * (c_ji - e_a*v_ji)
 * Input:
 *  v_hndl = global arrays handle for vectors, V
 *  c_hndl = global arrays handle for vectors, C=Hv
 *  r_hndl = global arrays handle for residual vector
 *  hevec  = eigenvectors of v.Hv
 *  heval  = eigenvalues of v.Hv
 *  ndets  = number of determinants (length of vectors V, C, & r)
 *  ckdim  = current dimension of krylov subspace
 *  croot  = current root being optimized.
 */
void generate_residual (int v_hndl, int c_hndl, int r_hndl, double **hevec,
                        double *heval, int ndets, int ckdim, int croot)
{
        int root_id;

        int r_lo[1] = {0};     /* starting indices of memory block */
        int r_hi[1] = {0};     /* ending indices of memory block */
        int r_ld[1] = {0};     /* Leading dimension of R local buffer */
        int v_lo[2] = {0, 0};  /* starting indices of memory block */
        int v_hi[2] = {0, 0};  /* ending indices of memory block */
        int v_ld[1] = {0};     /* Leading dimension of V local buffer */
        int c_lo[2] = {0, 0};  /* starting indices of memory block */
        int c_hi[2] = {0, 0};  /* ending indices of memory block */
        int c_ld[1] = {0};     /* Leading dimension of C local buffer */

        double **v_local = NULL;
        double *vdata = NULL;
        double **c_local = NULL;
        double *cdata = NULL;
        double *r_local = NULL;

        int i_start = 0, i_final = 0;
        int v_rows = 0;
        int v_cols = 0;

        int i, j;
        
        root_id = croot - 1;
        /* Determine which blocks of data for r is locally owned. Read
         * this data. Use these dimensions for getting chunks of V and C. */
        NGA_Zero(r_hndl);
        NGA_Distribution(r_hndl, mpi_proc_rank, r_lo, r_hi);
        // i_start = r_lo[0];
        // i_final = r_hi[0];
        v_rows = r_hi[0] - r_lo[0] + 1;
        v_cols = ckdim;
        if (mpi_proc_rank == mpi_root) {
                printf("i_start = %d\ni_final = %d\nv_rows = %d\n",
                       r_lo[0], r_hi[0], v_rows);
        }
        printf("vrows = %d\nvcols = %d\n", v_rows, v_cols);
        /* Allocate arrays */
        vdata = allocate_mem_double_cont(&v_local, v_rows, v_cols);
        cdata = allocate_mem_double_cont(&c_local, v_rows, v_cols);
        r_local = malloc(sizeof(double) * v_rows);
        if (r_local == NULL || cdata == NULL || vdata == NULL) {
                error_message(mpi_proc_rank, "Error allocating memory",
                              "generate_residual");
                return -100;
        }
        
        v_lo[0] = r_lo[0];
        v_lo[1] = 0; // Always grab first column.
        v_hi[0] = r_hi[0];
        v_hi[1] = v_cols; // Grab all computed vectors
        v_ld[0] = v_cols;
        printf(" (%d, %d); (%d, %d); %d\n", v_lo[0], v_lo[1],
               v_hi[0], v_hi[1], v_ld[0]);
        //r_ld[0] = 0;
        NGA_Get(v_hndl, v_lo, v_hi, vdata, v_ld);
        NGA_Get(c_hndl, v_lo, v_hi, cdata, v_ld);
        NGA_Get(r_hndl, r_lo, r_hi, r_local, r_ld);

        printf("vrows = %d ckdim = %d root_id = %d\n", v_rows, ckdim, root_id);
        init_dbl_array_0(r_local, v_rows);
        for (i = 0; i < v_cols; i++) {
                printf(" Eigenvector %d: %15.8lf ", i, heval[i]);
                for (j = 0; j < v_cols; j++) {
                        printf("%15.8lf ", hevec[i][j]);
                }
                printf("\n");
        }

        for (i = 0; i < v_rows; i++) {
                for (j = 0; j < v_cols; j++) {
                        r_local[i] += hevec[root_id][j] *
                                (c_local[j][i] - heval[root_id] * v_local[j][i]);
                }
        }

        NGA_Put(r_hndl, r_lo, r_hi, r_local, r_ld);

        //deallocate_mem_cont(&c_local, cdata);
        //deallocate_mem_cont(&v_local, vdata);
        //free(r_local);
        
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
 * get_upptri_element_position: get row and column indices of an n x n
 * upper triangle matrix list element.
 */
void get_upptri_element_position (int element, int n, int *i, int *j)
{
        int tmp = 0;
        *i = 0;
        *j = 0;
        for (*i = 0; *i < n; (*i)++) {
                for (*j = *i; *j < n; (*j)++) {
                        tmp++;
                        if (tmp == element) return;
                }
        }
        *i = 0;
        *j = 0;
        return;
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

/*
 * init_diag_H_subspace: generate reference vectors from diagonalization
 * of a subspace of Hij.
 */
void init_diag_H_subspace(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets, int dim,
                          double **refspace)
{
        struct det deti;
        struct det detj;
        double **hij = NULL;
        double *hij_data = NULL;
        double *rdata = NULL;
        double *rev = NULL;
        int **d_triplet = NULL;
        int *d_trip_dat = NULL;
        int pq_start = 0, pq_final = 0;
        int p_start = 0, p_final = 0;
        int q_start = 0, q_final = 0;
        int i = 0, j = 0, k = 0, l = 0;
        int ii;
        int error = 0;
        int cnti, cntj;
        int start, final;
        start = 0;
        final = dim - 1;
        /* Allocate h matrix subblock and refvec data */
        hij_data = allocate_mem_double_cont(&hij, dim, dim);
        rdata = malloc(sizeof(double) * dim * dim);
        rev = malloc(sizeof(double) * dim);
        /* Get starting/ending values for determinant-generating loop */
        determinant_string_info(start, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_start, &p_start, &q_start);
        determinant_string_info(final, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_final, &p_final, &q_final);
        /* Each determinant is associated with a triple: (pq, p, q).
         * Because of the loop structure this triple is dependent on, we
         * have to test for whether the iteration is the first/last
         * in an block to determine what values are set. */
        d_trip_dat = allocate_mem_int_cont(&d_triplet, 3, dim);
        generate_det_triples(dim, d_triplet, pq_start, p_start, q_start,
                             pq_final, p_final, q_final, pqs, num_pq, peosp,
                             pegrps, qeosp, qegrps);
        /* Loop through list of triplets for determinants <i| . */
        for (i = 0; i < dim; i++) {
                deti.astr = pstr[d_triplet[i][0]];
                deti.bstr = qstr[d_triplet[i][1]];
                deti.cas = d_triplet[i][2];
                /* Loop over determinants |j> */
                for (j = 0; j < dim; j++) {
                        detj.astr = pstr[d_triplet[j][0]];
                        detj.bstr = qstr[d_triplet[j][1]];
                        detj.cas = d_triplet[j][2];

                        hij[i][j] = hmatels(deti, detj, m1, m2, aelec,
                                            belec, intorb);
                }
        }

        /* Diagonalize hij */
        error = diagmat_dsyevr(hij_data, dim, rdata, rev);
        printf(" eigenvalues = %lf %lf %lf\n",
               rev[0] + total_core_e,
               rev[1] + total_core_e,
               rev[2] + total_core_e);
        /* Copy rdata into refspace */
        ii = 0;
        for (i = 0; i < dim; i++) {
                for (j = 0; j < dim; j++) {
                        refspace[i][j] = rdata[ii];
                        ii++;
                }
        }

        deallocate_mem_cont(&hij, hij_data);
        free(rdata);
        return;
}

/*
 * make_subspacehmat_ga build v.Hv matrix within GA toolkit
 */
void make_subspacehmat_ga (int v_hndl, int c_hndl, int ndets, int ckdim,
                        double **vhv)
{
        int alo[2] = {0, 0};
        int ahi[2] = {0, 0};
        int blo[2] = {0, 0};
        int bhi[2] = {0, 0};
        int i, j;
        ahi[1] = ndets - 1;
        bhi[1] = ndets - 1;
        for (i = 0; i < ckdim; i++) {
                alo[0] = i;
                ahi[0] = i;
                for (j = 0; j < ckdim; j++) {
                        blo[0] = j;
                        bhi[0] = j;
                        vhv[i][j] = NGA_Ddot_patch(v_hndl, 'n', alo, ahi,
                                                   c_hndl, 'n', blo, bhi);
                }
        }
        
        return;
}
/*
 * make_subspacehmat: build v.Hv matrix.
 */
void make_subspacehmat (int v_hndl, int c_hndl, int ndets, int ckdim,
                        double **vhv)
{
        int v_lo[2] = {0, 0};  /* starting indices of memory block */
        int v_hi[2] = {0, 0};  /* ending indices of memory block */
        int v_ld[1] = {0};     /* Leading dimension of V local buffer */
        int c_lo[2] = {0, 0};  /* starting indices of memory block */
        int c_hi[2] = {0, 0};  /* ending indices of memory block */
        int c_ld[1] = {0};     /* Leading dimension of C local buffer */

        double **v_local = NULL;   /* Local vectors v */
        double *vdata = NULL;      /* Local vector data */
        int vrows = 0, vcols = 0; /* Rows and columns of local V */
        double **c_local = NULL;   /* Local vectors c */
        double *cdata = NULL;      /* Local vector data */
        
        int i, j, k;
        
        /* Determine which block of data for V is locally owned. Read
         * in this data. Get corresponding data for C = HV. */
        NGA_Distribution(v_hndl, mpi_proc_rank, v_lo, v_hi);
        vrows = v_hi[0] - v_lo[0] + 1;
        vcols = v_hi[1] - v_lo[1] + 1;
        vcols = ckdim;
        v_hi[1] = ckdim - 1;
        v_ld[0] = ckdim;
        printf("Hello from, %d: %d %d %d %d %d %d\n",
               mpi_proc_rank, v_lo[0], v_lo[1], v_hi[0], v_hi[1], vrows, vcols);
        vdata = allocate_mem_double_cont(&v_local, vrows, vcols);
        cdata = allocate_mem_double_cont(&c_local, vrows, vcols);

        NGA_Get(v_hndl, v_lo, v_hi, vdata, v_ld);
        c_ld[0] = vcols;
        NGA_Get(c_hndl, v_lo, v_hi, cdata, c_ld);

        /* Perform v_i.c_j for i, j vectors */
        printf("ckdim = %d\n", ckdim);
        for (i = 0; i < ckdim; i++) {
                for (j = 0; j < ckdim; j++) {
                        vhv[i][j] = dot_product(&(v_local[i][0]), &(c_local[j][0]), vrows);
                }
        }
        
        GA_Sync();

        deallocate_mem_cont(&v_local, vdata);
        deallocate_mem_cont(&c_local, cdata);
        
        return;
}

/*
 * orthonormalize_newvector: orthogonalize new vector to rest of basis.
 * Normalize result.
 * Input:
 *  v_hndl = GA handle for basis vectors
 *  nvecs  = number of basis vectors
 *  ndets  = number of determinants
 *  n_hndl = GA handle for new vector
 */
void orthonormalize_newvector (int v_hndl, int nvecs, int ndets, int n_hndl)
{
        double *ovrlps = NULL;    /* Overlaps with basis space (local) */
        double *overlaps = NULL;  /* Overlaps with basis space (global) */
        double *nlocal = NULL;    /* Local new vector array */
        double **vlocal= NULL;    /* Local v vector array */
        double *vdata  = NULL;
        int nlen = 0;
        int vlo[2] = {0, 0};
        int vhi[2] = {0, 0};
        int vld[1] = {0};
        int nlo[1] = {0};
        int nhi[1] = {0};
        int nld[1] = {0};
        double alpha[1] = {0.0};
        double beta[1]  = {0.0};
        double nnrm_local  = 0.0;
        double nnrm_global = 0.0;
        int i, j;
        /* Get GA n vector held locally on process, and corresponding
         * basis vector patch, V.*/
        NGA_Distribution(n_hndl, mpi_proc_rank, nlo, nhi);
        nlen = nhi[0] - nlo[0] + 1;
        nlocal = malloc(sizeof(double) * nlen);
        NGA_Get(n_hndl, nlo, nhi, nlocal, nld);
        vhi[0] = nlen - 1;
        vhi[1] = nvecs - 1;
        vld[0] = nvecs;
        vdata = allocate_mem_double_cont(&vlocal, nlen, nvecs);
        NGA_Get(v_hndl, vlo, vhi, vdata, vld);
        /* Compute overlaps */
        ovrlps = malloc(sizeof(double) * nvecs);
        overlaps = malloc(sizeof(double) * nvecs);
        for (i = 0; i < nvecs; i++) {
                ovrlps[i] = dot_product(&(vlocal[i][0]), nlocal, nlen);
        }
        MPI_Allreduce(ovrlps, overlaps, nvecs, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        /* Subtract overlaps from vector N */
        for (i = 0; i < nvecs; i++) {
                for (j = 0; j < nlen; j++) {
                        nlocal[j] = nlocal[j] - overlaps[i]*vlocal[i][j];
                }
        }
        /* Normalize this vector */
        nnrm_local = dot_product(nlocal, nlocal, nlen);
        MPI_Allreduce(&nnrm_local, &nnrm_global, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        nnrm_global = sqrt(nnrm_global);
        for (i = 0; i < nlen; i++) {
                nlocal[i] = nlocal[i] / nnrm_global;
        }

        /* Check overlaps */
        for (i = 0; i < nvecs; i++) {
                ovrlps[i] = dot_product(&(vlocal[i][0]), nlocal, nlen);
                overlaps[i] = 0.0;
        }
        MPI_Allreduce(ovrlps, overlaps, nvecs, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        for (i = 0; i < nvecs; i++) {
                printf(" overlap %d = %14.8lf\n", i, overlaps[i]);
                if (overlaps[i] > 0.000001) {
                        error_message(mpi_proc_rank,
                                      " Warning! Non-zero overlap.",
                                      "orthonormalize_newvector");
                }
        }
        /* Check norm */
        nnrm_local = dot_product(nlocal, nlocal, nlen);
        MPI_Allreduce(&nnrm_local, &nnrm_global, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        nnrm_global = sqrt(nnrm_global);
        if ((nnrm_global - 1.0) > 0.000001) {
                error_message(mpi_proc_rank,
                              " Warning! New vector norm != 1.0",
                              "orthonormalize_newvector");
        }


        /* Update global array */
        NGA_Put(n_hndl, nlo, nhi, nlocal, nld);

        /* Deallocate arrays */
        deallocate_mem_cont(&vlocal, vdata);
        free(nlocal);
        free(overlaps);
        free(ovrlps);
        return;
}

/*
 * peform_hv_initspace: perform Hv=c on basis vectors v_i, i = 1, .., n.
 * Input:
 *  pstr  = alpha strings
 *  peosp = alpha electron spaces
 *  pegrps= number of alpha electron spaces
 *  qstr  = beta  strings
 *  qeosp = beta  electron spaces
 *  qegrps= number of beta  electron spaces
 *  pqs   = alpha and beta space pairs
 *  num_pq= number of alpha and beta space pairs
 *  m1    = 1-e integrals
 *  m2    = 2-e integrals
 *  aelec = alpha electrons
 *  belec = beta  electrons
 *  intorb= internal orbitals
 *  ndets = number of determinants
 *  core_e= core energy
 *  dim   = number of basis vectors
 *  mdim  = maximum size of krylov space
 *  v_hndl= (GLOBAL ARRAY HANDLE) basis vectors
 *  d_hndl= (GLOBAL ARRAY HANDLE) diagonal elements <i|H|i>
 * Output:
 *  c_hndl= (GLOBAL ARRAY HANDLE) Hv=c vectors
 */
void perform_hv_initspace(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets,
                          double core_e, int dim, int mdim, int v_hndl, int d_hndl,
                          int c_hndl)
{
        int error = 0; /* Error flag */
        
        double **c_local = NULL;  /* Local c array */
        double *cdata = NULL;    /* Local c array data */
        int c_rows = 0;          /* Local c rows */
        int c_cols = 0;          /* Local c columns */

        double **v_local = NULL;  /* Local v array */
        double *vdata = NULL;    /* Local v array data */
        int v_rows = 0;          /* Local v rows */
        int v_cols = 0;          /* LOcal v columns */

        /*
         * The following convention is used:
         *    H(i,j)*V(j,k)=C(i,k)
         */
        int start_det_i = 0;    /* Starting determinant index of i */
        int final_det_i = 0;    /* Ending determinant index of i   */
        int start_det_j = 0;    /* Starting determinant index of j */
        int final_det_j = 0;    /* Ending determinant index of j   */
        int ndetsi = 0;         /* Number of i determinants */
        
        int pq_start_i = 0;     /* Starting pq-pair index of i */
        int pstart_i   = 0;     /* Starting alpha string of i  */
        int qstart_i   = 0;     /* Starting beta  string of i  */
        int pq_final_i = 0;     /* Ending pq-pair index of i   */
        int pfinal_i   = 0;     /* Ending alpha string of i    */
        int qfinal_i   = 0;     /* Ending beta  string of i    */
        int pq_start_j = 0;     /* Starting pq-pair index of j */
        int pstart_j   = 0;     /* Starting alpha string of j  */
        int qstart_j   = 0;     /* Starting beta  string of j  */
        int pq_final_j = 0;     /* Ending pq-pair index of j   */
        int pfinal_j   = 0;     /* Ending alpha string of j    */
        int qfinal_j   = 0;     /* Ending beta  string of j    */

        
        int c_lo[2] = {0, 0}; /* starting indices for memory block */
        int c_hi[2] = {0, 0}; /* ending indices of memory block */
        int v_lo[2] = {0, 0}; /* starting indices of memory block */
        int v_hi[2] = {0, 0}; /* ending indices of memory block */

        int v_ld[1] = {0}; /* Leading dimensions of V local buffer */
        int c_ld[1] = {0}; /* Leading dimensions of C local buffer */

        double alpha[1] = {1.0}; /* Scale factor for c_local into c_global. */
        
        NGA_Zero(c_hndl);

        /* Determine which block of data is locally owned. And get the blocks of
         * V that are required to compute c. */
        NGA_Distribution(c_hndl, mpi_proc_rank, c_lo, c_hi);
        printf("c_lo = %d %d; c_hi = %d %d\n", c_lo[0], c_lo[1], c_hi[0], c_hi[1]);
        c_cols = c_hi[0] - c_lo[0] + 1;
        c_rows = c_hi[1] - c_lo[1] + 1;
        c_cols = int_min(c_cols, dim);
        c_hi[0] = c_cols - 1;
        c_ld[0] = c_rows;
        printf("C_LD[0] = %d\n", c_cols);
        cdata = allocate_mem_double_cont(&c_local, c_rows, c_cols);
        v_lo[0] = 0;
        v_lo[1] = 0;
        v_hi[0] = dim - 1;
        v_hi[1] = ndets - 1;
        v_cols = v_hi[0] - v_lo[0] + 1;
        v_rows = v_hi[1] - v_lo[1] + 1;
        v_ld[0] = v_rows;
        vdata = allocate_mem_double_cont(&v_local, v_rows, v_cols);

        NGA_Get(v_hndl, v_lo, v_hi, vdata, v_ld);
        NGA_Get(c_hndl, c_lo, c_hi, cdata, c_ld);

        printf("Perf_HV: \n");
        print_vector_space(v_hndl, v_cols, 21);
        print_vector_space(c_hndl, c_cols, 21);
        GA_Sync();
        
        /* Get values for determinant index at beginning/ending of block */
        start_det_i = c_lo[1];
        final_det_i = c_hi[1];
        start_det_j = v_lo[1];
        final_det_j = v_hi[1];
        determinant_string_info(start_det_i, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_start_i, &pstart_i, &qstart_i);
        determinant_string_info(final_det_i, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_final_i, &pfinal_i, &qfinal_i);
        determinant_string_info(start_det_j, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_start_j, &pstart_j, &qstart_j);
        determinant_string_info(final_det_j, peosp, pegrps, qeosp, qegrps, pqs,
                                num_pq, &pq_final_j, &pfinal_j, &qfinal_j);

        
        /* Evaluate block of H(i,j)*V(j,k)=C(i,k). This is done via the
         * space indexes. Each block is over all V vectors, k.*/
        evaluate_hdblock_ij(pq_start_i, pstart_i, qstart_i,
                            pq_final_i, pfinal_i, qfinal_i,
                            pq_start_j, pstart_j, qstart_j,
                            pq_final_j, pfinal_j, qfinal_j,
                            v_rows, v_cols, v_local,
                            c_rows, c_cols, c_local,
                            start_det_i, final_det_i,
                            start_det_j, final_det_j,
                            ndets, peosp, pegrps, pstr, qeosp, qegrps, qstr,
                            pqs, num_pq, m1, m2, aelec, belec, intorb);

        NGA_Acc(c_hndl, c_lo, c_hi, cdata, c_ld, alpha);
        //NGA_Put(c_hndl, c_lo, c_hi, cdata, c_ld);
        GA_Sync();

        printf("    1    \n---------\n");
        for (int w = 0; w < 21; w++) {
                printf("%d%14.5lf\n", (w + 1), v_local[0][w]);
        }
        printf("\n");
        deallocate_mem_cont(&v_local, vdata);
        deallocate_mem_cont(&c_local, cdata);
        
        return;
}

/*
 * print_subspace_eigeninfo: print diagonalization iformation for krylov
 * space.
 */
void print_subspace_eigeninfo(double **v, double *e, int kdim, double core_e)
{
        int i = 0, j = 0;
        if (mpi_proc_rank == mpi_root) {
                printf("\n");
                printf(" Eigenvalues:\n");
                for (i = 0; i < kdim; i++) {
                        printf("%14.7lf", (e[i] + core_e));
                }
                printf("\n");
                printf(" Eigenvectors: \n");
                for (i = 0; i < kdim; i++) {
                        for (j = 0; j < kdim; j++) {
                                printf("%14.7lf", v[j][i]);
                        }
                        printf("\n");
                }
        }
        return;
}

/*
 * print_subspacehmat: print the krylov space hmat, v.Hv
 */
void print_subspacehmat(double **vhv, int d)
{
        int i = 0, j = 0;
        if (mpi_proc_rank == 0) {
                printf(" Subspace v.Hv: (dimension: %d)\n", d);
                for (i = 0; i < d; i++) {
                        for (j = 0; j < d; j++) {
                                printf("%13.8lf", vhv[i][j]);
                        }
                        printf("\n");
                }
        }
        return;
}


#ifdef DEBUGGING
void print_vector_space(int v, int ckdim, int ndets)
{
        int v_lo[2] = {0, 0};
        int v_hi[2] = {0, 0};
        int v_ld[1] = {0};
        int pretty = 1;

        // rows
        v_lo[1] = 0;
        v_hi[1] = 20;
        // cols
        v_lo[0] = 0;
        v_hi[0] = ckdim - 1;

        NGA_Print_patch(v, v_lo, v_hi, pretty); 
        return;

}
void test_new_vector_space(int v, int ckdim, int ndets, int nv)
{

        int v_lo[2] = {0, 0};
        int v_hi[2] = {0, 0};
        int v_ld[1] = {0};
        int pretty = 1;

        // rows
        v_lo[0] = 0;
        v_hi[0] = 20;
        // cols
        v_lo[1] = 0;
        v_hi[1] = ckdim;

        NGA_Print_patch(v, v_lo, v_hi, pretty); 
        return;
}
#endif

