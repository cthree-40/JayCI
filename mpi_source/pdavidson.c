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

#include <mpi.h>
#include <ga.h>
#include <macdecls.h>

/* -- OpenMP options -- */
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif
/* -------------------- */

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
              int krymax, int nroots, int prediagr, int refdim, double restol,
              int ga_buffer_len)
{
        int v_hndl = 0;           /* GLOBAL basis vectors, V */
        int v_dims[2]  = {0, 0};  /* GLOBAL basis vectors dimensions */
        int v_chunk[2] = {0, 0};  /* GLOBAL basis vectors chunk sizes */
        int c_hndl = 0;           /* GLOBAL Hv=c  vectors, C */
        int n_hndl = 0;           /* GLOBAL new vector, N */
        int n_dims[1]  = {0};     /* GLOBAL new vector dimensions */
        int n_chunk[1] = {0};     /* GLOBAL new vector chunk size */
        int r_hndl = 0;           /* GLOBAL residual vector, R */
        int x_hndl = 0;           /* GLOBAL 1-D scratch array */
        int d_hndl = 0;           /* GLOBAL <i|H|i> vector, D */

        int w_hndl = 0;           /* GLOBAL |i> = |(pq, p, q)> array */
        int w_dims[2] = {0, 0};   /* GLOBAL |i> = |(pq, p, q)> array */
        int w_chunk[2] = {0, 0};  /* GLOBAL |i> = |(pq, p, q)> chunk sizes */
        
        double *d_local = NULL;   /* LOCAL <i|H|i> array. */

        double **vhv = NULL;      /* LOCAL v.Hv array */
        double *vhv_data = NULL;  /* LOCAL v.Hv memory block */
        double *vhv_scr  = NULL;  /* LOCAL v.Hv scratch array */
        
        double **hevec = NULL;    /* LOCAL v.Hv eigenvectors */
        double *hevec_data = NULL;/* LOCAL v.Hv eigenvectors memory block */
        double *heval = NULL;     /* LOCAL v.Hv eigenvalues */
        double *hevec_scr= NULL;  /* LOCAL v.Hv eigenvector scratch array */

        double rnorm = 0.0;       /* ||r|| */
        double nnorm = 0.0;       /* ||n|| */
        
        int citer = 0; /* current iteration */
        int croot = 0; /* current root */
        int ckdim = 0; /* current dimension of krylo space */
        int cflag = 0; /* convergence flag */
        
        int lo[2] = {0, 0}; 
        int hi[2] = {0, 0};
        int ld_ndets[1] = {0};

        double totcore_e = 0.0;
        double memusage = 0.0;  /* Estimated memory usage */
        int error = 0;
        
        totcore_e = nucrep_e + frzcore_e;
        ld_ndets[0] = ndets;

        /* Allocate GLOBAL arrays: V, Hv=c, N, R, and D */
        if (mpi_proc_rank == mpi_root) {
		printf("Creating global arrays...\n");
		memusage = krymax * ndets * 8 * 2; // C and V
		memusage += ndets * 8 * 3; // N, R, X, D
		memusage += ndets * 4 * 3; // W
		memusage = memusage / 1048576;
		printf(" Global Arrays memory usage: ");
		printf("  %10.2lf MB\n", memusage);
		memusage = ga_buffer_len * krymax * 8;
		memusage = memusage + ga_buffer_len * 4 * 3;
		memusage = memusage / 1048576;
		memusage = memusage * mpi_num_procs;
		printf(" Local buffer usage: ");
		printf("  %10.2lf MB\n", memusage);
		printf("  buflen = %d\n", ga_buffer_len);
		fflush(stdout);
	}
	GA_Sync();
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
        x_hndl = NGA_Duplicate(n_hndl, "Scratch 1-D vector");
        if (!x_hndl) GA_Error("Duplicate failed: Scratch 1-D vector", 1);
        d_hndl = NGA_Duplicate(n_hndl, "Diagonal vector");
        if (!d_hndl) GA_Error("Duplicate failed: Diagonal vectors", 1);
        /* Set wavefunction global array */
        w_dims[0] = ndets;
        w_dims[1] = 3;
        w_chunk[0] = -1; // Distribute evenly
        w_chunk[1] =  3;
        w_hndl = NGA_Create(C_INT, 2, w_dims, "Determinant Triples",w_chunk);
        if (!w_hndl) GA_Error("Create failed: Determinant Triples", 2);

        if (mpi_proc_rank == mpi_root) {
		printf("Global arrays created.\n\n");
		fflush(stdout);
	}
	
        /* Generate wavefunction list */
	if (mpi_proc_rank == mpi_root) {
		printf("Generating wlist...\n");
		fflush(stdout);
	}
	generate_wlist(w_hndl, ndets, pq_space_pairs, num_pq, peospace, pegrps,
                       qeospace, qegrps);
        
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
	if (mpi_proc_rank == mpi_root) {
		printf("Computing diagonal matrix elements...\n");
		fflush(stdout);
	}
	compute_diagonal_matrix_elements(d_local, lo[0], hi[0], moints1, moints2,
                                         aelec, belec, peospace, pegrps,
                                         qeospace, qegrps, pq_space_pairs,
                                         num_pq, pstrings, qstrings, intorb);
        NGA_Put(d_hndl, lo, hi, d_local, ld_ndets);
        if (printlvl > 0 && lo[0] == 0) {
                printf("<1|H|1> = ");
                printf("%lf\n", (d_local[0] + totcore_e));
        }

        if (mpi_proc_rank == mpi_root) {
                printf("\nBeginning Davidson algorithm...\n");
		fflush(stdout);
	}
        
        /* Build initial guess basis vectors. */
        build_init_guess_vectors(prediagr, v_hndl, refdim, krymin, ndets,
                                 pstrings, peospace, pegrps, qstrings,
                                 qeospace, qegrps, pq_space_pairs, num_pq,
                                 moints1, moints2, aelec, belec, intorb, w_hndl);
        if (mpi_proc_rank == mpi_root) {
                printf(" Initial guess vectors set.\n");
		fflush(stdout);
	}
        
        GA_Sync();

        /* .. MAIN LOOP .. */
        citer = 1; croot = 1; ckdim = krymin; cflag = 0;
        while (citer < maxiter && croot <= nroots) {

                GA_Sync();
                perform_hv_initspace(pstrings, peospace, pegrps, qstrings,
                                     qeospace, qegrps, pq_space_pairs, num_pq,
                                     moints1, moints2, aelec, belec, intorb,
                                     ndets, totcore_e, ckdim, krymax, v_hndl, d_hndl,
                                     c_hndl, w_hndl, ga_buffer_len);
                make_subspacehmat_ga(v_hndl, c_hndl, ndets, ckdim, vhv);
                print_subspacehmat(vhv, ckdim);
                error = diag_subspacehmat(vhv, hevec, heval, ckdim, krymax,
                                          vhv_scr, hevec_scr);
                if (error != 0)  return error;
                print_subspace_eigeninfo(hevec, heval, ckdim, totcore_e);

                GA_Sync();

                while (ckdim < krymax) {

                        GA_Sync();

                        generate_residual(v_hndl, c_hndl, r_hndl, hevec, heval,
                                          ndets, ckdim, croot, x_hndl);
                        compute_GA_norm(r_hndl, &rnorm);
                        print_iter_info(heval, ckdim, croot, rnorm, totcore_e);
                        
                        GA_Sync();
                        
                        cflag = test_convergence(rnorm, restol, croot, nroots);
                        if (cflag == 2) {
                                if (mpi_proc_rank == mpi_root) {
                                        printf("\n ** CI CONVERGED! **\n");
					fflush(stdout);
				}
                                break;
                        } else if (cflag == 1) {
                                if (mpi_proc_rank == mpi_root) {
                                        printf(" Root %d converged!\n", croot);
					fflush(stdout);
				}
                                croot++;
                                break;
                        } else {
                                if (mpi_proc_rank == mpi_root) {
                                        printf(" Root not yet converged.\n");
					fflush(stdout);
				}
                        }

                        GA_Sync();
                        
                        generate_newvector(r_hndl, d_hndl, heval[croot - 1],
                                           ndets, n_hndl, x_hndl);

                        compute_GA_norm(n_hndl, &nnorm);
                        if (mpi_proc_rank == mpi_root) {
                                printf("\n ||r|| = %12.8lf  ||n|| = %12.8lf\n",
                                       rnorm, nnorm);
				fflush(stdout);
                        }

                        orthonormalize_newvector(v_hndl, ckdim, ndets, n_hndl);

                        /* Increase current krylov space dimension */
                        ckdim++;

                        add_new_vector(v_hndl, ckdim, ndets, n_hndl);
                        compute_hv_newvector(v_hndl, c_hndl, ckdim, pstrings,
                                             peospace, pegrps, qstrings,
                                             qeospace, qegrps, pq_space_pairs,
                                             num_pq, moints1, moints2, aelec,
                                             belec, intorb, ndets, krymax,
                                             w_hndl, ga_buffer_len);
                        make_subspacehmat_ga(v_hndl, c_hndl, ndets, ckdim, vhv);
                        print_subspacehmat(vhv, ckdim);
                        error = diag_subspacehmat(vhv, hevec, heval, ckdim,
                                                  krymax, vhv_scr, hevec_scr);
                        if (error != 0)  return error;
                        print_subspace_eigeninfo(hevec, heval, ckdim, totcore_e);
                        citer++;
                }
                /* truncate the krylov space. Note: the Hv=c array is used
                 * as a scratch buffer for this routine. */
                truncate_krylov_space(v_hndl, ndets, krymin, krymax, croot,
                                      hevec, c_hndl);
                ckdim = krymin;
                GA_Zero(c_hndl);
                /* Check if CI has converged. If it has, leave loop. */
                if (cflag == 2) {
                        break;
                }
        }
        if (mpi_proc_rank == mpi_root) {
                printf(" Davidson algorithm finished. \n");
		fflush(stdout);
        }
        print_gavectors2file_dbl_ufmt(v_hndl, ndets, nroots,"civec");
        if (cflag != 2)
                print_gavectors2file_dbl_ufmt(c_hndl, ndets, nroots, "hvvec");
        return error;
}

/*
 * add_new_vector: add a new vector to basis space.
 * cdim = dimension INCLUDING new vector.
 */
void add_new_vector(int v_hndl, int cdim, int len, int n_hndl)
{
        int vlo[2] = {0, 0};
        int vhi[2] = {0, 0};
        int nlo[1] = {0};
        int nhi[1] = {0};
        vlo[0] = cdim - 1; vhi[0] = cdim - 1;
        vhi[1] = len - 1;
        nhi[0] = len - 1;
        NGA_Copy_patch('n', n_hndl, nlo, nhi, v_hndl, vlo, vhi);
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
                              int aelec, int belec, int intorb, int w_hndl)
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
                init_diag_H_subspace(w_hndl,
                                     pstr, qstr, m1, m2, aelec, belec, intorb,
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
 * Reads buffer of V, rather than entire array lengths.
 * Input:
 *  v_hndl = GA handle for basis vectors (ckdim is new vector)
 *  c_hndl = GA handle for Hv=c vectors
 *  ckdim  = current dimension of space
 *  w_hndl = wavefunction list (deteriminant triplets)
 */
void compute_hv_newvector(int v_hndl, int c_hndl, int ckdim, struct occstr *pstr,
                          struct eospace *peosp, int pegrps, struct occstr *qstr,
                          struct eospace *qeosp, int qegrps, int **pqs,
                          int num_pq, double *m1, double *m2, int aelec,
                          int belec, int intorb, int ndets, int kmax, int w_hndl,
                          int ga_buffer_len)
{
        
        int error = 0; /* Error flag */

        double *c_local = NULL;   /* Local c array */
        int c_rows = 0;           /* Local c rows  */
        int c_cols = 0;           /* Local c columns (=1) */

        double *v_local = NULL;   /* Local v array */
        int buflen = 0;
        int v_rows = 0;           /* Local v rows  */
        int v_cols = 0;           /* Local v columns (=1) */

        /*
         * The following convention is used:
         *    H(i,j)*V(j,k)=C(i,k)
         */
        int **wi       = NULL;     /* Local w array */
        int *widata    = NULL;     /* Local w array data (1-D) */
        int **wj       = NULL;     /* Local w array */
        int *wjdata    = NULL;     /* Local w array data (1-D) */
        int wjlen  = 0;            /* Triplets in w to evaluate */
        
        
        int c_lo[2]  = {0, 0};  /* starting indices for memory block */
        int c_hi[2]  = {0, 0};  /* ending indices of memory block */
        int v_lo[2]  = {0, 0};  /* starting indices of memory block */
        int v_hi[2]  = {0, 0};  /* ending indices of memory block */
        int wi_lo[2] = {0, 0};
        int wi_hi[2] = {0, 0};
        int wj_lo[2] = {0, 0};
        int wj_hi[2] = {0, 0};
        
        int v_ld[1]  = {0};  /* Leading dimensions of V local buffer */
        int c_ld[1]  = {0};  /* Leading dimensions of C local buffer */
        int wi_ld[1] = {0};
        int wj_ld[1] = {0};

        double alpha[1] = {1.0}; /* Scale factor for c_local into c_global. */

        int i, j, jmax;
        
        /*
         * Determine which block of data is locally owned. And get the blocks
         * of V that are required to compute c.
         * Recall: C is split by rows, and not columns, across processes.
         *         We only require the last columns of both V and C.
         */
        NGA_Distribution(c_hndl, mpi_proc_rank, c_lo, c_hi);
        c_rows  = c_hi[1] - c_lo[1] + 1;
        c_lo[0] = ckdim - 1;  // Last column
        c_hi[0] = ckdim - 1;  // Last column
        c_cols  = c_hi[0] - c_lo[0] + 1;
        c_ld[0] = 1;          // 1-D array
        /* Allocate local array and get C data */
        c_local = malloc(sizeof(double) * c_rows);
        NGA_Get(c_hndl, c_lo, c_hi, c_local, c_ld);
                
        /* Allocate local Wi array and get W data */
        widata = allocate_mem_int_cont(&wi, 3, c_rows);
        wi_lo[0] = c_lo[1];
        wi_lo[1] = 0;
        wi_hi[0] = c_hi[1];
        wi_hi[1] = 2;
        wi_ld[0] = 3;
        NGA_Get(w_hndl, wi_lo, wi_hi, widata, wi_ld); 

        v_lo[0] = ckdim - 1;
        v_hi[0] = ckdim - 1;
        v_lo[1] = 0;
        v_hi[1] = ndets - 1;
        v_rows  = v_hi[1] - v_lo[1] + 1;
        v_cols  = v_hi[0] - v_hi[0] + 1;
        v_ld[0] = 1;
        /* Allocate local array */
        buflen = ga_buffer_len;
        v_local = malloc(sizeof(double) * buflen);

        
        /* Allocate local Wj array  */
        wjdata = allocate_mem_int_cont(&wj, 3, buflen);
        wj_lo[1] = 0;
        wj_hi[1] = 2;
        wj_ld[0] = 3;

        GA_Sync();

        for (j = 0; j < ndets; j += buflen) {

                jmax = int_min((j + buflen - 1), (ndets - 1));

                v_lo[1] = j;
                v_hi[1] = jmax;
                v_rows = v_hi[1] - v_lo[1] + 1;
                NGA_Get(v_hndl, v_lo, v_hi, v_local, v_ld);

                wj_lo[0] = j;
                wj_hi[0] = jmax;
                wjlen = jmax - j + 1;
                NGA_Get(w_hndl, wj_lo, wj_hi, wjdata, wj_ld);
                
                
                /* Evaluate block of H(i,j)*V(j,k)=C(i,k). This is done via the   
                 * space indexes. Each block is over all V vectors, k.*/
                evaluate_hdblock_ij_1d2(wi, c_rows, wj, wjlen,
                                        v_rows, v_cols, v_local,
                                        c_rows, c_cols, c_local,
                                        ndets, peosp, pegrps, pstr, qeosp, qegrps,
                                        qstr, pqs, num_pq, m1, m2, aelec, belec,
                                        intorb);

        }
        NGA_Acc(c_hndl, c_lo, c_hi, c_local, c_ld, alpha);

        free(v_local);

        GA_Sync();
        
        return;
}


/*
 * compute_GA_norm: compute the norm of a global-array vector.
 */
void compute_GA_norm (int r_hndl,  double *norm)
{
        *norm = GA_Ddot(r_hndl, r_hndl);
        *norm = sqrt(*norm);
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

        /* OMP Section */
#pragma omp parallel                                                               \
        shared(ndeti,ndetj,mo1,mo2,aelec,belec,intorb,c,v,pstr,qstr,dtj,d_triplet) \
        private(deti,detj,i,j,l,hijval)
        {
#pragma omp for schedule(runtime)
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
#pragma omp atomic update
                                c[l][i] = c[l][i] + hijval * v[l][j];
                        }
                }
        }
        } /* End of OMP Section */
        return;
}

/*
 * evaluate_hdblock_ij2: evaluate a block of H:
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
void evaluate_hdblock_ij2(int **wi, int idets, int **wj, int jdets,
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
        struct det deti;        /* Determinant i */
        struct det detj;        /* Determinant j */

        double hijval = 0.0;    /* <i|H|j> value */
        
        int i, j, l;
        
        /* OMP Section */
#pragma omp parallel                                                               \
        shared(idets,jdets,mo1,mo2,aelec,belec,intorb,c,v,pstr,qstr,wi,wj,vcols) \
        private(deti,detj,i,j,l,hijval)
        {
#pragma omp for schedule(runtime)
        /* Loop through list of triplets for determinants |i>. */
        for (i = 0; i < idets; i++) {
                deti.astr = pstr[wi[i][0]];
                deti.bstr = qstr[wi[i][1]];
                deti.cas = wi[i][2];
                /* Loop over determinants |j> */
                for (j = 0; j < jdets; j++) {
                        detj.astr = pstr[wj[j][0]];
                        detj.bstr = qstr[wj[j][1]];
                        detj.cas  = wj[j][2];
                        
                        hijval = hmatels(deti, detj, mo1, mo2,
                                         aelec, belec, intorb);
                        /* H_ij*v_jl = c_il */
                        for (l = 0; l < vcols; l++) {
#pragma omp atomic update
                                c[l][i] = c[l][i] + hijval * v[l][j];
                        }
                }
        }
        } /* End of OMP Section */
        return;
}

/*
 * evaluate_hdblock_ij_1d: evaluate a block of H: **(For one column)**
 *  H(i,j)*V(j)=C(i)
 *  Input:
 *   vrows  = columns of C
 *   vcols  = rows of C
 *   crows  = columns of C
 *   ccols  = rows of C
 *   v      = V(j)
 *   c      = C(i)
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
 *   C      = C(i)
 */
void evaluate_hdblock_ij_1d(int pq_start_i, int pstart_i, int qstart_i,
                            int pq_final_i, int pfinal_i, int qfinal_i,
                            int pq_start_j, int pstart_j, int qstart_j,
                            int pq_final_j, int pfinal_j, int qfinal_j,
                            int vrows, int vcols, double *v,
                            int crows, int ccols, double *c,
                            int starti, int finali, int startj, int finalj,
                            int ndets, struct eospace *peosp, int pegrps,
                            struct occstr *pstr,
                            struct eospace *qeosp, int qegrps,
                            struct occstr *qstr,  int **pq,
                            int npq, double *mo1, double *mo2, int aelec,
                            int belec, int intorb)
{

        int **d_triplet = NULL; /* |i> = (p, q, CAS-flag) list*/
        int *d_trip_dat = NULL; /* d_triplet memory block */
        int **dtj = NULL;
        int *dtjdat=NULL;
        int ndetj;
        
        struct det deti;        /* Determinant i */
        struct det detj;        /* Determinant j */
        int ndeti;              /* Number of i determinants (chunk size) */

        double hijval = 0.0;    /* <i|H|j> value */
        
        int i = 0, j = 0, l = 0;
        
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

        /* OMP Section */
#pragma omp parallel                                                               \
        shared(ndeti,ndetj,mo1,mo2,aelec,belec,intorb,c,v,pstr,qstr,dtj,d_triplet) \
        private(deti,detj,i,j,l,hijval)
        {
#pragma omp for schedule(runtime)
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
                        /* H_ij*v_j = c_i */
#pragma omp atomic update
                        c[i] = c[i] + hijval * v[j];
                }
        }
        }
        return;
}

/*
 * evaluate_hdblock_ij_1d2: evaluate a block of H: **(For one column)**
 *  H(i,j)*V(j)=C(i)
 *  Input:
 *   vrows  = columns of C
 *   vcols  = rows of C
 *   crows  = columns of C
 *   ccols  = rows of C
 *   v      = V(j)
 *   c      = C(i)
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
 *   C      = C(i)
 */
void evaluate_hdblock_ij_1d2(int **wi, int idets, int **wj, int jdets,
                             int vrows, int vcols, double *v,
                             int crows, int ccols, double *c,
                             int ndets, struct eospace *peosp, int pegrps,
                             struct occstr *pstr,
                             struct eospace *qeosp, int qegrps,
                             struct occstr *qstr,  int **pq,
                             int npq, double *mo1, double *mo2, int aelec,
                             int belec, int intorb)
{
        struct det deti;        /* Determinant i */
        struct det detj;        /* Determinant j */

        double hijval = 0.0;    /* <i|H|j> value */
        
        int i, j;
        
        /* OMP Section */
#pragma omp parallel                                                               \
        shared(idets,jdets,mo1,mo2,aelec,belec,intorb,c,v,pstr,qstr,wi,wj) \
        private(deti,detj,i,j,hijval)
        {
#pragma omp for schedule(runtime)
        /* Loop through list of triplets for determinants |i>. */
        for (i = 0; i < idets; i++) {
                deti.astr = pstr[wi[i][0]];
                deti.bstr = qstr[wi[i][1]];
                deti.cas = wi[i][2];
                /* Loop over determinants |j> */
                for (j = 0; j < jdets; j++) {
                        detj.astr = pstr[wj[j][0]];
                        detj.bstr = qstr[wj[j][1]];
                        detj.cas  = wj[j][2];
                        
                        hijval = hmatels(deti, detj, mo1, mo2,
                                         aelec, belec, intorb);
                        /* H_ij*v_j = c_i */
#pragma omp atomic update
                        c[i] = c[i] + hijval * v[j];
                }
        }
        }
        return;
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
 *  x_hndl = GA handle for scratch array
 */
void generate_newvector (int r_hndl, int d_hndl, double eval, int ndets,
                         int n_hndl, int x_hndl)
{
        double e = 0.0;
        double neg1 = -1.0;
        
        NGA_Zero(x_hndl);
        GA_Copy(d_hndl, x_hndl);
        e = neg1 * eval;
        GA_Add_constant(x_hndl, &e);
        GA_Scale(x_hndl, &neg1);
        GA_Elem_divide(r_hndl, x_hndl, n_hndl);
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
 *  rscr_hndl = scratch vector global array
 */
void generate_residual (int v_hndl, int c_hndl, int r_hndl, double **hevec,
                        double *heval, int ndets, int ckdim, int croot,
                        int rscr_hndl)
{
        int root_id = 0; /* Index of root. (croot - 1) */
        /* For each column, j, alpha and beta equal
         * alpha = v_aj; beta = v_aj*e_a*(-1) */
        double alpha = 0.0, beta = 0.0;
        int alo[2] = {0, 0}, blo[2] = {0, 0}, clo[1] = {0};
        int ahi[2] = {0, 0}, bhi[2] = {0, 0}, chi[1] = {0};
        int i = 0;

        /* Set second indices for *hi for ndets (whole vector) */
        ahi[1] = ndets - 1;
        bhi[1] = ndets - 1;
        chi[0] = ndets - 1;
        
        root_id = croot - 1;
        NGA_Zero(r_hndl);
        for (i = 0; i < ckdim; i++) {
                NGA_Zero(rscr_hndl);
                alpha = hevec[root_id][i];
                beta  = hevec[root_id][i] * heval[root_id] * (-1.0);
                /* r = c*alpha + v*beta */
                alo[0] = i; ahi[0] = i;
                blo[0] = i; bhi[0] = i;
                NGA_Add_patch(&alpha, c_hndl, alo, ahi, &beta, v_hndl, blo, bhi,
                              rscr_hndl, clo, chi);
                /* Sum */
                alpha = 1.0;
                beta  = 1.0;
                GA_Add(&alpha, rscr_hndl, &beta, r_hndl, r_hndl);
        }
        return;
}

/*
 * generate_wlist: generate the wavefunction list of triplets
 * Input:
 *  hndl     = handle of global array of W
 *  ndets    = number of determinants
 *  pq       = alpha/beta pairs
 *  npq      = number of alpha/beta pairs
 *  peospace = alpha string electron orbital space
 *  pegrps   = number of alpha string orbital spaces
 *  qeospace = beta  string electron orbital space
 *  qegrps   = number of beta  string orbital spaces
 */
void generate_wlist(int hndl, int ndets, int **pq, int npq,
                    struct eospace *peosp, int pegrps,
                    struct eospace *qeosp, int qegrps)
{
        int lo[2] = {0, 0};
        int hi[2] = {0, 0};
        int ld[1] = {0};       /* Leading dimension of w buffer */
        int rows = 0, cols = 0;/* Dimensions of local w array */
        int **w = NULL;        /* Local w array */
        int *wdata = NULL;     /* Local w array data (1-d) */

        int pq_start = 0, pq_final = 0;
        int pstart = 0, pfinal = 0;
        int qstart = 0, qfinal = 0;
        
        /* Find distribution of W that is allocate on this process.
         * Allocate local W array */
        NGA_Distribution(hndl, mpi_proc_rank, lo, hi);
        cols = hi[0] - lo[0] + 1;
        rows = hi[1] - lo[1] + 1;
        if (rows != 3) {
                error_message(mpi_proc_rank, "rows != 3","generate_wlist");
                error_flag(mpi_proc_rank, rows, "generate_wlist");
                return;
        }
        wdata = allocate_mem_int_cont(&w, rows, cols);
        
        /* Get starting deteriminant index values, and generate the determinant
         * triplet list for this section. */
        determinant_string_info(lo[0], peosp, pegrps, qeosp, qegrps, pq, npq,
                                &pq_start, &pstart, &qstart);
        determinant_string_info(hi[0], peosp, pegrps, qeosp, qegrps, pq, npq,
                                &pq_final, &pfinal, &qfinal);
        generate_det_triples(cols, w, pq_start, pstart, qstart, pq_final,
                             pfinal, qfinal, pq, npq, peosp, pegrps,
                             qeosp, qegrps);

        /* Send this info to the global array W. */
        ld[0] = rows;
        NGA_Put(hndl, lo, hi, wdata, ld);

        deallocate_mem_cont_int(&w, wdata);

        GA_Sync();
                
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
void init_diag_H_subspace( int w_hndl, struct occstr *pstr, struct occstr *qstr,
                           double *m1, double *m2,
                           int aelec, int belec, int intorb, int ndets, int dim,
                           double **refspace)
{
        struct det deti;
        struct det detj;
        double **hij = NULL;
        double *hij_data = NULL;
        double *rdata = NULL;
        double *rev = NULL;

        int **w       = NULL;     /* Local w array */
        int *wdata    = NULL;     /* Local w array data (1-D) */
        int w_lo[2] = {0, 0};
        int w_hi[2] = {0, 0};
        int w_ld[1] = {0};
        
        int i, j, ii;
        int error = 0;

        /* Allocate h matrix subblock and refvec data */
        hij_data = allocate_mem_double_cont(&hij, dim, dim);
        rdata = malloc(sizeof(double) * dim * dim);
        rev = malloc(sizeof(double) * dim);

        /* Allocate w array and get wavefunction information */
        wdata = allocate_mem_int_cont(&w, 3, dim);
        w_lo[0] = 0;
        w_lo[1] = 0;
        w_hi[0] = dim - 1;
        w_hi[1] = 2;
        w_ld[0] = 3;
        NGA_Get(w_hndl, w_lo, w_hi, wdata, w_ld);

        /* Loop through list of triplets for determinants <i| . */
        for (i = 0; i < dim; i++) {
                deti.astr = pstr[w[i][0]];
                deti.bstr = qstr[w[i][1]];
                deti.cas = w[i][2];
                /* Loop over determinants |j> */
                for (j = 0; j < dim; j++) {
                        detj.astr = pstr[w[j][0]];
                        detj.bstr = qstr[w[j][1]];
                        detj.cas = w[j][2];

                        hij[i][j] = hmatels(deti, detj, m1, m2, aelec,
                                            belec, intorb);
                }
        }

        /* Diagonalize hij */
        error = diagmat_dsyevr(hij_data, dim, rdata, rev);
        if (error != 0) {
                error_message(mpi_proc_rank, "Error occured during DSYEVR",
                              "init_diag_H_subspace");
        }
        printf(" eigenvalues = %lf %lf %lf\n",
               rev[0] + total_core_e,
               rev[1] + total_core_e,
               rev[2] + total_core_e);
	fflush(stdout);
        /* Copy rdata into refspace */
        ii = 0;
        for (i = 0; i < dim; i++) {
                for (j = 0; j < dim; j++) {
                        refspace[i][j] = rdata[ii];
                        ii++;
                }
        }

        deallocate_mem_cont(&hij, hij_data);
        deallocate_mem_cont_int(&w, wdata);
        free(rdata);
        free(rev);
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
        double *overlaps = NULL;
        int vlo[2] = {0, 0}, vhi[2] = {0, 0};
        int nlo[1] = {0},    nhi[1] = {0};
        double alpha = 1.0;
        int i = 0;
        overlaps = malloc(sizeof(double) * nvecs);
        vhi[1] = ndets - 1;
        nhi[0] = ndets - 1;
        for (i = 0; i < nvecs; i++) {
                vlo[0] = i; vhi[0] = i;
                overlaps[i] = NGA_Ddot_patch(v_hndl, 'n', vlo, vhi,
                                             n_hndl, 'n', nlo, nhi);
        }
        for (i = 0; i < nvecs; i++) {
                overlaps[i] = overlaps[i] * (-1.0);
                vlo[0] = i; vhi[0] = i;
                NGA_Add_patch(&alpha, n_hndl, nlo, nhi, &overlaps[i], v_hndl,
                              vlo, vhi, n_hndl, nlo, nhi);
        }
        /* Get normalize the new, now orthogonal vector */
        alpha = GA_Ddot(n_hndl, n_hndl);
        alpha = sqrt(alpha);
        alpha = 1 / alpha;
        GA_Scale(n_hndl, &alpha);

        /* Check overlaps */
        for (i = 0; i < nvecs; i++) {
                vlo[0] = i; vhi[0] = i;
                overlaps[i] = NGA_Ddot_patch(v_hndl, 'n', vlo, vhi,
                                             n_hndl, 'n', nlo, nhi);
        }
        for (i = 0; i < nvecs; i++) {
                if (overlaps[i] > 0.000001) {
                        error_message(mpi_proc_rank,
                                      "Warning! Non-zero overlap.",
                                      "orthonormalize_newvector");
                }
        }
        /* Check norm */
        alpha = GA_Ddot(n_hndl, n_hndl);
        alpha = sqrt(alpha);
        if ((alpha - 1.0) > 0.000001) {
                error_message(mpi_proc_rank,
                              "Warning! New vector norm != 1.0",
                              "orthonormalize_newvector");
        }
        free(overlaps);
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
 *  w_hndl= (GLOBAL ARRAY HANDLE) wavefunction
 * Output:
 *  c_hndl= (GLOBAL ARRAY HANDLE) Hv=c vectors
 */
void perform_hv_initspace(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets,
                          double core_e, int dim, int mdim, int v_hndl, int d_hndl,
                          int c_hndl, int w_hndl, int ga_buffer_len)
{
        
        double **c_local = NULL;  /* Local c array */
        double *cdata = NULL;     /* Local c array data */
        int c_rows = 0;           /* Local c rows */
        int c_cols = 0;           /* Local c columns */

        double **v_local = NULL;  /* Local v array */
        double *vdata = NULL;     /* Local v array data */
        double buflen = 0;        /* Length of ONE COLUMN of v buffer */
        int v_rows = 0;           /* Local v rows */
        int v_cols = 0;           /* LOcal v columns */

        /*
         * The following convention is used:
         *    H(i,j)*V(j,k)=C(i,k)
         */
        int **wi       = NULL;     /* Local w array */
        int *widata    = NULL;     /* Local w array data (1-D) */
        int **wj       = NULL;     /* Local w array */
        int *wjdata    = NULL;     /* Local w array data (1-D) */
        int wjlen  = 0;            /* Triplets in w to evaluate */

        int start_det_i = 0;    /* Starting determinant index of i */
        int final_det_i = 0;    /* Ending determinant index of i   */
        int start_det_j = 0;    /* Starting determinant index of j */
        int final_det_j = 0;    /* Ending determinant index of j   */
        
        int c_lo[2] = {0, 0}; /* starting indices for memory block */
        int c_hi[2] = {0, 0}; /* ending indices of memory block */
        int v_lo[2] = {0, 0}; /* starting indices of memory block */
        int v_hi[2] = {0, 0}; /* ending indices of memory block */
        int wi_lo[2] = {0, 0};
        int wi_hi[2] = {0, 0};
        int wj_lo[2] = {0, 0};
        int wj_hi[2] = {0, 0};

        int v_ld[1] = {0}; /* Leading dimensions of V local buffer */
        int c_ld[1] = {0}; /* Leading dimensions of C local buffer */
        int wi_ld[1] = {0};
        int wj_ld[1] = {0};
        
        double alpha[1] = {1.0}; /* Scale factor for c_local into c_global. */

        int j, jmax;

	if (mpi_proc_rank == mpi_root) {
		printf(" Performing Hv=c on initial vector space...\n");
		fflush(stdout);
	}
	
        NGA_Zero(c_hndl);

        /* Determine which block of data is locally owned. And get the blocks of
         * V that are required to compute c. */
        NGA_Distribution(c_hndl, mpi_proc_rank, c_lo, c_hi);
        c_cols = c_hi[0] - c_lo[0] + 1;
        c_rows = c_hi[1] - c_lo[1] + 1;
        c_cols = int_min(c_cols, dim);
        c_hi[0] = c_cols - 1;
        c_ld[0] = c_rows;
        /* Allocate local array and get C data */
        cdata = allocate_mem_double_cont(&c_local, c_rows, c_cols);
        NGA_Get(c_hndl, c_lo, c_hi, cdata, c_ld);

        /* Allocate local Wi array and get W data */
        widata = allocate_mem_int_cont(&wi, 3, c_rows);
        wi_lo[0] = c_lo[1];
        wi_lo[1] = 0;
        wi_hi[0] = c_hi[1];
        wi_hi[1] = 2;
        wi_ld[0] = 3;
        NGA_Get(w_hndl, wi_lo, wi_hi, widata, wi_ld); 

        
        buflen = ga_buffer_len; /* Set buffer size for each column */
        v_lo[0] = 0;
        v_lo[1] = 0;
        v_hi[0] = dim - 1;
        v_hi[1] = ndets - 1;
        v_cols = v_hi[0] - v_lo[0] + 1;
        //v_rows = v_hi[1] - v_lo[1] + 1;
        v_rows = buflen;
        v_ld[0] = v_rows;
        /* Allocate local array (buflen x v_cols) */
        vdata = allocate_mem_double_cont(&v_local, v_rows, v_cols);

        /* allocate local Wj array */
        wjdata = allocate_mem_int_cont(&wj, 3, buflen);
        wj_lo[1] = 0;
        wj_hi[1] = 2;
        wj_ld[0] = 3;

        GA_Sync();

        for (j = 0; j < ndets; j += buflen) {

                jmax = int_min((j + buflen - 1), (ndets - 1));

                v_lo[1] = j;
                v_hi[1] = jmax;
                v_rows = v_hi[1] - v_lo[1] + 1;
                NGA_Get(v_hndl, v_lo, v_hi, vdata, v_ld);

                wj_lo[0] = j;
                wj_hi[0] = jmax;
                wjlen = jmax - j + 1;
                NGA_Get(w_hndl, wj_lo, wj_hi, wjdata, wj_ld);

                /* Evaluate block of H(i,j)*V(j,k)=C(i,k). This is done via the
                 * space indexes. Each block is over all V vectors, k.*/
                evaluate_hdblock_ij2(wi, c_rows, wj, wjlen,
                                     v_rows, v_cols, v_local,
                                     c_rows, c_cols, c_local,
                                     start_det_i, final_det_i,
                                     start_det_j, final_det_j,
                                     ndets, peosp, pegrps, pstr, qeosp,
                                     qegrps, qstr, pqs, num_pq, m1, m2,
                                     aelec, belec, intorb);
        }
        NGA_Acc(c_hndl, c_lo, c_hi, cdata, c_ld, alpha);

        deallocate_mem_cont(&v_local, vdata);
        deallocate_mem_cont(&c_local, cdata);

        GA_Sync();
                
        return;
}

/*
 * print_iter_info: print iteration information.
 */
void print_iter_info(double *heval, int ckdim, int croot, double rnorm,
                     double totfrze)
{
	int i;
        if (mpi_proc_rank == mpi_root) {
                for (i = 0; i < 70; i++) {
                        printf("-");
                }
                printf("\n  Eigenvalues:\n");
                for (i = 0; i < ckdim; i++) {
                        if ((i + 1) < croot) {
                                printf("   Root #%2d  %15.8lf *CONVERGED*\n",
                                       (i + 1), (heval[i] + totfrze));
                        } else if ((i + 1) == croot) {
                                printf("   Root #%2d  %15.8lf ||r||: %15.8lf\n",
                                       (i + 1), (heval[i] + totfrze), rnorm);
                        } else {
                                printf("   Root #%2d  %15.8lf\n",
                                       (i + 1), (heval[i] + totfrze));
                        }
                }
                for (i = 0; i < 70; i++) {
                        printf("-");
                }
                fprintf(stdout,"\n");
        }
	fflush(stdout);
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
	fflush(stdout);
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
	fflush(stdout);
        return;
}

/*
 * test_convergence: test convergence of davidson algorithm. returns
 * convergence flag.
 * Input:
 *  rnorm  = norm of residual, ||r||
 *  restol = convergence tolerance
 *  croot  = current root being optimized
 *  nroot  = total number of roots
 *  f      = flag
 */
int test_convergence(double rnorm, double restol, int croot, int nroot)
{
        /*
         * f = 0, no convergence
         * f = 1, convergence is reached. on to croot + 1
         * f = 2, convergence is reached. calc finished. */
        int f = 0;
        if (rnorm >= restol) f = 0;
        if (rnorm < restol) {
                if (croot < nroot) {
                        f = 1;
                } else {
                        f = 2;
                }
        }
        return f;
}

/*
 * truncate_krylov_space: truncate the krylov space from krymax to krymin.
 * Input:
 *  v_hndl = GA handle of basis vectors array
 *  ndets  = number of determinants
 *  krymin = minimum dimension of krylov space
 *  krymax = maximum dimension of krylov space
 *  croot  = current root
 *  hevec  = eigenvectors
 *  x_hndl = scratch array
 */
void truncate_krylov_space(int v_hndl, int ndets, int krymin, int krymax,
                           int croot, double **hevec, int x_hndl)
{
        int vlo[2] = {0, 0}, vhi[2] = {0, 0};
        int xlo[2] = {0, 0}, xhi[2] = {0, 0};
        double alpha = 0.0, beta = 0.0;
        int i = 0, j = 0;
        if (mpi_proc_rank == mpi_root) {
                printf(" Truncating krylov space...\n");
		fflush(stdout);
        }
        NGA_Zero(x_hndl);
        vhi[1] = ndets - 1;
        xhi[1] = ndets - 1;
        /* b_i = e_ik * v_k */
        for (i = 0; i < krymin; i++) {
                xlo[0] = i;
                xhi[0] = i;
                for (j = 0; j < krymax; j++) {
                        vlo[0] = j;
                        vhi[0] = j;
                        alpha = hevec[i][j];
                        beta  = 1.0;
                        NGA_Add_patch(&alpha, v_hndl, vlo, vhi, &beta, x_hndl,
                                      xlo, xhi, x_hndl, xlo, xhi);
                }
#ifdef DEBUGGING
                print_vector_space(x_hndl, krymin, 10);
#endif
        }
        NGA_Zero(v_hndl);
        GA_Copy(x_hndl, v_hndl);
        if (mpi_proc_rank == mpi_root) {
                printf(" Completed krylov space truncation.\n");
		fflush(stdout);
        }
        return;
}

#ifdef DEBUGGING
void print_vector_space(int v, int ckdim, int ndets)
{
        int v_lo[2] = {0, 0};
        int v_hi[2] = {0, 0};
        int pretty = 1;

        // rows
        v_lo[1] = 0;
        v_hi[1] = ndets-1;
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

