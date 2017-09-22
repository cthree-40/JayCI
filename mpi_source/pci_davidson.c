// File: pci_davidson.c
/*
 * Davidson algorithm for pjayci.c
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "ioutil.h"
#include "arrayutil.h"
#include "allocate_mem.h"
#include "initguess_sbd.h"
#include "pci_davidson.h"

/*
 * pdvdalg: parallel davidson algorithm. Hv=c is parallelized.
 */
int pdvdalg(struct det *dlist, int ndets, double *moints1, double *moints2,
	    int aelec, int belec, int ninto, double totfrze, int maxiter,
	    int krymin, int krymax, int nroots, double restol,
	    struct rowmap *hmap, int lwrbnd, int uppbnd, int predr, int plvl,
	    int mpi_proc_rank, int nrows)
{
	int error = 0; /* Error flag */
	int const mpi_root = 0; /* Root process is always 0 */
	
	double *hdgls = NULL; /* diagonal matrix elements <i|H|i> */

        double **hscr = NULL; /* krylov space hamiltonian v.Hv array */
        double *hscr_1d = NULL; /* one-dimensional form of v.Hv for interface
                                 * with F90 */
        double **cscr = NULL; /* Hv=c vectors scratch array */
        double **vscr = NULL; /* v vectors scratch array */
        double *vscr_1d = NULL; /* one-dimensional form of basis vectors */
        double *tscr_1d = NULL; /* two-dimensional scratch for truncated basis
                                 * vectors */
        double **hevec = NULL; /* krylov eigenvectors array */
        double *hevec_1d = NULL; /* one-dimensional form of eigenvectors */
        double *heval = NULL; /* krylov eigenvalues array */
        double *rvec = NULL; /* residual vector: Hv - (v.Hv)c = r */
        double *nvec = NULL; /* new vector */
        double rnorm = 0.0;

        double **mpi_cvec = NULL; /* hv=c (mpi portion) */
        double **mpi_vvec = NULL; /* hv=c (mpi portion) */
        
        int croot = 0;
        int ckdim = 0;
        int citer = 0;
        int i;

        /* Allocate scratch arrays. Dimensions are kept constant throughout
         * routine. These will only be allocated, accessed, and freed on the
         * root process. */
        if (mpi_proc_rank == mpi_root) {
                if (plvl > 0) printf("Allocating memory.\n");
                error = allocate_mem_double(&hscr, krymax, krymax);
                error = error + allocate_mem_double(&cscr, ndets, krymax);
                error = error + allocate_mem_double(&vscr, ndets, krymax);
                error = error + allocate_mem_double(&hevec, krymax, krymax);
                if (error != 0) {
                        error_flag(error,"mem allocation: dvdalg");
                        return error;
                }
                init_dbl_2darray_0(hscr, krymax, krymax);
                init_dbl_2darray_0(cscr, ndets, krymax);
                init_dbl_2darray_0(vscr, ndets, krymax);
                init_dbl_2darray_0(hevec, krymax, krymax);
     
                heval = (double *) malloc(krymax * sizeof(double));
                init_dbl_array_0(heval, krymax);
                rvec = (double *) malloc(ndets * sizeof(double));
                init_dbl_array_0(rvec, ndets);
                nvec = (double *) malloc(ndets * sizeof(double));
                init_dbl_array_0(nvec, ndets);
                hscr_1d = (double *) malloc(krymax * krymax * sizeof(double));
                init_dbl_array_0(hscr_1d, (krymax * krymax));
                hevec_1d = (double *) malloc(krymax * krymax * sizeof(double));
                init_dbl_array_0(hevec_1d, (krymax * krymax));
                vscr_1d = (double *) malloc(krymax * ndets * sizeof(double));
                init_dbl_array_0(vscr_1d, (krymax * krymax));
                tscr_1d = (double *) malloc(krymax * ndets * sizeof(double));
                init_dbl_array_0(tscr_1d, (krymax * krymax));
        }

        /* Allocate arrays to compute Hv=c for each mpi process. */
        error = allocate_mem_double(&mpi_cvec, ndets, krymax);
        error = error + allocate_mem_double(&mpi_vvec, ndets, krymax);
        if (error != 0) {
                error_flag(error, "mem allocation: dvdalg");
                return error;
        }
        
        /* Compute diagonal matrix elements <i|H|i>. This is done on the root
	 * process. It is needed for the Davidson algorithm. */
	if (mpi_proc_rank == mpi_root) {
		if (plvl > 0) printf("Computing diagonal elements.\n");
		hdgls = (double *) malloc(ndets * sizeof(double));
		compute_hdgls(detlist, ndets, moints1, moints2, aelec, belec,
			      hdgls, ninto);
		if (plvl > 0) printf("Hartree Fock energy = %15.8lf\n",
				     (hdgls[0] + totfrze));
	}

        /* Generate initial guess vectors by diagonalizing an upper block of
         * the Hamiltonian. Send the initial guess vectors to the slave nodes */
        if (mpi_proc_rank == mpi_root) {
                if (plvl > 0) printf("Computing initial guess.\n");
                error = initguess_sbd(dlist, ndets, moints1, moints2, aelec,
                                      belec, ninto, krymin, mpi_vvec);
                if (error != 0) {
                        error_flag(error, "pdvdalg");
                        return error;
                }
                
        }
        for (i = 0; i < krymin; i++) {
                MPI_Bcast(mpi_vvec[i], ndets, MPI_DOUBLE, mpi_root,
                          MPI_COMM_WORLD);
        }

        /* Perform H[n,n+rows;ndets]v[ndets,i]=c[ndets,i] for row slice of H */
        for (i = 0; i < krymin; i++) {
                compute_hv(detlist, ndets, moints1, moints2, aelec, belec,
                           mpi_vvec[i], mpi_cvec[i], ninto, hmap);
                /* Collect c contributions, and add them for each process */
                
        }
        
}
