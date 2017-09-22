// File: execute_pci_calculation.c
/*
 * Main driver for parallel CI calculation.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "timestamp.h"
#include "errorlib.h"
#include "allocate_mem.h"
#include "binary.h"
#include "binarystr.h"
#include "ioutil.h"
#include "cimapping.h"
#include "genbindet.h"
#include "action_util.h"
#include "execute_pci_calculation.h"

/*
 * execute_pci_calculation: main driver for parallel CI calculation.
 */
void execute_pci_calculation(int aelec, int belec, int orbs, int nastr,
			     int nbstr, int ndets, int ndocc, int nactv,
			     double *moints1, double *moints2, double nucrep,
			     double frzcore, int plvl)
{
	/* .. MPI Proccess scalars .. */
	int mpi_num_procs; /* MPI: Number of mpi processes */
	int mpi_proc_rank; /* MPI: Process rank */
	const int mpi_root = 0;

	/* .. MPI Data types .. */
	int blocklengths[4] = {1,1,1,1};
	MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Datatype mpi_occstr_type;
	MPI_Datatype mpi_det_type;
	MPI_Aint offsets[4] = {0,0,0,0};
	
	int error = 0; /* error flag */
	double totfrze = 0.0; /* nucrep + frzcore */
	int ninto = 0; /* internal orbitals (docc + active) */
	
	int maxiter = 0; /* maximum number of davidson algorithm iterations */
	int krymin = 0; /* minimum size of the krylov space */
	int krymax = 0; /* maximum size of the krylov space */
	int nroots = 0; /* number of roots to solve for */
	int prediag_routine = 0; /* prediagonalization routine */
	int refdim = 0; /* prediagonalization reference space */
	double restol = 0.0; /* residual norm converegence tolerance */
	
	struct det *detlist = NULL; /* determinant list */
	double *hdgls = NULL; /* diagonal elements of hamiltonian */
	double **civec = NULL; /* final CI eigenvectors */
	double *cival = NULL; /* final ci eigenvalues */
	struct rowmap *hmap = NULL; /* valid <i|H|j> combinations */
	int totelm; /* total matrix elements to evaluate */\
		
	int map_chunk = 0; /* Map chunk size */
	int map_lwrbnd = 0, map_uppbnd = 0; /* Map chunk lower/upper bounds */
	int nrows = 0; /* number of rows in hmap */
	clock_t curr_time, prev_time;

	/* Create MPI Type for BCast */
	blocklengths[2] = 2;
	types[0] = MPI_LONG_LONG_INT;
	types[1] = MPI_LONG_LONG_INT;
	offsets[0] = offsetof(struct occstr, byte1);
	offsets[1] = offsetof(struct occstr, byte2);
	offsets[2] = offsetof(struct occstr, virtx);
	offsets[3] = offsetof(struct occstr, nvrtx);
	MPI_Type_create_struct(4, blocklengths, offsets, types, &mpi_occstr_type);
	MPI_Type_commit(&mpi_occstr_type);
	blocklengths[2] = 1;
	types[0] = mpi_occstr_type;
	types[1] = mpi_occstr_type;
	offsets[0] = offsetof(struct det, astr);
	offsets[1] = offsetof(struct det, bstr);
	offsets[2] = offsetof(struct det, cas);
	MPI_Type_create_struct(3, blocklengths, offsets, types, &mpi_det_type);
	MPI_Type_commit(&mpi_det_type);
	/* */
	
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_proc_rank);

	/* Read davidson algorithm input. compute total frozen core enrgy
	 * and number of internal orbitals. broadcast the total number of
	 * internal orbitals to all processes. */
	if (mpi_proc_rank == mpi_root) {
		fprintf(stdout, " Start time: ");
		timestamp();
		totfrze = nucrep + frzcore;
		ninto = ndocc + nactv;
		/* Read &dgalinfo namelist. */
		readdaiinput(&maxiter, &krymin, &krymax, &nroots,
			     &prediag_routine, &refdim, &restol, &error);
		if (error != 0) {
			error_flag(error, "execute_pci_calculation");
			return;
		}
	}
	MPI_Bcast(&ninto, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);

	/* Read in det.list on master process. Broadcast the list to slaves. */
	detlist = (struct det *) malloc(ndets * sizeof(struct det));
	init_detlist(detlist, ndets);
	if (mpi_proc_rank == mpi_root) {
		error = genbinarydetlist(detlist, aelec, belec, orbs, ndocc,
					 nactv, ndets);
		if (error != 0) {
			error_flag(error, "execute_pci_calculation");
			return;
		}
	}
	MPI_Bcast(detlist, ndets, mpi_det_type, mpi_root, MPI_COMM_WORLD);

	/* Generate map of valid CI matrix elements */
	/* Divide <i|H|j> matrix up into sections for each process to handle. */
	if (plvl > 0 && mpi_proc_rank == mpi_root) printf("Generating cimap.\n");
	map_chunk = get_upptri_size(ndets);
	map_chunk = map_chunk / mpi_num_procs;
	map_lwrbnd = get_upptri_row((mpi_proc_rank * map_chunk), ndets);
	map_uppbnd = get_upptri_row(((mpi_proc_rank + 1) * map_chunk), ndets);
	if (plvl > 1 && mpi_proc_rank == mpi_root)
		printf("Chunksize = %15d\n", map_chunk);
	if (plvl > 1) printf("proc %2d: lower bound = %10d, upper bound = %10d\n",
			     mpi_proc_rank, map_lwrbnd, map_uppbnd);
	/* Generate partial hmap for selected rows. */
	nrows = map_uppbnd - map_lwrbnd;
	hmap = (struct rowmap *) malloc(sizeof(struct rowmap));
	error = pci_generate_cimap(detlist, ndets, nactv, map_lwrbnd, map_uppbnd,
				   hmap, nrows);
	/* Get total "on" elements for each process. Sort sections of hmap to
	 * balance the computational load. */
	error = get_total_nelm(hmap, nrows, &totelm);
	printf("%d: %d \n", mpi_proc_rank, totelm);

        /* Call davidson algorithm subroutine. */
        
	return;
}
		
