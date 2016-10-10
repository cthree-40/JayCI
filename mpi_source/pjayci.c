// File: pjayci.c
/*
 * Parallel implementation of jayci.x determinant-based CI algorithm.
 *
 * Written by Christopher L Malbon
 * Yarkony Group, Dept of Chemistry
 * Johns Hopkins University
 * 2016-10-10
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "ioutil.h"
#include "progheader.h"
#include "moindex.h"
#include "execute_pci_calculation.h"
/*
 * pjayci: parallel implementation of determinant-based CI algorithm.
 */
int main(int argc, char *argv[])
{
	int error = 0;

	/* .. MPI Proccess scalars .. */
	int mpi_num_procs; /* MPI: Number of mpi processes */
	int mpi_proc_rank; /* MPI: Process rank */
	const int mpi_root = 0;

	int general_input[8]; /* &general namelist input buffer */
	int injayci_input[6]; /* input.jayci input buffer */

	int m1len = 0, m2len = 0; /* 1-e, 2-e integral lengths */
	double *moints1 = NULL;
	double *moints2 = NULL;
	char moflname[FLNMSIZE] = {""}; /* SIFS integral filename */
	int itype = 1; /* integral type. Always 1. */
	
	double nucrep_energy = 0.0;  /* nuclear repulsion energy */
	double frzcore_energy = 0.0; /* frozen core orbital energy */
	
	int i, j;
	
	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_proc_rank);

	/* Print header. Read scalar inputs. Broadcast values. */
	if (mpi_proc_rank == 0) {
		print_pjayciheader();
		readgeninput(&general_input[0], &general_input[1],
			     &general_input[2], &general_input[3],
			     &general_input[4], &general_input[5],
			     &general_input[6], &general_input[7], &error);
		error = readinputjayci(&injayci_input[0], &injayci_input[1],
				       &injayci_input[2], &injayci_input[3],
				       &injayci_input[4], &injayci_input[5]);
	}
	MPI_Bcast(general_input, 8, MPI_INT, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast(injayci_input, 6, MPI_INT, mpi_root, MPI_COMM_WORLD);

	/* Allocate molecular integral arrays. Read integrals on
	 * master process and broadcast them. */
	m1len = index1e(general_input[1], general_input[1]);
	m2len = index2e(general_input[1], general_input[1],
			general_input[1], general_input[1]);
	moints1 = (double *) malloc(m1len * sizeof(double));
	moints2 = (double *) malloc(m2len * sizeof(double));
	init_dbl_array_0(moints1, m1len);
	init_dbl_array_0(moints2, m1len);
	if (mpi_proc_rank == 0) {
		strncpy(moflname, "moints", FLNMSIZE);
		readmointegrals(moints1, moints2, itype, injayci_input[2],
				moflname, m1len, m2len, &nucrep_energy,
				&frzcore_energy);
	}
	MPI_Bcast(moints1, m1len, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast(moints2, m2len, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);

	/* Execute parallel CI calculation */
	execute_pci_calculation(injayci_input[0], injayci_input[1],
				injayci_input[2], injayci_input[3],
				injayci_input[4], injayci_input[5],
				general_input[3], general_input[4],
				moints1, moints2, nucrep_energy, frzcore_energy,
				general_input[7]);
	MPI_Finalize();
	return 0;
}
