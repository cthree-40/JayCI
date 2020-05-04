// File: pjayci_global.h
/*
 * Global variables for pjayci execution.
 */
#ifndef pjayci_global_h
#define pjayci_global_h

int mpi_num_procs;          /* MPI: Number of mpi processes */
int mpi_proc_rank;          /* MPI: Processor rank */

extern const int mpi_root;  /* MPI: Root process is always 0. */

double nuc_rep_e;           /* Nuclear repulsion energy */
double total_core_e;        /* Total core energy */

int ga_buffer_len;          /* Length of GA buffer to read during Hv=c */

#endif
