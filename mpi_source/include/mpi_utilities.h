// File: mpi_utilities.h
/*
 * Utilities for MPI implementations of pjayci. The routines contained
 * herein modify global pjayci.x variables, and are responsible for
 * inter-node communication. Mostly interfaces for MPI Library routines.
 */

#ifndef mpi_utilities_h
#define mpi_utilities_h

int mpi_num_procs;          /* MPI: Number of mpi processes */
int mpi_proc_rank;          /* MPI: Processor rank */

extern const int mpi_root;  /* MPI: Root process is always 0. */

/*
 * mpi_end_program: end program if error occurs.
 */
void mpi_end_program ();

/*
 * mpi_error_check_msg: check for error. Print message if error has
 * occured.
 * Input:
 *  error    = error flag
 *  fcn_name = calling function name
 *  message  = error message to print
 */
void mpi_error_check_msg (int error, char *fcn_name, char *message);

/*
 * mpi_split_work_array_1d: get first and last elements for partitioning
 * a 1d array amongst work processes.
 * Input:
 *  len  = length of vector (total)
 * Output:
 *  chunk = chunk size
 *  lo    = first element
 *  hi    = last  element
 */
void mpi_split_work_array_1d (int len, int *chunk, int *lo, int *hi);

/*
 * print_gavectors2file_int: print a set of GA vectors to a file.
 * Input:
 *  hndl  = global arrays handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void print_gavectors2file_int(int hndl, int len, int dim, char *fname);

/*
 * print_gavectors2file_int_ufmt: print a set of GA vectors to a file.
 * UNFORMATTED
 * Input:
 *  hndl  = global arrays handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void print_gavectors2file_int_ufmt(int hndl, int len, int dim, char *fname);

/*
 * print_gavectors2file_dbl: print a set of GA vectors to a file.
 */
void print_gavectors2file_dbl(int hndl, int len, int dim, char *fname);

/*
 * print_gavectors2file_dbl_ufmt: print a set of GA vectors to a file.
 * UNFORMATTED
 * Input:
 *  hndl  = global arrays handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void print_gavectors2file_dbl_ufmt(int hndl, int len, int dim, char *fname);

/*
 * read_gavectorsfile_dbl_ufmt: read vectors into GA arrays
 * UNFORMATTED
 * Input:
 *  hndl  = GA handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void read_gavectorsfile_dbl_ufmt(int hndl, int len, int dim, char *fname);

/*
 * set_mpi_process_number_and_rank: set global variables $mpi_num_procs
 * and $mpi_proc_rank.
 */
void set_mpi_process_number_and_rank ();

/*
 * set_ga_process_number_and_rank: set global variables $mpi_num_procs
 * and $mpi_proc_rank with GA wrappers.
 */
void set_ga_process_number_and_rank ();

#endif
