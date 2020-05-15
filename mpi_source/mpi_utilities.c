// File: mpi_utilities.c
/*
 * Utilities for MPI implementations of pjayci. The routines contained
 * herein modify global pjayci.x variables, and are responsible for
 * inter-node communication. Mostly interfaces for MPI Library routines.
 */

#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "iminmax.h"
#include "mpi_utilities.h"
#include <mpi.h>
#include <ga.h>
#include <macdecls.h>

/*
 * mpi_error_check_msg: check for error. Print message if error has
 * occured.
 * Input:
 *  error    = error flag
 *  fcn_name = calling function name
 *  message  = error message to print
 */
void mpi_error_check_msg (int error, char *fcn_name, char *message)
{
        int eflag;
        MPI_Allreduce(&error, &eflag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (eflag != 0) {
                MPI_Allreduce(&error, &eflag, 1, MPI_INT, MPI_MINLOC,
                              MPI_COMM_WORLD);
                if (mpi_proc_rank == mpi_root) {
                        error_message(eflag, message, fcn_name);
                }
                MPI_Finalize();
                exit(-1);
        }
}

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
void mpi_split_work_array_1d (int len, int *chunk, int *lo, int *hi)
{
        *chunk = 0;
        *chunk = len / mpi_num_procs;
        *lo = mpi_proc_rank * (*chunk);
        *hi = (mpi_proc_rank + 1) * (*chunk) - 1;
        if ((*chunk * mpi_num_procs) != len) {
                if (mpi_proc_rank == (mpi_num_procs - 1)) {
                        *hi = len - 1; // Last element of array
                }
        }
        return;
}

/*
 * print_gavectors2file_int: print a set of GA vectors to a file.
 * Input:
 *  hndl  = global arrays handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void print_gavectors2file_int(int hndl, int len, int dim, char *fname)
{
#define BUFSIZE 6
        int vdata[BUFSIZE];          /* GA Vector buffer */
        int ld = 1;                  /* Local vector buffer leading dimension */
        FILE *fptr = NULL;           /* File pointer */
        int hi[2] = {0, 0};
        int lo[2] = {0, 0};
        int ga_type    = 0;          /* Data type of GA */
        int ga_ndim    = 0;          /* Number of dimensions of GA */
        int ga_dims[2] = {0, 0};     /* Dimensions of GA */
        int buflen = BUFSIZE;
        int i, j, jj;
        int jmax;
        
        if (mpi_proc_rank == mpi_root) {
                
                NGA_Inquire(hndl, &ga_type, &ga_ndim, ga_dims);
                if (dim > ga_dims[0] || len > ga_dims[1]) {
                        error_message(mpi_proc_rank, "Dimension out of range.",
                                      "print_gavectors2file");
                        fprintf(stderr, "(%d, %d) != GA[%d, %d]\n", dim, len,
                                ga_dims[0], ga_dims[1]);
                        return;
                }

                fptr = fopen(fname, "w");
                fprintf(fptr, " %10d %10d\n", dim, len);

                for (i = 0; i < dim; i++) {
                        for (j = 0; j < len; j += buflen) {
                                jmax = int_min((j + buflen - 1), (len - 1));
                                lo[0] = i;
                                lo[1] = j;
                                hi[0] = i;
                                hi[1] = jmax;

                                NGA_Get(hndl, lo, hi, vdata, &ld);

                                for (jj = 0; jj < (jmax - j + 1); jj++) {
                                        fprintf(fptr, " %20d", vdata[jj]);
                                }
                        }
                        fprintf(fptr, "\n");
                }

                fflush(fptr);

                fclose(fptr);

        }

        GA_Sync();
        return;
}
/*
 * print_gavectors2file_int_ufmt: print a set of GA vectors to a file.
 * UNFORMATTED
 * Input:
 *  hndl  = global arrays handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void print_gavectors2file_int_ufmt(int hndl, int len, int dim, char *fname)
{
#define BUFSIZE 6
        int vdata[BUFSIZE];          /* GA Vector buffer */
        int ld = 1;                  /* Local vector buffer leading dimension */
        FILE *fptr = NULL;           /* File pointer */
        int hi[2] = {0, 0};
        int lo[2] = {0, 0};
        int ga_type    = 0;          /* Data type of GA */
        int ga_ndim    = 0;          /* Number of dimensions of GA */
        int ga_dims[2] = {0, 0};     /* Dimensions of GA */
        int buflen = BUFSIZE;
        int i, j, jj;
        int jmax;
        
        if (mpi_proc_rank == mpi_root) {
                
                NGA_Inquire(hndl, &ga_type, &ga_ndim, ga_dims);
                if (dim > ga_dims[0] || len > ga_dims[1]) {
                        error_message(mpi_proc_rank, "Dimension out of range.",
                                      "print_gavectors2file");
                        fprintf(stderr, "(%d, %d) != GA[%d, %d]\n", dim, len,
                                ga_dims[0], ga_dims[1]);
                        return;
                }

                fptr = fopen(fname, "w");
                if (fwrite(&dim, sizeof(int), 1, fptr) != 1) {
                        fprintf(stderr, "Error occured writing to file.\n");
                        fclose(fptr);
                        return;
                }
                if (fwrite(&len, sizeof(int), 1, fptr) != 1) {
                        fprintf(stderr, "Error occured writing to file.\n");
                        fclose(fptr);
                        return;
                }

                for (i = 0; i < dim; i++) {
                        for (j = 0; j < len; j += buflen) {
                                jmax = int_min((j + buflen - 1), (len - 1));
                                lo[0] = i;
                                lo[1] = j;
                                hi[0] = i;
                                hi[1] = jmax;

                                NGA_Get(hndl, lo, hi, vdata, &ld);

                                if (fwrite(vdata, sizeof(double), (jmax - j + 1),
                                           fptr) != (jmax - j + 1)) {
                                        fprintf(stderr, "Error ");
                                        fprintf(stderr, "writing vector!\n");
                                        fclose(fptr);
                                        return;
                                }
                        }
                }

                fflush(fptr);

                fclose(fptr);

        }

        GA_Sync();
        return;
}

/*
 * print_gavectors2file_dbl: print a set of GA vectors to a file.
 * Input:
 *  hndl  = global arrays handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void print_gavectors2file_dbl(int hndl, int len, int dim, char *fname)
{
#define BUFSIZE 6
        double vdata[BUFSIZE];       /* GA Vector buffer */
        int ld = 1;                  /* Local vector buffer leading dimension */
        FILE *fptr = NULL;           /* File pointer */
        int hi[2] = {0, 0};
        int lo[2] = {0, 0};
        int ga_type    = 0;          /* Data type of GA */
        int ga_ndim    = 0;          /* Number of dimensions of GA */
        int ga_dims[2] = {0, 0};     /* Dimensions of GA */
        int buflen = BUFSIZE;
        int i, j, jj;
        int jmax;
        
        if (mpi_proc_rank == mpi_root) {
                
                NGA_Inquire(hndl, &ga_type, &ga_ndim, ga_dims);
                if (dim > ga_dims[0] || len > ga_dims[1]) {
                        error_message(mpi_proc_rank, "Dimension out of range.",
                                      "print_gavectors2file");
                        fprintf(stderr, "(%d, %d) != GA[%d, %d]\n", dim, len,
                                ga_dims[0], ga_dims[1]);
                        return;
                }

                fptr = fopen(fname, "w");
                fprintf(fptr, " %10d %10d\n", dim, len);

                for (i = 0; i < dim; i++) {
                        for (j = 0; j < len; j += buflen) {
                                jmax = int_min((j + buflen - 1), (len - 1));
                                lo[0] = i;
                                lo[1] = j;
                                hi[0] = i;
                                hi[1] = jmax;

                                NGA_Get(hndl, lo, hi, vdata, &ld);

                                for (jj = 0; jj < (jmax - j + 1); jj++) {
                                        fprintf(fptr, " %20.15lf", vdata[jj]);
                                }
                        }
                        fprintf(fptr, "\n");
                }

                fflush(fptr);

                fclose(fptr);

        }

        GA_Sync();
        return;
}

/*
 * print_gavectors2file_dbl_ufmt: print a set of GA vectors to a file.
 * UNFORMATTED
 * Input:
 *  hndl  = global arrays handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void print_gavectors2file_dbl_ufmt(int hndl, int len, int dim, char *fname)
{
#define BUFSIZE 6
        double vdata[BUFSIZE];       /* GA Vector buffer */
        int ld = 1;                  /* Local vector buffer leading dimension */
        FILE *fptr = NULL;           /* File pointer */
        int hi[2] = {0, 0};
        int lo[2] = {0, 0};
        int ga_type    = 0;          /* Data type of GA */
        int ga_ndim    = 0;          /* Number of dimensions of GA */
        int ga_dims[2] = {0, 0};     /* Dimensions of GA */
        int buflen = BUFSIZE;
        int i, j, jj;
        int jmax;
        int eflag = 0;
        
        if (mpi_proc_rank == mpi_root) {
                
                NGA_Inquire(hndl, &ga_type, &ga_ndim, ga_dims);
                if (dim > ga_dims[0] || len > ga_dims[1]) {
                        error_message(mpi_proc_rank, "Dimension out of range.",
                                      "print_gavectors2file");
                        fprintf(stderr, "(%d, %d) != GA[%d, %d]\n", dim, len,
                                ga_dims[0], ga_dims[1]);
                        return;
                }

                fptr = fopen(fname, "w");
                if (fwrite(&dim, sizeof(int), 1, fptr) != 1) {
                        fprintf(stderr, "Error occured writing to file.\n");
                        fclose(fptr);
                        return;
                }
                if (fwrite(&len, sizeof(int), 1, fptr) != 1) {
                        fprintf(stderr, "Error occured writing to file.\n");
                        fclose(fptr);
                        return;
                }
                
                for (i = 0; i < dim; i++) {
                        for (j = 0; j < len; j += buflen) {
                                jmax = int_min((j + buflen - 1), (len - 1));
                                lo[0] = i;
                                lo[1] = j;
                                hi[0] = i;
                                hi[1] = jmax;

                                NGA_Get(hndl, lo, hi, vdata, &ld);

                                if (fwrite(vdata, sizeof(double), (jmax - j + 1),
                                           fptr) != (jmax - j + 1)) {
                                        fprintf(stderr, "Error ");
                                        fprintf(stderr, "writing vector!\n");
                                        fclose(fptr);
                                        return;
                                }
                        }
                }

                fflush(fptr);

                fclose(fptr);

        }

        GA_Sync();
        return;
}

/*
 * read_gavectorsfile_dbl_ufmt: read vectors into GA arrays
 * UNFORMATTED
 * Input:
 *  hndl  = GA handle
 *  len   = length of vectors
 *  dim   = number of vectors
 *  fname = file name
 */
void read_gavectorsfile_dbl_ufmt(int hndl, int len, int dim, char *fname)
{
#define BUFSIZE 6
        double vdata[BUFSIZE];     /* GA Vector buffer */
        int ld = 1;                /* Local vector buffer leading dimension */
        FILE *fptr = NULL;         /* File pointer */
        int hi[2] = {0, 0};
        int lo[2] = {0, 0};
        int ga_dims[2] = {0, 0};   /* Dimensions of GA */
        int ga_ndim    = {0};      /* Number of dimensions of GA */
        int ga_type    = 0;        /* Data type of GA */
        
        int buflen = BUFSIZE;
        int i, j, jj;
        int jmax;

        if (mpi_proc_rank == mpi_root) {

                NGA_Inquire(hndl, &ga_type, &ga_ndim, ga_dims);
                if (dim > ga_dims[0] || len > ga_dims[1]) {
                        error_message(mpi_proc_rank, "Dimension of of range.",
                                      "read_gavectorsfile_dbl");
                        fprintf(stderr, "(%d, %d) != GA[%d, %d]\n", dim, len,
                                ga_dims[0], ga_dims[1]);
                        return;
                }

                fptr = fopen(fname, "r");
                

                
        }
        return;
}
/*
 * set_mpi_process_number_and_rank: set global variables $mpi_num_procs
 * and $mpi_proc_rank.
 */
void set_mpi_process_number_and_rank ()
{
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_proc_rank);
}

/*
 * set_ga_process_number_and_rank: set global variables $mpi_num_procs
 * and $mpi_proc_rank with GA wrappers.
 */
void set_ga_process_number_and_rank ()
{
        mpi_proc_rank = GA_Nodeid();
        mpi_num_procs = GA_Nnodes();
}
