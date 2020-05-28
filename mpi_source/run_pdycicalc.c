// File: run_pdycical.c
#include <stdio.h>
#include <stdlib.h>
#include "mpi_utilities.h"
#include "errorlib.h"
#include "allocate_mem.h"
#include "ioutil.h"
#include "combinatorial.h"
#include "abecalc.h"
#include "binarystr.h"
#include "iminmax.h"
#include "citruncate.h"
#include "run_pdycicalc.h"

#include <ga.h>
#include <macdecls.h>

#include <mpi.h>

/*
 * run_pdycicalc: Execute CI dyson orbital calculation.
 */
int run_pdycicalc ()
{
        int error = 0; /* Error flag */

        int w0_hndl = 0;               /* GA wavefunction 0 */
        int w0_dims[2] = {0, 0};       /* GA wavefunction 0 dimension */
        int w0_chunk[2]= {0, 0};       /* GA wavefunction 0 chunk sizes */
        int w1_hndl = 0;               /* GA wavefunction 1 */
        int w1_dims[2] = {0, 0};       /* GA wavefunction 1 dimension */
        int w1_chunk[2]= {0, 0};       /* GA wavefunction 1 chunk sizes */

        /*
         * Wavefunction information.
         * 0 refers to the N+1 electron wavefunction, and 1 refers to the
         * N electron wavefunction.
         */
        int nstates0;  /* Number of states */
        int nelecs0;   /* Number of electrons */
        int naelec0;   /* Number of alpha electrons */
        int nbelec0;   /* Number of beta  electrons */
        int nfrzc0;    /* Number of frozen core orbitals */
        int ndocc0;    /* Number of DOCC orbitals */
        int nactv0;    /* Number of CAS  orbitals */
        int nfrzv0;    /* Number of frozen virtual orbitals */
        int xlvl0;     /* CI excitation level */
        int norbs0;    /* Number of orbitals */
        int ciorbs0;   /* NUmber of ci orbitals */
        int ciaelec0;  /* Number of CI alpha electrons */
        int cibelec0;  /* Number of CI beta  electrons */
        int pstr0_len; /* Number of alpha strings */
        int qstr0_len; /* Number of beta  strings */

        struct occstr  *pstrings0;
        struct occstr  *qstrings0;
        struct eospace *peospace0;
        struct eospace *qeospace0;
        int **pq_space_pairs0 = NULL;
        int pegrps0 = 0;
        int qegrps0 = 0;
        int num_pq0 = 0;
        int dtrm0_len = 0;
        
        int nstates1;  /* Number of states */
        int nelecs1;   /* Number of electrons */
        int naelec1;   /* Number of alpha electrons */
        int nbelec1;   /* Number of beta  electrons */
        int nfrzc1;    /* Number of frozen core orbitals */
        int ndocc1;    /* Number of DOCC orbitals */
        int nactv1;    /* Number of CAS  orbitals */
        int nfrzv1;    /* Number of frozen virtual orbitals */
        int xlvl1;     /* CI excitation level */
        int norbs1;    /* Number of orbitals */
        int ciorbs1;   /* NUmber of ci orbitals */
        int ciaelec1;  /* Number of CI alpha electrons */
        int cibelec1;  /* Number of CI beta  electrons */
        int pstr1_len; /* Number of alpha strings */
        int qstr1_len; /* Number of beta  strings */

        struct occstr  *pstrings1;
        struct occstr  *qstrings1;
        struct eospace *peospace1;
        struct eospace *qeospace1;
        int **pq_space_pairs1 = NULL;
        int pegrps1 = 0;
        int qegrps1 = 0;
        int num_pq1 = 0;
        int dtrm1_len = 0;

        int v0_hndl = 0;          /* GA 0 (Anion) CI vectors */
        int v0_dims[2] = {0, 0};  /* GA 0 (Anion) CI vector dimensions */
        int v0_chunk[2]= {0, 0};  /* GA 0 (Anion) CI vector chunk sizes */
        int v1_hndl = 0;          /* GA 1 (Neutral) CI vectors */
        int v1_dims[2] = {0, 0};  /* GA 1 (Neutral) CI vector dimensions */
        int v1_chunk[2]= {0, 0};  /* GA 1 (Neutral) CI vector chunk sizes */
        
        int ndyorbs = 0;              /* Number of dyson orbitals to compute */
        int *dysnst0 = NULL;          /* Anion states of dyson orbital */
        int *dysnst1 = NULL;          /* Neutral states of dyson orbital */
        int maxstates = 5;            /* Max Anion/Neutral states */
        int ndyst0 = 0;               /* Number of anion staets in dyson orb. */
        int ndyst1 = 0;               /* Number of neutral states in dyson orb.*/
        double **dyorb = NULL;        /* LOCAL dyson orbitals */
        double *dyorb_data = NULL;    /* LOCAL dyson orbital memory block */

        double memusage = 0.0;  /* Estimated memory usage. */

        int i = 0;
        
        /* Read wavefunction input. */
        if (mpi_proc_rank == mpi_root) {
                readwf0input(&nelecs0, &norbs0, &nfrzc0, &ndocc0, &nactv0,
                             &xlvl0, &nfrzv0, &nstates0, &error);
        }
        mpi_error_check_msg(error, "run_pdycicalc", "Error reading wf0 input!");
        if (mpi_proc_rank == mpi_root) {
                readwf1input(&nelecs1, &norbs1, &nfrzc1, &ndocc1, &nactv1,
                             &xlvl1, &nfrzv1, &nstates1, &error);
        }
        mpi_error_check_msg(error, "run_pdycicalc", "Error reading wf0 input!");
        if (mpi_proc_rank == mpi_root) {
                print_wavefunction_info("Wavefunction 0", nelecs0, norbs0,
                                        nfrzc0, ndocc0, nactv0, nfrzv0,
                                        xlvl0, nstates0);
                print_wavefunction_info("Wavefunction 1", nelecs1, norbs1,
                                        nfrzc1, ndocc1, nactv1, nfrzv1,
                                        xlvl1, nstates1);
        }


        MPI_Bcast(&nelecs0,  1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&norbs0,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzc0,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&ndocc0,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nactv0,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzv0,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&xlvl0,    1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nstates0, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);

        MPI_Bcast(&nelecs1,  1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&norbs1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzc1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&ndocc1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nactv1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzv1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&xlvl1,    1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nstates1, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);

        error = check_wavefunction_input(nelecs0, norbs0, nfrzc0, ndocc0,
                                         nactv0, nfrzv0, nelecs1, norbs1,
                                         nfrzc1, ndocc1, nactv1, nfrzv1);
        mpi_error_check_msg(error, "run_pdycicalc",
                            "Error in wavefunction input");

        /* Read dysonorbital input */
        dysnst0 = malloc(sizeof(int) * maxstates);
        dysnst1 = malloc(sizeof(int) * maxstates);
        if (mpi_proc_rank == mpi_root) {
                readdysoninput(dysnst0, dysnst1, maxstates, &ndyst0, &ndyst1,
                               &error);
                printf("\nComputing dyson orbitals between:\n");
                printf("Anion states:   ");
                for (i = 0; i < ndyst0; i++) {
                        printf(" %d", dysnst0[i]);
                }
                printf("\n");
                printf("Neutral states: ");
                for (i = 0; i < ndyst1; i++) {
                        printf(" %d", dysnst1[i]);
                }
                printf("\n\n");
        }
        mpi_error_check_msg(error, "run_dycicalc", "Error reading dyson input.");
        MPI_Bcast(&dysnst0,  maxstates, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&ndyst0,           1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&dysnst1,  maxstates, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&ndyst1,           1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        
        /* Set up wavefunctions */
        /* W0 (Anion) */
        abecalc(nelecs0, &naelec0, &nbelec0);
        compute_ci_elecs_and_orbitals(naelec0, nbelec0, norbs0, nfrzc0, nfrzv0,
                                      &ciaelec0, &cibelec0, &ciorbs0);
        pstr0_len = compute_stringnum(ciorbs0, ciaelec0, ndocc0, nactv0, xlvl0);
        qstr0_len = compute_stringnum(ciorbs0, cibelec0, ndocc0, nactv0, xlvl0);
        pstrings0 = allocate_occstr_arrays(pstr0_len);
        qstrings0 = allocate_occstr_arrays(qstr0_len);
        peospace0 = allocate_eospace_array(ciaelec0, ciorbs0, ndocc0, nactv0,
                                           xlvl0, &pegrps0);
        qeospace0 = allocate_eospace_array(cibelec0, ciorbs0, ndocc0, nactv0,
                                           xlvl0, &qegrps0);
        num_pq0 = pegrps0 * qegrps0;
        error = allocate_mem_int(&pq_space_pairs0, 2, num_pq0);
        error = citrunc(naelec0, nbelec0, norbs0, nfrzc0, ndocc0, nactv0, nfrzv0,
                        xlvl0, pstrings0, pstr0_len, qstrings0, qstr0_len,
                        peospace0, pegrps0, qeospace0, qegrps0, &dtrm0_len,
                        pq_space_pairs0, &num_pq0);
        mpi_error_check_msg(error, "run_pdycicalc",
                            "Error during wavefunction generation.");
        GA_Sync();
        if (mpi_proc_rank == mpi_root) {
                printf("Wavefunction 0:\n");
                printf(" Determinants   = %15d\n", dtrm0_len);
                printf("  Alpha strings = %15d\n", pstr0_len);
                printf("  Beta  strings = %15d\n", qstr0_len);
		memusage = ((pstr0_len + qstr0_len) * (8 + 4 + 4 + 4) +
			    (pegrps0 * qegrps0) * 2 * 4 +
			    (pegrps0 + qegrps0) * (4 + 4 + 4 + 4 + 4))
			/ 1048576;
	}
        
        /* W1 (Neutral) */
        abecalc(nelecs1, &naelec1, &nbelec1);
        compute_ci_elecs_and_orbitals(naelec1, nbelec1, norbs1, nfrzc1, nfrzv1,
                                      &ciaelec1, &cibelec1, &ciorbs1);
        pstr1_len = compute_stringnum(ciorbs1, ciaelec1, ndocc1, nactv1, xlvl1);
        qstr1_len = compute_stringnum(ciorbs1, cibelec1, ndocc1, nactv1, xlvl1);
        pstrings1 = allocate_occstr_arrays(pstr1_len);
        qstrings1 = allocate_occstr_arrays(qstr1_len);
        peospace1 = allocate_eospace_array(ciaelec1, ciorbs1, ndocc1, nactv1,
                                           xlvl1, &pegrps1);
        qeospace1 = allocate_eospace_array(cibelec1, ciorbs1, ndocc1, nactv1,
                                           xlvl1, &qegrps1);
        num_pq1 = pegrps1 * qegrps1;
        error = allocate_mem_int(&pq_space_pairs1, 2, num_pq1);
        error = citrunc(naelec1, nbelec1, norbs1, nfrzc1, ndocc1, nactv1, nfrzv1,
                        xlvl1, pstrings1, pstr1_len, qstrings1, qstr1_len,
                        peospace1, pegrps1, qeospace1, qegrps1, &dtrm1_len,
                        pq_space_pairs1, &num_pq1);
        mpi_error_check_msg(error, "run_pdycicalc",
                            "Error during wavefunction generation.");
        GA_Sync();
        if (mpi_proc_rank == mpi_root) {
                printf("Wavefunction 1:\n");
                printf(" Determinants   = %15d\n", dtrm1_len);
                printf("  Alpha strings = %15d\n", pstr1_len);
                printf("  Beta  strings = %15d\n", qstr1_len);
		memusage += ((pstr1_len + qstr1_len) * (8 + 4 + 4 + 4) +
                             (pegrps1 * qegrps1) * 2 * 4 +
                             (pegrps1 + qegrps1) * (4 + 4 + 4 + 4 + 4))
			/ 1048576;
	}

        if (mpi_proc_rank == mpi_root) {
                if (memusage < 1.0) {
                        printf("\nEstimated local memory usage: %10.2f KB\n",
                               memusage * 1024);
                } else {
                        printf("\nEstimated local memory usage: %10.2f MB\n",
                               memusage);
                }
        }
        
        /* Allocate GLOBAL arrays: W0, W1, V0, V1 */
        /* WX arrays are 3 x ndet arrays with (p, q, cas) where p is the
         * alpha string index, q is the beta string index, and cas is the
         * CAS flag. */
        if (mpi_proc_rank == mpi_root) printf("\nCreating global arrays...\n");
        w0_dims[0] = dtrm0_len;
        w0_dims[1] = 3;
        w0_chunk[0] = -1; // Distribute evenly
        w0_chunk[1] = 3;
        w0_hndl = NGA_Create(C_INT, 2, w0_dims, "W0: (p, q, cas)", w0_chunk);
        if (!w0_hndl) GA_Error("Create failed: W0: (p, q, cas)", 2);
        w1_dims[0] = dtrm1_len;
        w1_dims[1] = 3;
        w1_chunk[0] = -1;
        w1_chunk[1] = 3;
        w1_hndl = NGA_Create(C_INT, 2, w1_dims, "W1: (p, q, cas)", w1_chunk);
        if (!w1_hndl) GA_Error("Create failed: W1: (p, q, cas)", 2);
        /* Vx arrays are nstate x ndet arrays. They are the CI vectors */
        v0_dims[0] = nstates0;
        v0_dims[1] = dtrm0_len;
        v0_chunk[0] = nstates0;
        v0_chunk[1] = -1; // Distribute evenly
        v0_hndl = NGA_Create(C_DBL, 2, v0_dims, "V0: CI Vector", v0_chunk);
        if (!v0_hndl) GA_Error("Create failed: V0: CI Vector", 2);
        v1_dims[0] = nstates1;
        v1_dims[1] = dtrm1_len;
        v1_chunk[0] = nstates1;
        v1_chunk[1] = -1; // Distribute evenly
        v1_hndl = NGA_Create(C_DBL, 2, v1_dims, "V1: CI Vector", v1_chunk);
        if (!v1_hndl) GA_Error("Create failed: V1: CI Vector", 2);
        if (mpi_proc_rank == mpi_root) printf("Global arrays created.\n");
        
        /* Allocate LOCAL arrays: dyorb */
        ndyorbs = ndyst0 * ndyst1;
        dyorb_data = allocate_mem_double_cont(&dyorb, norbs0, ndyorbs);
        GA_Sync();
        
        /* Read civectors. */
        read_gavectorsfile_dbl_ufmt(v0_hndl, dtrm0_len, nstates0, "anion.ci",
                                    &error);
        mpi_error_check_msg(error, "run_dycicalc", "Error reading ci vectors.");
        read_gavectorsfile_dbl_ufmt(v1_hndl, dtrm1_len, nstates1, "neutral.ci",
                                    &error);
        mpi_error_check_msg(error, "run_dycicalc", "Error reading ci vectors.");
        
        /* Compute dyson orbital */
        
        return error;
}

/*
 * check_wavefunction_input: check user input of global wavefunction variables.
 * This is executed on all processes.
 */
int check_wavefunction_input(int nelecs0, int norbs0, int nfrzc0,  int ndocc0,
                             int nactv0,  int nfrzv0, int nelecs1, int norbs1,
                             int nfrzc1,  int ndocc1, int nactv1,  int nfrzv1)
{
        int error = 0;

        if (nelecs1 >= nelecs0) {
                printf("Error! nelecs1 >= nelecs0: %d >= %d\n",
                       nelecs1, nelecs0);
                error_message(mpi_proc_rank, "nelecs1 >= nelecs0",
                              "check_wavefunction_input");
                error = -1;
                return error;
        }
        if ((nelecs0 - nelecs1) > 1) {
                printf("Error! nelecs0 - nelecs1 > 1: %d - %d > 1\n",
                       nelecs0, nelecs1);
                error_message(mpi_proc_rank, "nelecs0 - nelecs1 > 1",
                              "check_wavefunction_input");
                error = -2;
                return error;
        }
        if (norbs1 != norbs0 || nfrzc1 != nfrzc0 || ndocc1 != ndocc0 ||
            nactv1 != nactv0 || nfrzv1 != nfrzv0) {
                printf("Error: Wavefunction active space differs!\n");
                error_message(mpi_proc_rank, "Wavefunction input incorrect.",
                              "check_wavefunction_input");
                error = -3;
                return error;
        }
        return error;
}
