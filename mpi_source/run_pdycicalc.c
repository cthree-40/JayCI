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
#include "dysoncomp.h"
#include "run_pdycicalc.h"

#include <ga.h>
#include <macdecls.h>

#include <mpi.h>

/*
 * run_pdycicalc: Execute CI dyson orbital calculation.
 */
int run_pdycicalc ()
{
#define MAXSTATES_CI 5
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
        int ninto0;    /* Number of internal (DOCC+CAS) orbitals */
	int nvirt0;    /* Number of external orbitals */
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
        int *pqsp0data = NULL;
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
        int ninto1;    /* Number of internal (DOCC+CAS) orbitals */
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
        int *pqsp1data = NULL;
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
        int dysnst0[MAXSTATES_CI];    /* Anion states of dyson orbital */
        int dysnst1[MAXSTATES_CI];    /* Neutral states of dyson orbital */
        int maxstates = 0;            /* Max Anion/Neutral states */
        int ndyst0 = 0;               /* Number of anion staets in dyson orb. */
        int ndyst1 = 0;               /* Number of neutral states in dyson orb.*/
        double **dyorb_lc = NULL;     /* LOCAL dyson orbitals */
        double *dyorb_lc_data = NULL; /* LOCAL dyson orbital memory block */
        double **dyorb_gl = NULL;     /* GLOBAL dyson orbitals */
        double *dyorb_gl_data = NULL; /* GLOBAL dyson orbitals memory block */

	int **strcont   = NULL;  /* String contributions to dyson orbital */
	int *strcont1d  = NULL;  /* 1-d data */
	
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
        ninto0 = ndocc0 + nactv0;
	nvirt0 = norbs0 - nfrzc0 - nfrzv0 - ninto0;
        
        MPI_Bcast(&nelecs1,  1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&norbs1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzc1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&ndocc1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nactv1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzv1,   1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&xlvl1,    1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nstates1, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        ninto1 = ndocc1 + nactv1;
        
        error = check_wavefunction_input(nelecs0, norbs0, nfrzc0, ndocc0,
                                         nactv0, nfrzv0, nelecs1, norbs1,
                                         nfrzc1, ndocc1, nactv1, nfrzv1);
        mpi_error_check_msg(error, "run_pdycicalc",
                            "Error in wavefunction input");

        /* Read dysonorbital input */
        maxstates = MAXSTATES_CI;
        if (mpi_proc_rank == mpi_root) {
                readdysoninput(dysnst0, dysnst1, maxstates, &ndyst0, &ndyst1,
                               &error);

                printf("\nComputing dyson orbitals between:\n");
                printf(" Anion states (%d):   ", ndyst0);
                for (i = 0; i < ndyst0; i++) {
                        printf(" %d", dysnst0[i]);
                        /* Decrement value for C array indexing. */
                        dysnst0[i] = dysnst0[i] - 1;
                }
                printf("\n");
                printf(" Neutral states (%d): ", ndyst1);
                for (i = 0; i < ndyst1; i++) {
                        /* Decrement value for C array indexing. */
                        printf(" %d", dysnst1[i]);
                        dysnst1[i] = dysnst1[i] - 1;
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
        pqsp0data = allocate_mem_int_cont(&pq_space_pairs0, 2, num_pq0);
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
        pqsp1data = allocate_mem_int_cont(&pq_space_pairs1, 2, num_pq1);
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
        if (mpi_proc_rank == mpi_root) fflush(stdout);
	
        /* Allocate LOCAL arrays: dyorb_lc, dyorb_gl */
        ndyorbs = ndyst0 * ndyst1;
        dyorb_lc_data = allocate_mem_double_cont(&dyorb_lc, norbs0, ndyorbs);
        dyorb_gl_data = allocate_mem_double_cont(&dyorb_gl, norbs0, ndyorbs);
        GA_Sync();
        
        /* Read civectors. */
        read_gavectorsfile_dbl_ufmt(v0_hndl, dtrm0_len, nstates0, "anion.ci");
        read_gavectorsfile_dbl_ufmt(v1_hndl, dtrm1_len, nstates1, "neutral.ci");

        /* Generate W arrays for each wavefunction. */
        generate_wlist(w0_hndl, dtrm0_len, pq_space_pairs0, num_pq0, peospace0,
                       pegrps0, qeospace0, qegrps0);
        generate_wlist(w1_hndl, dtrm1_len, pq_space_pairs1, num_pq1, peospace1,
                       pegrps1, qeospace1, qegrps1);

	/* Generate list of N-electron strings that pair with N+1-electron
	 * strings */
	strcont1d = allocate_mem_int_cont(&strcont, cibelec0, qstr0_len);
	generate_strcontlist(qstrings0, qstr0_len, qeospace0, qegrps0, ndocc0,
			     nactv0, nvirt0, strcont, cibelec1, qeospace1,
			     qegrps1);
	   
        /* Compute dyson orbitals */
        /* If number of alpha/beta electrons are equal in both (N+1) and (N)
         * electron wavefunctions, then alpha/beta strings are IDENTICAL, and
         * they do not need to be compared. (Contribution is 1.)
         * An important point must be made: the abecalc() routine will always
         * make the number of alpha electrons greater than the number of beta
         * electrons.
         * THEREFORE:
         *  1) Singlet/Triplet (alpha = beta) -> Doublet (alpha > beta)
         *         :: compare beta  strings only
         *  2) Doublet (alpha > beta) -> Singlet/Triplet (alpha = beta)
         *         :: compare alpha strings only
         */
        if (nelecs0 % 2 == 0) {
                if (mpi_proc_rank == mpi_root) {
                        printf("N+1 wavefunction is singlet/triplet.\n");
                        printf(" Comparing beta electron strings...\n");
			fflush(stdout);
		}
                /* S/T (N+1) wavefunction, compare beta  strings. */
		compute_dyson_orbital(v0_hndl, v1_hndl, w0_hndl, w1_hndl,
				      pstrings0, peospace0, pegrps0,
				      qstrings0, qeospace0, qegrps0,
				      pstrings1, peospace1, pegrps1,
				      qstrings1, qeospace1, qegrps1,
				      ciorbs1, ndocc1, nactv1, ndyst0,
				      dysnst0, ndyst1, dysnst1, dtrm0_len,
				      dtrm1_len, 1, 0);
		
        } else {
                if (mpi_proc_rank == mpi_root) {
                        printf("N+1 wavefunction is doublet.\n");
                        printf(" Comparing alpha electron strings...\n");
			fflush(stdout);
		}
                /* D   (N+1) wavefunction, compare alpha strings. */
		compute_dyson_orbital(v0_hndl, v1_hndl, w0_hndl, w1_hndl,
				      pstrings0, peospace0, pegrps0,
				      qstrings0, qeospace0, qegrps0,
				      pstrings1, peospace1, pegrps1,
				      qstrings1, qeospace1, qegrps1,
				      ciorbs1, ndocc1, nactv1, ndyst0,
				      dysnst0, ndyst1, dysnst1, dtrm0_len,
				      dtrm1_len, 0, 1);
        }

        /* Accumulate dyson orbitals from each process. */
        GA_Sync();
        MPI_Allreduce(dyorb_lc_data, dyorb_gl_data, (ndyorbs * norbs0),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (mpi_proc_rank == mpi_root) {
                printf("Finished computing %d dyson orbitals.\n", ndyorbs);
        }
        /* Print dyson orbitals to file */
        if (mpi_proc_rank == mpi_root) {
                print_dysonorbitals_to_file("dysonorb.dat", ndyorbs, norbs0,
                                            dyorb_gl, dysnst0, ndyst0, dysnst1,
                                            ndyst1);
        }

        /* Deallocate pstrings and qstrings */
        free(pstrings0);
        free(pstrings1);
        free(qstrings0);
        free(qstrings1);
        free(peospace0);
        free(peospace1);
        free(qeospace0);
        free(qeospace1);
        /* Deallocate dyorb_*_data */
        deallocate_mem_cont(&dyorb_lc, dyorb_lc_data);
        deallocate_mem_cont(&dyorb_gl, dyorb_gl_data);
        /* Deallocate pqspace pairs */
        deallocate_mem_cont_int(&pq_space_pairs0, pqsp0data);
        deallocate_mem_cont_int(&pq_space_pairs1, pqsp1data);
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
 * print_dysonorbitals_to_file: print the compute dyson orbitals to file.
 */
void print_dysonorbitals_to_file(char *filename, int ndyorbs, int orbs,
                                 double **dyorbs, int *dysnst0, int ndyst0,
                                 int *dysnst1, int ndyst1)
{
        FILE* fptr = NULL; /* File pointer */
        int i, j;

        fptr = fopen(filename, "w");
        if (fptr == NULL) {
                printf("Could not open file: %s\n", filename);
                return;
        }
        /* Print dyson orbital information first */
        /*
         * 1st line: # of... dyson orbitals, anion states, neutral states
         * 2nd line: # of molecular orbitals
         * 3rd line: Anion states
         * 4th line: Neutral states
         */
        fprintf(fptr, " %d %d %d\n", ndyorbs, ndyst0, ndyst1);
        fprintf(fptr, " %d\n", orbs);
        for (i = 0; i < ndyst0; i++) {
                fprintf(fptr, " %d", (dysnst0[i] + 1));
        }
        fprintf(fptr, "\n");
        for (i = 0; i < ndyst1; i++) {
                fprintf(fptr, " %d", (dysnst1[i] + 1));
        }
        fprintf(fptr, "\n");
        for (i = 0; i < ndyorbs; i++) {
                for (j = 0; j < orbs; j++) {
                        fprintf(fptr, "%15.8lf\n", dyorbs[i][j]);
                }
                fprintf(fptr,"\n");
        }

        fclose(fptr);
        
        return;
}
