// File: execute_pjayci.c
/*
 * Execute parallel, determinant-based CI algorithm.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pjayci_global.h"
#include "errorlib.h"
#include "arrayutil.h"
#include "ioutil.h"
#include "combinatorial.h"
#include "moindex.h"
#include "abecalc.h"
#include "mpi_utilities.h"
#include "binarystr.h"
#include "citruncate.h"
#include "execute_pjayci.h"
#include <mpi.h>

/*
 * execute_pjayci: execute parallel CI.
 */
int execute_pjayci ()
{
        int error = 0;     /* Error flag */

        int electrons = 0; /* total number of electrons */
        int orbitals = 0;  /* total number of orbitals */
        int nfrzc = 0;     /* number of frozen core orbitals */
        int nfrzv = 0;     /* number of frozen virtual orbitals */
        int nactv = 0;     /* number of CAS active orbitals */
        int ndocc = 0;     /* number of CAS doubly-occupied orbitals */
        int xlvl = 0;      /* CI expansion excitation level. must be <= 2 */
        int printlvl = 0;  /* print level */
        int printwvf = 0;  /* print wavefunction info */

        int aelec = 0;    /* Alpha electrons */
        int belec = 0;    /* Beta  eletrons */
        int ci_aelec = 0; /* CI alpha electrons */
        int ci_belec = 0; /* CI beta  electrons */
        int ci_orbs = 0;  /* active CI orbitals */

        struct occstr  *pstrings;  /* Alpha electron strings */
        struct occstr  *qstrings;  /* Beta  electron strings */
        struct eospace *peospace;  /* Alpha electron string spaces */
        struct eospace *qeospace;  /* Beta  electron string spaces */
        int pegrps = 0;            /* Alpha electron occupation groups */
        int qegrps = 0;            /* Beta  electron occupation groups */
        int pstr_len = 0;          /* Number of alpha strings */
        int qstr_len = 0;          /* Number of beta  strings */
        int dtrm_len = 0;          /* Total number of determinants */

        int m1len = 0;                  /* Number of 1-e integrals */
        int m2len = 0;                  /* Number of 2-e integrals */
        double *moints1 = NULL;         /* 1-e integrals */
        double *moints2 = NULL;         /* 2-e integrals */
        double frzcore_e = 0.0;         /* Frozen core energy */
        double nucrep_e  = 0.0;         /* Nuclear repulsion energy */
        char moflname[FLNMSIZE] = {""}; /* SIFS integral filename */
        int itype = 1;                  /* Integral type. Always 1 */
        
        /* Read in the &general namelist. Ensure that the expansion's
         * CAS space is not greater than 64 orbitals. */
        if (mpi_proc_rank == mpi_root) {
                printf("Reading &general input\n");
                readgeninput(&electrons, &orbitals, &nfrzc, &ndocc, &nactv,
                             &xlvl, &nfrzv, &printlvl, &printwvf, &error);
                if (error != 0) {
                        error_flag(mpi_proc_rank, error, "execute_pjayci");
                        return error;
                }
                if ((ndocc + nactv) > 64) {
                        error = ndocc + nactv;
                        error_flag(mpi_proc_rank, error, "execute_pjayci");
                        error_message(mpi_proc_rank,
                                      "CI-active orbital number must be < 64.\n",
                                      "execute_pjayci");
                        return error;
                }
        }
        /* Broadcast &general values. These are needed on all processes. */
        MPI_Bcast(&electrons, 1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&orbitals,  1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzc,     1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&ndocc,     1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nactv,     1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&xlvl,      1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&nfrzv,     1, MPI_INT, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(&printlvl,  1, MPI_INT, mpi_root, MPI_COMM_WORLD);

        /* Get number of alpha/beta electrons */
        abecalc(electrons, &aelec, &belec);
        
        /* Generate pstrings and qstrings. */
        compute_ci_elecs_and_orbitals(aelec, belec, orbitals, nfrzc, nfrzv,
                                      &ci_aelec, &ci_belec, &ci_orbs);
        pstr_len = compute_stringnum(ci_orbs, ci_aelec, ndocc, nactv, xlvl);
        qstr_len = compute_stringnum(ci_orbs, ci_belec, ndocc, nactv, xlvl);
        pstrings = allocate_occstr_arrays(pstr_len);
        qstrings = allocate_occstr_arrays(qstr_len);
        peospace = allocate_eospace_array(ci_aelec, ci_orbs, ndocc, nactv, xlvl,
                                       &pegrps);
        qeospace = allocate_eospace_array(ci_belec, ci_orbs, ndocc, nactv, xlvl,
                                       &qegrps);
        error = citrunc(aelec, belec, orbitals, nfrzc, ndocc, nactv, nfrzv,
                        xlvl, pstrings, pstr_len, qstrings, qstr_len,
                        peospace, pegrps, qeospace, qegrps, &dtrm_len);
        MPI_Barrier(MPI_COMM_WORLD);
       
        /* Read the molecular orbitals */
        m1len = index1e(orbitals, orbitals);
        m2len = index2e(orbitals, orbitals, orbitals, orbitals);
        strncpy(moflname, "moints", FLNMSIZE);
        moints1 = malloc(sizeof(double) * m1len);
        moints2 = malloc(sizeof(double) * m2len);
        init_dbl_array_0(moints1, m1len);
        init_dbl_array_0(moints2, m2len);
        if (mpi_proc_rank == mpi_root) {
                readmointegrals(moints1, moints2, itype, ci_orbs, moflname,
                                m1len, m2len, &nucrep_e, &frzcore_e);
        }
        MPI_Bcast(moints1, m1len, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
        MPI_Bcast(moints2, m2len, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);


        free(moints1);
        free(moints2);
        return error;
}
                
