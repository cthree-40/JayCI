/* File: execute_dycicalc.c
 *
 * execute_dycicalc: execute dyson orbital calculation from determinant-based
 * CI calculations.
 */
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "abecalc.h"
#include "allocate_mem.h"
#include "binarystr.h"
#include "ioutil.h"
#include "read_wavefunction.h"
#include "dysoncomp.h"
#include "buildao.h"
#include "execute_dycicalc.h"

/*
 * execute_dycicalc: execute dyson orbital calculation from determinant-based
 * CI calculation. (JayCI).
 * Input:
 *  nstates0 = N+1 electron wavefunction CI roots
 *  nelec0   = N+1 electron wavefunction electrons
 *  ndets0   = N+1 electron wavefunction determinants
 *  nstates1 = N   electron wavefunction CI roots
 *  nelec1   = N   electron wavefunction electrons
 *  ndets1   = N   electron wavefunction determinants
 *  orbitals = Total molecular orbitals
 */
int execute_dycicalc(char *wvfcn_file0, char *wvfcn_file1,
                     int nstates0, int nstates1, int nelecs0, int nelecs1,
                     int orbitals, int ndets0, int ndets1)
{
        int error = 0; /* Error flag */
        double **dysonorb = NULL; /* Dyson orbital one-electron function */
        int ndysonorbs = 0; /* number of dyson orbitals to compute */
        
        int naelec0= 0; /* N+1 electron wavefunction alpha electrons */
        int nbelec0= 0; /* n+1 electron wavefunction beta  electrons */
        struct det *detlist0 = NULL; /* N+1 electron wavefunction */
        double **civec0 = NULL; /* N+1 electron wavefunction CI vector(s) */
        double *cival0 = NULL; /* N+1 electron eigenvalues */
        int ninto0 = 0; /* N+1 electron wavefunction internal orbitals */
        int nfrzc0 = 0; /* Number of frozen core orbitals */
        int ndocc0 = 0; /* Number of docc orbitals */
        int nactv0 = 0; /* Number of active orbitals */
        int nfrzv0 = 0; /* NUmber of frozen virtual orbitals */
        int orbs0 = 0; /* Total number of orbitals */
        int elec0 = 0; /* Total number of electrons */

        int xlvl = 0; /* Excitation level (not used) */
        
        int naelec1= 0; /* N electron wavefunction alpha electrons */
        int nbelec1= 0; /* N electron wavefunction beta  electrons */
        struct det *detlist1 = NULL; /* N electron wavefunction */
        double **civec1 = NULL; /* N electron wavefunction CI vector(s) */
        double *cival1 = NULL; /* N electron eigenvalues */
        int ninto1 = 0; /* N electron wavefunction internal orbitals */
        int nfrzc1 = 0; /* Number of frozen core orbitals */
        int ndocc1 = 0; /* Number of docc orbitals */
        int nactv1 = 0; /* Number of active orbitals */
        int nfrzv1 = 0; /* NUmber of frozen virtual orbitals */
        int orbs1 = 0; /* Total number of orbitals */
        int elec1 = 0; /* Total number of electrons */
        
        double *mocoeff = NULL; /* Molecular orbital coefficients */
        int moclen = 1046529; /* Molecular orbital coefficient list size */

        struct ao_basisfunc *aobasis = NULL; /* Atomic orbital basis */
        
        int i = 0;

        /* Check input values. Allocate wavefunctions and civector structures. */
        if (nstates0 != 1) {
                error_message("nstates0 != 1", "execute_dycicalc");
                error = 100;
                return error;
        }
        error = allocate_mem_double(&civec0, ndets0, nstates0);
        error = error + allocate_mem_double(&civec1, ndets1, nstates1);
        if (error != 0) {
                error_message("error allocating civectors.", "execute_dycicalc");
                error_flag(error, "execute_dycicalc");
                return error;
        }
        cival0 = (double *) malloc(nstates0 * sizeof(double));
        cival1 = (double *) malloc(nstates1 * sizeof(double));
        detlist0 = (struct det *) malloc(ndets0 * sizeof(struct det));
        init_detlist(detlist0, ndets0);
        detlist1 = (struct det *) malloc(ndets1 * sizeof(struct det));
        init_detlist(detlist1, ndets1);

        /* Set alpha/beta electron numbers for N+1 and N electron 
         * wavefunctions. */
        abecalc(nelecs0, &naelec0, &nbelec0);
        abecalc(nelecs1, &naelec1, &nbelec1);
        
        /* Read in N+1 electron wavefunction from file. */
        error = read_wavefunction(detlist0, ndets0, nstates0, civec0, cival0,
                                  orbitals, nelecs0, ninto0, wvfcn_file0,
                                  &orbs0, &elec0, &nfrzc0, &ndocc0, &nactv0,
                                  &nfrzv0, &xlvl);
        if (error != 0) {
                error_message("Error reading N+1 wavefunction",
                              "execute_dycicalc");
                return error;
        }
        /* Read in N electron wavefunction from file. */
        error = read_wavefunction(detlist1, ndets1, nstates1, civec1, cival1,
                                  orbitals, nelecs1, ninto1, wvfcn_file1,
                                  &orbs1, &elec1, &nfrzc1, &ndocc1, &nactv1,
                                  &nfrzv1, &xlvl);
        if (error != 0) {
                error_message("Error reading N wavefunction",
                              "execute_dycicalc");
                return error;
        }
        
        /* Compare orbital bassis sets */
        if (orbs1 != orbs0) {
                error++;
                printf("orbs1 = %d\norbs2 = %d\n", orbs1, orbs0);
                error_message("MO Bases are different!", "execute_dycicalc");
                return error;
        }
        
        /* Compute dyson orbital */
        ndysonorbs = nstates1;
        error = allocate_mem_double(&dysonorb, orbitals, ndysonorbs);
        if (error != 0) {
                error_message("error allocating dyson orbital vector.",
                              "execute_dycicalc");
                return error;
        }
        /* If number of alpha/beta electrons are equal between (N+1) and (N)
         * wavefunctions, then alpha/beta strings must be IDENTICAL. */
        if (nelecs0 % 2) {
                /* NELECS0 is ODD */
                for (i = 0; i < ndysonorbs; i++) {
                        compute_dyson_orbital_beq(detlist0, civec0[0], ndets0,
                                                  detlist1, civec1[i], ndets1,
                                                  orbitals, dysonorb[i]);
                }
        } else {
                for (i = 0; i < ndysonorbs; i++) {
                        compute_dyson_orbital_aeq(detlist0, civec0[0], ndets0,
                                                  detlist1, civec1[i], ndets1,
                                                  orbitals, dysonorb[i]);
                }
        }

        /* Read in molecular coefficients */
        mocoeff = (double *) malloc(moclen * sizeof(double));
        readmocoeffs(mocoeff, moclen);

        /* Build atomic orbitals. Need full set:
         * (orbitals + nfrzc + nfrzv) = orbs{0/1}*/
        aobasis = ao_initialize_atomic_orbital_set((orbs0));
        if (aobasis == NULL) {
                error_message("Error initializing atomic orbital set!",
                              "execute_dycicalc");
                error++;
                return error;
        }
        error = ao_buildao(aobasis);
        for (i = 0; i < orbitals; i++) {
                printf(" --- Orbital #%3d ---\n", (i + 1));
                ao_print_orbital_information(aobasis[i]);
        }
        return error;
}
