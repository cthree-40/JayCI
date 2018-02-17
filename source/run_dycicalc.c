// File: run_dycicalc.c
/*
 * run_dycicalc: Execute CI dyson orbital calculation.
 */
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "ioutil.h"
//#include "execute_dycicalc.h"
#include "citruncate.h"
#include "run_dycicalc.h"

/*
 * run_dycicalc: Execute CI dyson orbital calculation.
 * (N+1)electron wavefunction = 0
 * (N)  electron wavefunction = 1
 */
int run_dycicalc ()
{
        /* input wavefunction information */
        int nstates0 = 0, nstates1 = 0; /* Number of states */
        int nelecs0  = 0,  nelecs1 = 0; /* Number of electrons */
        int naelec0  = 0,  naelec1 = 0; /* Number of alpha elctrons */
        int nbelec0  = 0,  nbelec1 = 0; /* Number of beta  electrons */
        int nfrzc0   = 0,   nfrzc1 = 0; /* Number of frozen core orbitals */
        int ndocc0   = 0,   ndocc1 = 0; /* Number of DOCC orbitals */
        int nactv0   = 0,   nactv1 = 0; /* Number of CAS  orbitals */
        int nfrzv0   = 0,   nfrzv1 = 0; /* Number of frozen virtual orbitals */
        int xlvl0    = 0,    xlvl1 = 0; /* CI excitation level */
        int norbs0   = 0,   norbs1 = 0; /* Number of orbitals */
        int ciorbs0  = 0,  ciorbs1 = 0; /* Number of CI orbitals */
        int ciaelec0 = 0, ciaelec1 = 0; /* CI alpha/beta electrons */
        int cibelec0 = 0, cibelec1 = 0; /* CI alpha/beta electrons */
        
        /* wavefunctions */
        int ndets0 = 0;              /* Number of (anion)   determinants */
        int ndets1 = 0;              /* Number of (neutral) determinants */
        struct det *detlist0 = NULL; /* (Anion)   determinants */
        struct det *detlist1 = NULL; /* (Neutral) determinants */
        double **civec0 = NULL;      /* (Anion)   CI vectors */
        double **civec1 = NULL;      /* (Neutral) CI vectors */
        double *cival0  = NULL;      /* (Anion)   CI eigenvalues */
        double *cival1  = NULL;      /* (Neutral  CI eigenvalues */
        
        int astrlen = 0, bstrlen = 0; /* Number of a/b strings (not used) */
        
        int error = 0; /* Error flag */

        /* Read the &wavefcn0 and &wavefcn1 namelists. Build wavefunctions. */
        readwf0input(&nelecs0, &norbs0, &nfrzc0, &ndocc0, &nactv0, &xlvl0,
                     &nfrzv0,  &nstates0, &error);
        if (error != 0) {
                error_message("Error reading &wavefcn0!","run_dycicalc");
                return error;
        }
        readwf1input(&nelecs1, &norbs1, &nfrzc1, &ndocc1, &nactv1, &xlvl1,
                     &nfrzv1,  &nstates1, &error);
        if (error != 0) {
                error_message("Error reading &wavefcn1!","run_dycicalc");
                return error;
        }
        abecalc(nelecs0, &naelec0, &nbelec0);
        abecalc(nelecs1, &naelec1, &nbelec1);
        printf("(N+1): a=%d  b=%d\n", naelec0, nbelec0);
        printf("(N):   a=%d  b=%d\n", naelec1, nbelec1);
        detlist0 = citrunc_rtnlist(naelec0, nbelec0, norbs0, nfrzc0, ndocc0,
                                   nactv0, nfrzv0, xlvl0, &astrlen, &bstrlen,
                                   &ndets0, &ciorbs0, &ciaelec0, &cibelec0);
        detlist1 = citrunc_rtnlist(naelec1, nbelec1, norbs1, nfrzc1, ndocc1,
                                   nactv1, nfrzv1, xlvl1, &astrlen, &bstrlen,
                                   &ndets1, &ciorbs1, &ciaelec1, &cibelec1);

        return error;
}
        
