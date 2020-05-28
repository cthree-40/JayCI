// File: run_pdycicalc.h

#ifndef run_pdycicalc_h
#define run_pdycicalc_h

/*
 * run_pdycicalc: Execute CI dyson orbital calculation.
 */
int run_pdycicalc ();

/*
 * check_wavefunction_input: check user input of global wavefunction variables.
 * This is executed on all processes.
 */
int check_wavefunction_input(int nelecs0, int norbs0, int nfrzc0,  int ndocc0,
                             int nactv0,  int nfrzv0, int nelecs1, int norbs1,
                             int nfrzc1,  int ndocc1, int nactv1,  int nfrzv1);
        
#endif
