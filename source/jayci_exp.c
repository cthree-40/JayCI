// FILE: jayci_exp.c
/*********************************************************************
 * jayci_exp
 * ---------
 * Program to perform truncation of CI wavefunction.
 * 
 * Needs:
 *  jayci.in
 *
 * -----------------------
 * By Christopher L Malbon
 * Dept of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "ioutil.h"
#include "abecalc.h"
#include "citruncate.h"

int main()
/* main driver
 * -----------
 * Execution outline:
 *  Checks for all necessary files.
 *  Reads input.(&general)
 *  Computes number of alpha/beta electrons.
 *  Truncates space generating the following files:
 *   det.list (list of determinants)
 *   str.list (alpha/beta string lists)
 *  Generates binary string lists:
 *   bin.det.list ([alpha string] [beta string] binary digits)
 *   bin.str.list (alpha/beta binary digit lists)
 */
{
     /* &general scalars 
      * elec   = total electrons in system 
      * orbs   = MO's in system
      * nfrzc  = number of frozen core orbitals
      * ndocc  = number of doubly-occupied orbitals
      * nactv  = number of active orbitals
      * nfrzv  = number of frozen virtual orbitals
      * xlevel = excitation level
      * prntlvl= print level 
      * printwvf = print wavefunction */
     int elec, orbs;
     int nfrzc, ndocc, nactv, nfrzv;
     int xlevel;
     int prntlvl;
     int printwvf;

     /* local scalars
      * err    = error handling 
      * aelec  = alpha electrons
      * belec  = beta  electrons 
      * astr_len = number of alpha strings
      * bstr_len = number of beta  strings 
      * dtrm_len = total number of truncated determiants */ 
     int err;
     int aelec, belec;
     int astr_len, bstr_len, dtrm_len;
     
     err = 0;
     
     /* check for presence of input files */
     err = checkinputfiles();
     if (err == 1) {
	  /* missing jayci.in, terminate. */
	  fprintf(stderr,"*** ERROR! Missing input file: jayci.in ***\n");
	  exit(err);
     }
     if (err == 2 || err == 3) {
	  /* missing input.jayci (OK, we are making this file now.) */
	  /*  or missing moints (OK, we do not need integrals yet.) */
	  fprintf(stderr,"WARNING! Missing integral file: moints\n");
     }

     /* read the &general namelist */
     readgeninput(&elec, &orbs, &nfrzc, &ndocc, &nactv, &xlevel,
		  &nfrzv, &prntlvl, &printwvf, &err);
     if (err != 0) {
	  /* error reading namelist */
	  fprintf(stderr,
		  "*** ERROR! Error reading &general namelist. ***\n");
	  exit(err);
     }

     /* compute alpha/beta electron numbers */
     abecalc(elec, &aelec, &belec);
     printf("Alpha electrons = %3d\nBeta electrons  = %3d\n",
            aelec, belec);
     /* truncate ci-space, generating expansion */
     err = citrunc(aelec, belec, orbs, nfrzc, ndocc, nactv, nfrzv, xlevel,
		   &astr_len, &bstr_len, &dtrm_len);
     if (err != 0) {
	  // ERROR HANDLING //
	  exit(err);
     }
     
     fprintf(stdout, "Alpha strings = %d\n", astr_len); 
     fprintf(stdout, "Beta  strings = %d\n", bstr_len); 
     fprintf(stdout, "Determinants  = %d\n", dtrm_len);

     /* generate input for jayci.x */
     fprintf(stdout, "Generating input for jayci.x\n");
     err = geninput(dtrm_len, astr_len, bstr_len, aelec, belec, orbs, nfrzc,
		    nfrzv);
     
     return err;
}
