/* File: dysoncomp.c
 *
 * Functions to compute dyson orbital.
 * comparedets_dyson
 * compute_dyson_orbitals
 */
#include <stdio.h>
#include <stdlib.h>
#include "errorlib.h"
#include "binary.h"
#include "bitutil.h"
#include "binarystr.h"
#include "cimapping.h"
#include "action_util.h"
#include "dysoncomp.h"

/*
 * compute_dyson_orbital: compute the dyson orbital between one N+1 electron
 * wavefunction and on N electron wavefunction.
 */
void compute_dyson_orbital(struct det *dlist0, double *civec0, int ndets0,
                           struct det *dlist1, double *civec1, int ndets1,
                           int orbs, double *dysonorb)
{
        int i = 0, j = 0;
        int ndiff = 0; /* number of differences */
        int orbindx = 0; /* orbital index of contribution */

#ifdef DEBUGGING
        char str1[65];
        char str2[65];
#endif
        
        /* Loop over determinants. N electrons must be in the same slots
         * for both determinants. The N+1 electron, in a unique slot, is
         * the MO index for which the matrix element contributes to the
         * dyson orbital. */
        for (i = 0; i < ndets0; i++) {
                for (j = 0; j < ndets1; j++) {
                        /* Compare N and N+1 electron determinants */
                        orbindx = comparedets_dyson(dlist0[i],dlist1[j]);
                        if (orbindx == 0) continue;
#ifdef DEBUGGING
                        printf(" -----------------------------------\n");
                        printf(" Matrix element = < %d | %d >\n", (i+1), (j+1));
                        llint2bin(dlist0[i].astr.byte1, str1);
                        llint2bin(dlist1[j].astr.byte1, str2);
                        printf(" Alpha strings: \n");
                        printf("   %s\n", str1);
                        printf("   %s\n", str2);
                        printf("  Virtual orbitals: \n");
                        printf("   %d  %d\n", dlist0[i].astr.virtx[0],
                               dlist0[i].astr.virtx[1]);
                        printf("   %d  %d\n", dlist1[j].astr.virtx[0],
                               dlist1[j].astr.virtx[1]);
                        llint2bin(dlist0[i].bstr.byte1, str1);
                        llint2bin(dlist1[j].bstr.byte1, str2);
                        printf(" Beta  strings: \n");
                        printf("   %s\n", str1);
                        printf("   %s\n", str2);
                        printf("  Virtual orbitals: \n");
                        printf("   %d  %d\n", dlist0[i].bstr.virtx[0],
                               dlist0[i].bstr.virtx[1]);
                        printf("   %d  %d\n", dlist1[j].bstr.virtx[0],
                               dlist1[j].bstr.virtx[1]);
                        printf("\n");
                        printf(" Orbital index = %d\n", orbindx);
                        printf(" -----------------------------------\n");
#endif

                        dysonorb[(orbindx - 1)] = civec0[i]*civec1[j];
                }
        }
        printf(" Dyson orbital:\n");
        for (i = 0; i < orbs; i++) {
                printf(" %14.8lf\n", dysonorb[i]);
        }
        
                
        return;
}

/*
 * compute_dyson_orbital_beq: compute the dyson orbital between one N+1 electron
 * wavefunction and on N electron wavefunction when both wavefunctions have an
 * equal number of beta electrons.
 */
void compute_dyson_orbital_beq(struct det *dlist0, double *civec0, int ndets0,
                               struct det *dlist1, double *civec1, int ndets1,
                               int orbs, double *dysonorb)
{
        int i = 0, j = 0;
        int ndiff = 0; /* number of differences */
        int orbindx = 0; /* orbital index of contribution */

#ifdef DEBUGGING
        char str1[65];
        char str2[65];
#endif
        
        /* Loop over determinants. N electrons must be in the same slots
         * for both determinants. The N+1 electron, in a unique slot, is
         * the MO index for which the matrix element contributes to the
         * dyson orbital. */
        for (i = 0; i < ndets0; i++) {
                for (j = 0; j < ndets1; j++) {
                        /* Compare N and N+1 electron determinants */
                        orbindx = comparedets_dyson_beq(dlist0[i],dlist1[j]);
                        if (orbindx == 0) continue;
#ifdef DEBUGGING
                        printf(" ----------- POST EVAL -------------\n");
                        printf(" Matrix element = < %d | %d >\n", (i+1), (j+1));
                        llint2bin(dlist0[i].astr.byte1, str1);
                        llint2bin(dlist1[j].astr.byte1, str2);
                        printf(" Alpha strings: \n");
                        printf("   %s\n", str1);
                        printf("   %s\n", str2);
                        printf("  Virtual orbitals: \n");
                        printf("   %d  %d\n", dlist0[i].astr.virtx[0],
                               dlist0[i].astr.virtx[1]);
                        printf("   %d  %d\n", dlist1[j].astr.virtx[0],
                               dlist1[j].astr.virtx[1]);
                        llint2bin(dlist0[i].bstr.byte1, str1);
                        llint2bin(dlist1[j].bstr.byte1, str2);
                        printf(" Beta  strings: \n");
                        printf("   %s\n", str1);
                        printf("   %s\n", str2);
                        printf("  Virtual orbitals: \n");
                        printf("   %d  %d\n", dlist0[i].bstr.virtx[0],
                               dlist0[i].bstr.virtx[1]);
                        printf("   %d  %d\n", dlist1[j].bstr.virtx[0],
                               dlist1[j].bstr.virtx[1]);
                        printf("\n");
                        printf(" Orbital index = %d\n", orbindx);
                        printf(" -----------------------------------\n");
#endif

                        dysonorb[(orbindx - 1)] = civec0[i]*civec1[j];
                }
        }
        printf(" Dyson orbital:\n");
        for (i = 0; i < orbs; i++) {
                printf(" %14.8lf\n", dysonorb[i]);
        }
        
                
        return;
}

/*
 * compute_dyson_orbital_aeq: compute the dyson orbital between one N+1 electron
 * wavefunction and on N electron wavefunction when both wavefunctions have an
 * equal number of alpha electrons.
 */
void compute_dyson_orbital_aeq(struct det *dlist0, double *civec0, int ndets0,
                               struct det *dlist1, double *civec1, int ndets1,
                               int orbs, double *dysonorb)
{
        int i = 0, j = 0;
        int ndiff = 0; /* number of differences */
        int orbindx = 0; /* orbital index of contribution */

#ifdef DEBUGGING
        char str1[65];
        char str2[65];
#endif
        
        /* Loop over determinants. N electrons must be in the same slots
         * for both determinants. The N+1 electron, in a unique slot, is
         * the MO index for which the matrix element contributes to the
         * dyson orbital. */
        for (i = 0; i < ndets0; i++) {
                for (j = 0; j < ndets1; j++) {
                        /* Compare N and N+1 electron determinants */
                        orbindx = comparedets_dyson_aeq(dlist0[i],dlist1[j]);
                        if (orbindx == 0) continue;
#ifdef DEBUGGING
                        printf(" -----------------------------------\n");
                        printf(" Matrix element = < %d | %d >\n", (i+1), (j+1));
                        llint2bin(dlist0[i].astr.byte1, str1);
                        llint2bin(dlist1[j].astr.byte1, str2);
                        printf(" Alpha strings: \n");
                        printf("   %s\n", str1);
                        printf("   %s\n", str2);
                        printf("  Virtual orbitals: \n");
                        printf("   %d  %d\n", dlist0[i].astr.virtx[0],
                               dlist0[i].astr.virtx[1]);
                        printf("   %d  %d\n", dlist1[j].astr.virtx[0],
                               dlist1[j].astr.virtx[1]);
                        llint2bin(dlist0[i].bstr.byte1, str1);
                        llint2bin(dlist1[j].bstr.byte1, str2);
                        printf(" Beta  strings: \n");
                        printf("   %s\n", str1);
                        printf("   %s\n", str2);
                        printf("  Virtual orbitals: \n");
                        printf("   %d  %d\n", dlist0[i].bstr.virtx[0],
                               dlist0[i].bstr.virtx[1]);
                        printf("   %d  %d\n", dlist1[j].bstr.virtx[0],
                               dlist1[j].bstr.virtx[1]);
                        printf("\n");
                        printf(" Orbital index = %d\n", orbindx);
                        printf(" -----------------------------------\n");
#endif

                        dysonorb[(orbindx - 1)] = civec0[i]*civec1[j];
                }
        }
        printf(" Dyson orbital:\n");
        for (i = 0; i < orbs; i++) {
                printf(" %14.8lf\n", dysonorb[i]);
        }
        
                
        return;
}

/*
 * comparedets_dyson: compare two determinants for computation of dyson
 * orbital.
 */
int comparedets_dyson(struct det det0, struct det det1)
{
        int oindex = 0; /* Orbital index of difference */
        int numaxv = 0; /* alpha virtual orbital diffs  */
        int numbxv = 0; /* beta  virtual orbital diffs  */
        int numaxc = 0; /* alpha cas byte orbital diffs */
        int numbxc = 0; /* beta  cas byte orbital diffs */
        int ifo[2] = {0}; /* initial, 'final' orbital */
        long long int diffsa = 0; /* difference byte */
        long long int diffsb = 0; /* difference byte */
        
        int ninto = 64; /* Number of internal orbitals */
        
        /* Test virtual orbital blocks for differences. There can be no
         * more than one difference. */
        numaxv = compute_virt_diffs(det0.astr, det1.astr);
        numbxv = compute_virt_diffs(det0.bstr, det1.bstr);
        if ((numaxv + numbxv) > 1) return oindex;

        /* If there is one virtual difference, than there can be no 
         * differences in the CAS byte. */
        if ((numaxv + numbxv) == 1) {
                if ((det0.astr.byte1 ^ det1.astr.byte1) != 0x00) {
                        return oindex;
                }
                if ((det0.bstr.byte1 ^ det1.bstr.byte1) != 0x00) {
                        return oindex;
                }
                if (numaxv == 1) {
                        virtdiffs_single_rep(det0.astr.virtx, det1.astr.virtx,
                                             ifo);
                } else {
                        virtdiffs_single_rep(det0.bstr.virtx, det1.bstr.virtx,
                                             ifo);
                }
                oindex = ifo[0];
                return oindex;
        }

        /* If here there are no differences in the virtual orbital blocks.
         * The difference must be in the active space bytes. Again, there
         * can only be one difference. */
        numaxc = ndiffbytes(det0.astr.byte1, det1.astr.byte1, ninto, &diffsa);
        numbxc = ndiffbytes(det0.bstr.byte1, det1.bstr.byte1, ninto, &diffsb);
        if ((numaxc + numbxc) > 1) return oindex;

        /* There must be one difference. */
        if (numaxc == 1) {
                nonzerobits(diffsa, ninto, &oindex);
        } else {
                nonzerobits(diffsb, ninto, &oindex);
        }
        
        return oindex;
}

/*
 * comparedets_dyson_beq: compare two determinants for computation of dyson
 * orbital with equal beta electrons.
 */
int comparedets_dyson_beq(struct det det0, struct det det1)
{
        int oindex = 0; /* Orbital index of difference */
        int numaxv = 0; /* alpha virtual orbital diffs  */
        int numbxv = 0; /* beta  virtual orbital diffs  */
        int numaxc = 0; /* alpha cas byte orbital diffs */
        int numbxc = 0; /* beta  cas byte orbital diffs */
        int ifo[2] = {0}; /* initial, 'final' orbital */
        long long int diffsa = 0; /* difference byte */
        
        int ninto = 64; /* Number of internal orbitals */

        /* Beta strings must be equal. */
        numbxv = compute_virt_diffs(det0.bstr, det1.bstr);
        if (numbxv != 0) return oindex;
        if ((det0.bstr.byte1 ^ det1.bstr.byte1) != 0x00) return oindex;
        
        
        /* Test virtual orbital blocks for differences. There can be no
         * more than one difference. */
        numaxv = compute_virt_diffs(det0.astr, det1.astr);
        if (numaxv > 1) return oindex;

        /* If there is one virtual difference, than there can be no 
         * differences in the CAS byte. */
        if (numaxv == 1) {
                if ((det0.astr.byte1 ^ det1.astr.byte1) != 0x00) {
                        return oindex;
                }
                virtdiffs_single_rep(det0.astr.virtx, det1.astr.virtx,
                                     ifo);
                oindex = ifo[0];
                return oindex;
        }

        /* If here there are no differences in the virtual orbital blocks.
         * The difference must be in the active space bytes. Again, there
         * can only be one difference. */
        numaxc = ndiffbytes(det0.astr.byte1, det1.astr.byte1, ninto, &diffsa);
        if (numaxc > 1) return oindex;

        /* There must be one difference. */
        if (numaxc == 1) {
                nonzerobits(diffsa, ninto, &oindex);
        }
        
        return oindex;
}

/*
 * comparedets_dyson_aeq: compare two determinants for computation of dyson
 * orbital with equal alpha electrons.
 */
int comparedets_dyson_aeq(struct det det0, struct det det1)
{
        int oindex = 0; /* Orbital index of difference */
        int numaxv = 0; /* alpha virtual orbital diffs  */
        int numbxv = 0; /* beta  virtual orbital diffs  */
        int numaxc = 0; /* alpha cas byte orbital diffs */
        int numbxc = 0; /* beta  cas byte orbital diffs */
        int ifo[2] = {0}; /* initial, 'final' orbital */
        long long int diffsb = 0; /* difference byte */
        
        int ninto = 64; /* Number of internal orbitals */

        /* Alpha strings must be equal. */
        numaxv = compute_virt_diffs(det0.astr, det1.astr);
        if (numaxv != 0) return oindex;
        if ((det0.astr.byte1 ^ det1.astr.byte1) != 0x00) return oindex;
        
        
        /* Test virtual orbital blocks for differences. There can be no
         * more than one difference. */
        numbxv = compute_virt_diffs(det0.bstr, det1.bstr);
        if (numbxv > 1) return oindex;

        /* If there is one virtual difference, than there can be no 
         * differences in the CAS byte. */
        if (numbxv == 1) {
                if ((det0.bstr.byte1 ^ det1.bstr.byte1) != 0x00) {
                        return oindex;
                }
                virtdiffs_single_rep(det0.bstr.virtx, det1.bstr.virtx,
                                     ifo);
                oindex = ifo[0];
                return oindex;
        }

        /* If here there are no differences in the virtual orbital blocks.
         * The difference must be in the active space bytes. Again, there
         * can only be one difference. */
        numbxc = ndiffbytes(det0.bstr.byte1, det1.bstr.byte1, ninto, &diffsb);
        if (numbxc > 1) return oindex;

        /* There must be one difference. */
        if (numbxc == 1) {
                nonzerobits(diffsb, ninto, &oindex);
        }
        
        return oindex;
}

