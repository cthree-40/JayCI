// File: straddress.c
/********************************************************************
 * straddress
 * ----------
 * Contains string addressing utilities
 *
 * Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "combinatorial.h"
#include "straddress.h"

/* str_adr2str: Computes orbital index string from address
 * -------------------------------------------------------------------
 * [citation]
 * 
 * Input:
 *  index = index of string
 *  elec  = number of electrons
 *  orbs  = number of orbitals 
 *  scr   = scratch array
 * Output:
 *  str = orbital index string */
void str_adr2str(int index, int *scr, int elec, int orbs, int *str)
{
    int i;
    if (index == 0) return;
    for (i = 0; i < elec; i++) {
	scr[i] = i + 1;
    }
    str_strfind2(scr, 1, elec, orbs, index, str);
    return;

}    
/* str_adrfind: Computes address of orbital index string
 * -------------------------------------------------------------------
 * [citation]
 *
 * Input:
 *  str  = orbital index string
 *  elec = electron number
 *  orbs = MO's in system */
int str_adrfind(int *str, int elec, int orbs)
{
     int address;
     int i, j;

     address = 1;
     for (i = 2; i <= elec; i++) {
	  for (j = (str[i - 2] + 1); j <= (str[i - 1] - 1); j++ ) {
	       address = address +
		    binomial_coef2((orbs - j),(elec - i));
	  }
     }
     for (i = 1; i < str[0]; i++) {
	  address = address +
	       binomial_coef2((orbs - j),(elec - 1));
     }
     return address;
}
/* str_strfind1: 
 *   Compute orbital occupation string from preceding string 
 * -------------------------------------------------------------------
 * [citation]
 *
 * Input:
 *  str1 = orbital index string 1  
 *  elec = electrons
 *  orbs = orbitals
 * Output:
 *  str2 = orbital index string 2 */
void str_strfind1(int *str1, int elec, int orbs, int *str2)
{
     int i, j;

     /* loop over electrons */
     for (i = 0; i < elec; i++) {
	  if (str1[(elec - i) - 1] == (orbs - i)) {
	       continue; //i->i+1
	  } else {
	       /* determine new string */
	       for (j = 0; j < elec - i - 1; j++) {
		    str2[j] = str1[j];
	       }
	       for (j = (elec - i); j <= elec; j++) {
		    str2[j-1] =
			 str1[(elec - i) - 1] + 1 + j - (elec - i);
	       }
	  }
	  break;
     }

     return;

}
/* str_strfind2:
 *  Compute orbital occupation string for jth string from ith
 * -------------------------------------------------------------------
 * [citation]
 *
 * Input:
 *  str1  = ith string (this string becomes j-1)
 *  indx1 = ith index
 *  elec  = electrons
 *  orbs  = orbitals
 *  indx2 = jth index
 * Output:
 *  str2  = jth string */
void str_strfind2(int *str1, int indx1, int elec, int orbs, int indx2,
		  int *str2)
{
    int i, j;

    for (i = 0; i < elec; i++) {
            str2[i] = str1[i];
    }
    for (i = 0; i < indx2 - indx1; i++) {
	str_strfind1(str1, elec, orbs, str2);
	for (j = 0; j < elec; j++) {
	    str1[j] = str2[j];
	}
    }
    
    return;
}
	  
