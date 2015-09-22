#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "combinatorial.h"
#include "straddress.h"
#include "binarystr.h"
#include "binary.h"
#include "moindex.h"
#include "ioutil.h"

void main()
{
     int orb = 23;
     int m1len, m2len;

     double *moints1, *moints2;
     double nrep, fce;
     unsigned char moflname[FLNMSIZE];
     int type;

     int elec, orbs;
     int nfrzc, nfrzv, ndocc, nactv;
     int xlvl, prtlvl;

     int *astring, *bstring;
     int *astring2, *bstring2;
     int aelec, belec;
     int aindx, bindx;

     struct occstr string1;
     char binstr1[65];
     
     int err;
     int i;
     
     m1len = index1e(orb, orb);
     m2len = index2e(orb, orb, orb, orb);

     moints1 = malloc(m1len * sizeof(double));
     moints2 = malloc(m2len * sizeof(double));

     type = 1;
     strncpy(moflname,"moints", FLNMSIZE);

     printf(" Reading molecular integrals from file: %s\n", moflname);
     printf("  %5d 1-e integrals\n  %5d 2-e integrals\n", m1len, m2len);

     readmointegrals(moints1, moints2, type, orb,
		     moflname, m1len, m2len, &nrep, &fce);

     for (i=0; i<m1len; i++) {
	  printf(" %15.8lf\n", moints1[i]);
     }
     printf("Nuc Rep   = %15.8lf\n", nrep);
     printf("FC Energy = %15.8lf\n", fce);
     
     printf(" Reading namelist.\n");

     readgeninput(&elec, &orbs, &nfrzc, &ndocc, &nactv,
		  &xlvl, &nfrzv, &prtlvl, &err);
     printf("  electrons = %5d\n", elec);
     printf("  orbitals  = %5d\n", orbs);
     printf("  nfrzc     = %5d\n", nfrzc);
     printf("  ndocc     = %5d\n", ndocc);
     printf("  nactv     = %5d\n", nactv);
     printf("  nfrzv     = %5d\n", nfrzv);

     /* test addressing */
     aelec = 3;
     astring = malloc(aelec * sizeof(int));
     astring2= malloc(aelec * sizeof(int));
     astring[0] = 1;
     astring[1] = 3;
     astring[2] = 4;
     aindx = str_adrfind(astring, aelec, 6);
     printf("Address = %d\n", aindx);
     
     str_strfind1(astring, aelec, 6, astring2);
     printf("String1 = %d %d %d\n", astring[0], astring[1],
	    astring[2]);
     printf("String2 = %d %d %d\n", astring2[0], astring2[1],
	    astring2[2]);

     string1 = str2occstr(astring2, aelec, 1, 3);
     printf("%d %d %d %d\n", string1.byte1, string1.byte2,
	    string1.virtx[0], string1.virtx[1]);

     llint2bin(string1.byte1, binstr1);
     binstr1[64] = '\0';
     printf(" %.*s\n", 3, binstr1);
     printf(" ---\n");
     free(moints1);
     free(moints2);

}

	
