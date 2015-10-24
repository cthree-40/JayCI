/* test.c */
/*
 * Test driver for subfunctions of jayci.x
 */

#include <stdio.h>
#include <stdlib.h>
#include "arrayutil.h"
#include "binarystr.h"
#include "action_util.h"

int main()
{
        long long int stri = 0;
        long long int xi = 0, xf = 0;
        int ninto = 7;
        int istr[5];
        int i;
        int pindx;
        /* stri = 1,3,4,5,6
         * strj = 1,2,3,4,6
         * xi = 5
         * xj = 2 */
        istr[0] = 1;
        istr[1] = 3;
        istr[2] = 4;
        istr[3] = 5;
        istr[4] = 6;

        for (i = 0; i < 5; i++) {
                stri = stri + (1 << (istr[i] - 1));
        }
        xi = (1 << 4);
        xf = (1 << 1);
        
        pindx = pindex_single_rep_cas(stri, xi, xf, ninto);

        printf("pindx = %d\n", pindx);
}
