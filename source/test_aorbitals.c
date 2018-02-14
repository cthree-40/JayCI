// test_aorbitals.c
/*
 * Test the functions evaluating atomic orbitals.
 */
#include <stdio.h>
#include <stdlib.h>
#include "buildao.h"
#include "mathutil.h"
#include "atomic_orbitals.h"

main ( )
{
        double *origin, *alphac, *ccoef;
        double *evalpt;
        int ng;
        double gscale;

        double result;

        /* One gaussian */
        ng = 1;
        
        origin = (double *) malloc(sizeof(double) * 3);
        alphac = (double *) malloc(sizeof(double) * ng);
        ccoef  = (double *) malloc(sizeof(double) * ng);

        evalpt = (double *) malloc(sizeof(double) * 3);
   
        gscale = 1.0;

        /* Origin placed at (0,0,0) */
        origin[0] = 0.0; 
        origin[1] = 0.0; 
        origin[2] = 0.0;
        /* Alpha coefficient */
        alphac[0] = 1.0;
        /* Contraction coefficient */
        ccoef[0]  = 1.0;

        /* Evaluation point at (0.5,1.0,0) */
        evalpt[0] = 0.5;
        evalpt[1] = 1.0;
        evalpt[2] = 0.0;

        result = evaluate_s(evalpt, alphac, ccoef, ng, gscale, origin);

        printf(" Result = %lf\n", result);

}
        
        
