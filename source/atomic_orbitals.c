// File: atomic_orbitals.c
/*
 * atomic_orbitals: contains atomic orbital evaluation functions.
 * --------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "buildao.h"
#include "mathutil.h"
#include "atomic_orbitals.h" 

/*
 * compute_d_normconst: Compute the d orbital normalization constant.
 * (2048a^7/Pi^3)^(1/4).
 */
double compute_d_normconst(double alpha)
{
        double val;
        val = 2048 * pow(alpha, 7) / M_PI_3;
        val = pow(val, 1.0/4.0);
        return val;
}

/*
 * compute_p_normconst: Compute the p orbital normalization constant.
 * (128a^5/Pi^3)^(1/4).
 */
double compute_p_normconst(double alpha)
{
        double val;
        val = 128 * pow(alpha, 5) / M_PI_3;
        val = pow(val, 1.0/4.0);
        return val;
}

/*
 * compute_s_normconst: Compute the s orbital normalization constant.
 * (8a^3/PI^3)^(1/4).
 */
double compute_s_normconst(double alpha)
{
        double val;
        val = 8 * pow(alpha, 3) / M_PI_3;
        val = pow(val, 1.0/4.0);
        return val;
}

/*
 * evaluate_orbital: evaluate an atomic orbital at some position (x, y, z).
 */
double evaluate_orbital(struct ao_basisfunc orb, double *pos)
{
        double value = 0.0; /* Return value */

        switch (orb.type) {
        case 1:
                value = evaluate_s(pos, orb.alpha, orb.ccoef, orb.ugaus,
                                   orb.gscale, orb.geom);
                break;
        case 2: case 3: case 4:
                value = evaluate_p(pos, orb.alpha, orb.ccoef, orb.ugaus,
                                   orb.gscale, orb.geom, orb.type);
                break;
        case 5: case 6: case 7: case 8: case 9:
                value = 0.0;
                break;
        default:
                error_message("Unknown orbital type", "evaluate_orbital");
                break;
        }

        return value;
}

/*
 * evaluate_p: evaluate a contracted p orbital at the position (x, y, z).
 *
 * Input:
 *  pos     = x, y, z to evalute orbital at
 *  alphac  = alpha coefficients for contracted gaussians
 *  ccoef   = contraction coefficients
 *  ng      = number of gaussians in contraction
 *  gscale  = gaussian scaling coefficient
 *  atompos = x, y, z of nuclear center
 *  typ     = px=2, py=3, pz=4
 */
double evaluate_p(double *pos, double *alphac, double *ccoef, int ng,
                  double gscale, double *atompos, int typ)
{
        double value = 0.0; /* Return value */
        double rval;
        double pcmpnt; /* p_{type} component {type} (x, y, or z) */
        double nrmcnst;
        int i;
        switch (typ) {
        case 2:
                pcmpnt = pos[0] - atompos[0];
                break;
        case 3:
                pcmpnt = pos[1] - atompos[1];
                break;
        case 4:
                pcmpnt = pos[2] - atompos[2];
                break;
        default:
                error_message("Unknown type!", "evaluate_p");
                printf(" type = %d\n", typ);
                return value;
        }
                
        rval = compute_3d_distance(atompos, pos);
        for (i = 0; i < ng; i++) {
                nrmcnst = compute_p_normconst(alphac[i]);
                value = value + (gscale * ccoef[i] * nrmcnst *
                                 pcmpnt * gauss_fcn(alphac[i], rval));
        }
        return value;
}

/*
 * evaluate_s: evaluate a contracted s orbital at the position (x, y, z).
 * Input:
 *  pos     = x, y, z to evalute orbital at
 *  alphac  = alpha coefficients for contracted gaussians
 *  ccoef   = contraction coefficients
 *  ng      = number of gaussians in contraction
 *  gscale  = gaussian scaling coefficient
 *  atompos = x, y, z of nuclear center
 */
double evaluate_s(double *pos, double *alphac, double *ccoef, int ng,
                  double gscale, double *atompos)
{
        double value = 0.0; /* Return value */
        double rval;  /* Sqrt(pos - atompos) */
        double nrmcnst;
        int i;
        rval = compute_3d_distance(atompos, pos);
        for (i = 0; i < ng; i++) {
                nrmcnst = compute_s_normconst(alphac[i]);
                value = value + (gscale * ccoef[i] *
                                 nrmcnst * gauss_fcn(alphac[i], rval));
        }
        return value;
}
        
        
