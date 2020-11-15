// File: atomic_orbitals.c
/*
 * atomic_orbitals: contains atomic orbital evaluation functions.
 * --------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "errorlib.h"
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
 * compute_f_normconst: Compute the d orbital normalization constant.
 * (32768a^9/Pi^3)^(1/4).
 */
double compute_f_normconst(double alpha)
{
    double val;
    val = 32768 * pow(alpha, 9) / M_PI_3;
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
 *
 * Input:
 *  orb = orbital basis function
 *  pos = (x, y, z)
 *  scr = scratch array to hold |R - r| (3 * sizeof(double))
 */
double evaluate_orbital(struct ao_basisfunc orb, double *pos, double *scr)
{
    double value = 0.0; /* Return value */
    
    switch (orb.type) {
    case 1:
        value = evaluate_s(pos, orb.alpha, orb.ccoef, orb.ugaus,
                           orb.gscale, orb.geom);
        break;
    case 2: case 3: case 4:
        value = evaluate_p(pos, orb.alpha, orb.ccoef, orb.ugaus,
                           orb.gscale, orb.geom, orb.type, scr);
        break;
    case 5: case 6: case 7: case 8: case 9:
        value = evaluate_d(pos, orb.alpha, orb.ccoef, orb.ugaus,
                           orb.gscale, orb.geom, orb.type, scr);
        break;
    case 10: case 11: case 12: case 13: case 14: case 15: case 16:
        value = evaluate_f(pos, orb.alpha, orb.ccoef, orb.ugaus,
                           orb.gscale, orb.geom, orb.type, scr);
        break;
    default:
        error_message(0, "Unknown orbital type", "evaluate_orbital");
        break;
    }
    
    return value;
}

/*
 * evaluate_d: evaluate a contracted d orbital at the position (x, y, z).
 *
 * Input:
 *  pos     = x, y, z to evalute orbital at
 *  alphac  = alpha coefficients for contracted gaussians
 *  ccoef   = contraction coefficients
 *  ng      = number of gaussians in contraction
 *  gscale  = gaussian scaling coefficient
 *  atompos = x, y, z of nuclear center
 *  typ     = 5=d2-(xy), 6=d1-(xz), 7=d0 (z2), 8=d1+(yz), 9=d2+(x2-y2)
 *  scr     = scratch array to hold |R - r|
 */
double evaluate_d(double *pos, double *alphac, double *ccoef, int ng,
                  double gscale, double *atompos, int typ, double *scr)
{
    double value = 0.0; /* Return value */
    double rval;
    double dcmpnt; /* d_{type} component {type} (xy, xz, z2, yz, x2-y2) */
    double nrmcnst;
    int i;
    /* Compute {X-x,Y-y,Z-z} and ||{X-x,Y-y,Z-z}||*/
    vector_difference(pos, atompos, 3, scr);
    rval = compute_vector_norm(scr, 3);
    switch (typ) {
    case 5:
        /* xy */
        dcmpnt = scr[0] * scr[1];
        break;
    case 6:
        /* xz */
        dcmpnt = scr[0] * scr[2];
        break;
    case 7:
        /* z2 = 2z^2 - x^2 - y^2*/
        dcmpnt =  2 * scr[2] * scr[2];
        dcmpnt -= scr[0] * scr[0];
        dcmpnt -= scr[1] * scr[1];
        break;
    case 8:
        /* yz */
        dcmpnt = scr[1] * scr[2];
        break;
    case 9:
        /* x2 - y2 */
        dcmpnt =  scr[0] * scr[0];
        dcmpnt -= scr[1] * scr[1];
        break;
    default:
        error_message(0, "Unknown type!", "evaluate_d");
        printf(" type = %d\n", typ);
        return value;
    }
    for (i = 0; i < ng; i++) {
        nrmcnst = compute_d_normconst(alphac[i]);
        value = value + (gscale * ccoef[i] * nrmcnst *
                         dcmpnt * gauss_fcn(alphac[i], rval));
    }
    return value;
}

/*
 * evaluate_f: evaluate a contracted f orbital at the position (x, y, z).
 *
 * Input:
 *  pos     = x, y, z to evalute orbital at
 *  alphac  = alpha coefficients for contracted gaussians
 *  ccoef   = contraction coefficients
 *  ng      = number of gaussians in contraction
 *  gscale  = gaussian scaling coefficient
 *  atompos = x, y, z of nuclear center
 *  typ     = 10=f3-, 11=f2-, 12=f1-, 13=f0,
 *            14=f1+, 15=f2+, 16=f3+
 *  scr     = scratch array to hold |R - r|
 */
double evaluate_f(double *pos, double *alphac, double *ccoef, int ng,
                  double gscale, double *atompos, int typ, double *scr)
{
    double value = 0.0; /* Return value */
    double rval;
    double fcmpnt; /* f_{type} component {type} () */
    double nrmcnst;
    int i;
    /* Compute {X-x,Y-y,Z-z} and ||{X-x,Y-y,Z-z}||*/
    vector_difference(pos, atompos, 3, scr);
    rval = compute_vector_norm(scr, 3);
    switch (typ) {
    case 10:
        /* m=-3, xxy - yyy */
        fcmpnt =  pow(scr[0], 2) * scr[1];
        fcmpnt -= pow(scr[1], 3);
        break;
    case 11:
        /* m=-2, xyz */
        fcmpnt = scr[0] * scr[1] * scr[2];
        break;
    case 12:
        /* m=-1, 4*yzz - 3*yyy - xxy */
        fcmpnt =  4 * scr[1] * scr[2] * scr[2];
        fcmpnt -= 3 * pow(scr[1], 3);
        fcmpnt -= pow(scr[0], 2) * scr[1];
        break;
    case 13:
        /* m=0, xxz + yyz - 2*zzz */
        fcmpnt =  pow(scr[0], 2);
        fcmpnt += pow(scr[1], 2) * scr[2];
        fcmpnt -= 2 * pow(scr[2], 3);
        break;
    case 14:
        /* m=+1, 4*xzz - 3*xxx - xyy */
        fcmpnt =  4 * scr[0] * pow(scr[2], 2);
        fcmpnt -= 3 * pow(scr[0], 3);
        fcmpnt -= scr[0] * pow(scr[1], 2);
        break;
    case 15:
        /* m=+2, xxz - yyz */
        fcmpnt =  pow(scr[0], 2) * scr[2];
        fcmpnt -= pow(scr[1], 2) * scr[2];
        break;
    case 16:
        /* m=+3, xyy - xxx */
        fcmpnt =  scr[0] * pow(scr[1], 2);
        fcmpnt -= pow(scr[0], 3);
        break;
    default:
        error_message(0, "Unknown type!", "evaluate_f");
        printf(" type = %d\n", typ);
        return value;
    }
    for (i = 0; i < ng; i++) {
        nrmcnst = compute_f_normconst(alphac[i]);
        value = value + (gscale * ccoef[i] * nrmcnst *
                         fcmpnt * gauss_fcn(alphac[i], rval));
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
 *  scr     = scratch array to hold |R - r|
 */
double evaluate_p(double *pos, double *alphac, double *ccoef, int ng,
                  double gscale, double *atompos, int typ, double *scr)
{
    double value = 0.0; /* Return value */
    double rval;
    double pcmpnt; /* p_{type} component {type} (x, y, or z) */
    double nrmcnst;
    int i;
    /* Compute {X-x,Y-y,Z-z} and ||{X-x,Y-y,Z-z}||*/
    vector_difference(pos, atompos, 3, scr);
    rval = compute_vector_norm(scr, 3);
    switch (typ) {
    case 2:
        pcmpnt = scr[0];
        break;
    case 3:
        pcmpnt = scr[1];
        break;
    case 4:
        pcmpnt = scr[2];
        break;
    default:
        error_message(0, "Unknown type!", "evaluate_p");
        printf(" type = %d\n", typ);
        return value;
    }
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
    printf("%lf %lf %lf\n", pos[0], pos[1], pos[2]);
    printf("%lf %lf %lf\n", atompos[0], atompos[1], atompos[2]);
    printf("rval = %lf\n", rval);
    for (i = 0; i < ng; i++) {
        nrmcnst = compute_s_normconst(alphac[i]);
        value = value + (gscale * ccoef[i] *
                         nrmcnst * gauss_fcn(alphac[i], rval));
    }
    return value;
}


