// File: atomic_orbitals.h

#ifndef atomic_orbitals_h
#define atomic_orbitals_h

/*
 * compute_d_normconst: Compute the d orbital normalization constant.
 * (2048a^7/Pi^3)^(1/4).
 */
double compute_d_normconst(double alpha);

/*
 * compute_f_normconst: Compute the d orbital normalization constant.
 * (32768a^9/Pi^3)^(1/4).
 */
double compute_f_normconst(double alpha);

/*
 * compute_p_normconst: Compute the p orbital normalization constant.
 * (128a^5/Pi^3)^(1/4).
 */
double compute_p_normconst(double alpha);

/*
 * compute_s_normconst: Compute the s orbital normalization constant.
 * (8a^3/PI^3)^(1/4).
 */
double compute_s_normconst(double alpha);

/*
 * evaluate_orbital: evaluate an atomic orbital at some position (x, y, z).
 *
 * Input:
 *  orb = orbital basis function
 *  pos = (x, y, z)
 *  scr = scratch array to hold |R - r| (3 * sizeof(double))
 */
double evaluate_orbital(struct ao_basisfunc orb, double *pos, double *scr);

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
 *  typ     = 5=d2-, 6=d1-, 7=d0, 8=d1+, 9=d2+
 */
double evaluate_d(double *pos, double *alphac, double *ccoef, int ng,
                  double gscale, double *atompos, int typ, double *scr);

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
                  double gscale, double *atompos, int typ, double *scr);

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
                  double gscale, double *atompos, int typ, double *scr);

/*
 * evaluate_s: evaluate a contracted s orbital at the position (x, y, z).
 * Input:
 *  pos    = x, y, z to evalute orbital at
 *  alphac = alpha coefficients for contracted gaussians
 *  ccoef  = contraction coefficients
 *  ng     = number of gaussians in contraction
 *  gscale = gaussian scaling coefficient
 *  atompos= x, y, z of nuclear center
 */
double evaluate_s(double *pos, double *alphac, double *ccoef, int ng,
                  double gscale, double *atompos);

#endif
