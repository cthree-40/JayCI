// File: transmo.c
/*
 * Routines to compute transition moments of AO with free
 * electrons.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "transmo.h"

/*
 * compute_k: compute the wavevector k of a free electron given
 * the kinetic energy. (Atomic units: hbar = 1, m_e = 1)
 */
double compute_k(double e)
{
    double kval = 0.0;
    kval = 2.0 * e;
    kval = sqrt(kval);
    return kval;
}
