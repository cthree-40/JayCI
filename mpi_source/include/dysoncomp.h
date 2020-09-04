// File: dysoncomp.h

#ifndef dysoncomp_h
#define dysoncomp_h

/*
 * build_ppo_triples: build (p0, p1, o) triples for dyson evaluation.
 */
int build_ppo_triples(int **ppo, struct occstr *str0, int nstr0,
		      struct occstr *str1, int nstr1, int ninto);

/*
 * comparestrings_dyson: compare the strings of N and N+1 determinants
 * to get the index of orbital contribution.
 */
int comparestrings_dyson(struct occstr str0, struct occstr str1, int ninto);

/*
 * compute_dyson_orbital: compute the dyson orbital between electronic
 * states of N+1 and N electron wavefuntions by comparing alpha/beta strings.
 * Input:
 *  w0_hndl = N+1 wavefunction handle
 *  dlen0   = number of N+1 determinants
 *  w1_hndl = N   wvaefunction handle
 *  dlen1   = number of N determinants
 *  v0_hndl = N+1 CI vectors handle
 *  v1_hndl = N   CI vectors handle
 *  str0    = N+1 alpha/beta strings
 *  str1    = N   alpha/beta strings
 *  spindx  = alpha: 0,  beta: 1
 *  spindx2 = beta:  1, alpha: 0
 *  ninto0  = N+1 internal orbitals
 *  ninto1  = N   internal orbitals
 *  norbs   = total number of orbitals
 *  dysnst0 = dyson orbital N+1 states
 *  ndyst0  = number of N+1 states for dyson orbitals
 *  dysnst1 = dyson orbital N   states
 *  ndyst1  = number of N   states for dyson orbitals
 *  ndyorbs = number of dyson orbitals
 *  dyorb   = dyson orbital
 */
void compute_dyson_orbital(int w0_hndl, int dlen0, int w1_hndl, int dlen1,
                           int v0_hndl, int v1_hndl,
                           struct occstr *str0, struct occstr *str1,
                           int spindx1, int spindx2, int ninto0, int ninto1,
                           int norbs, int *dysnst0, int ndyst0, int *dysnst1,
                           int ndyst1, int ndyorbs, double **dyorb);

/*
 * compute_det_contributions: compute determinant contributions to dyson
 * orbitals.
 * Input:
 *  w0      = p, q, CAS
 *  v0      = CI vectors
 *  v0_rows = rows of CI vectors
 *  v0_cols = columns of CI vectors
 *  w1      = p, q, CAS
 *  v1      = CI vectors
 *  v1_rows = rows of CI vectors
 *  v1_cols = columns of CI vectors
 *  sp      = spin index (0 = alpha, 1 = beta)
 *  spx     = opposite spin index (1 = beta, 0 = alpha)
 *  str0    = p/q strings
 *  str1    = p/q strings
 *  dyst0   = dyson orbital N+1 electron states
 *  ndyst0  = number of dyson orbital N+1 electron states
 *  dyst1   = dyson orbital N electron states
 *  ndyst1  = number of dyson orbital N electron states
 *  ninto   = number of internal orbitals
 * Output:
 *  dyorb   = dyson orbitals
 */
void compute_det_contributions(int **w0, double **v0, int v0_rows, int v0_cols,
                               int **w1, double **v1, int v1_rows, int v1_cols,
                               int sp, int spx,
                               struct occstr *str0, struct occstr *str1,
                               int *dyst0, int ndyst0, int *dyst1, int ndyst1,
                               double **dyorb, int ninto);

#endif
