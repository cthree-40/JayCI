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
 * compute_dyson_orbital_a: compute the dyson orbital between electronic
 * states of N+1 and N electron wavefuntions by comparing alpha strings.
 * Input:
 *
 */
void compute_dyson_orbital_a(int v0_hndl, int v1_hndl, int w0_hndl, int w1_hndl,
                             struct occstr *pstr0, struct eospace *peosp0, int npe0,
                             struct occstr *qstr0, struct eospace *qeosp0, int nqe0,
                             struct occstr *pstr1, struct eospace *peosp1, int npe1,
                             struct occstr *qstr1, struct eospace *qeosp1, int nqe1,
                             int **pq1, int npq1,
                             int norbs, int ndocc, int nactv, int ndyst0,
                             int *dysnst0, int ndyst1, int *dysnst1, int ndets0,
                             int ndets1, int **strcont, int naelec0,
                             double **dyorb);

/*
 * compute_dyson_orbital_b: compute the dyson orbital between electronic
 * states of N+1 and N electron wavefuntions by comparing beta strings.
 * Input:
 *
 */
void compute_dyson_orbital_b(int v0_hndl, int v1_hndl, int w0_hndl, int w1_hndl,
                             struct occstr *pstr0, struct eospace *peosp0, int npe0,
                             struct occstr *qstr0, struct eospace *qeosp0, int nqe0,
                             struct occstr *pstr1, struct eospace *peosp1, int npe1,
                             struct occstr *qstr1, struct eospace *qeosp1, int nqe1,
                             int **pq1, int npq1,
                             int norbs, int ndocc, int nactv, int dyst0,
                             int *dysnst0, int dyst1, int *dysnst1, int ndets0,
                             int ndets1, int **strcont, int nbelec0,
                             double **dyorb);

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

/*
 * generate_strcontlist: generate contribution list for each string.
 * Input:
 *  str    = string list
 *  nstr   = number of alpha/beta strings
 *  eosp0  = electron orbital space list (N+1 electrons)
 *  ne0    = number of electron orbital spaces (N+1 electrons)
 *  ndocc  = number of docc orbitals
 *  nactv  = number of actv orbitals
 *  nvirt  = number of virt orbitals
 *  nelec1 = number of N-electron alpha/beta electrons
 *  eosp1  = electron orbital space list (N-electron)
 *  ne1    = number of electron orbital spaces (N-electron)
 * Output:
 *  strcont = string contribution list
 */
void generate_strcontlist(struct occstr *str, int nstr, struct eospace *eosp0,
			  int ne0, int ndocc, int nactv, int nvirt,
			  int **strcont, int nelec1, struct eospace *eosp1,
			  int ne1);

/*
 * string_info_to_determinant: compute the determinant index given
 * the p and q string information.
 * Input:
 *  pval    = p string
 *  qval    = q string
 *  peosp   = alpha electron orbital spaces
 *  pegrps  = number of alpha electron orbital spaces
 *  qeosp   = beta  electron orbital spaces
 *  qegrps  = number of beta  electron orbital spaces
 *  pq      = (p,q)-space pairings
 *  npq     = number of (p,q)-space pairings
 * Output:
 *  detindx = determinant index in expansion.
 */
int string_info_to_determinant(int pval, int qval, struct eospace *peosp,
                               int pegrps, struct eospace *qeosp, int qegrps,
                               int **pq, int npq);

/*
 * set_ga_det_indexes_spec: set the array of indices to gather from global array.
 * Input:
 *  jindx   = list of row numbers in vector V
 *  num     = number of row numbers
 *  cols    = number of columns of V
 *  colindx = column indices
 * Output:
 *  vindx   = indices of global array V to gather
 */
void set_ga_det_indexes_spec(int *jindx, int num, int cols, int *colindx,
                             int **vindx);

#endif
