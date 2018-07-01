// File: execute_pjayci.h

#ifndef execute_pjayci_h
#define execute_pjayci_h

/*
 * execute_pjayci: execute parallel CI.
 */
int execute_pjayci ();

/*
 * compute_cimat_chunks: compute chunksize of bounds of H for evaluation.
 * Input:
 *  dlen = number of determinants
 */
void compute_cimat_chunks (int dlen, int *chunk, int *lwrbnd, int *uppbnd);

/*
 * compute_diagonal_matrix_elements: compute <i|H|i> elements.
 * Input:
 *  hdgls   = digonal elements <i|H|i> for (i = start,...,final)
 *  start   = starting determinant index
 *  final   = final determinant index
 *  mo1     = 1-e integrals
 *  mo2     = 2-e integrals
 *  aelec   = CI alpha electrons
 *  belec   = CI beta  electrons
 *  peosp   = alpha electron orbital spaces
 *  pegrps  = number of alpha electron orbital spaces
 *  qeosp   = beta  electron orbital spaces
 *  qegrps  = number of beta  electron orbital spaces
 *  pq      = (p,q)-space pairings
 *  npq     = number of (p,q)-space pairings
 *  pstrings= alpha strings
 *  qstrings= beta  strings
 *  intorb  = internal orbitals (DOCC + CAS)
 */
void compute_diagonal_matrix_elements(double *hdgls, int start, int final,
                                      double *mo1, double *mo2, int aelec,
                                      int belec, struct eospace *peosp,
                                      int pegrps,struct eospace *qeosp,
                                      int qegrps, int **pq, int npq,
                                      struct occstr *pstrings,
                                      struct occstr *qstrings,
                                      int intorb);

/*
 * determinant_string_info: compute string information given
 * a determinant.
 * Input:
 *  detindx = determinant index (position in list)
 *  peosp   = alpha electron orbital spaces
 *  pegrps  = number of alpha electron orbital spaces
 *  qeosp   = beta  electron orbital spaces
 *  qegrps  = number of beta  electron orbital spaces
 *  pq      = (p,q)-space pairings
 *  npq     = number of (p,q)-space pairings
 * Output:
 *  pqval = (p,q) pairing start
 *  pval  = p string of (p,q)-space pairing
 *  qval  = q string of (p,q)-space pairing
 */
void determinant_string_info(int detindx, struct eospace *peosp, int pegrps,
                             struct eospace *qeosp, int qegrps, int **pq,
                             int npq, int *pqval, int *pval, int *qval);

/*
 * get_upptri_element_index: get index of an element (i,j) in
 * list of elements in upper triangle of H matrix (n x n).
 */
int get_upptri_element_index (int i, int j, int n);

/*
 * get_upptri_element_rownumber: get row number of upper triangle
 * matrix list element.
 */
int get_upptri_element_rownumber (int element, int n);

/*
 * get_upptri_size: compute the size of H matrix upper triangle.
 */
int get_upptri_size (int n);



#endif
