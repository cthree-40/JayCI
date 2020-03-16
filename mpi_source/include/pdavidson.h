// File: pdavidson.h

#ifndef pdavidson_h
#define pdavidson_h

/*
 * pdavidson: parallel implementation of davidson algorithm.
 */
int pdavidson (struct occstr *pstrings, struct eospace *peospace, int pegrps,
               struct occstr *qstrings, struct eospace *qeospace, int qegrps,
               int **pq_space_pairs, int num_pq, double *moints1, double *moints2,
               int aelec, int belec, int intorb, int ndets, double nucrep_e,
               double frzcore_e, int printlvl, int maxiter, int krymin,
               int krymax, int nroots, int prediagr, int refdim, double restol);

/*
 * add_new_vector: add a new vector to basis space.
 */
void add_new_vector(int v_hndl, int cdim, int len, int n_hndl);


/*
 * build_init_guess_vectors: build initial guess basis vectors for davidson
 * procedure.
 * Input:
 *  n     = input routine: 2 = subblock diagonalization, 1 = unit vectors
 *  v     = basis vectors handle (GLOBAL ARRAY)
 *  dim   = dimension of subspace
 *  kmin  = krylov space minimum dimension
 *  ndets = number of determinants
 */
void build_init_guess_vectors(int n, int v, int dim, int kmin, int ndets,
                              struct occstr *pstr, struct eospace *peosp,
                              int pegrps, struct occstr *qstr,
                              struct eospace *qeosp, int qegrps,
                              int **pqs, int num_pq, double *m1, double *m2,
                              int aelec, int belec, int intorb);

/*
 * compute_cimat_chunks: compute chunksize of bounds of H for evaluation.
 * Input:
 *  dlen = number of determinants
 */
void compute_cimat_chunks (int dlen, int *chunk, int *lwrbnd, int *uppbnd);

/*
 * compute_cimat_chunks_2: compute chunksize of bounds of H for evaluation.
 * Only elements in the upper triangle are computed.
 * Input:
 *  dlen = number of determinants
 */
void compute_cimat_chunks_2 (int dlen, int *chunk, int *row_lo, int *row_hi,
                             int *col_lo, int *col_hi);

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
 * compute_hv_newvector: compute Hv=c for newest vector in basis space.
 * Input:
 *  v_hndl = GA handle for basis vectors (ckdim is new vector)
 *  c_hndl = GA handle for Hv=c vectors
 *  ckdim  = current dimension of space
 */
void compute_hv_newvector(int v_hndl, int c_hndl, int ckdim, struct occstr *pstr,
                          struct eospace *peosp, int pegrps, struct occstr *qstr,
                          struct eospace *qeosp, int qegrps, int **pqs,
                          int num_pq, double *m1, double *m2, int aelec,
                          int belec, int intorb, int ndets, int kmax);


/*
 * compute_GA_norm: compute the norm of a GA vector.
 */
void compute_GA_norm (int r_hndl,  double *norm);

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
 * diag_subspacehmat: diagonalize the v.Hv subspace hamiltonian, returning
 * eigenvectors and eigenvalues.
 * Input:
 *  vhv = v.Hv
 *  hevec = Eigenvectors of v.Hv
 *  heval = Eigenvalues of v.Hv
 *  cdim  = current dimension of v.Hv
 *  maxdim = maximum dimension of v.Hv
 * Output:
 *  error = error flag
 */
int diag_subspacehmat(double **vhv, double **hevec, double *heval, int cdim,
                      int maxdim, double *vhv_scr, double *hevec_scr);

/*
 * evaluate_hdblock_ij: evaluate a block of H:
 *  H(i,j)*V(j,k)=C(i,k)
 *  Input:
 *   vrows  = columns of C
 *   vcols  = rows of C
 *   crows  = columns of C
 *   ccols  = rows of C
 *   v      = V(j,k)
 *   c      = C(i,k)
 *   starti = starting index i
 *   finali = ending index i
 *   startj = starting index j
 *   finalj = ending index j
 *   ndets  = number of determinants
 *   peosp  = alpha electron orbital spaces
 *   pegrps = number of alpha string orbital spaces
 *   pstr   = alpha electron strings
 *   qeosp  = beta electron orbital spaces
 *   qegrps = number of beta string orbital spaces
 *   qstr   = beta electron strings
 *   pq     = (p,q)-space pairings
 *   npq    = number of (p, q)-space pairings
 *   mo1    = 1-e integrals
 *   mo2    = 2-e integrals
 *   aelec  = alpha electrons
 *   belec  = beta  electrons
 *   intorb = internal orbitals (DOCC + CAS)
 *  Output:
 *   C      = C(i,k)
 */
void evaluate_hdblock_ij(int pq_start_i, int pstart_i, int qstart_i,
                         int pq_final_i, int pfinal_i, int qfinal_i,
                         int pq_start_j, int pstart_j, int qstart_j,
                         int pq_final_j, int pfinal_j, int qfinal_j,
                         int vrows, int vcols, double **v,
                         int crows, int ccols, double **c,
                         int starti, int finali, int startj, int finalj,
                         int ndets, struct eospace *peosp, int pegrps,
                         struct occstr *pstr,
                         struct eospace *qeosp, int qegrps, 
                         struct occstr *qstr,int **pq,
                         int npq, double *mo1, double *mo2, int aelec,
                         int belec, int intorb);

/*
 * evaluate_hdblock_ij_1d: evaluate a block of H: **(For one column)**
 *  H(i,j)*V(j)=C(i)
 *  Input:
 *   vrows  = columns of C
 *   vcols  = rows of C
 *   crows  = columns of C
 *   ccols  = rows of C
 *   v      = V(j)
 *   c      = C(i)
 *   starti = starting index i
 *   finali = ending index i
 *   startj = starting index j
 *   finalj = ending index j
 *   ndets  = number of determinants
 *   peosp  = alpha electron orbital spaces
 *   pegrps = number of alpha string orbital spaces
 *   pstr   = alpha electron strings
 *   qeosp  = beta electron orbital spaces
 *   qegrps = number of beta string orbital spaces
 *   qstr   = beta electron strings
 *   pq     = (p,q)-space pairings
 *   npq    = number of (p, q)-space pairings
 *   mo1    = 1-e integrals
 *   mo2    = 2-e integrals
 *   aelec  = alpha electrons
 *   belec  = beta  electrons
 *   intorb = internal orbitals (DOCC + CAS)
 *  Output:
 *   C      = C(i)
 */
void evaluate_hdblock_ij_1d(int pq_start_i, int pstart_i, int qstart_i,
                            int pq_final_i, int pfinal_i, int qfinal_i,
                            int pq_start_j, int pstart_j, int qstart_j,
                            int pq_final_j, int pfinal_j, int qfinal_j,
                            int vrows, int vcols, double *v,
                            int crows, int ccols, double *c,
                            int starti, int finali, int startj, int finalj,
                            int ndets, struct eospace *peosp, int pegrps,
                            struct occstr *pstr,
                            struct eospace *qeosp, int qegrps,
                            struct occstr *qstr,  int **pq,
                            int npq, double *mo1, double *mo2, int aelec,
                            int belec, int intorb);

/*
 * generate_det_triples: generate list of triplets for each determinant:
 *  |i> = |(pq, p, q)>
 */
void generate_det_triples (int ndeti, int **d_triplet, int pq_start,
                           int pstart, int qstart, int pq_final,
                           int pfinal, int qfinal, int **pq, int npq,
                           struct eospace *peosp, int pegrps,
                           struct eospace *qeosp, int qegrps);

/*
 * generate_newvector: generate new vector to be added to space.
 *  n = - r / (d - E)
 * Input:
 *  r_hndl = GA handle for residual vector
 *  d_hndl = GA handle for diagonal elements
 *  n_hndl = GA handle for new vector
 *  x_hndl = GA handle for scratch array
 */
void generate_newvector (int r_hndl, int d_hndl, double eval, int ndets,
                         int n_hndl, int x_hndl);


/*
 * generate_residual: generate residual vector.
 *  r_i = v_aj * (c_ji - e_a*v_ji)
 * Input:
 *  v_hndl = global arrays handle for vectors, V
 *  c_hndl = global arrays handle for vectors, C=Hv
 *  r_hndl = global arrays handle for residual vector
 *  hevec  = eigenvectors of v.Hv
 *  heval  = eigenvalues of v.Hv
 *  ndets  = number of determinants (length of vectors V, C, & r)
 *  ckdim  = current dimension of krylov subspace
 *  croot  = current root being optimized.
 *  rscr_hndl = GA handle for residual vector scratch array
 */
void generate_residual (int v_hndl, int c_hndl, int r_hndl, double **hevec,
                        double *heval, int ndets, int ckdim, int croot,
                        int rscr_hndl);


/*
 * get_upptri_element_index: get index of an element (i,j) in
 * list of elements in upper triangle of H matrix (n x n).
 */
int get_upptri_element_index (int i, int j, int n);

/*
 * get_upptri_element_position: get row and column indices of an n x n
 * upper triangle matrix list element.
 */
void get_upptri_element_position (int element, int n, int *i, int *j);

/*
 * get_upptri_element_rownumber: get row number of upper triangle
 * matrix list element.
 */
int get_upptri_element_rownumber (int element, int n);

/*
 * get_upptri_size: compute the size of H matrix upper triangle.
 */
int get_upptri_size (int n);

/*
 * init_diag_H_subspace: generate reference vectors from diagonalization
 * of a subspace of Hij.
 */
void init_diag_H_subspace(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets, int dim,
                          double **refspace);

/*
 * make_subspacehmat_ga build v.Hv matrix within GA toolkit
 */
void make_subspacehmat_ga (int v_hndl, int c_hndl, int ndets, int ckdim,
                           double **vhv);


/*
 * make_subspacehmat: build v.Hv matrix.
 */
void make_subspacehmat (int v_hndl, int c_hndl, int ndets, int ckdim,
                        double **vhv);

/*
 * orthonormalize_newvector: orthogonalize new vector to rest of basis.
 * Normalize result.
 * Input:
 *  v_hndl = GA handle for basis vectors
 *  nvecs  = number of basis vectors
 *  ndets  = number of determinants
 *  n_hndl = GA handle for new vector
 */
void orthonormalize_newvector (int v_hndl, int nvecs, int ndets, int n_hndl);


/*
 * peform_hv_initspace: perform Hv=c on basis vectors v_i, i = 1, .., n.
 * Input:
 *  pstr  = alpha strings
 *  peosp = alpha electron spaces
 *  pegrps= number of alpha electron spaces
 *  qstr  = beta  strings
 *  qeosp = beta  electron spaces
 *  qegrps= number of beta  electron spaces
 *  pqs   = alpha and beta space pairs
 *  num_pq= number of alpha and beta space pairs
 *  m1    = 1-e integrals
 *  m2    = 2-e integrals
 *  aelec = alpha electrons
 *  belec = beta  electrons
 *  intorb= internal orbitals
 *  ndets = number of determinants
 *  core_e= core energy
 *  dim   = number of basis vectors
 *  mdim  = maximum size of krylov space
 *  v_hndl= (GLOBAL ARRAY HANDLE) basis vectors
 *  d_hndl= (GLOBAL ARRAY HANDLE) diagonal elements <i|H|i>
 * Output:
 *  c_hndl= (GLOBAL ARRAY HANDLE) Hv=c vectors
 */
void perform_hv_initspace(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets,
                          double core_e, int dim, int mdim, int v_hndl, int d_hndl,
                          int c_hndl);

/*
 * print_subspace_eigeninfo: print diagonalization information for krylov
 * space.
 */
void print_subspace_eigeninfo(double **v, double *e, int kdim, double core_e);

/*
 * print_subspacehmat: print the krylov space hmat, v.Hv
 */
void print_subspacehmat(double **vhv, int d);

/*
 * test_convergence: test convergence of davidson algorithm. returns
 * convergence flag.
 * Input:
 *  rnorm  = norm of residual, ||r||
 *  restol = convergence tolerance
 *  croot  = current root being optimized
 *  nroot  = total number of roots
 *  f      = flag
 */
int test_convergence(double rnorm, double restol, int croot, int nroot);

/*
 * truncate_krylov_space: truncate the krylov space from krymax to krymin.
 * Input:
 *  v_hndl = GA handle of basis vectors array
 *  ndets  = number of determinants
 *  krymin = minimum dimension of krylov space
 *  krymax = maximum dimension of krylov space
 *  croot  = current root
 *  hevec  = eigenvectors
 *  x_hndl = scratch array
 */
void truncate_krylov_space(int v_hndl, int ndets, int krymin, int krymax,
                           int croot, double **hevec, int x_hndl);


#ifdef DEBUGGING
void test_new_vector_space(int v, int ckdim, int ndets, int nv);
void print_vector_space(int v, int ckdim, int ndets);
#endif

#endif
