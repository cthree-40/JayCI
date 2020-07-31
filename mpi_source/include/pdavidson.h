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
               int krymax, int nroots, int prediagr, int refdim, double restol,
               int ga_buffer_len, int nmos, int ndocc, int nactv);

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
                              int aelec, int belec, int intorb, int w);

/*
 * compute_cblock_H: compute values for a block from the vectors, C.
 * Input:
 *  c      = local c array
 *  ccols  = columns of c array
 *  crows  = rows of c array
 *  wi     = p, q, cas triples for c elements
 *  w_hndl = GA handle of wavefunction info
 *  v_hndl = GA handle for vectors, V
 *  d_hndl = GA handle for diagonal vectors, D
 *  buflen = length of buffer (set by user during input)
 *  pstr   = alpha strings
 *  peosp  = alpha electron occupation spaces
 *  pegrps = number of alpha electron occupation spaces
 *  qstr   = beta  strings
 *  qeosp  = beta  electron occupation spaces
 *  qegrps = number of beta  electron occupation spaces
 *  pq     = valid elec-occupation spaces of expansion
 *  npq    = number of valid elec-occupation spaces of expansion
 *  m1     = 1-e integrals
 *  m2     = 2-e integrals
 *  aelec  = alpha electrons
 *  belec  = beta  electrons
 *  intorb = internal orbitals (DOCC + ACTV)
 *  ndets  = total number of determinants
 *  nmos   = total number of molecular orbitals
 *  ndocc  = DOCC orbitals
 *  nactv  = active orbitals
 */
void compute_cblock_H(double **c, int ccols, int crows, int **wi, int w_hndl,
                      int v_hndl, int d_hndl, int buflen, struct occstr *pstr,
                      struct eospace *peosp, int pegrps,  struct occstr *qstr,
                      struct eospace *qeosp, int qegrps, int **pq, int npq,
                      double *m1, double *m2, int aelec, int belec, int intorb,
                      int ndets, int nmos, int ndocc, int nactv);

/*
 * compute_cblock_Hfast: compute values for a block from the vectors, C.
 * Input:
 *  c      = local c array
 *  ccols  = columns of c array
 *  crows  = rows of c array
 *  wi     = p, q, cas triples for c elements
 *  w_hndl = GA handle of wavefunction info
 *  v_hndl = GA handle for vectors, V
 *  d_hndl = GA handle for diagonal vectors, D
 *  buflen = length of buffer (set by user during input)
 *  pstr   = alpha strings
 *  peosp  = alpha electron occupation spaces
 *  pegrps = number of alpha electron occupation spaces
 *  qstr   = beta  strings
 *  qeosp  = beta  electron occupation spaces
 *  qegrps = number of beta  electron occupation spaces
 *  pq     = valid elec-occupation spaces of expansion
 *  npq    = number of valid elec-occupation spaces of expansion
 *  m1     = 1-e integrals
 *  m2     = 2-e integrals
 *  aelec  = alpha electrons
 *  belec  = beta  electrons
 *  intorb = internal orbitals (DOCC + ACTV)
 *  ndets  = total number of determinants
 *  nmos   = total number of molecular orbitals
 */
void compute_cblock_Hfast(double **c, int ccols, int crows, int **wi, int w_hndl,
                          int v_hndl, int d_hndl, int buflen,struct occstr *pstr,
                          struct eospace *peosp, int pegrps, struct occstr *qstr,
                          struct eospace *qeosp, int qegrps, int **pq, int npq,
                          double *m1, double *m2, int aelec,int belec,int intorb,
                          int ndets, int nmos, int ndocc, int nactv, int c_hndl,
                          int cstep);

/*
 * compute_cblock_Hfaster: compute values for a block from the vectors, C.
 * Input:
 *  c      = local c array
 *  ccols  = columns of c array
 *  crows  = rows of c array
 *  wi     = p, q, cas triples for c elements
 *  w_hndl = GA handle of wavefunction info
 *  v_hndl = GA handle for vectors, V
 *  d_hndl = GA handle for diagonal vectors, D
 *  buflen = length of buffer (set by user during input)
 *  pstr   = alpha strings
 *  peosp  = alpha electron occupation spaces
 *  pegrps = number of alpha electron occupation spaces
 *  qstr   = beta  strings
 *  qeosp  = beta  electron occupation spaces
 *  qegrps = number of beta  electron occupation spaces
 *  pq     = valid elec-occupation spaces of expansion
 *  npq    = number of valid elec-occupation spaces of expansion
 *  m1     = 1-e integrals
 *  m2     = 2-e integrals
 *  aelec  = alpha electrons
 *  belec  = beta  electrons
 *  intorb = internal orbitals (DOCC + ACTV)
 *  ndets  = total number of determinants
 *  nmos   = total number of molecular orbitals
 */
void compute_cblock_Hfaster(double *c, int ccols, int crows, int **wi, int w_hndl,
                          int v_hndl, int d_hndl, int buflen,struct occstr *pstr,
                          struct eospace *peosp, int pegrps, struct occstr *qstr,
                          struct eospace *qeosp, int qegrps, int **pq, int npq,
                          double *m1, double *m2, int aelec,int belec,int intorb,
                          int ndets, int nmos, int ndocc, int nactv, int c_hndl,
                            int cstep, int *colnums);

/*
 * compute_cblock_Hfastest: compute values for a block from the vectors, C.
 * Computes upper diagonal.
 *
 *  H(i,j)*V(j,k)=C(i,k)
 *        *V(i,k)=C(j,k)
 *
 * Input:
 *  c      = local c array
 *  ccols  = columns of c array
 *  crows  = rows of c array
 *  wi     = p, q, cas triples for c elements
 *  w_hndl = GA handle of wavefunction info
 *  v_hndl = GA handle for vectors, V
 *  d_hndl = GA handle for diagonal vectors, D
 *  buflen = length of buffer (set by user during input)
 *  pstr   = alpha strings
 *  peosp  = alpha electron occupation spaces
 *  pegrps = number of alpha electron occupation spaces
 *  qstr   = beta  strings
 *  qeosp  = beta  electron occupation spaces
 *  qegrps = number of beta  electron occupation spaces
 *  pq     = valid elec-occupation spaces of expansion
 *  npq    = number of valid elec-occupation spaces of expansion
 *  m1     = 1-e integrals
 *  m2     = 2-e integrals
 *  aelec  = alpha electrons
 *  belec  = beta  electrons
 *  intorb = internal orbitals (DOCC + ACTV)
 *  ndets  = total number of determinants
 *  nmos   = total number of molecular orbitals
 *  ndocc  = number of docc orbitals
 *  nactv  = number of active orbitals
 *  c_hndl = GLOBAL ARRAY handle for Hv=c vectors
 *  cstep  = first row index in block
 *  cmax   = last  row index in block
 *  colnums= indices of C_i to evaluate Hv_i=c_i 
 */
void compute_cblock_Hfastest(double *c1d, int ccols, int crows, int **wi, int w_hndl,
                             int v_hndl, int d_hndl, int buflen,struct occstr *pstr,
                             struct eospace *peosp, int pegrps, struct occstr *qstr,
                             struct eospace *qeosp, int qegrps, int **pq, int npq,
                             double *m1, double *m2, int aelec,int belec,int intorb,
                             int ndets, int nmos, int ndocc, int nactv, int c_hndl,
                             int cstep, int cmax, int *colnums);

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
 * compute_diagonal_iHi: compute <i|H|i> elements.
 * Input:
 *  hdgls   = digonal elements <i|H|i> for (i = start,...,final)
 *  start   = starting determinant index
 *  final   = final determinant index
 *  mo1     = 1-e integrals
 *  mo2     = 2-e integrals
 *  aelec   = CI alpha electrons
 *  belec   = CI beta  electrons
 *  intorb  = internal orbitals (DOCC + CAS)
 *  w_hndl  = global array handle for wavefunction
 */
void compute_diagonal_iHi(double *hdgls, int start, int final,
			  double *mo1, double *mo2, int aelec,
			  int belec, int intorb, int w_hndl,
			  struct occstr *pstr, struct occstr *qstr);


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
                          int belec, int intorb, int ndets, int kmax, int w_hndl,
                          int ga_buffer_len);

/*
 * compute_hv_newvectorfaster: compute Hv=c for newest vector in basis space.
 * Input:
 *  v_hndl = GA handle for basis vectors (ckdim is new vector)
 *  c_hndl = GA handle for Hv=c vectors
 *  ckdim  = current dimension of space
 *  w_hndl = wavefunction list (deteriminant triplets)
 */
void compute_hv_newvectorfaster(struct occstr *pstr, struct eospace *peosp, int pegrps,
                                struct occstr *qstr, struct eospace *qeosp, int qegrps,
                                int **pqs, int num_pq, double *m1, double *m2,
                                int aelec, int belec, int intorb, int ndets,
                                double core_e, int ckdim, int mdim, int v_hndl,
                                int d_hndl, int c_hndl, int w_hndl, int ga_buffer_len,
                                int nmo, int ndocc, int nactv);

/*
 * compute_hij_eosp: compute hij for an electron-occupation space.
 * Uses OpenMP.
 * Compute only upper triangle.
 */
void compute_hij_eosp(double *ci, int ccols, int crows, int **wi,
                      struct occstr *pstr, struct eospace *peosp, int pegrps,
                      struct occstr *qstr, struct eospace *qeosp, int qegrps,
                      int **pq, int npq, double *m1, double *m2, int aelec,
                      int belec, int intorb,
                      int nmos, int ndocc, int nactv, int cstep, int *cnums,
                      int jstart, int jmax, int jstartp, int jstartq,
                      int jfinalp, int jfinalq, int *jpair, double *vj,
                      double *vi, double *cj);

/*
 * compute_hvc_diagonal_ga: compute <i|H|i>*v(i,j)=c(i,j) using global arrays.
 * Subscript 1 is column. Subscript 2 is row.
 */
void compute_hvc_diagonal_ga(int c_hndl, int v_hndl, int d_hndl, int start,
                             int final, int ndets);

/*
 * compute_hvnewvecfast: perform Hv=c on basis vector v_n+1.
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
 *  w_hndl= (GLOBAL ARRAY HANDLE) wavefunction
 * Output:
 *  c_hndl= (GLOBAL ARRAY HANDLE) Hv=c vectors
 */
void compute_hvnewvecfast(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets,
                          double core_e, int dim, int mdim, int v_hndl,
                          int c_hndl, int w_hndl, int d_hndl, int ga_buffer_len,
                          int nmo, int ndocc, int nactv);

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
void evaluate_hdblock_ij_1d2(int **wi, int idets, int **wj, int jdets,
                             int vrows, int vcols, double *v,
                             int crows, int ccols, double *c,
                             int ndets, struct eospace *peosp, int pegrps,
                             struct occstr *pstr,
                             struct eospace *qeosp, int qegrps,
                             struct occstr *qstr,  int **pq,
                             int npq, double *mo1, double *mo2, int aelec,
                             int belec, int intorb);

/*
 * evaluate_hij_pxlist1x: evaluate hij for single replacements in alpha strings.
 */
void evaluate_hij_pxlist1x(struct det deti, struct xstr *pxlist, int npx,
                           int qindx,
                           int nqx, struct occstr *pstr, struct eospace *peosp,
                           int npe, struct occstr *qstr, struct eospace *qeosp,
                           int nqe, int **pq, int npq, double *m1, double *m2,
                           int aelec, int belec, int intorb, double *c,
                           int vrows, int vcols, int **vindx, int **windx,
                           int *jindx, double **v, double *v1d, int **w,
                           int *w1d, double *hijval, int w_hndl, int v_hndl);
/*
 * evaluate_hij_pxlist1x_ut: evaluate hij for single replacements in alpha
 * strings. Only upper triangle is computed.
 */
void evaluate_hij_pxlist1x_ut(struct det deti, struct xstr *pxlist, int npx,
                              int qindx,
                              int nqx, struct occstr *pstr, struct eospace *peosp,
                              int npe, struct occstr *qstr, struct eospace *qeosp,
                              int nqe, int **pq, int npq, double *m1, double *m2,
                              int aelec, int belec, int intorb, double *c,
                              int vrows, int vcols, int **vindx, int **windx,
                              int *jindx, double **v, double *v1d, int **w,
                              int *w1d, double *hijval, int w_hndl, int v_hndl,
                              int c_hndl, int cindx, double *vik, double *cjk,
                              int **vx2, int *cnums);

/*
 * evaluate_hij_pxlist1x_ut2: evaluate hij for single replacements in alpha
 * strings.
 * Input:
 *  deti   = <i| determinant
 *  pxlist = list of single replacements
 *  npx    = number of single replacments
 *  qindx  = index of beta string of |j>
 *  nqx    = always 1
 *  pstr   = alpha occupation strings
 *  peosp  = alpha electron occupation spaces
 *  npe    = number of alpha electron occupation spaces
 *  qstr   = beta  occupation strings
 *  qeosp  = beta  electron occupation spaces
 *  nqe    = number of beta  electron occupation spaces
 *  pq     = p,q space pairings
 *  npq    = number of p,q space pairings
 *  m1     = 1-e integrals
 *  m2     = 2-e integrals
 *  aelec  = alpha electrons
 *  belec  = beta  electrons
 *  intorb = internal orbitals
 *  vrows  = number of rows j
 *  vcols  = number of columns k
 *  jstep  = index in wavefunction of first determinant in this buffer j
 *  cik    = C(i,k)
 *  vjk    = V(j,k)
 *  vik    = V(i,k)
 *  cjk    = C(j,k)
 *  hijval = <i|H|j> values
 *  jindx  = array for determinant indices
 */
void evaluate_hij_pxlist1x_ut2(struct det deti, struct xstr *pxlist, int npx,
                               int qindx, int nqx,
                               struct occstr *pstr, struct eospace *peosp, int npe,
                               struct occstr *qstr, struct eospace *qeosp, int nqe,
                               int **pq, int npq, double *m1, double *m2, int aelec,
                               int belec, int intorb, int vrows, int vcols,
                               int jstep, double *cik, double *vjk, double *vik,
                               double *cjk, double *hijval, int *jindx);

/*
 * evaluate_hij_pxlist2x: evaluate hij for double replacements in alpha strings.
 */
void evaluate_hij_pxlist2x(struct det deti, struct xstr *pxlist, int npx,
                           int qindx,
                           int nqx, struct occstr *pstr, struct eospace *peosp,
                           int npe, struct occstr *qstr, struct eospace *qeosp,
                           int nqe, int **pq, int npq, double *m1, double *m2,
                           int aelec, int belec, int intorb, double *c,
                           int vrows, int vcols, int **vindx, int **windx,
                           int *jindx, double **v, double *v1d, int **w,
                           int *w1d, double *hijval, int w_hndl, int v_hndl);

/*
 * evaluate_hij_pxlist2x_ut: evaluate hij for double replacements in alpha
 * strings. Only upper triangle is computed.
 */
void evaluate_hij_pxlist2x_ut(struct det deti, struct xstr *pxlist, int npx,
                              int qindx,
                              int nqx, struct occstr *pstr, struct eospace *peosp,
                              int npe, struct occstr *qstr, struct eospace *qeosp,
                              int nqe, int **pq, int npq, double *m1, double *m2,
                              int aelec, int belec, int intorb, double *c,
                              int vrows, int vcols, int **vindx, int **windx,
                              int *jindx, double **v, double *v1d, int **w,
                              int *w1d, double *hijval, int w_hndl, int v_hndl,
                              int c_hndl, int cindx, double *vik, double *cjk,
                              int **vx2, int *cnums);

/*
 * evaluate_hij_pxqxlist2x: evaluate hij for single replacements in alpha
 * and beta strings.
 */
void evaluate_hij_pxqxlist2x(struct det deti, struct xstr *pxlist, int npx,
                             struct xstr *qxlist,
                             int nqx, struct occstr *pstr, struct eospace *peosp,
                             int npe, struct occstr *qstr, struct eospace *qeosp,
                             int nqe, int **pq, int npq, double *m1, double *m2,
                             int aelec, int belec, int intorb, double *c,
                             int vrows, int vcols, int **vindx, int **windx,
                             int *jindx, double **v, double *v1d, int **w,
                             int *w1d, double *hijval, int w_hndl, int v_hndl);

/*
 * evaluate_hij_pxqxlist2x_ut: evaluate hij for single replacements in alpha
 * and beta strings. Upper triangle only.
 */
void evaluate_hij_pxqxlist2x_ut(struct det deti, struct xstr *pxlist, int npx,
                                struct xstr *qxlist,
                                int nqx, struct occstr *pstr, struct eospace *peosp,
                                int npe, struct occstr *qstr, struct eospace *qeosp,
                                int nqe, int **pq, int npq, double *m1, double *m2,
                                int aelec, int belec, int intorb, double *c,
                                int vrows, int vcols, int **vindx, int **windx,
                                int *jindx, double **v, double *v1d, int **w,
                                int *w1d, double *hijval, int w_hndl, int v_hndl,
                                int c_hndl, int cindx, double *vik, double *cjk,
                                int **vx2, int *cnums);

/*
 * evaluate_hij_pxlist2x_ut2: evaluate hij for double replacements in alpha
 * strings. Only upper triangle is computed.
 */
void evaluate_hij_pxlist2x_ut2(struct det deti, struct xstr *pxlist, int npx,
                               int qindx, int nqx,
                               struct occstr *pstr, struct eospace *peosp, int npe,
                               struct occstr *qstr, struct eospace *qeosp, int nqe,
                               int **pq, int npq, double *m1, double *m2, int aelec,
                               int belec, int intorb, int vrows, int vcols,
                               int jstep, double *cik, double *vjk, double *vik,
                               double *cjk, double *hijval, int *jindx);

/*
 * evaluate_hij_pxqxlist2x_ut2: evaluate hij for single replacements in alpha
 * and beta strings. Upper triangle only.
 */
void evaluate_hij_pxqxlist2x_ut2(struct det deti, struct xstr *pxlist, int npx,
                                 struct xstr *qxlist, int nqx,
                                 struct occstr *pstr, struct eospace *peosp, int npe,
                                 struct occstr *qstr, struct eospace *qeosp, int nqe,
                                 int **pq, int npq, double *m1, double *m2, int aelec,
                                 int belec, int intorb, int vrows, int vcols,
                                 int jstep, double *cik, double *vjk, double *vik,
                                 double *cjk, double *hijval, int *jindx);
/*
 * evaluate_hij_qxlist1x: evaluate hij for single replacements in alpha strings.
 */
void evaluate_hij_qxlist1x(struct det deti, int pindx, int npx,
                           struct xstr *qxlist,
                           int nqx, struct occstr *pstr, struct eospace *peosp,
                           int npe, struct occstr *qstr, struct eospace *qeosp,
                           int nqe, int **pq, int npq, double *m1, double *m2,
                           int aelec, int belec, int intorb, double *c,
                           int vrows, int vcols, int **vindx, int **windx,
                           int *jindx, double **v, double *v1d, int **w,
                           int *w1d, double *hijval, int w_hndl, int v_hndl);

/*
 * evaluate_hij_qxlist1x_ut: evaluate hij for single replacements in alpha
 * strings. Only upper triangle is computed.
 */
void evaluate_hij_qxlist1x_ut(struct det deti, int pindx, int npx,
                              struct xstr *qxlist,
                              int nqx, struct occstr *pstr, struct eospace *peosp,
                              int npe, struct occstr *qstr, struct eospace *qeosp,
                              int nqe, int **pq, int npq, double *m1, double *m2,
                              int aelec, int belec, int intorb, double *c,
                              int vrows, int vcols, int **vindx, int **windx,
                              int *jindx, double **v, double *v1d, int **w,
                              int *w1d, double *hijval, int w_hndl, int v_hndl,
                              int c_hndl, int cindx, double *vik, double *cjk,
                              int **vx2, int *cnums);

/*
 * evaluate_hij_qxlist1x_ut2: evaluate hij for single replacements in beta
 * strings.
 * Input:
 *  deti   = <i| determinant
 *  pxlist = list of single replacements
 *  npx    = number of single replacments
 *  qindx  = index of beta string of |j>
 *  nqx    = always 1
 *  pstr   = alpha occupation strings
 *  peosp  = alpha electron occupation spaces
 *  npe    = number of alpha electron occupation spaces
 *  qstr   = beta  occupation strings
 *  qeosp  = beta  electron occupation spaces
 *  nqe    = number of beta  electron occupation spaces
 *  pq     = p,q space pairings
 *  npq    = number of p,q space pairings
 *  m1     = 1-e integrals
 *  m2     = 2-e integrals
 *  aelec  = alpha electrons
 *  belec  = beta  electrons
 *  intorb = internal orbitals
 *  vrows  = number of rows j
 *  vcols  = number of columns k
 *  jstep  = index in wavefunction of first determinant in this buffer j
 *  cik    = C(i,k)
 *  vjk    = V(j,k)
 *  vik    = V(i,k)
 *  cjk    = C(j,k)
 *  hijval = <i|H|j> values
 *  jindx  = array for determinant indices
 */
void evaluate_hij_qxlist1x_ut2(struct det deti, int pindx, int npx,
                               struct xstr *qxlist, int nqx,
                               struct occstr *pstr, struct eospace *peosp, int npe,
                               struct occstr *qstr, struct eospace *qeosp, int nqe,
                               int **pq, int npq, double *m1, double *m2, int aelec,
                               int belec, int intorb, int vrows, int vcols,
                               int jstep, double *cik, double *vjk, double *vik,
                               double *cjk, double *hijval, int *jindx);

/*
 * evaluate_hij_qxlist2x: evaluate hij for double replacements in beta strings.
 */
void evaluate_hij_qxlist2x(struct det deti, int pindx, int npx,
                           struct xstr *qxlist,
                           int nqx, struct occstr *pstr, struct eospace *peosp,
                           int npe, struct occstr *qstr, struct eospace *qeosp,
                           int nqe, int **pq, int npq, double *m1, double *m2,
                           int aelec, int belec, int intorb, double *c,
                           int vrows, int vcols, int **vindx, int **windx,
                           int *jindx, double **v, double *v1d, int **w,
                           int *w1d, double *hijval, int w_hndl, int v_hndl);

/*
 * evaluate_hij_qxlist2x_ut: evaluate hij for double replacements in alpha
 * strings. Only upper triangle is computed.
 */
void evaluate_hij_qxlist2x_ut(struct det deti, int pindx, int npx,
                              struct xstr *qxlist,
                              int nqx, struct occstr *pstr, struct eospace *peosp,
                              int npe, struct occstr *qstr, struct eospace *qeosp,
                              int nqe, int **pq, int npq, double *m1, double *m2,
                              int aelec, int belec, int intorb, double *c,
                              int vrows, int vcols, int **vindx, int **windx,
                              int *jindx, double **v, double *v1d, int **w,
                              int *w1d, double *hijval, int w_hndl, int v_hndl,
                              int c_hndl, int cindx, double *vik, double *cjk,
                              int **vx2, int *cnums);

/*
 * evaluate_hij_qxlist2x_ut2: evaluate hij for double replacements in alpha
 * strings. Only upper triangle is computed.
 */
void evaluate_hij_qxlist2x_ut2(struct det deti, int pindx, int npx,
                               struct xstr *qxlist, int nqx,
                               struct occstr *pstr, struct eospace *peosp, int npe,
                               struct occstr *qstr, struct eospace *qeosp, int nqe,
                               int **pq, int npq, double *m1, double *m2, int aelec,
                               int belec, int intorb, int vrows, int vcols,
                               int jstep, double *cik, double *vjk, double *vik,
                               double *cjk, double *hijval, int *jindx);


/*
 * evaluate_hij_jindx: evaluate hij given a determinant |i> and a list
 * of excitations: |r,s> = |p',q>, |p,q'>, |p',q'>, |p",q>, |p,q">
 */
void evaluate_hij_pxqxlist(struct det deti, struct xstr *pxlist, int npx,
                           struct xstr *qxlist,
                           int nqx, struct occstr *pstr, struct eospace *peosp,
                           int npe, struct occstr *qstr, struct eospace *qeosp,
                           int nqe, int **pq, int npq, double *m1, double *m2,
                           int aelec, int belec, int intorb, double *c,
                           int vrows, int vcols, int **vindx, int **windx,
                           int *jindx, double **v, double *v1d,int **w,
                           int *w1d, double *hijval, int w_hndl, int v_hndl);

/*
 * evaluate_hij_jindx_1d: evaluate hij given a determinant |i> and a list
 * of excitations: |r,s> = |p',q>, |p,q'>, |p',q'>, |p",q>, |p,q">, for
 * the vector V_i to make the vector C_i.
 */
void evaluate_hij_pxqxlist_1d(struct det deti, int *pxlist, int npx, int *qxlist,
                              int nqx, struct occstr *pstr, struct eospace *peosp,
                              int npe, struct occstr *qstr, struct eospace *qeosp,
                              int nqe, int **pq, int npq, double *m1, double *m2,
                              int aelec, int belec, int intorb, double *c,
                              int vrows, int vvecx, int **vindx, int **windx,
                              int *jindx, double *v, int **w,
                              int *w1d, double *hijval, int w_hndl, int v_hndl);

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
 * generate_wlist: generate the wavefunction list of triplets
 * Input:
 *  hndl     = handle of global array of W
 *  ndets    = number of determinants
 *  pq       = alpha/beta pairs
 *  npq      = number of alpha/beta pairs
 *  peospace = alpha string electron orbital space
 *  pegrps   = number of alpha string orbital spaces
 *  qeospace = beta  string electron orbital space
 *  qegrps   = number of beta  string orbital spaces
 */
void generate_wlist(int hndl, int ndets, int **pq, int npq,
                    struct eospace *peosp, int pegrps,
                    struct eospace *qeosp, int qegrps);


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
int get_upptri_element_rownumber (long long int element, int n);

/*
 * get_upptri_size: compute the size of H matrix upper triangle.
 */
long long int get_upptri_size (int n);

/*
 * init_diag_H_subspace: generate reference vectors from diagonalization
 * of a subspace of Hij.
 */
void init_diag_H_subspace(int w_hndl,
                          struct occstr *pstr,// struct eospace *peosp, int pegrps,
                          struct occstr *qstr,// struct eospace *qeosp, int qegrps,
                          double *m1, double *m2,
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
 *  w_hndl= (GLOBAL ARRAY HANDLE) wavefunction
 * Output:
 *  c_hndl= (GLOBAL ARRAY HANDLE) Hv=c vectors
 */
void perform_hv_initspace(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets,
                          double core_e, int dim, int mdim, int v_hndl, int d_hndl,
                          int c_hndl, int w_hndl, int ga_buffer_len);

/*
 * peform_hvispacefast: perform Hv=c on basis vectors v_i, i = 1, .., n.
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
 *  w_hndl= (GLOBAL ARRAY HANDLE) wavefunction
 * Output:
 *  c_hndl= (GLOBAL ARRAY HANDLE) Hv=c vectors
 */
void perform_hvispacefast(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets,
                          double core_e, int dim, int mdim, int v_hndl,
                          int d_hndl, int c_hndl, int w_hndl, int ga_buffer_len,
                          int nmo, int ndocc, int nactv);

void perform_hvispacefast_debug(struct occstr *pstr, struct eospace *peosp, int pegrps,
                          struct occstr *qstr, struct eospace *qeosp, int qegrps,
                          int **pqs, int num_pq, double *m1, double *m2,
                          int aelec, int belec, int intorb, int ndets,
                          double core_e, int dim, int mdim, int v_hndl,
                          int d_hndl, int c_hndl, int w_hndl, int ga_buffer_len,
                                int nmo, int ndocc, int nactv);

/*
 * print_iter_info: print iteration information.
 */
void print_iter_info(double *heval, int ckdim, int croot, double rnorm,
                     double totfrze);

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
 * set_ga_det_indexes: set the array of indices to gather from global array.
 * Input:
 *  jindx = list of row numbers in vector V
 *  num   = number of row numbers
 *  cols  = number of columns of V
 * Output:
 *  vindx = indices of global array V to gather
 */
void set_ga_det_indexes(int *jindx, int num, int cols, int **vindx);

/*
 * set_ga_det_indexes: set the array of indices to gather from global array.
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

/*
 * set_ga_det_indexes_1D: set the array of indices to gather from global array.
 * Input:
 *  jindx = list of row numbers in vector V
 *  num   = number of row numbers
 *  icol  = column of V
 * Output:
 *  vindx = indices of global array V to gather
 */
void set_ga_det_indexes_1D(int *jindx, int num, int icol, int **vindx);

/*
 * set_ga_det_indexes_trans: set the array of indices to gather from
 * global array. Transpose of above.
 * Input:
 *  jindx = list of row numbers in vector V
 *  num   = number of row numbers
 *  cols  = number of columns of V
 * Output:
 *  vindx = indices of global array V to gather
 */
void set_ga_det_indexes_trans(int *jindx, int num, int cols, int **vindx);

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
