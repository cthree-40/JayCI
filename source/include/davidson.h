// File: davidson.h

/* Subfunctions for davidson algorithm. */

#ifndef davidson_h
#define davidson_h

/*
 * append_vector_to_space: append new orthonormal basis vector to space.
 */
void append_vector_to_space(
    double **bvecs, /* basis vectors */
    int spsize, /* number of basis vectors in current space minus */
    int vlen, /* length of vectors */
    double *newvec /* vector being added to space */
    );

/*
 * diag_subspacehmat: diagonalize subspace matrix, returning eigenvectors
 */
int diag_subspacehmat(
    double **hmat, /* 2-d array v.Hv */
    double **hevec, /* 2-d array for eigenvectors of v.Hv */
    double *heval, /* eigenvalues of v.Hv */
    int dim, /* dimension of v.Hv */
    double *hmat_scr, /* 1-d scratch array for v.Hv */
    double *hevec_scr /* 1-d scratch array for eigenvectors of v.Hv */
    );
    
/*
 * dvdalg: perform davidson algorithm on CI space to find lowest nroots
 * roots.
 *
 * dlist = determinant list
 * ndets = number of determinants
 * moints1 = 1-e integrals
 * moints2 = 2-e integrals
 * aelec = alpha electrons
 * belec = beta electrons
 * hdgls = <i|H|i>
 * ninto = (ndocc + nactv) CAS DOCC + CAS ACTIVE orbitals
 * totgrze = nucrep + frozen core
 * maxiter = maximum davidson iterations to exectute
 * krymin = minimum size of krylov space
 * krymax = maximum size of krylov space
 * nroots = number of roots to compute
 * restol = convergence tolerance of residual vector
 * hmap = map of nonzero matrix elements <i|H|j>
 * civec = civectors; civec[nroots][ndets] double array
 * cival = eigenvalues; cival[nroots] double array
 * predr = prediagonalization routine
 * plvl = print level
 */
int dvdalg(struct det *dlist, int ndets, double *moints1, double *moints2,
	   int aelec, int belec, double *hdgls, int ninto, double totfrze,
	   int maxiter, int krymin, int krymax, int nroots, double restol,
	   struct rowmap *hmap, double **civec, double *cival, int predr,
	   int plvl);

/*
 * generate_newvector: build correction vector.
 */
void generate_newvector(
    double *rvec, /* residual vector */
    double *dgls, /* diagonals, <i|H|i> */
    double eval, /* current root's eigenvalue */
    int ndets, /* dimension of residual vector */
    double *nvec /* new vector */
    );

/*
 * generate_residual: generate residual vector for davidson algorithm.
 */
void generate_residual(
    double **vvecs, /* basis vectors v */
    double **cvecs, /* Hv=c vectors c */
    int ndets, /* length of vectors */
    int nvec, /* size of vector basis */
    double **hevec, /* eigenvectors of v.Hv */
    double *heval, /* eigenvalues of v.Hv */
    int croot, /* current root being solved for */
    double *rvec /* residual vector */
    );

/*
 * make_subspacehmat: build krylov space hamiltonian v.Hv.
 */
void make_subspacehmat(
    double **v, /* basis vectors, v */
    double **hv,/* Hv=c vectors, c */
    int ndets, /* dimension of basis vectors (number of determinants) */
    int ckdim, /* current dimension of krylov space */
    double **hmat /* krylov space hamiltonian */
    );

/*
 * orthonormalize_vector: orthogonalize and normalize new vector to space
 * of basis vectors.
 */
void orthonormalize_vector(
    double **bvecs, /* basis vectors */
    int nbv, /* number of basis vectors */
    int vlen, /* length of each vector */
    double *nvec /* new vector; vector to be normalized */
    );

/*
 * perform_hv_initspace: perform hv on inital vectors.
 */
void perform_hv_initspace(
    struct det *dlist, /* determinant list */
    int ndets, /* number of determinants */
    double *moints1, /* 1-e integrals */
    double *moints2, /* 2-e integrals */
    int aelec, /* alpha electrons */
    int belec, /* beta electrons */
    int ninto, /* docc + active orbitals */
    struct rowmap *hmap, /* map of nonzero <i|H|j> */
    int nvec, /* number of vectors to compute (krymin) */
    double **vecs, /* vectors v */
    double **hvecs /* vectors c */
    );

/*
 * print_iteration_info: print davidson iteration information.
 */
void print_iteration_info(
    double *heval, /* eigenvalues */
    int ckdim, /* current dimension of krylov space */
    int croot, /* current root being solved for */
    double rnorm, /* norm of residual vector */
    double totfrze /* nuc rep + frozen core energy contribution */
    );

/*
 * truncate_kspace: truncate the krylov space from krymax to krymin.
 */
int truncate_kspace(
    double **vscr, /* basis vectors (2-d) */
    int vlen, /* length of basis vectors */
    int krymin, /* minimum dimension of krylov space */
    int krymax, /* maximum dimension of krylov space */
    int croot, /* current root being solved for */
    double **hevec, /* krylov space eigenvectors (2-d) */
    double *vscr1d, /* 1-d scratch array for basis vectors */
    double *hevec1d, /* 1-d scratch array for eigenvectors */
    double *scr /* scratch array for new basis vectors */
    );

#endif
