// File: mathutil.h
/*
 * Subfunctions to perform basic functions. Some functions interface
 * with LAPACK.
 *
 * Subfunctions
 * ------------
 * compute_vector_norm: compute norm of a vector
 * diagmat_dsyevr:      diagonalizes a square matrix using dsyevr
 * dot_product:         computes dot product of two vectors.
 * gauss_fcn:           gaussian function evaluated at x.
 * matmul_dgemm:        matrix multiplication with dgemm.  
 * orthonormalize_vector: orthonormalize a vector to a space.
 */
#ifndef mathutil_h
#define mathutil_h

/*
 * compute_3d_distance: compute eucliden distance of two points.
 */
double compute_3d_distance(double *v, double *u);

/*
 * compute_vector_norm: compute euclidean norm of vector.
 */
double compute_vector_norm(
        double   *vector, 
        int    dimension
        );

/*
 * diagmat_dsyevr: diagonalizes a square, symmetric matrix using DSYEVR
 * LAPACK subroutine.
 */
int diagmat_dsyevr(
        double           *square_matrix,
        int         dimension_of_matrix,
        double           *eigen_vectors,
        double            *eigen_values
        );

/*
 * dot_product: compute the dot product of two vectors
 *  U.V=d
 */
double dot_product(
        double            *vector_u, 
        double            *vector_v, 
        int    dimension_of_vectors 
        );


/*
 * gauss_fcn: gaussian function. exp(-1*alpha*x^2)
 */
double gauss_fcn(
        double alpha,
        double x);

/*
 * matmul_dgemm: matrix multiply two arrays via DGEMM.
 */
int matmul_dgemm(
    double *mata,
    int row_a,
    int col_a,
    double *matb,
    int row_b,
    int col_b,
    double *matc,
    int row_c,
    int col_c
    );

/*
 * orthonormalize_vector: orthogonalize and normalize new vector to space
 * of basis vectors.
 */
void orthonormalize_vector(
    double **bvecs, /* basis vectors */
    int nbv,        /* number of basis vectors */
    int vlen,       /* length of each vector */
    double *nvec    /* new vector; vector to be normalized */
    );

/*
 * vector_difference: compute the difference between two vectors of length n.
 */
void vector_difference(
        double *v1,
        double *v2,
        int n,
        double *v3);
#endif
