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
 *  
 */
#ifndef mathutil_h
#define mathutil_h

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

#endif
