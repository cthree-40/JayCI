// File: allocate_mem.h

/*
 * Allocate memory.
 */

#ifndef allocate_m_h
#define allocate_m_h

/*
 * allocate_mem_char: allocated a 2d character array.
 */
int allocate_mem_char(
        char ***array,
        int n,
        int m
        );

/*
 * allocate_mem_double: allocate a 2d array.
 */
int allocate_mem_double(
    double ***array, /* array */
    int n, /* rows */
    int m  /* columns */
    );
/*
 * allocate_mem_int: allocate a 2d integer array
 * Array = M[n x m]
 */
int allocate_mem_int(
        int ***array,
        int n,
        int m
        );

/*
 * deallocate_mem: deallocate 2d array
 */
int deallocate_mem(
    double ***array, /* array */
    int n, /* rows */
    int m  /* columns */
    );

#endif
