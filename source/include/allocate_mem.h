// File: allocate_mem.h

/*
 * Allocate memory.
 */

#ifndef allocate_m_h
#define allocate_m_h

/*
 * allocate_mem_double: allocate a 2d array.
 */
int allocate_mem_double(
    double ***array, /* array */
    int n, /* rows */
    int m  /* columns */
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
