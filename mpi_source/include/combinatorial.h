/* File: combinatorial.h */
#ifndef combinatorial_h
#define combinatorial_h

/* Size of binomial coefficient array */
#define MAX_N 301
#define MAX_K  31

int binom_data[MAX_N][MAX_K];
int binom_init;


int binomial_coef(int m, int n);
/*
 * binomial_coef2: return binomial coefficient of n and k
 */
int binomial_coef2(int n, int k);
int fact_product(int m, int n);
/*
 * initialize_binom_coef: initialize the binomial coefficient array.
 */
void initialize_binom_coef();

int minimum(int m, int n);

#endif
