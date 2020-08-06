// File: combinatorial.c
/********************************************************************* 
 * combinatorial
 * ------------
 * Contains subfunctions for performing binomial coefficients
 *
 * By Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "iminmax.h"
#include "combinatorial.h"

/* Size of binomial coefficient array */
#define MAX_N 301
#define MAX_K  31

int binom_data[MAX_N][MAX_K];
int binom_init;

/* subfunctions
 * ------------
 * int binomial_coef(int m, int n);
 * int binomial_coef2(int m, int n);
 * int fact_product(int m, int n);
 * int minimum(int m, int n);
 */

int binomial_coef(int n, int k)
/* binomial_coef
 * -------------
 * Returns binomial coefficient of n and k. */
{
     int* c = (int*)calloc(k + 1, sizeof(int));
     int i, j, result;
     
     c[0] = 1;
     
     for (i = 1; i <= n; i++) {
	  for (j = minimum(i, k); j > 0; j--)
	       c[j] = c[j] + c[j - 1];
     }
     
     result = c[k];
     
     free(c);
     
     return result;
}

/*
 * binomial_coef2: return binomial coefficient of n and k
 */
int binomial_coef2(int n, int k)
{
        int result;
        if (binom_init == 0) {
                initialize_binom_coef();
                binom_init++;
        }
        result = binom_data[n][k];
        return result;
};

/*
 * binomial_coef3: return binomial coefficient of n and k.
 * WARNING: Does not check for initialization of binom_data[][].
 */
int binomial_coef3(int n, int k)
{
        int result;
        result = binom_data[n][k];
        return result;
};


int fact_product(int m, int n)
/* fact_product
 * ------------
 * Computes the product: m * (m+1) * ... * (n - 1) * (n) */
{
     int i, result;
  
     if (n < 0) {
	  result = 0;
     }
     else if (n == 0) {
	  result = 1;
     }
     else {
	  result = m;
	  for (i = m + 1; i <= n; i++) {
	       result = result * i;
	  }
     }
     
  return result;
}

/*
 * initialize_binom_coef: initialize the binomial coefficient array.
 */
void initialize_binom_coef()
{
        int i, j;
        for (i = 0; i <= MAX_N; i++) {
                for (j = 0; j <= int_min(i, MAX_K); j++) {
                        if (j == 0 || j == i) {
                                binom_data[i][j] = 1;
                        } else {
                                binom_data[i][j] = binom_data[i - 1][j - 1] +
                                        binom_data[i - 1][j];
                        }
                }
        }
        return;
}

int minimum(int m, int n)
/* minimum
 * -------
 * Returns minimum of m and n */
{
     int result;
     
     if (m <= n) {
	  result = m;
     }
     else{
	  result = n;
     }
     
     return result;
}
