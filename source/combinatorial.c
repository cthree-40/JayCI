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
#include "combinatorial.h"

/* subfunctions
 * ------------
 * int binomial_coef(int m, int n);
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