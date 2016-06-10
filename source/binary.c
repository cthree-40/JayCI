// File: binary.c
/********************************************************************
 * binary.c
 * --------
 * Operations for binary digit computation.
 *
 * By Christopher L Malbon
 * Dept. of Chemistry, The Johns Hopkins University
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "binary.h"

/* llint2bin: convert long long int to binary representation
 * -------------------------------------------------------------------
 * Input: 
 *  n = 64 bit integer
 *  buffer = buffer to print binary representation to
 *  buf_size = size of buffer
 *
 * *Note: LSB printed far left */
char *llint2bin(long long int n, char *buffer)
{
     int i;

     for (i = 63; i>= 0; i--) {
	  *buffer++ = (n & 1) + '0';
	  n >>= 1;
     }

     buffer[64]='\0';
     
     return buffer;
}

