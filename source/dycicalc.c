/* File: dycicalc.c
 *
 * Main driver for dycicalc program.
 * Computes dyson orbital between N+1 and N electron CI wavefunctions.
 *
 * By Christopher Malbon
 * Yarkony Group, Department of Chemistry
 * The Johns Hopkins University
 */
#include <stdio.h>
#include <stdlib.h>
#include "progheader.h"
#include "run_dycicalc.h"
/*
 * Main routine.
 */
int main(void)
{
        int error = 0;
        
        print_dyciheader();
        error = run_dycicalc();
        return error;
}

        
