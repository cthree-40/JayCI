// File: dysonio.c
/*
 * Input/Output routines for Dyson orbitals.
 */
#include <stdio.h>
#include <stdlib.h>
#include "dysonio.h"

/*
 * print_dysonorbital_info: print dyson orbital information
 */
void print_dysonorbital_info(int num_dyo, int ndyst0, int ndyst1, int norbs,
			     int *dyst0, int *dyst1)
{
    int i = 0;

    printf("Number of dyson orbitals: %d\n", num_dyo);
    printf(" N+1 states: ");
    for (i = 0; i < ndyst0; i++) {
	printf("%d ", dyst0[i]);
    }
    printf("\n");
    printf(" N   states: ");
    for (i = 0; i < ndyst1; i++) {
	printf("%d ", dyst1[i]);
    }
    printf("\n");
    printf("Molecular orbitals: %d\n", norbs);

    return;
}


/*
 * read_dysonorbital_info: read dyson orbital information
 */
int read_dysonorbital_info(char *filename, int *num_dyo, int *ndyst0,
			   int *ndyst1, int *norbs, int *dyst0, int *dyst1)
{
    int error = 0;
    FILE* fptr = NULL; /* dyson orbital file */
    int i;
    
    /* open file and read:
     * 1st line: # of... dyson orbitals, anion states, neutral states
     * 2nd line: # of molecular orbitals
     * 3rd line: Anion states
     * 4th line: Neutral states
     */
    fptr = fopen(filename, "r");
    if (fptr == NULL) {
	printf("Could not open file: %s\n", filename);
	return -1;
    }
    fscanf(fptr, " %d %d %d\n", num_dyo, ndyst0, ndyst1);
    fscanf(fptr, " %d\n", norbs);
    for (i = 0; i < *ndyst0; i++) {
	fscanf(fptr, " %d", &dyst0[i]);
    }
    fscanf(fptr, "\n");
    for (i = 0; i < *ndyst1; i++) {
	fscanf(fptr, " %d", &dyst1[i]);
    }
    fscanf(fptr, "\n");
    fclose(fptr);
    return error;
}

/*
 * read_dysonorbitals_from_file: read dyson orbitals from file
 */
int read_dysonorbitals_from_file(char *filename, int num_dyo, int norbs,
				 double **dyson)
{
    int error = 0;
    char buff[255];
    FILE* fptr = NULL;
    int i, j;
    
    fptr = fopen(filename, "r");
    if (fptr == NULL) {
	return -1;
    }
    /* Skip first four lines */
    for (i = 0; i < 4; i++) {
        fgets(buff, 255, fptr);
    }
    /* Read in dyson orbitals */
    for (i = 0; i < num_dyo; i++) {
	for (j = 0; j < norbs; j++) {
	    fscanf(fptr, "%lf\n", &(dyson[i][j]));
	}
	fscanf(fptr, "\n");
    }
    return error;
}
