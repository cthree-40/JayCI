// File: dysonio.h
/* Header file for Input/Ouput Dyson orbital routines. */
#ifndef dysonio_h_
#define dysonio_h_

/*
 * print_dysonorbital_info: print dyson orbital information
 */
void print_dysonorbital_info(int num_dyo, int ndyst0, int ndyst1, int norbs,
			     int *dyst0, int *dyst1);


/*
 * read_dysonorbital_info: read dyson orbital information
 */
int read_dysonorbital_info(char *filename, int *num_dyo, int *ndyst0,
			   int *ndyst1, int *norbs, int *dyst0, int *dyst1);

/*
 * read_dysonorbitals_from_file: read dyson orbitals from file
 */
int read_dysonorbitals_from_file(char *filename, int num_dyo, int norbs,
				 double **dyson);

#endif
