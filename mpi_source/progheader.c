// File: progheader.c
/* Print header for jayci.x run. */

#include <stdio.h>
#include "progheader.h"

/*
 * print_bars: prints '='.
 */
void print_bars(int num)
{
	int i;
	for (i = 0; i < num; i++) {
		fprintf(stdout, "=");
	}
	fprintf(stdout, "\n");
	return;
}

/*
 * print_jayciheader: print header for jayci.x program to stdout.
 */
void print_pjayciheader()
{
	char *program = "** pjayci.x **";
	char *name = "Christopher L Malbon (cthree-40)";
	char *group = "Yarkony Group";
	char *institution="Department of Chemistry, The Johns Hopkins University";
	
	print_bars(70);
	fprintf(stdout,"%30s\n", program);
	fprintf(stdout,"%10s\n", name);
	fprintf(stdout,"%10s\n", group);
	fprintf(stdout,"%10s\n", institution);
	print_minibars(70);
}


/*
 * print_minibars: prints '-'.
 */
void print_minibars(int num)
{
	int i;
	for (i = 0; i < num; i++) {
		fprintf(stdout, "-");
	}
	fprintf(stdout, "\n");
	return;
}
