// File: cimapping.h

#ifndef cimapping_h
#define cimapping_h

/* section: section of valid matrix elements */
struct section {
	int first; /* first j of section */
	int last;  /* last j of section */
};

/* sectionll: linked list of STRUCT SECTION. */
struct sectionll {
	int first; /* first j of section */
	int last;  /* last j of section */
	struct sectionll *next;
	struct sectionll *prev;
};

/* rowmap: map of nonzero sections of row */
struct rowmap {
	int row;  /* row index */
	int nsec; /* number of sections to read in */
	int nelm; /* number of nonzero elements */
	struct section *sec; /* sections */
	struct rowmap *next; /* next row mapping */
};

/*
 * convertrow2map: convert row of bits into rowmap structure.
 */
void convertrow2map(
	int *row, /* H(i,:): 1 or 0 */
	int n2v, /* Last column needed to be evaluated */
	struct rowmap *map /* row map */
	);

/*
 * generate_cimap: generate map for row of H(i,:).
 */
int pci_generate_cimap(
	struct det *dlist, /* determinant list */
	int ndets, /* number of determinants */
	int nactv, /* number of active orbitals */
	int lwrbnd, /* least row index */
	int upprbnd, /* greatest row index */
	struct rowmap *hrow, /* map of row being evaluated */
	int nrows /* upprbnd - lwrbnd */
	);

/*
 * get_upptri_row: get row of upper triangle matrix.
 */
int get_upptri_row(
	int element, /* i of element M(i,j) */
	int matdim /* Dimension of matrix */
	);

/*
 * get_upptri_size: get size of upper triangle of square matrix.
 */
int get_upptri_size(
	int matdim /* Dimension of matrix */
	);

/*
 * upptri_rowpack_index: return index of uppertriangle element packed by rows.
 */
int upptri_rowpack_index(
	int i, /* Row number */
	int j, /* Column number */
	int matdim /* Dimension of matrix */
	);


#endif
