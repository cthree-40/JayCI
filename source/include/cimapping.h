// File: cimapping.h

#ifndef cimapping_h
#define cimapping_h

/* section: section of valid matrix elements */
struct section {
    int first; /* first j of section */
    int last;  /* final j of section */
    struct section *next; /* next section */
    struct section *prev; /* previous section */
};

struct rowmap {
    int nsec; /* number of sections to read in */
    struct section *sec; /* sections */
};

/*
 * convertrow2map: convert a row of bits into rowmap structure
 */
void convertrow2map(
    int *row, /* array of bits indicating nonzero matrix elements */
    int n2v,  /* length of row */
    struct rowmap *v2v /* rowmap structure of above array */
    );


/*
 * deallocate_cimap: deallocate a cimap linked list array.
 */
int deallocate_cimap(
	struct rowmap *hmap,
	int ndets
	);

	
/*
 * generate_cimap: generate map for evaluating nonzero matrix elements of
 * the CI Hamiltonian.
 */
int generate_cimap(
    struct det *dlist, /* determinant list */
    int ndets, /* number of determinants */
    int nactv, /* active orbitals */
    struct rowmap *hmap /* linked list of nonzero elements of H */
    );

int get_detdiffs(struct det d1, struct det d2, int nactv);

#endif
