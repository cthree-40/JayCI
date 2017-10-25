// File: buildao.h

#ifndef buildao_h
#define buildao_h

/* Pi^3 (From mathematica) */
#define M_PI_3 31.0062766802998201

/*
 * ao_basisfunc: an atomic orbital basis function
 */
struct ao_basisfunc {
        int type; /* s=1,px=2,py=3,pz=4,d2-=5,d1-=6,d0=7,
                   * d1+=8,d2+=9,f3-=10,f2-=11,f1-=12,... */
        int lval;
        int mval;
        int ugaus;      /* Uncontracted gaussians */
        double *ccoef;  /* Contraction coefficients */
        int atom;       /* Center index */
        double *geom;   /* Center geometry */
        double gscale;  /* Gaussian scaling */
};


/*
 * ao_basis: atomic orbital basis information.
 */
struct ao_basis {
        int ugaus;      /* uncontracted gaussians */
        int cgaus;      /* contracted gaussians */
        double *alpha;  /* alpha values */
        double **ccoef; /* contraction coefficients. Stored: M[cgaus][ugaus] */      
};
        
/*
 * ao_atomdata: atom data for atomic orbitals 
 */
struct ao_atomdata {
        char name[2];           /* Atom name */
        double nucchrg;         /* Nuclear charge */
        int natm;               /* Number of atoms in type. */
        int maxlval;            /* Maximum l value */
        int *naob;              /* Number of aoblocks per l value */
        double **geom;          /* atom geometries */
        int totblk;             /* total number of ao blocks */
        int *lval;              /* lvalues for each block: s = 0, p = 1, etc */
        struct ao_basis *basis; /* Basis information for each block */
};

/*
 * ao_buildao: build atomic orbitals.
 */
int ao_buildao();

/*
 * ao_check_for_daltaoin: Check for daltaoin file.
 */
int ao_check_for_daltaoin();

/*
 * ao_check_for_soinfodat: Check for soinfo.dat file.
 */
int ao_check_for_soinfodat();

/*
 * ao_get_orbitalnum: get orbital number of soinfo.dat file.
 */
int ao_get_orbitalnum(
        FILE *soinfo,
        int *norbs);

/*
 * ao_open_daltonfile: open dalton file returing file stream.
 */
FILE *ao_open_daltonfile(
        FILE *daltin); /* Dalton input file stream */

/*
 * ao_open_soinfodat: open soinfo.dat file, returning file stream.
 */
FILE *ao_open_soinfodat(
        FILE *soinfo);

/*
 * ao_read_daltontitle: read first three lines of daltaoin file.
 * INTGRL
 * X
 * X
 */
int ao_read_daltontitle(
        FILE *daltin); /* Dalton input file stream */

/*
 * ao_read_dalton_atombasis1: read first line of atom basis set iformation.
 * (Reads in: nuclear charge, number of atoms, max L value, blocks per L value
 */
int ao_read_dalton_atombasis1(
        FILE *daltin,
        struct ao_atomdata *adata);

/*
 * ao_read_dalton_atombasis2: read atomic orbital basis blocks from daltaoin
 * file into the AO_ATOMDATA.AO_ABASIS type.
 */
int ao_read_dalton_atombasis2(
        FILE *daltin,
        struct ao_atomdata *adata);

/*
 * ao_print_dalton_basisinfo1: print first line of values read from daltaoin.
 */
void ao_print_dalton_basisinfo1(
        char      *crt, /* spherical harmonics flag */
        int     atypes, /* Unique atom types */
        int     molchg, /* Molecular charge */
        char   *symtxt, /* symmetry operations test */
        double    thrs);/* integralthreshold */

/*
 * ao_print_dalton_atombasis1(adata)
 */
void ao_print_dalton_atombasis1(
        struct ao_atomdata adata);


/*
 * ao_read_dalton_basisinfo1: read first line of basis set information in
 * dalton file.
 */
int ao_read_dalton_basisinfo1(
        FILE *daltin, 
        char    *crt,
        int  *atypes,
        int  *molchg,
        char *symtxt,
        double *thrs,
        char    *id3);
/*
 * ao_read_dalton_contractionvals: read contraction values for a uncontracted
 * gaussian.
 */
int ao_read_dalton_contractionvals(
        FILE *daltin,
        struct ao_basis *basis,
        int ug);

#endif
