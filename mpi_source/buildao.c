// File: buildao.c
/*
 * Build atomic orbitals
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "errorlib.h"
#include "allocate_mem.h"
#include "arrayutil.h"
#include "ioutil.h"
#include "mathutil.h"
#include "buildao.h"

/*
 * ao_buildao: build atomic orbitals
 */
int ao_buildao(struct ao_basisfunc *aobasis)
{
    int error = 0;
    int natoms = 0;       /* Number of atoms */
    FILE *daltin = NULL;  /* Dalton input file */
    char crt[2]; /* Spherical or cartesian functions test*/
    int atypes = 0; /* Atom types */
    int molchg = 0; /* Molecular charge */
    char symtxt[2]; /* sym operations test */
    double thrs = 0.0; /* Integral threshold */
    char id3[2];
    int norbs = 0; /* Numer of orbitals */
    struct ao_atomdata *adata = NULL; /* Atom basis data */
    //struct ao_basisfunc *aobasis = NULL; /* Atomic orbital basis */
    int i = 0;
    
    /* Check for daltaoin file. Open the file and read the first
     * three lines. The ao_read_dalton_title subroutine ensures
     * the file is of the correct form. */
    error = ao_check_for_daltaoin();
    if (error != 0) return error;
    daltin = ao_open_daltonfile(daltin);
    if (daltin == NULL) {
        error++;
        return error;
    }
    error = ao_read_daltontitle(daltin);
    if (error != 0) return error;
    
    /* Read the first line of basis set information. */
    error = ao_read_dalton_basisinfo1(daltin, crt, &atypes, &molchg,
                                      symtxt, &thrs, id3);
    if (error != 0) return error;
    ao_print_dalton_basisinfo1(crt, atypes, molchg, symtxt, thrs);
    
    
    /* Allocate atomdata structure type. Read in basis set information
     * for each atom type. */
    adata = (struct ao_atomdata *)
        malloc(sizeof(struct ao_atomdata) * atypes);
    for (i = 0; i < atypes; i++) {
        error = ao_read_dalton_atombasis1(daltin, &adata[i]);
        if (error != 0) return error;
        ao_print_dalton_atombasis1(adata[i]);
        error = ao_read_dalton_atombasis2(daltin, &adata[i]);
    }
    
    /* Process the basis set information and construct the atomic
     * orbital basis. */
    error = ao_build_atomic_orbitalset(adata, aobasis, atypes, &norbs);
    return error;
}

/*
 * ao_build_atomic_orbitalset: build atomic orbital basis set from adata.
 * This routine will build the array (struct ao_basisfunc) aobasis. AOBASIS
 * will hold l, m, contraction information, and geometry information. This
 * routine requires the soinfo.dat generated by DALTON.
 */
int ao_build_atomic_orbitalset(struct ao_atomdata *adata,
                               struct ao_basisfunc *aobasis,
                               int atypes, int *norbs)
{
    int error = 0;        /* Error flag */
    FILE *soinfo  = NULL; /* soinfo.dat file pointer */
    double *gscal = NULL; /* Gaussian scalings */
    int ngaus = 0;        /* Number of gaussians to scale */
    char **ainfo = NULL;  /* Atom info array [natom][MAX_LINE_SIZE] */
    int natoms = 0;       /* Number of atoms */
    int i = 0;
    
    /* Check for soinfo.dat file. Get number of atomic orbitals for
     * the molecule from the soinfo.dat file. Allocate the aobasis
     * array to hold each atomic orbital. */
    error = ao_check_for_soinfodat();
    if (error != 0) return error;
    soinfo = ao_open_soinfodat(soinfo);
    if (soinfo == NULL) {
        error++;
        return error;
    }
    error = ao_get_orbitalnum(soinfo, norbs);
    if (error != 0) return error;
    printf("Atomic orbital number: %d\n", *norbs);
//        aobasis = (struct ao_basisfunc *)
//                malloc(*norbs * sizeof(struct ao_basisfunc));
//        for (i = 0; i < *norbs; i++) {
//                aobasis[i].geom = (double *) malloc(3 * sizeof(double));
//        }
    /* Get gaussian scaling information */
    error = ao_get_gscalings(soinfo, &gscal, &ngaus);
    printf("Gaussian scalings:\n");
    for (i = 0; i < ngaus; i++) {
        printf(" %4d %16.10lf\n", (i + 1), gscal[i]);
    }
    /* Get atominfo. This will be used to assign basis functions */
    error = ao_get_atominfo(soinfo, &ainfo, &natoms);
    for (i = 0; i < natoms; i++) {
        printf("%s", ainfo[i]);
    }
    error = ao_process_orbitaldata(adata, aobasis, soinfo, *norbs, gscal,
                                   ngaus, ainfo, atypes);
    for (i = 0; i < *norbs; i++) {
        printf(" --- Orbital #%3d ---\n", (i + 1));
        ao_print_orbital_information(aobasis[i]);
    }
    
    return error;
}

/*
 * ao_check_for_daltaoin: Check for daltaoin file.
 */
int ao_check_for_daltaoin()
{
    int error = 0;
    error = check_for_file("daltaoin","r");
    return error;
}

/*
 * ao_check_for_soinfodat: Check for soinfo.dat file.
 */
int ao_check_for_soinfodat()
{
    int error = 0;
    error = check_for_file("soinfo.dat","r");
    return error;
}

/*
 * ao_cp_alphavals_ccoefs: copy alpha values and contraction coefficients for
 * ao block and contracted gaussian.
 */
void ao_cp_alphavals_ccoefs(int block, int gauss, struct ao_atomdata adata,
                            struct ao_basisfunc *aobasis)
{
    /* Copy alpha values */
    cparray_1d1d(adata.basis[block].alpha, adata.basis[block].ugaus,
                 aobasis->alpha, aobasis->ugaus);
    cparray_1d1d(adata.basis[block].ccoef[gauss], aobasis->ugaus,
                 aobasis->ccoef, aobasis->ugaus);
    return;
}

/*
 * ao_cp_gaussian_scaling: copy gaussian scaling constant into new orbital.
 */
void ao_cp_gaussian_scaling(struct ao_basisfunc *bas, double *gscal,
                            int *norb_per_l, int lval)
{
    int gaussn = 0;
    int i = 0;
    
    for (i = 0; i <= lval; i++) {
        gaussn = gaussn + norb_per_l[i]; 
    }
    gaussn = gaussn - 1;
    bas->gscale = gscal[gaussn];
    return;
}

/*
 * ao_get_atominfo: get atominfo character arrays from soinfo.dat file.
 */
int ao_get_atominfo(FILE *soinfo, char ***ainfo, int *natoms)
{
    int error = 0; /* Error flag */
    char line[MAX_LINE_SIZE];
    int i = 0;
    
    *natoms = find_str_count_in_file("ATOMINFO", soinfo);
    error = allocate_mem_char(ainfo, MAX_LINE_SIZE, *natoms);
    soinfo = find_str_line("ATOMINFO", soinfo);
    for (i = 0; i < *natoms; i++) {
        fgets((*ainfo)[i], MAX_LINE_SIZE, soinfo);
    }
    return error;
}

/*
 * ao_get_lvalue_from_type: get lvalue from the orbital type.
 * 0 = s(1), 1 = p(2,3,4), 2 = d(5,6,7,8,9), 3 = f(10,11,12,13,14,15,16),...
 */
int ao_get_lvalue_from_type(int oindex)
{
    int lvalue = 0; /* Block index */
    
    switch (oindex) {
    case 1: 
        lvalue = 0; /* s orbital */
        break;
    case 2: case 3: case 4:
        lvalue = 1; /* p orbital */
        break;
    case 5: case 6: case 7: case 8: case 9:
        lvalue = 2; /* d orbital */
        break;
    case 10: case 11: case 12: case 13: case 14: case 15: case 16:
        lvalue = 3; /* f orbital */
        break;
    default:
        lvalue = -1; /* basis orbital not implemented */
        error_message(0, ">f orbitals not yet implemented.",
                      "ao_get_block_index_from_type");
        break;
    }
    
    return lvalue;
}

/*
 * ao_get_gscalings: get gaussian scalings from soinfo.dat file.
 */
int ao_get_gscalings(FILE *soinfo, double **gptr, int *ngaus)
{
    int error = 0;              /* Error flag */
    char line[MAX_LINE_SIZE];   /* Scratch line */
    char *lnptr;                /* Line pointer */
    double *gscal;              /* Gaussian scalings */
    int i = 0;
    
    /* Get the total number of gaussians and read in each scaling. */
    *ngaus = find_str_count_in_file("GTOSCALE", soinfo);
    if (*ngaus <= 0) {
        error_message(0, "GTOSCALE not found in soinfo.dat",
                      "ao_build_atomic_orbitalset");
        error++;
        return error;
    }
    gscal = (double *) malloc(sizeof(double) * (*ngaus));
    
    /* Find gaussian scalings in soinfo.dat file. These lines start
     * with 'GTOSCALE'. The format is '(A10,16.10)' */
    soinfo = find_str_line("GTOSCALE", soinfo);
    for (i = 0; i < *ngaus; i++) {
        fgets(line, MAX_LINE_SIZE, soinfo);
        lnptr = substring(line, 10, 16);
        sscanf(lnptr, "%lf", &gscal[i]);
    }
    *gptr = gscal; /* Set pointer to allocated array and return. */ 
    return error;
}

/*
 * ao_get_orbitalnum: get orbital number of soinfo.dat file.
 */
int ao_get_orbitalnum(FILE *soinfo, int *norbs)
{
    int error = 0; /* Error flag */
    
    *norbs = find_str_count_in_file("CAOINFO", soinfo);
    if (*norbs <= 0) {
        error_message(0, "CAOINFO not found in soinfo.dat",
                      "ao_get_orbitalnum");
        error++;
    }
    return error;
}

/*
 * ao_get_block_and_gauss_number: get block number and contracted gaussian
 * number from the current orbital.
 */
void ao_get_block_and_gauss_number(int *norbpl, struct ao_atomdata adata,
                                   int lvalue, int *blockn, int *gaussn)
{
    int minblock = 0; /* Lowest block number for l value */
    int count = 0;
    int i = 0;
    
    /* Get lowest block for l value */
    for (i = 0; i < adata.totblk; i++) {
        if (lvalue == adata.lval[i]) break;
    }
    minblock = i;
    
    /* Count gaussians through l-value blocks */
    for (i = minblock; i < minblock + adata.naob[lvalue]; i++) {
        count = count + adata.basis[i].cgaus;
        if (count >= norbpl[lvalue]) break;
    }
    *blockn = i;
    *gaussn = adata.basis[*blockn].cgaus - (count - norbpl[lvalue]);
    *gaussn = *gaussn - 1;
    return;
}

/*
 * ao_increment_norb_per_l: increment number of orbital per l value if
 * px, d2-, etc.
 */
void ao_increment_norb_per_l(int *nopl, int oindex)
{
    switch (oindex) {
    case 1: case 2: case 5: case 10: case 17:
        (*nopl)++;
        break;
    default:
        break;
    }
    return;
}

/*
 * ao_initialize_atomic_orbital_set: initialize atomic orbital array.
 */
struct ao_basisfunc *ao_initialize_atomic_orbital_set(int norbs)
{
    struct ao_basisfunc *ptr; /* Return value */
    int error = 0; /* Error flag */
    int i = 0;
    
    ptr = (struct ao_basisfunc *)
        malloc(sizeof(struct ao_basisfunc) * norbs);
    if (ptr == NULL) {
        error_message(0, "Error allocating atomic orbital array!",
                      "ao_initialize_atomic_orbital_set");
        return ptr;
    }
    for (i = 0; i < norbs; i++) {
        ptr[i].geom = (double *) malloc(3 * sizeof(double));
    }
    return ptr;
}

/*
 * ao_open_daltonfile: open dalton file returing file stream.
 */
FILE *ao_open_daltonfile(FILE *daltin)
{
    daltin = fopen("daltaoin", "r");
    if (daltin == NULL) {
        error_message(0, "Cannot open daltaoin!", "ao_open_daltonfile");
    }
    return daltin;
}

/*
 * ao_open_soinfodat: open soinfo.dat file, returning file stream.
 */
FILE *ao_open_soinfodat(FILE *soinfo)
{
    soinfo = fopen("soinfo.dat", "r");
    if (soinfo == NULL) {
        error_message(0, "Cannot open soinfo.dat!", "ao_open_soinfodat");
    }
    return soinfo;
}

/*
 * ao_print_dalton_atombasis1(adata)
 */
void ao_print_dalton_atombasis1(struct ao_atomdata adata)
{
    int i = 0;
    printf("Atom name: %s\n", adata.name);
    printf("Nuclear charge: %lf\n", adata.nucchrg);
    printf("Number of atoms: %d\n", adata.natm);
    printf("Largest L value: %d\n", adata.maxlval);
    printf("Atom blocks per l value: ");
    for (i = 0; i <= adata.maxlval; i++) {
        printf(" %d", adata.naob[i]);
    }
    printf("\n");
    printf("Geometries: \n");
    for (i = 0; i < adata.natm; i++) {
        printf("%d %20.15lf%20.15lf%20.15lf\n",
               (i + 1), adata.geom[i][0], adata.geom[i][1],
               adata.geom[i][2]);
    }
    return;
}

/*
 * ao_print_dalton_basisinfo1: print first line of values read from daltaoin.
 */
void ao_print_dalton_basisinfo1(char *crt, int atypes, int molchg,
                                char *symtxt, double thrs)
{
    printf("Spherical harmonics = %s\n", crt);
    printf("Unique atom types =   %d\n", atypes);
    printf("Molecular charge =    %d\n", molchg);
    printf("Symm ops to perform = %s\n", symtxt);
    printf("\nIntegral threshold = %8.2e\n", thrs);
    
    return;
}

/*
 * ao_print_orbital_information: print final atomic orbital information.
 */
void ao_print_orbital_information(struct ao_basisfunc aobasis)
{
    int i = 0;
    printf(" Orbital type:           %10d\n", aobasis.type);
    printf(" Orbital lval:           %10d\n", aobasis.lval);
    printf(" Uncontracted gaussians: %10d\n", aobasis.ugaus);
    printf(" Atomic center index:    %10d\n", aobasis.atom);
    printf(" Geometry:\n");
    printf(" %15.8lf %15.8lf %15.8lf\n",
           aobasis.geom[0], aobasis.geom[1], aobasis.geom[3]);
    printf(" Gaussian scaling:       %10.8lf\n", aobasis.gscale);
    printf(" Alpha & contraction coeffs: \n");
    for (i = 0; i < aobasis.ugaus; i++) {
        printf(" %15.10lf %15.10lf\n", aobasis.alpha[i],
               aobasis.ccoef[i]);
    }
    return;
}

/*
 * ao_process_aorbital_dataline: process atomic orbital data line from
 * soinfo.dat. This line describes the orbital type in the atomic orbital
 * ordering.
 */
void ao_process_aorbital_dataline(char *line, int *anum, char *atom,
                                  int *atype_anum, int *onum, char *ao_desc,
                                  int *oindex, char *atom_desc)
{
    char *lnptr = NULL; /* Line pointer */
    
    /* Format: (A8,I3,A1,I3,I2,A5,I3)
     * "CAOINFO: [atom number][atom][atom number w/in type][num of orbital]
     *           [orbital name desc.][orbital index] */
    lnptr = substring(line, 8, 3);
    sscanf(lnptr, "%d", anum);
    
    /* Read in atom description as character array to compare with ATOMINFO
     * arrays. This will help to switch appropriate geometry information */
    lnptr = substring(line,11, 4);
    sscanf(lnptr, "%s", atom_desc);
    
    lnptr = substring(line,11, 1);
    sscanf(lnptr, "%s", atom);
    lnptr = substring(line,12, 3);
    sscanf(lnptr, "%d", atype_anum);
    lnptr = substring(line,15, 2);
    sscanf(lnptr, "%d", onum);
    lnptr = substring(line,17, 5);
    sscanf(lnptr, "%s", ao_desc);
    lnptr = substring(line,22, 3);
    sscanf(lnptr, "%d", oindex);
    
    return;
}

/*
 * ao_process_orbitaldata: process atomic orbital data, generating an array of
 * atomic orbitals.
 */
int ao_process_orbitaldata(struct ao_atomdata *adata,
                           struct ao_basisfunc *aobasis,
                           FILE *soinfo, int norbs, double *gscal,
                           int ngaus, char **ainfo, int atypes)
{
    int error = 0;            /* Error flag */
    char line[MAX_LINE_SIZE]; /* Scratch line array */
    char *lnptr = NULL;       /* Line pointer */
    int anum = 0;             /* Atom number */
    int old_anum = 0;         /* New atom check. */
    char atom[2];             /* Atom name */
    char atom_desc[5];        /* [Atom name][Atom number within type] */
    int atype_anum = 0;       /* Atom number within type */
    int onum = 0;             /* Orbital number */
    char ao_desc[6];          /* Atomic orbital description */
    int oindex = 0;           /* Orbital index s=1,px=2,py=3,etc. */
    int ocnt = 0;             /* Orbital counter */
    int atyp = 0;             /* Atom type. */
    int block_index = 0;      /* AO Basis block index */
    int lvalue = 0;           /* Orbital lvalue */
    int norb_per_l[10];       /* Number of orbitals per l value */
    int blockn = 0;           /* Block number of new orbital */
    int gaussn = 0;           /* Gaussian number of new orbital */
    int i = 0, j = 0;        
    
    /* Locate 'CAOINFO' in soinfo.dat file. */
    soinfo = find_str_line("CAOINFO", soinfo);
    
    /* Read over each atomic orbital and enter its basis set data.
     * Format: (A8,I3,A1,I3,I2,A5,I3)
     * "CAOINFO: [atom number][atom][atom number w/in type][num of orbital]
     *           [orbital name desc.][orbital index] */
    old_anum = 1;
    init_int_array_0(norb_per_l, 10);
    for (i = 0; i < norbs; i++) {
        fgets(line, MAX_LINE_SIZE, soinfo);
        ao_process_aorbital_dataline(line, &anum, atom, &atype_anum,
                                     &onum, ao_desc, &oindex, atom_desc);
        /* Zero out the l value array for a new atom */
        if (old_anum != anum) {
            init_int_array_0(norb_per_l, 10);
            old_anum = anum;
        }
        /* Find matching atom type */
        for (atyp = 0; atyp <= atypes; atyp++) {
            if (strstr(atom,adata[atyp].name) != NULL) {
                break;
            }
        }
        if (atyp == atypes) return atyp;
        /* Copy geometry. Write orbital index and atom number. */
        cparray_1d1d(adata[atyp].geom[atype_anum - 1], 3,
                     aobasis[i].geom, 3);
        aobasis[i].type = oindex;
        aobasis[i].atom = anum;
        /* Get l value. Increment orbital type if necessary, and
         * get orbital block & gaussian information */
        lvalue = ao_get_lvalue_from_type(aobasis[i].type);
        if (lvalue < 0) return lvalue;
        aobasis[i].lval = lvalue;
        ao_increment_norb_per_l(&norb_per_l[lvalue], oindex);
        ao_get_block_and_gauss_number(norb_per_l, adata[atyp],
                                      lvalue, &blockn, &gaussn);
        aobasis[i].ugaus = adata[atyp].basis[blockn].ugaus;
        aobasis[i].ccoef = (double *)
            malloc(sizeof(double) * aobasis[i].ugaus);
        aobasis[i].alpha = (double *)
            malloc(sizeof(double) * aobasis[i].ugaus);
        /* Copy alpha values */
        ao_cp_alphavals_ccoefs(blockn, gaussn, adata[atyp], &aobasis[i]);
        /* Copy gaussian scaling */
        ao_cp_gaussian_scaling(&aobasis[i], gscal, norb_per_l, lvalue);
        
    }
    
    return error;
}

/*
 * ao_read_dalton_aoblock: read am atomic orbital basis block from daltaoin
 */
int ao_read_dalton_aoblock(FILE *daltin, struct ao_basis *basis)
{
    int error = 0;              /* Error flag */
    char line[MAX_LINE_SIZE];   /* scratch line */
    char *lnptr = NULL;         /* Line pointer */
    char cscr[2];               /* Single character scratch */
    double dscr = 0.0;          /* Single double scratch */
    int newg = 0;               /* New gaussian switch */
    int i = 0, j = 0;
    fpos_t ptr;    /* File position pointer */
    int pos = 0;   /* position */
    int cnt = 0;   /* counter */
    int cclns = 0; /* contraction coefficient lines per gaussian */
    
    /* Read first line of block. This line contains:
     * '(A1,i5,i5)' H   [uncontracted]   [contracted] */
    fgets(line, MAX_LINE_SIZE, daltin);
    lnptr = substring(line, 0, 1);
    sscanf(lnptr, "%s", cscr);
    if ((strstr(cscr,"H")) == NULL) {
        error++;
        return error;
    }
    lnptr = substring(line, 1, 5);
    sscanf(lnptr, "%d", &basis->ugaus);
    lnptr = substring(line, 6, 5);
    sscanf(lnptr, "%d", &basis->cgaus);
    if ((basis->ugaus == 0) || (basis->cgaus == 0)) {
        error++;
        return error;
    }
    basis->alpha = (double *) malloc(sizeof(double) * basis->ugaus);
    error = allocate_mem_double(&basis->ccoef, basis->cgaus, basis->ugaus);
    
    /* Get number of lines per contraction. Dalton output prints 3 
     * coefficients per line, so the number of lines to read must be 
     * computed. Read in alpha and contraction values for block */
    cclns = basis->cgaus / 3;
    if ((basis->cgaus % 3) != 0) cclns++;
    for (i = 0; i < basis->ugaus; i++) {
        fgetpos(daltin, &ptr); // Get position for first line read
        fgets(line, MAX_LINE_SIZE, daltin);
        lnptr = substring(line, 0, 20);
        sscanf(lnptr, "%lf", &basis->alpha[i]);
        fsetpos(daltin, &ptr); // Reset position
        /* Read in contraction values */
        error = ao_read_dalton_contractionvals(daltin, basis, i);
    }
    for (i = 0; i < basis->ugaus; i++) {
        printf("%20.8lf", basis->alpha[i]);
        for (j = 0; j < basis->cgaus; j++) {
            printf("%20.8lf", basis->ccoef[j][i]);
        }
        printf("\n");
    }
    return error; 
}

/*
 * ao_read_dalton_atombasis1: read first line of atom basis set iformation.
 * (Reads in: nuclear charge, number of atoms, max L value, blocks per L value
 */
int ao_read_dalton_atombasis1(FILE *daltin, struct ao_atomdata *adata)
{
    int error = 0;            /* Error flag */
    char line[MAX_LINE_SIZE]; /* scratch line */
    char *lnptr = NULL;       /* line pointer */
    char scr[2];
    int i = 0, pos = 0;
    
    fgets(line, MAX_LINE_SIZE, daltin);
    
    /* Read in values according to the following format: 
     * '(BN,6x,F4.0,I5,24I5)' */
    lnptr = substring(line, 5, 4);
    sscanf(lnptr, "%lf", &adata->nucchrg);
    lnptr = substring(line, 10, 5);
    sscanf(lnptr, "%d", &adata->natm);
    lnptr = substring(line, 15, 5);
    sscanf(lnptr, "%d", &adata->maxlval);
    adata->naob = (int *) malloc(sizeof(int) * adata->maxlval);
    for (i = 0; i < adata->maxlval; i++) {
        pos = 20 + i*5;
        lnptr = substring(line, pos, 5);
        sscanf(lnptr, "%d", &adata->naob[i]);
    }
    
    /* Read in geometries */
    error = allocate_mem_double(&adata->geom, 3, adata->natm);
    if (error != 0) {
        error_message(0, "*** Error allocating data.geom!",
                      "ao_read_dalton_atombasis1");
        return error;
    }
    for (i = 0; i < adata->natm; i++) {
        fgets(line, MAX_LINE_SIZE, daltin);
        lnptr = substring(line, 0, 1);
        sscanf(lnptr, "%s", scr);
        strcpy(adata->name, scr);
        lnptr = substring(line, 4, 60);
        sscanf(lnptr, "%lf %lf %lf", &adata->geom[i][0],
               &adata->geom[i][1], &adata->geom[i][2]);
    }
    
    /* Adjust adata->maxlval to read l value: s=0, p=1, d=2, ... */
    adata->maxlval--;
    
    return error;
};

/*
 * ao_read_dalton_atombasis2: read atomic orbital basis blocks from daltaoin
 * file into the AO_ATOMDATA.AO_ABASIS type.
 */
int ao_read_dalton_atombasis2(FILE *daltin, struct ao_atomdata *adata)
{
    int error = 0;            /* Error flag */
    int i = 0, j = 0;
    int count = 0;
    
    /* Get total number of blocks. Then loop over each block reading in
     * the AO basis definition. */
    adata->totblk = 0;
    for (i = 0; i <= adata->maxlval; i++) {
        adata->totblk = adata->totblk + adata->naob[i];
    }
    adata->lval = (int *) malloc(sizeof(int) * adata->totblk);
    adata->basis = (struct ao_basis *)
        malloc(sizeof(struct ao_basis) * adata->totblk);
    for (i = 0; i <= adata->maxlval; i++) {
        for (j = 0; j < adata->naob[i]; j++) {
            adata->lval[count] = i;
            count++;
        }
    }
    printf("Total blocks to read: %d\n", adata->totblk);
    for (i = 0; i < adata->totblk; i++) {
        printf("Block %d: l=%d\n", (i+1), adata->lval[i]);
        error = ao_read_dalton_aoblock(daltin, &adata->basis[i]);
    }
    return error;
}

/*
 * ao_read_dalton_basisinfo1: read first line of basis set information in
 * dalton file.
 * Input:
 *  daltin = dalton input file file stream.
 * Output:
 *  crt = slpherical harmonics flag
 *  atypes = atom types
 *  molchg = molecular charge
 *  symtxt = symmetry operations test
 *  thrs  = integral threshold
 *  id3 = not referenced
 */
int ao_read_dalton_basisinfo1(FILE *daltin, char *crt, int *atypes,
                              int *molchg, char *symtxt,
                              double *thrs, char *id3)
{
    int error = 0; /* error flag */
    char line[MAX_LINE_SIZE];
    char *lnptr = NULL; /* Line pointer */
    
    fgets(line, MAX_LINE_SIZE, daltin);
    
    /* Read in values according to the following format:
     * '(BN,A1,I4,I3,A2,10A1,D10.2,6I5)  */
    lnptr = substring(line, 0, 1);
    sscanf(lnptr, "%s", crt);
    lnptr = substring(line, 1, 4);
    sscanf(lnptr, "%d", atypes);
    lnptr = substring(line, 5, 3);
    sscanf(lnptr, "%d", molchg);
    lnptr = substring(line, 8, 2);
    sscanf(lnptr, "%s", symtxt);
    lnptr = substring(line, 19, 1);
    sscanf(lnptr, "%s", id3);
    
    /* Set threshold. */
    *thrs = pow(10,-15);
    
    /* Check results. */
    if ((strstr(crt,"s")) == NULL) {
        error++;
        error_message(0, "crt != 's'\nMust use spherical harmonic basis!",
                      "ao_read_dalton_basisinfo1");
        return error;
    }
    if ((strstr(symtxt,"0")) == NULL) {
        error++;
        error_message(0, "symtxt != '0'\nThis is not implemented yet.",
                      "ao_read_dalton_basisinfo1");
        return error;
    }
    return error;
}

/*
 * ao_read_dalton_contractionvals: read contraction values for a uncontracted
 * gaussian.
 */
int ao_read_dalton_contractionvals(FILE *daltin, struct ao_basis *basis, int ug)
{
    int error = 0; /* Error flag */
    int cclns = 0; /* Contraction coefficients line to read in */
    char line[MAX_LINE_SIZE];
    int i = 0, j = 0;
    int cnt = 0; /* Counter */
    int pos = 0; /* stream position */
    char *lnptr = NULL;
    
    cclns = basis->cgaus / 3;
    if ((basis->cgaus % 3) != 0) cclns++;
    /* Loop over contraction coeffcients */
    cnt = 0;
    for (i = 0; i < cclns; i++) {
        fgets(line, MAX_LINE_SIZE, daltin);
        for (j = 0; j < 3; j++) {
            pos = 20 + j * 20;
            lnptr = substring(line, pos, 20);
            sscanf(lnptr, "%lf", &(basis->ccoef[cnt][ug]));
            cnt++;
            if (cnt == basis->cgaus) break;
        }
    }
    
    return error;
}


/*
 * ao_read_daltontitle: read first three lines of daltaoin file.
 * INTGRL
 * X
 * X
 */
int ao_read_daltontitle(FILE *daltin)
{
    int error = 0;
    char firstline[MAX_LINE_SIZE];
    char str[MAX_LINE_SIZE];
    
    sprintf(firstline,"INTGRL");
    fgets(str, MAX_LINE_SIZE, daltin);
    if ((strstr(str,firstline)) == NULL) {
        error++;
        error_message(0, "Error reading daltaoin", "ao_read_daltontitle");
        return error;
    }
    /* Skip next two lines. (They are empty.) */
    fgets(str, MAX_LINE_SIZE, daltin);
    fgets(str, MAX_LINE_SIZE, daltin);
    
    return error;
}
