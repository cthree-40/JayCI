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
#include "ioutil.h"
#include "mathutil.h"
#include "buildao.h"

/*
 * ao_buildao: build atomic orbitals
 */
int ao_buildao()
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
        struct ao_atomdata *adata = NULL; /* Atom basis data */
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
 * ao_open_daltonfile: open dalton file returing file stream.
 */
FILE *ao_open_daltonfile(FILE *daltin)
{
        daltin = fopen("daltaoin", "r");
        if (daltin == NULL) {
                error_message("Cannot open daltaoin!", "ao_open_daltonfile");
        }
        return daltin;
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
        printf("%s\n", lnptr);
        sscanf(lnptr, "%d", &adata->natm);
        lnptr = substring(line, 15, 5);
        sscanf(lnptr, "%d", &adata->maxlval);
        adata->naob = (int *) malloc(sizeof(int) * adata->maxlval);
        for (i = 0; i < adata->maxlval; i++) {
                pos = 20 + i;
                lnptr = substring(line, pos, 5);
                sscanf(lnptr, "%d", &adata->naob[i]);
        }

        /* Read in geometries */
        error = allocate_mem_double(&adata->geom, 3, adata->natm);
        if (error != 0) {
                error_message("*** Error allocating data.geom!",
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
                printf("Block %d\n", i);
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
                error_message("crt != 's'\nMust use spherical harmonic basis!",
                              "ao_read_dalton_basisinfo1");
                return error;
        }
        if ((strstr(symtxt,"0")) == NULL) {
                error++;
                error_message("symtxt != '0'\nThis is not implemented yet.",
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
                error_message("Error reading daltaoin", "ao_read_daltontitle");
                return error;
        }
        /* Skip next two lines. (They are empty.) */
        fgets(str, MAX_LINE_SIZE, daltin);
        fgets(str, MAX_LINE_SIZE, daltin);

        return error;
}
