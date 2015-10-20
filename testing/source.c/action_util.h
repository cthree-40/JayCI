/* File: action_util.h */
/* 
 * note: the subfunctions in this file are listed alphabetically,
 *       not as they appear in action_util.c 
 */
#ifndef action_util_h
#define action_util_h

/* 
 * cas_to_virt_replacements: compute excitations for cas<->virt replacements
 */
void cas_to_virt_replacements(
	int nreps,           /* number of cas->virt replacements */
	int ncr,             /* number of cas->cas replacements */
	int nvr,             /* number of virt->virt replacements */
	long long int xi,    /* cas byte of initial orbitals */
	long long int xf,    /* cas byte of final orbitals */
	int *restrict vxi,   /* inital virtual orbitals of excitation */
	int *restrict vxj,   /* final virtual orbitals of excitation */
	int *restrict reps); /* replacement arrays */

/*
 * eval0_cas: evaluate diagonal matrix elements between cas-flagged
 * determinants
 */
double eval0_cas(
	struct det    deti,  /* <i|H|i> */ 
	int          aelec,  /* alpha electrons */
	int          belec,  /* beta  electrons */
	double    *moints1,  /* 1-e integrals   */
	double    *moints2); /* 2-e integrals   */

/* 
 * eval0_ncas: evaluate diagonal matrix elements between non-cas-flagged
 * determinants 
 */
double eval0_ncas(
	struct det    deti,  /* <i|H|i> */ 
	int          aelec,  /* alpha electrons */
	int          belec,  /* beta  electrons */
	double    *moints1,  /* 1-e integrals   */
	double    *moints2); /* 2-e integrals   */

/*
 * eval1_10_cas: evaluate the matrix element of a single replacement between
 * two cas-flagged determinants.
 */
double eval1_10_cas(
	struct occstr ostr1, /* (alpha/beta) occupation string */ 
	long long int    xi, /* initial orbitals of excitation */
	long long int    xf, /* final orbitals of excitation   */
	struct occstr ostr2, /* (beta/alpha) occupation string */
	int             ne1, /* (alpha/beta) electrons */
	int             ne2, /* (beta/alpha) electrons */
	double     *moints1, /* 1-e integrals */
	double     *moints2);/* 2-e integrals */

/*
 * eval1_ncas_c0cv0v1: evaluate the matrix element of a single virtual 
 * virtual replacement between two non-cas-flagged determinants. 
 */
double eval1_ncas_c0cv0v1(
	struct occstr ostr1i, /* (alpha/beta) occupation string of det i */ 
	struct occstr ostr1j, /* (alpha/beta) occupation string of det j */
	int ne1,              /* (alpha/beta) electrons */
	struct occstr ostr2i, /* (beta/alpha) occupation string of det i */
	int ne2,              /* (beta/alpha) electrons */
	double *moints1,      /* 1-e integrals */
	double *moints2);     /* 2-e integrals */

/* 
 * eval1_ncas_c0cv1v0: evaluate the matrix element of a single cas->virtual
 * replacement between two non-cas-flagged determinants.
 */
double eval1_ncas_c0cv1v0(
	long long int xi,    /* initial CAS orbitals of excitation */
	long long int xf,    /* final CAS orbitals of excitation   */
	struct occstr str1i, /* deti (alpha/beta) occupation string */
	struct occstr str1j, /* detj (alpha/beta) occupation string */
	struct occstr str2i, /* deti (beta/alpha) occupation string */
	int ne1,             /* (alpha/beta) electrons */
	int ne2,             /* (beta/alpha) electrons */
	double *moints1,     /* 1-e integrals          */
	double *moints2);    /* 2-e integrals          */

/*
 * eval1_ncas_c1cv0v0: evaluate the matrix element of a single cas->cas 
 * replacement between two non-cas-flagged determinants.
 */
double eval1_ncas_c1cv0v0(
	struct occstr       occ_str1,  /* (alpha/beta) occupation string */ 
	long long int  init_orbs_cas,  /* inital orbitals of excitation  */
	long long int  finl_orbs_cas,  /* final orbitals of excitation   */
	struct occstr       occ_str2,  /* (beta/alpha) occupation string */
	int                   nelec1,  /* (alpha/beta) electrons */
	int                   nelec2,  /* (beta/alpha) electrons */
	double              *moints1,  /* 1-e integrals */
	double              *moints2); /* 2-e integrals */

/*
 * eval2_ncas_c00cv00v11: evaluate single virtual replacements in both
 * alpha and beta strings of non-cas-flagged determinants
 */
double eval2_ncas_c00cv00v11(
	int *avirti,       /* alpha virtual orbitals of determinant i */
	int *avirtj,       /* alpha virtual orbitals of determinant j */
	int *bvirti,       /* beta virtual orbitals of determinant i */
	int *bvirtj,       /* beta virtual orbitals of determinant j */
	int navrtx,       /* number of alpha virtual orbitals */
	int nbvrtx,       /* number of beta virtual orbitals */
	double *moints2);  /* 2-e integrals */

/* 
 * eval2_ncas_c00cv10v01: evaluate double replacements for non-cas-flagged
 * determinants with one cas->virt replacement in one string and one
 * virt->virt replacement in the other string.
 */
double eval2_ncas_c00cv10v01(
	long long int xi,    /* (alpha/beta) initial orbitals of excitation */ 
	long long int xf,    /* (alpha/beta) final orbitals of excitation */
	int *vx1i,           /* (alpha/beta) virtual orbitals of det i */
	int *vx1j,           /* (alpha/beta) virtual orbitals of det j */
	int *vx2i,           /* (beta/alpha) virtual orbitals of det i */
	int *vx2j,           /* (beta/alpha) virtual orbitals of det j */
	double *moints2);    /* 2-e integrals */


double eval2_ncas_c00cv10v10(
	long long int xi, 
	long long int xf,
	struct occstr stri, 
	struct occstr strj,
	double *moints2);

/* 
 * eval2_ncas_c00cv11v00: evaluate double replacements for non-cas-flagged
 * determinants with one cas->virt replacement in each string. 
 */
double eval2_ncas_c00cv11v00(
	long long int xi1, /* (alpha/beta) initial orbitals of excitation */ 
	long long int xf1, /* (alpha/beta) final orbitals of excitation */
	int *vxi1,         /* (alpha/beta) virtual orbitals of det i */
	int *vxj1,         /* (alpha/beta) virtual orbitals of det j */
	long long int xi2, /* (beta/alpha) initial orbitals of excitation */
	long long int xf2, /* (beta/alpha) final orbitals of excitation */
	int *vxi2,         /* (beta/alpha) virtual orbitals of det i */
	int *vxj2,         /* (beta/alpha) virtual orbitals of det j */
	double *moints2);  /* 2-e integrals */

/* 
 * eval2_ncas_c01cv10v00: evaluate double replcaements for non-cas-flagged
 * determinants with one cas->virt replacement in one string and one
 * cas->cas replacement in the other string.
 */
double eval2_ncas_c01cv10v00(
	struct occstr str1, /* (alpha/beta) determinant string */
	struct occstr str2, /* (beta/alpha) determinant string */
	long long int xi1,  /* (alpha/beta) initial orbitals of excitation */ 
	long long int xf1,  /* (alpha/beta) final orbitals of excitation */
	int *vx1i,          /* (alpha/beta) virtual orbital occupations */
	int *vx1j,          /* (alpha/beta) virtual orbital occupations */
	long long int xi2,  /* (beta/alpha) initial orbitals of excitation */
	long long int xf2,  /* (beta/alpha) final orbitals of excitation */
	int ne1,            /* (alpha/beta) electrons */
	int ne2,            /* (beta/alpha) electrons */
	double *moints2);   /* 2-e integrals */

/*
 * eval2_ncas_c10cv10v00: evaluate double replacements for non-cas-flagged
 * determinants with one cas->virt replacement and one cas->cas replacement
 * int the same string.
 */
double eval2_ncas_c10cv10v00(
	struct occstr str, /* (alpha/beta) determinant string */
	long long int xi,  /* (alpha/beta) initial orbitals of excitation */ 
	long long int xf,  /* (alpha/beta) final orbitals of excitaiton */
	int *vxi,          /* (alpha/beta) virtual orbitals of det i */
	int *vxj,          /* (alpha/beta) virtual orbitals of det j */
	int ne,            /* (alpha/beta) electrons */
	double *moints2);  /* 2-e integrals */

double eval2_ncas_c0cv0v2(
	int *vxi, 
	int *vxj, 
	double *moints2);

double eval2_ncas_c0cv2v0(
	long long int xi, 
	long long int xf, 
	int *vxi, 
	int *vxf, 
	double *moints2);

/*
 * eval2_ncas_c1cv0v1: evaluate cas + virtual replacements for non-cas-flagged
 * determinants.
 */
double eval2_ncas_c1cv0v1(
	struct occstr str, /* string containing CAS excitation */
	long long int xi,  /* initial orbitals of CAS excitation */
	long long int xf,  /* final orbitals of CAS excitation */
	int *vxi,          /* inital virtual orbitals */
	int *vxj,          /* final virtual orbitals */
	int ne,            /* number of electrons in CAS excitation string */
	int nvx,           /* number of virtual orbitals in vxi */
	double *moints2);  /* 2-e integrals */

double eval2_11_cas(long long int axi, 
		    long long int axf, 
		    long long int bxi,
		    long long int bxf, 
		    double   *moints2);

double eval2_20_cas(long long int xi, 
		    long long int xf, 
		    double  *moints2);

double evaluate_dets_cas(int         ndiff,
			 struct det   deti,
			 struct det   detj,
			 int         numax,
			 int         numbx,
			 long long int axi,
			 long long int axf,
			 long long int bxi,
			 long long int bxf,
			 int         aelec,
			 int         belec, 
			 double   *moints1,
			 double   *moints2);

double evaluate_dets_ncas(
	int          ndiff,
	struct det     deti,
	struct det     detj,
	int          numaxc,
	int          numbxc,
	int         numaxcv,
	int         numbxcv,
	int          numaxv,
	int          numbxv,
	long long int   axi,
	long long int   axf,
	long long int   bxi,
	long long int   bxf,
	int           aelec,
	int           belec, 
	double     *moints1,
	double     *moints2);

double evaluate_dets_virt(
	int detdiff, 
	struct det deti,
	struct det detj,
	int numaxc, 
	int numbxc,
	int numaxv,
	int numbxv,
	long long int axi,
	long long int axf, 
	long long int bxi,
	long long int bxf,
	double *moints1,
	double *moints2,
	int alec,
	int belec);

double evaluate_dets_virtcas(
	int ddiff, 
	struct det deti, 
	struct det detj,
	int naxc, 
	int nbxc, 
	int naxv, 
	int nbxv,
	long long int axi, 
	long long int axf,
	long long int bxi, 
	long long int bxf,
	double *moints1, 
	double *moints2, 
	int aelec, 
	int belec);

/*
 * compute_hv: perform Hv=c
 */
void compute_hv(
	struct det *dlist, /* list of determinants */ 
	int ndets,         /* number of determinants */
	double *moints1,   /* 1-e integrals */
	double *moints2,   /* 2-e integrals */
	int aelec,         /* alpha electrons */
	int belec,         /* beta  electrons */
	double *restrict v,         /* input vector */
	double *restrict c);        /* output vector */

double hmatels(
	struct det deti,
	struct det detj,
	double *moints1,
	double *moints2,
	int aelec,
	int belec);

void make_orbital_strings_virt(
	struct occstr ostr1,
	int *eostr1,
	int nelec1);

/* 
 * pindes_single_rep: compute permutational index of single replacement
 */
int pindex_single_rep(
	int *str,  /* electron orbital string */ 
	int io,    /* initial orbital */
	int fo,    /* final orbital */
	int lstr);  /* number of electrons */

double single_rep_2e_contribution(
	int *eocc_str1, 
	int init_orb, 
	int finl_orb,
	int perm_index, 
	int *eocc_str2, 
	int nelec1,
	int nelec2, 
	double *moints2);

void virtdiffs_single_cas_to_virt(
	int *vxi, 
	int *vxj, 
	int *repo);

void virtdiffs_single_rep(
	int *vxi, 
	int *vxj, 
	int *ifo);

#endif
