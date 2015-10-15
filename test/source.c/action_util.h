// File: action_util.h
/* 
 * note: the subfunctions in this file are listed alphabetically,
 *       not as they appear in action_util.c 
 */
#ifndef action_util_h
#define action_util_h

double eval0_cas(struct det    deti, 
		 int          aelec, 
		 int          belec, 
		 double    *moints1, 
		 double    *moints2);

double eval0_ncas(struct det    deti, 
 	 	  int          aelec, 
		  int          belec, 
		  double    *moints1, 
		  double    *moints2);

double eval1_10_cas(struct occstr ostr1, 
		    long long int    xi, 
		    long long int    xf,
		    struct occstr ostr2, 
		    int             ne1, 
		    int             ne2, 
		    double     *moints1,
		    double     *moints2);
double eval1_ncas_c0cv0v1(
	struct occstr ostr1i, 
	struct occstr ostr1j,
	int ne1, 
	struct occstr ostr2i, 
	int ne2,
	double *moints1, 
	double *moints2);

double eval1_ncas_c1cv0v0(struct occstr       occ_str1, 
			  long long int  init_orbs_cas,
			  long long int  finl_orbs_cas, 
			  struct occstr       occ_str2,
			  int                   nelec1, 
			  int                   nelec2, 
			  double              *moints1, 
			  double              *moints2);

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

double hmatels(
	struct det deti,
	struct det detj,
	double *moints1,
	double *moints2,
	int m1len,
	int m2len,
	int aelec,
	int belec);

void make_orbital_strings_virt(
	struct occstr ostr1,
	int *eostr1,
	int nelec1);

double single_rep_2e_contribution(
	int *eocc_str1, 
	int init_orb, 
	int finl_orb,
	int perm_index, 
	int *eocc_str2, 
	int nelec1,
	int nelec2, 
	double *moints2);


void virtdiffs_single_rep(
	int *vxi, 
	int *vxj, 
	int *ifo);

#endif
