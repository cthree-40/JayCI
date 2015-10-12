// File: action_util.h
#ifndef action_util_h
#define action_util_h

double eval0_cas(struct det    deti, 
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

double evaluate_dets_virt(int detdiff, 
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

double evaluate_dets_virtcas(int ddiff, 
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

double eval1_1000_vc(struct occstr ostr1i, 
		     struct occstr ostr2i,
		     struct occstr ostr1j,
		     long long int xi1, 
		     long long int xf1, 
		     int nelec1,
		     int nelec2, 
		     double *moints1, 
		     double *moints2);

double eval2_0110_vc(struct occstr ostr1i, 
		     struct occstr ostr1j,
		     int nelec1, 
		     int nelec2, 
		     long long int xi2,
		     long long int xf2, 
		     double *moints2);

double eval0_00_virt(struct det        deti, 
		     int              aelec, 
		     int              belec, 
		     double        *moints1, 
		     double        *moints2);


double eval1_10_virt(struct occstr ostr1i,
		     struct occstr ostr1j,
		     struct occstr ostr2i,
		     int nelec1,
		     int nelec2,
		     double *moints1,
		     double *moints2);

double eval2_11_virt(int *avxi, 
		     int *bvxi, 
		     int *avxj, 
		     int *bvxj, 
		     double *moints2);

double eval2_20_virt(int *vxi, 
		     int *vxj, 
		     double *moints2);

double hmatels(struct det deti,
	       struct det detj,
	       double *moints1,
	       double *moints2,
	       int m1len,
	       int m2len,
	       int aelec,
	       int belec);

void make_orbital_strings_virt(struct occstr ostr1,
			       int *eostr1,
			       int nelec1);

void virtdiffs_1(int *vxi, 
		 int *vxj, 
		 int *ifo);

#endif
