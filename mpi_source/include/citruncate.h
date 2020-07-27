// FILE: citruncate.h
#ifndef citruncate_h
#define citruncate_h

/*
 * struct eospace: electron numbers in a space making
 * a valid determinant.
 * eospace.docc + eospace.actv + eospace.virt = total elec (alpha/beta) */
struct eospace {
    int docc;
    int actv;
    int virt;
    int nstr;          /* Number of strings in grouping */
    int start;         /* Starting string index in list */
    int pairs[20];     /* Spaces of opposite spin that pair with this space
                        * to form a determinant in the expansion. */
    int npairs;        /* Number of spaces that can be paired */
};

/*
 * struct eostring: a valid alpha/beta electron string.
 */
struct eostring {
        int doccx;
        int actvx;
        int index;
        int *string;
};

/*
 * struct xstr: excitation string information.
 */
struct xstr {
    int index;   /* Position of string in pstrings/qstrings */
    int io[2];   /* Initial orbitals */
    int fo[2];   /* Final orbitals */
    int permx;   /* Permutation coefficient -1 or 1 */
};

/*
 * struct xstrmap: starting string index and ending string index in
 * struct eostring list of valid strings for excitation.
 */
struct xstrmap {
        int start;
        int finish;
};

/*
 * citrunc: truncate the CI expansion.
 * Input:
 *  aelec = alpha electrons
 *  belec = beta  electrons
 *  orbs  = molecular orbitals
 *  nfrzc = number of frozen core orbitals
 *  ndocc = number of doubly-occupied reference orbitals
 *  nactv = number of active reference orbitals
 *  nfrzv = number of frozen virtual orbitals
 *  xlvl  = excitation level (DOCC->ACTV & (DOCC + ACTV)->VIRT)
 *  astrings = alpha electron strings
 *  astr_len = alpha string number
 *  bstrings = beta  electron strings
 *  bstr_len = beta  string number
 *  peosp    = alpha electron space array
 *  pegrps   = alpha electron groupings
 *  qeosp    = beta electron space array
 *  qegrps   = beta electron groupings
 *  astr_len = alpha string number
 *  bstr_len = beta  string number
 *  dtrm_len = determinant number
 *  pq_spaces = (p,q) space pairings
 *  num_pq    = number of (p,q) pairings
 */
int citrunc(int aelec, int belec, int orbs, int nfrzc, int ndocc,
	    int nactv, int nfrzv, int xlvl, struct occstr *astrings,
            int astr_len, struct occstr *bstrings, int bstr_len,
            struct eospace *peosp, int pegrps,
            struct eospace *qeosp, int qegrps, int *dtrm_len,
            int **pq_spaces, int *num_pq);

/*
 * allocate_eospace_array: allocate the electron number space array.
 */
struct eospace *allocate_eospace_array(int nelec, int norbs, int ndocc,
                                       int nactv, int xlvl,  int *ngrps);

/*
 * allocate_strings_array: allocate the strings array.
 */
struct eostring *allocate_strings_array(
        int nstr,
        int nelecs
        );
/*
 * allocate_occstr_arrays: allocate binary electron orbital occupation
 * string arrays.
 */
struct occstr *allocate_occstr_arrays(int nstr);

/*
 * allocate_xmap: allocate 2d excitation map.
 */
struct xstrmap **allocate_xmap(int xlvl);

/*
 * construct_xlist: construct list of X excitations for electron strings
 */
int **construct_xlist(struct occstr *strlist, int nstr, int intorb, int x, int *max);

/*
 * compute_ci_elecs_and_orbitals: compute number of ci electrons and orbitals
 */
void compute_ci_elecs_and_orbitals(int aelec, int belec, int orbitals, int nfrzc,
                                   int nfrzv, int *ci_aelec, int *ci_belec,
                                   int *ci_orbitals);

/*
 * compute_detnum: compute the number of determinants in expansion.
 */
int compute_detnum(struct eospace *peosp, int pegrps, struct eospace *qeosp,
                   int qegrps, int ndocc, int nactv, int xlvl,
                   int **pq_spaces, int *num_pq);

/*
 * compute_eostrings_sr: compute electron occupation strings for all string
 * combinations of DOCC, and VIRT electrons.
 */
void compute_eostrings_sr(struct eostring *strlist, int *pos, int ci_orbs,
                          int ndocc, int docc_elec, int virt_elec, int *docc,
                          int *virt, int *doccscr, int *virtscr);

/*
 * compute_stringnum: compute number of valid strings in expansion.
 */
int compute_stringnum(int orbs, int elecs, int ndocc, int nactv, int xlvl);

/*
 * deallocate_eostrings_array: deallocate a struct eostring *array.
 */
void deallocate_eostrings_array(struct eostring *array, int nstr);

/*
 * generate_binstring_list: convert the integer string lists to struct occstr
 */
void generate_binstring_list(struct eostring *str, int nstr, int elec,
                             int ndocc, int nactv, struct occstr *binstr);


/*
 * generate_determinant_list: generate determinant list.
 */
void generate_determinant_list(struct eostring *pstrlist, int npstr, int aelec,
                               struct eostring *qstrlist, int nqstr, int belec,
                               struct eospace *peosp, int pegrps,
                               struct eospace *qeosp, int qegrps,
                               int ndocc, int nactv, int xlvl, int *dcnt);

/*
 * generate_determinant_list_rtnlist: generate determinant list, returning
 * the list.
 */
void generate_determinant_list_rtnlist(struct eostring *pstrlist, int npstr,
                                       int aelec, struct eostring *qstrlist,
                                       int nqstr, int belec,
                                       struct eospace *peosp, int pegrps,
                                       struct eospace *qeosp, int qegrps,
                                       int ndocc, int nactv, int xlvl,
                                       int dtrm_len, struct det *dtlist);

/*
 * generate_actvx: generate excitations within the ACTV space given a
 * string.
 */
void generate_actvx(int nrep, struct occstr str, int str_docc, int str_actv,
                    int str_virt, struct eospace eosp, int ndocc, int nactv,
                    int nvirt, int elec, int *scr, struct xstr *xlist, int *numx);

/*
 * generate_actv2virtx: generate replacements between DOCC and ACTV spaces.
 */
void generate_actv2virtx(int nrep, struct occstr str, int str_docc, int str_actv,
                         int str_virt, struct eospace eosp, int ndocc, int nactv,
                         int nvirt, int elec, int *scr, struct xstr *xlist, int *numx);

/*
 * generate_actv2virtx_actv1: generate replacements between ACTV and VIRT spaces.
 * This is a special case where one replacement occurs within ACTV, in addition
 * to the ACTV <-> VIRT replacement.
 */
void generate_actv2virtx_actv1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);

/*
 * generate_actv2virtx_docc1: generate replacements between DOCC, ACTV and VIRT
 * spaces. This is a special case where one replacement occurs within DOCC, in addition
 * to the ACTV <-> VIRT replacement.
 */
void generate_actv2virtx_docc1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);


/*
 * generate_actv2virtx_actv1: generate replacements between ACTV and VIRT spaces.
 * This is a special case where one replacement occurs within VIRT, in addition
 * to the ACTV <-> VIRT replacement.
 */
void generate_actv2virtx_virt1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);

/*
 * generate_docc2actvx_docc1: generate excitations from DOCC -> ACTV with
 * a replacement within the DOCC.
 */
void generate_docc2actvx_docc1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);

/*
 * generate_docc2actvx_actv1: generate excitations from DOCC -> ACTV with
 * a replacement within the ACTV.
 */
void generate_docc2actvx_actv1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);

/*
 * generate_docc2actvx_virt1: generate excitations from DOCC -> ACTV with
 * a replacement within the VIRT.
 */
void generate_docc2actvx_virt1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);
/*
 * generate_docc2virtx_docc1: generate excitations from DOCC -> VIRT with
 * a replacement within the DOCC.
 */
void generate_docc2virtx_docc1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);
/*
 * generate_docc2virtx_actv1: generate excitations from DOCC -> VIRT with
 * a replacement within the ACTV.
 */
void generate_docc2virtx_actv1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);

/*
 * generate_docc2virtx_actv1: generate excitations from DOCC -> VIRT with
 * a replacement within the virt.
 */
void generate_docc2virtx_virt1(int nrep, struct occstr str, int str_docc,
                               int str_actv, int str_virt, struct eospace eosp,
                               int ndocc, int nactv, int nvirt, int elec,
                               int *scr, struct xstr *xlist, int *numx);

/*
 * generate_doccx_actvx: generate replacements  DOCC + ACTV.
 */
void generate_doccx_actvx(int nrep, struct occstr str, int str_docc,
                          int str_actv, int str_virt, struct eospace eosp,
                          int ndocc, int nactv, int nvirt, int elec,
                          int *scr, struct xstr *xlist, int *numx);

/*
 * generate_doccx_virtx: generate replacements  DOCC + VIRT.
 */
void generate_doccx_virtx(int nrep, struct occstr str, int str_docc,
                          int str_actv, int str_virt, struct eospace eosp,
                          int ndocc, int nactv, int nvirt, int elec,
                          int *scr, struct xstr *xlist, int *numx);

/*
 * generate_actv_virtx: generate replacements  ACTV + VIRT.
 */
void generate_actvx_virtx(int nrep, struct occstr str, int str_docc,
                          int str_actv, int str_virt, struct eospace eosp,
                          int ndocc, int nactv, int nvirt, int elec,
                          int *scr, struct xstr *xlist, int *numx);

/*
 * generate_actv2doccvirtx: generate replacements between ACTV and (DOCC,VIRT)
 * spaces.
 */
void generate_actv2doccvirtx(int nrep, struct occstr str, int str_docc,
                             int str_actv, int str_virt, struct eospace eosp,
                             int ndocc, int nactv, int nvirt, int elec, int *scr,
                             struct xstr *xlist, int *numx);

/*
 * generate_doccx: generate excitations within the DOCC space given a
 * string.
 */
void generate_doccx(int nrep, struct occstr str, int str_docc, int str_actv,
                    int str_virt, struct eospace eosp, int ndocc, int nactv,
                    int nvirt, int elec, int *scr, struct xstr *xlist, int *numx);

/*
 * generate_docc2actvx: generate replacements between DOCC and ACTV spaces.
 */
void generate_docc2actvx(int nrep, struct occstr str, int str_docc, int str_actv,
                         int str_virt, struct eospace eosp, int ndocc, int nactv,
                         int nvirt, int elec, int *scr, struct xstr *xlist, int *numx);

/*
 * generate_docc2actvvirtx: generate excitations from DOCC -> (ACTV, VIRT)
 * for a given string.
 */
void generate_docc2actvvirtx(int nrep, struct occstr str, int str_docc,
                             int str_actv, int str_virt, struct eospace eosp,
                             int ndocc, int nactv, int nvirt, int elec, int *scr,
                             struct xstr *xlist, int *numx);

/*
 * generate_docc2virtx: generate excitations from DOCC -> VIRT spaces for
 * a given string.
 */
void generate_docc2virtx(int nrep, struct occstr str, int str_docc, int str_actv,
                         int str_virt, struct eospace eosp, int ndocc, int nactv,
                         int nvirt, int elec, int *scr, struct xstr *xlist, int *numx);

/*
 * generate_virtx: generate excitations within the VIRT space given a
 * string.
 */
void generate_virtx(int nrep, struct occstr str, int str_docc, int str_actv,
                    int str_virt, struct eospace eosp, int ndocc, int nactv,
                    int nvirt, int elec, int *scr, struct xstr *xlist, int *numx);

/*
 * generate_virt2doccactvx: generate replacements between VIRT and (DOCC,ACTV)
 * spaces.
 */
void generate_virt2doccactvx(int nrep, struct occstr str, int str_docc,
                             int str_actv, int str_virt, struct eospace eosp,
                             int ndocc, int nactv, int nvirt, int elec, int *scr,
                             struct xstr *xlist, int *numx);

/*
 * generate_single_excitations: generate the single excitations in
 * EOSPACE for an input string.
 * Note: DOCC has no replacements within an EOSPACE.
 */
int generate_single_excitations(struct occstr str, struct eospace eosp,
                                 int elec, int ndocc, int nactv,
                                 int intorb, int vorbs,
                                 struct xstr *singlex,
                                 int *elecs, int *orbsx);
/*
 * generate_double_excitations: generate the double excitations in
 * EOSPACE for an input string.
 */
int generate_double_excitations(struct occstr str, struct eospace eosp,
                                int nelec, int ndocc, int nactv,
                                int intorb, int vorbs,
                                struct xstr *doublex,
                                int *elecs, int *orbsx);


/*
 * generate_string_list: generate full *valid* alpha/beta string lists.
 */
void generate_string_list(struct eostring *strlist, int nstr, int orbs,
                          int elecs, int ndocc, int nactv, int xlvl,
                          struct eospace *eosp, int egrps);

/*
 * get_eospace_detrange: get starting and final determinant indexes
 * for an eospace pairing.
 */
void get_eospace_detrange(int **pq, int npq, int pair, struct eospace *peosp,
                          int pegrps, struct eospace *qeosp, int qegrps,
                          int *start, int *final);

/*
 * get_string_eospace: get eospace of input string.
 */
int get_string_eospace(struct occstr str, int ndocc, int nactv,
                       struct eospace *esp, int egrps);

/*
 * occstr2address: compute the string index of given an occupation string.
 */
int occstr2address(struct occstr str, struct eospace eosp, int ndocc, int nactv,
                   int nvirt, int nelec, int *elecs);

/*
 * remove_grt_xstr: remove replacements with indices greater than
 * val. (Modification of remove_leq_int in arrayutil.c)
 */
void remove_grt_xstr(int val, struct xstr *list, int *nx, struct xstr *scr);

/*
 * remove_leq_xstr: remove replacements with indices less than, or equal to
 * val. (Modification of remove_leq_int in arrayutil.c)
 */
void remove_leq_xstr(int val, struct xstr *list, int *nx, struct xstr *scr);

/*
 * setup_eostrings_compute: set up a electron, orbital, strings arrays for
 * compuation of the electron occupation strings.
 */
void setup_eostrings_compute(int *elecs, int *orbs, int *nstr, int *pegs,
                             int *nspcs,
                             int docc_elec, int actv_elec, int virt_elec,
                             int dstr, int vstr, int astr, int ndocc,
                             int nactv, int vorbs);

/*
 * string_number: generate number of possible *valid* strings given the number
 * of DOCC electrons, ACTV electrons, and VIRT electrons.
 */
int string_number(int nde, int nae, int nve, int ndocc, int nactv, int orbs,
                  int elecs);

/*
 * write_determinant_strpairs: write determinant alpha/beta string pairs for
 * valid determinants to file.
 */
void write_determinant_strpairs(FILE *fptr, int pstart, int pnstr, int qstart,
                                int qnstr, int *cnt,
                                struct eostring *pstr, struct eostring *qstr);

/*
 * write_determinant_strpairs_dtlist: write determinant alpha/beta string
 * pairs for valid determinants to dtlist.
 */
void write_determinant_strpairs_dtlist(int pstart, int pnstr, int qstart,
                                       int qnstr, struct eostring *pstr,
                                       struct eostring *qstr, int aelec,
                                       int belec, int ndocc, int nactv,
                                       struct det *dtlist, int *dptr,
                                       int cflag);

#endif

