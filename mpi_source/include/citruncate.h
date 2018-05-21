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
        int nstr;  /* Number of strings in grouping */
        int start; /* Starting string index in list */
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
 * Output:
 *  astr_len = alpha string number
 *  bstr_len = beta  string number
 *  dstr_len = determinant number
 */
int citrunc(int aelec, int belec, int orbs, int nfrzc, int ndocc,
	    int nactv, int nfrzv, int xlvl, struct occstr *astrings,
            int astr_len, struct occstr *bstrings, int bstr_len,
            struct eospace *peosp, int pegrps,
            struct eospace *qeosp, int qegrps, int *dtrm_len);

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
 * compute_ci_elecs_and_orbitals: compute number of ci electrons and orbitals
 */
void compute_ci_elecs_and_orbitals(int aelec, int belec, int orbitals, int nfrzc,
                                   int nfrzv, int *ci_aelec, int *ci_belec,
                                   int *ci_orbitals);

/*
 * compute_detnum: compute the number of determinants in expansion.
 */
int compute_detnum(struct eospace *peosp, int pegrps, struct eospace *qeosp,
                   int qegrps, int ndocc, int nactv, int xlvl);

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
 * generate_string_list: generate full *valid* alpha/beta string lists.
 */
void generate_string_list(struct eostring *strlist, int nstr, int orbs,
                          int elecs, int ndocc, int nactv, int xlvl,
                          struct eospace *eosp, int egrps);

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

