// File: citruncate_2.h

#ifndef citruncate_2_h
#define citruncate_2_h

struct eostring {
        int index;
        int *string;
};

/*
 * allocate_strings_array: allocate the strings array.
 */
struct eostring *allocate_strings_array(
        int nstr,
        int nelecs
        );

/*
 * check_ci_input_values: check the validity of truncation values before
 * execution.
 */
int check_ci_input_values(int cia, int cib, int cio, int nd);

/*
 * citrunc: truncate the ci expansion.
 */
int citrunc(int aelec, int belec, int orbs, int nfrzc, int ndocc,
	    int nactv, int nfrzv, int xlvl, int *astr_len,
	    int *bstr_len, int *dtrm_len);

/*
 * compute_eostrings: compute electron occupation strings for all string
 * combinations of DOCC, ACTV, and VIRT electrons.
 */
void compute_eostrings(struct eostring *strlist, int *pos, int ci_orbs,
                       int ndocc, int nactv, int docc_elec, int actv_elec,
                       int virt_elec, int *docc, int *actv, int *virt,
                       int *doccscr, int *actvscr, int *virtscr);

/*
 * compute_stringnum: compute number of valid strings in expansion.
 */
int compute_stringnum(int orbs, int elecs, int ndocc, int nactv, int xlvl);


/*
 * generate_string_list: generate full *valid* alph/beta string lists.
 */
void generate_string_list(struct eostring *strlist, int nstr, int orbs,
                          int elecs, int ndocc, int nactv, int xlvl);

/*
 * generate_actv_strings: generate ACTV orbital strings array.
 */
int generate_actv_strings(int **qactv_array, int nstr, int nactv, int elecs,
                          int ndocc, int xlvl);

/*
 * generate_docc_strings: generate DOCC orbital strings array.
 */
int generate_docc_strings(int **qdocc_array, int nstr, int orbs, int elec,
                          int xlvl);
/*
 * generate_virt_strings: generate VIRT orbital strings array.
 */
int generate_virt_strings(int **virt_array, int nstr, int orbs, int elec,
                          int ndocc, int nactv, int xlvl);

/*
 * max_actv_space_strings: return the maximum number of ACTV space strings.
 */
int max_actv_space_strings(int orbs, int elecs, int ndocc, int xlvl);

/*
 * max_docc_space_strings: return the maximum number of DOCC space strings.
 */
int max_docc_space_strings(int orbs, int elecs, int xlvl);

/*
 * string_number: generate number of possible *valid* strings given the number
 * of DOCC electrons, ACTV electrons, and VIRT electrons.
 */
int string_number(int nde, int nae, int nve, int ndocc, int nactv, int orbs,
                  int elecs);

#endif

