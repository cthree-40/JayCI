#ifndef straddress_h_
#define straddress_h_

void str_adr2str(int index,
		 int *scr,
		 int elec,
		 int orbs,
		 int *str);

int str_adrfind(int *str,
		int elec,
		int orbs);

/* str_adrfind_fast: Computes address of orbital index string
 * -------------------------------------------------------------------
 * WARNING: binomial_coef3 does not check for initialization of binom data.
 *
 * Input:
 *  str  = orbital index string
 *  elec = electron number
 *  orbs = MO's in system */
int str_adrfind_fast(int *str, int elec, int orbs);

void str_strfind1(int *str1,
		 int elec,
		 int orbs,
		 int *str2);

void str_strfind2(int *str1,
		 int indx1,
		 int elec,
		 int orbs,
		 int indx2,
		 int *str2);

#endif
