#ifndef straddress_h_
#define straddress_h_

int str_adrfind(int *str,
		int elec,
		int orbs);

int str_strfind1(int *str1,
		 int elec,
		 int orbs,
		 int *str2);

int str_strfind2(int *str1,
		 int indx1,
		 int elec,
		 int orbs,
		 int indx2,
		 int *str2);

#endif
