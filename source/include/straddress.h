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
