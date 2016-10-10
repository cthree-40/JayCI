// File: execute_pci_calculation.h

#ifndef execute_pci_calculation_h
#define execute_pci_calculation_h

/*
 * execute_pci_calculation: main driver for parallel CI calculation.
 */
void execute_pci_calculation(int aelec, int belec, int orbs, int nastr,
			     int nbstr, int ndets, int ndocc, int nactv,
			     double *moints1, double *moints2, double nucrep,
			     double frzcore, int plvl);


#endif
