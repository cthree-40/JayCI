// File: run_pdycicalc.h

#ifndef run_pdycicalc_h
#define run_pdycicalc_h

/*
 * run_pdycicalc: Execute CI dyson orbital calculation.
 */
int run_pdycicalc ();

/*
 * check_wavefunction_input: check user input of global wavefunction variables.
 * This is executed on all processes.
 */
int check_wavefunction_input(int nelecs0, int norbs0, int nfrzc0,  int ndocc0,
                             int nactv0,  int nfrzv0, int nelecs1, int norbs1,
                             int nfrzc1,  int ndocc1, int nactv1,  int nfrzv1);

/*
 * determinant_string_info: compute string information given
 * a determinant.
 * Input:
 *  detindx = determinant index (position in list)
 *  peosp   = alpha electron orbital spaces
 *  pegrps  = number of alpha electron orbital spaces
 *  qeosp   = beta  electron orbital spaces
 *  qegrps  = number of beta  electron orbital spaces
 *  pq      = (p,q)-space pairings
 *  npq     = number of (p,q)-space pairings
 * Output:
 *  pqval = (p,q) pairing start
 *  pval  = p string of (p,q)-space pairing
 *  qval  = q string of (p,q)-space pairing
 */
void determinant_string_info(int detindx, struct eospace *peosp, int pegrps,
                             struct eospace *qeosp, int qegrps, int **pq,
                             int npq, int *pqval, int *pval, int *qval);

/*
 * generate_det_triples: generate list of triplets for each determinant:
 *  |i> = |(pq, p, q)>
 */
void generate_det_triples (int ndeti, int **d_triplet, int pq_start,
                           int p_start, int q_start, int pq_final,
                           int p_final, int q_final, int **pq, int npq,
                           struct eospace *peosp, int pegrps,
                           struct eospace *qeosp, int qegrps);

/*
 * generate_wlist: generate the wavefunction list of triplets
 * Input:
 *  hndl     = handle of global array of W
 *  ndets    = number of determinants
 *  pq       = alpha/beta pairs
 *  npq      = number of alpha/beta pairs
 *  peospace = alpha string electron orbital space
 *  pegrps   = number of alpha string orbital spaces
 *  qeospace = beta  string electron orbital space
 *  qegrps   = number of beta  string orbital spaces
 */
void generate_wlist(int hndl, int ndets, int **pq, int npq,
                    struct eospace *peosp, int pegrps,
                    struct eospace *qeosp, int qegrps);

        
#endif
