/*
 *  Solve.h
 *  
 *  Created by Peter Dunning on 26/02/2007.
 *
 */

/*Function that sets up and runs the MA57 Multifrontal Solver*/
void solve(int n, int ne, int *irn, int *jcn, double *A, int n_rhs, double *rhs, int pinfo);

/*Function that uses MA44 to solve a overdtermined set of equations by the least squares method*/
void leastsq(int m, int n, int m1, double *A, double *B, double *Sol, int pinfo);

/*Function that sets up and runs the MA57 Multifrontal Solver*/
void solve2(int n, int s, int *irn, int *jcn, double *A, double *rhs);
