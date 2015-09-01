/*
 *  Solve.c
 *  
 *  Created by Peter Dunning on 26/02/2007.
 *
 */

#include "Solve.h"
#include <stdlib.h>
#include <stdio.h>

/*Function that sets up and runs the MA57 Multifrontal Solver*/
void solve(int n, int ne, int *irn, int *jcn, double *A, int n_rhs, double *rhs, int pinfo)
{
int go;
int inc; /*incrementor*/

/*printf("\nMatrixOrder= %i, NumEnt=%i, n_rhs=%i",n,ne,n_rhs);
for(inc=0;inc<ne;inc++)
{
	printf("\nirn[%i]=%i",inc,irn[inc]);
	printf("\tjcn[%i]=%i",inc,jcn[inc]);
	printf("\tA[%i]=%fl",inc,A[inc]);
}
for(inc=0;inc<(n * n_rhs);inc++)
{
	printf("\n rhs[%i]=%fl",inc,rhs[inc]);
}*/

/*Lots of variables for MA57 fortran multi-frontal solver*/
int lkeep,lfact,lifact,licntl,lcntl,linfo,lrinfo,lrhs,lwork,liwork;
int *keep,*icntl,*info,*ifact,*iwork;
double *rinfo,*fact,*cntl,*work;

/*set arrays for control data*/
lcntl = 5;
licntl = 20;
icntl = malloc(licntl * sizeof(int));
cntl = malloc(lcntl * sizeof(double));

/*initalise  control data arrays*/
ma57id_(cntl,icntl);
/*for(inc=0;inc<licntl;inc++)
{
	printf("\n icntl[%i]=%i",inc,icntl[inc]);
}*/

icntl[4] = pinfo; /*info printout on screen, 0 = none, 4 = Full*/
/*icntl[8] = 10; /*3 runs of IR*/
icntl[5] = 5; /*Metis ordering*/
icntl[7] = 1; /*don't restart if space for factorization runs out*/

/*setup infomation arrays*/
linfo = 40;
info = malloc(linfo * sizeof(int));

lrinfo = 20;
rinfo = malloc(lrinfo * sizeof(double));

int job = 0; /*solve for 1 array on rhs i.e. Ax=b*/

/*setup workspace arrays*/
lkeep = 7*n+2*ne+42;
keep = malloc(lkeep * sizeof(int));

liwork = 5 * n;
iwork = malloc(liwork * sizeof(int));

/*Do Analysis*/
ma57ad_(&n,&ne,irn,jcn,&lkeep,keep,iwork,icntl,info,rinfo);

/*re-size integer work array for next phase*/
free(iwork);
iwork = malloc(n * sizeof(int));

/*Setup Factorization arrays*/
lfact = 2 * info[8]; 
fact = malloc(lfact * sizeof(double));

lifact = 2 * info[9]; 
ifact = malloc(lifact * sizeof(int));

/* Do factorization */
ma57bd_(&n,&ne,A,fact,&lfact,ifact,&lifact,&lkeep,keep,iwork,icntl,cntl,info,rinfo);

/*Setup doubleing point workspace*/
lwork = n * n_rhs;
work  = malloc(lwork * sizeof(double));

 /* do solve. NB: lrhs = n*/
lrhs = n;
ma57cd_(&job,&n,fact,&lfact,ifact,&lifact,&n_rhs,rhs,&lrhs,work,&lwork,iwork,icntl,info);

/*if(pinfo > 0)
{
	printf("\nLargest scaled residule (1) = %fl", rinfo[5]);
	printf("\nLargest scaled residule (2) = %fl", rinfo[6]);
}*/

/*ma57d_(&job,&n,&ne,A,irn,jcn,fact,&lfact,ifact,&lifact,rhs,xxx,resid,work,iwork,icntl,cntl,info,rinfo);*/

/* free workspace */
	free(keep);
	free(rinfo);
	free(icntl);
	free(cntl);
	free(iwork);
	free(work);
	free(info);
	free(fact);
	free(ifact);
}

/*Function that uses MA44 to solve a overdtermined set of equations by the least squares method*/
void leastsq(int m, int n, int m1, double *A, double *B, double *Sol, int pinfo)
{
	/*printf("\nleastsq called");*/
	/*setup workspace*/
	int *iw;
	iw = malloc(n * sizeof(int));
	
	double *w;
	w = malloc(2 * n * sizeof(double));
	
	double err; /*error message varible*/
	
	/*NB: Modified MA44 source code to remove boolean variables, replaced with simple integer*/
	int res = 0;
	
	/*printf("\ntime for MA44...");*/
	
	ma44ad_(&m,&n,&m1,A,&m,B,&res,Sol,iw,w,&err,&pinfo,&pinfo);
	
	if(err < 0.0)
	{
		printf("\nError in MA44: %lf",err);
	}
	
	free(iw);
	free(w);
}
