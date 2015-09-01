/*
 *  KMatrix.c
 *  
 *  Created by Peter Dunning on 25/02/2007.
 *
 */
 
#include "KMatrix.h"
#include "ls_types.h"
#include <stdio.h>

/*Function that assemebles the element stiffness matricies into the Index arrays used by the MA57 solver*/
void Assemble2(int begin, int nNod, int *tnodes, double *KE, int *irn, int *jcn, double *A)
{
int start,start2,row,col,n,m,rn,cm,num,p,q,i,j,o,temp; /*Incrementors*/
num = begin;
int dof = 6 * nNod; /*dof for 2D shell element*/

	for(i=0;i<nNod;i++)
	{
		for(j=i;j<nNod;j++)
		{
			p = (tnodes[i] <= tnodes[j]) ? i : j;
			q = (tnodes[j] <= tnodes[i]) ? i : j;	/*Ensure upper triangle entries are read in*/
			rn = 6 * p;
			cm = 6 * q; /*element matrix row and column numbers -> for upper triangle*/
			row = tnodes[p] * 6;
			col = tnodes[q] * 6; /*Global matrix row and column numbers -> for upper triangle*/
			/*printf("\np=%i, q=%i, rn=%i, cm=%i, row=%i, col=%i",p,q,rn,cm,row,col);*/
			start = 0; 
			for(n=0;n<6;n++)
			{
				for(m=start;m<6;m++)
				{
					irn[num] = row + n + 1;
					jcn[num] = col + m + 1;
					temp = (dof * (n + rn)) + m + cm;
					A[num++] = KE[temp];
					/*printf("\nirn=%i, jcn=%i, A=%lf",irn[num-1],jcn[num-1],A[num-1]);*/
				}
				start2 = (i == j) ? 1 : 0;	/*Ensures only upper triangle for diagonal 6x6 blocks are read in*/
				start += start2;
				/*printf("\t Start = %i", start);*/
			}
		}
	}
/*printf("\n Start = %i", start);
printf("\n Start = %i", start);
double X,Y;
	X = 7/0;*/
}
