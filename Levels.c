/*
 *  Levels.c
 *  
 *
 *  Created by Peter Dunning on 17/02/2009.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ls_types.h"
#include "Levels.h"

/*Function to calculate active and mine nodes for the narrow band method*/
void NarBand(int NumNodes, double h, double lBand, int *NumAct, int *NumMine, double *lsf, int *active, int *mine, int NumFix, int *fixed, int *activeHoles, int *NumActHoles)
{
	*NumAct = 0;
	*NumMine = 0;
	*NumActHoles = 0; /*re-set these counts*/
	
	int i,j,k,count; /*Incrementors*/
	int temp;
	double ftemp; /*temporary variables*/
	/*double lBand = 4.0 * h; /*use a fixed bandwidth*/
	double lMine = lBand - h;
	
	for(i=0;i<NumNodes;i++)
	{
		ftemp = fabs(lsf[i]); /*absolute value of the signed distance function*/
		
		/*search fixed array to see if node should remain unchanged*/
		temp = InSet(i,fixed,NumFix);
		
		/*If node within narrow band then it is active - unless it is fixed!*/
		if( ((ftemp - lBand) < 0.000001) && (temp == 0))
		{
			j = *NumAct;
			active[j] = i;
			*NumAct = ++j;
			/*printf("\nNode %i",i+1);*/

			/*If node is also outside lMine then define it as a mine*/
			if((ftemp - lMine) > -0.000001)
			{
				j = *NumMine;
				mine[j] = i;
				*NumMine = ++j;
			}
		}
		if( ((ftemp - 1.74) < 0.000001) && (temp == 0))
		{
			j = *NumActHoles;
			activeHoles[j] = i;
			*NumActHoles = ++j;
			/*printf("\nNode %i",i+1);*/
		}
	}
	
	printf("\nNumAct = %i, NumMine = %i",*NumAct,*NumMine);
}

/*Function to calculates boundary length and strain enegy associated with boundary nodes*/
void LamCalc3a(double h, int elemX, int elemY, Elem Number[elemX][elemY], short *ElemStat, short *NodeStat, 
				Aux *auxNode, double *Nenergy, int NumNodes, Coord NodeCoord[NumNodes], int Ntot, Aux *Lbound, double hxz[elemX])
{
	int i,j,k,k2,count,num; /*incrementors*/
	int n1,n2; /*boundary point node numbers*/
	int *nodes; /*Grid Element node numbers*/
	nodes = malloc(4 * sizeof(int));
	double x1,x2,y1,y2; /*Co-ordinates*/
	int temp;
	double ftemp,len; /*tempary variables*/
	Coord *Ltemp;	/*tempory array to store boundary data*/
	Ltemp= calloc(Ntot,sizeof(Coord));
	
	count = 0; /*Initalise NIO element counter to zero*/
	
	for(j=0;j<elemY;j++)
	{
		for(i=0;i<elemX;i++)
		{
			num = Number[i][j].n -1; /*Read in number of current element*/
			
			/*read in node numbers*/
			nodes[0] = Number[i][j].a -1;
			nodes[1] = Number[i][j].b -1;
			nodes[2] = Number[i][j].c -1;
			nodes[3] = Number[i][j].d -1;
			
			/*---------------If element is NIO-------------*/
			if((ElemStat[num] != 4) && (ElemStat[num] != 0))
			{
				/*printf("\nElement %i on boundary",num+1);*/
				/*----------If element is Pentagonal----------*/
				if(ElemStat[num] == 3)
				{
					/*Read in auxillary node numbers and co-ordinates*/
					temp = 2 * count;
					n1 = auxNode[temp].n -1;
					n2 = auxNode[temp+1].n -1;
					x1 = auxNode[temp].xz;
					x2 = auxNode[temp+1].xz;
					y1 = auxNode[temp].y;
					y2 = auxNode[temp+1].y;
				}
				
				/*----------If element is Quadrilateral----------*/
				else if(ElemStat[num] == 2)
				{
					/*Read in auxillary node number and co-ordinates*/
					temp = 2 * count;
					n1 = auxNode[temp].n -1;
					x1 = auxNode[temp].xz;
					y1 = auxNode[temp].y;
					/*count number of OUT nodes*/
					k2 = 0;
					for(k=0;k<4;k++)
					{
						k2 += (NodeStat[nodes[k]] == 0) ? 1 : 0;
					}
					
					if(k2 == 2) /*if there are 2 OUT nodes then there should be two aux nodes logged*/
					{
						n2 = auxNode[temp+1].n -1;
						x2 = auxNode[temp+1].xz;
						y2 = auxNode[temp+1].y;
					}
					
					else /*otherwise the other aux node is one of the boundary nodes*/
					{
						for(k=0;k<4;k++)
						{
							k2 = (k > 1) ? k-2 : k+2; /*Find opposite node*/
							if((NodeStat[nodes[k]] == 2) && (NodeStat[nodes[k2]] == 1)) /*If current node on boundary and is opposite an IN node*/
							{
								n2 = nodes[k];
								x2 = NodeCoord[n2].xz;
								y2 = NodeCoord[n2].y;
							}
						}
					}
				}
				
				/*----------If element is Triangular----------*/
				else if(ElemStat[num] == 1)
				{
					/*count number of Boundary nodes (excluding node opposite IN node)*/
					k2 = 0;
					for(k=0;k<4;k++)
					{
						temp = (k < 2) ? (k+2) : (k-2); /*locate opposite node*/
						k2 += ((NodeStat[nodes[k]] == 2) && (NodeStat[nodes[temp]] != 1)) ? 1 : 0; /*increase count if opposite node isn't IN*/
					}
					
					/*If there are two boundary nodes then bounary length is sqrt(2)xh */
					if(k2 == 2)
					{
						/*Find node numbers of boundary nodes*/
						for(temp=0;temp<4;temp++)
						{
							if(NodeStat[nodes[temp]] == 2)
							{
								n1 = nodes[temp];
								break;
							}
						}
						for(k=temp+1;k<4;k++)
						{
							if(NodeStat[nodes[k]] == 2)
							{
								n2 = nodes[k];
								break;
							}
						}
						
						/*values set so that length is calculated as: sqrt(2) x h */
						x1 = hxz[i];
						y1 = h;
						x2 = 0.0;
						y2 = 0.0;
					}
					
					else if(k2 == 1)
					{
						/*Read in auxillary node numbers and co-ordinates*/
						temp = 2 * count;
						n1 = auxNode[temp].n -1;
						x1 = auxNode[temp].xz;
						y1 = auxNode[temp].y;
						
						/*Now find the boundary node and read in data*/
						for(k=0;k<4;k++)
						{
							if(NodeStat[nodes[k]] == 2)
							{
								/*printf("\nIn tri Node %i is on Boundary -- %i",nodes[k]+1,NodeStat[nodes[k]]);*/
								n2 = nodes[k];
								x2 = NodeCoord[n2].xz;
								y2 = NodeCoord[n2].y;
								break;
							}
						}
					}
					
					else if(k2 == 0)
					{
						/*Read in auxillary node numbers and co-ordinates*/
						temp = 2 * count;
						n1 = auxNode[temp].n -1;
						n2 = auxNode[temp+1].n -1;
						x1 = auxNode[temp].xz;
						x2 = auxNode[temp+1].xz;
						y1 = auxNode[temp].y;
						y2 = auxNode[temp+1].y;
					}
				}
					
				/*---------Now update the boundary data array-------*/

				/*printf("\nn1=%i, n2=%i, x1=%lf, x2=%lf, y1=%lf, y2=%lf",n1,n2,x1,x2,y1,y2);*/
				/*Calculate distance between auxillary nodes (using pythag)*/
				ftemp = x1 - x2;
				len = ftemp * ftemp;
				ftemp = y1 - y2;
				len += ftemp * ftemp;
				ftemp = len;
				len = 0.5 * sqrt(ftemp); /*only need half length for each update*/
				
				Ltemp[n1].x += len;
				Ltemp[n2].x += len; /*update bounary lengths*/
				/*printf("\nStn Add = %lf + %lf",Nenergy[n1], Nenergy[n2]);*/
				ftemp = (Nenergy[n1] + Nenergy[n2]) * len * 0.5;
				Ltemp[n1].y += ftemp;
				Ltemp[n2].y += ftemp; /*Update boundary strain energy*/
				
				count++; /*update NIO element counter*/
			}
			
			/*If element isn't out check all edges to see if any lie completely on the boundary*/
			if(ElemStat[num] != 0)
			{				
				for(k=0;k<4;k++)
				{
					k2 = (k == 3) ? 0 : k+1;
					if((NodeStat[nodes[k]] == 2) && (NodeStat[nodes[k2]] == 2)) /*if current edge lies on the boundary*/
					{
						/*printf("\nElement %i, Edge %i on boundary",num+1,k+1);*/		
						if((k==0)||(k==2))
						{
							Ltemp[nodes[k]].x += hxz[i] * 0.5;
							Ltemp[nodes[k2]].x += hxz[i] * 0.5; /*update bounary lengths*/
							ftemp = (Nenergy[nodes[k]] + Nenergy[nodes[k2]]) * hxz[i] * 0.25;
						}
						else
						{
							Ltemp[nodes[k]].x += h * 0.5;
							Ltemp[nodes[k2]].x += h * 0.5; /*update bounary lengths*/
							ftemp = (Nenergy[nodes[k]] + Nenergy[nodes[k2]]) * h * 0.25;
						}
						/*printf("\nStn Add = %lf + %lf",Nenergy[n1], Nenergy[n2]);*/
						Ltemp[nodes[k]].y += ftemp;
						Ltemp[nodes[k2]].y += ftemp; /*Update boundary strain energy*/
					}
				}
			}
		}
	}
	
	/*Now reduce Lbound array to required size*/
	count = 1;
	for(i=0;i<Ntot;i++)
	{
		if(Ltemp[i].x > 0.000001) /*if associated boundary length is > than small value then include*/
		{
			Lbound[count].n = i;
			Lbound[count].x = Ltemp[i].x;
			Lbound[count++].y = Ltemp[i].y;
		}
		
		else if( (i < NumNodes) && (NodeStat[i] == 2) ) /*if node designated a boundary node, then include anyway */
		{
			Lbound[count].n = i;
			Lbound[count].x = Ltemp[i].x;
			Lbound[count++].y = Ltemp[i].y;
		}
		
		else if(fabs(Ltemp[i].x) > 0.0000001) /*otherwise print length that is ignored*/
		{
			printf("\nBoundary length at node %i = %lf and is ignored", i+1, Ltemp[i].x);
		}
	}
	
	Lbound[0].n = count; /*store length of Lbound array as first entry of array*/
	Lbound = realloc(Lbound, count * sizeof(Aux));
	
	free(nodes);
	free(Ltemp);
	
	/* Write Boundary data information file*/
	FILE *outfile;	/*File varible for output files*/
	outfile = fopen("Boundary.txt", "w");
	if(outfile == NULL){
		printf("\nFailed to open Boundary writefile\n");
		}
	else{
	fprintf(outfile,"node\tUn\n"); /*column headings*/
		for(i=1;i<count;i++)
		{
			fprintf(outfile,"%i\t",Lbound[i].n+1);
			fprintf(outfile,"%lf\t",Lbound[i].x);
			fprintf(outfile,"%lf\t",Lbound[i].y);

			fprintf(outfile,"\n");
		}
	}

	fclose(outfile);

	printf("\nBoundary data File written\n");
}

void ArrayCopy(int length, double *original, double *copy)
{
	int i;
	for(i=0;i<length;i++)
	{
		copy[i] = original[i];
	}
}

double LamCalc4(double lambda, double dt, int NumNodes, double *Nenergy, double *lsf, double *lsfHole, int NumFix, int *fixed,
					int NumAct, int *active, double h, double hbar, double HoleLim, double Voltot, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], short *NodeStat, int elemX, double hxz[elemX])
{
	double vol1,vol2,lam1,lam2;
	double *lsf_temp;
	lsf_temp = malloc(NumNodes * sizeof(double)); /*temp array for lsfHole*/
	ArrayCopy(NumNodes, lsfHole, lsf_temp);
	
	lam1 = lambda;
	/*first clalcualte hole insertion volume reduction using bounary lambda value*/
	vol1 = LamCalc4a(lam1, dt, NumNodes, Nenergy, lsf, lsf_temp, NumFix, fixed, NumAct, active, h, hbar, NodeX, NodeY, Nodes2, NodeStat, elemX, hxz);
	
	double Volin = HoleLim * Voltot; /*set a hole creation limit of 1% current total volume*/
	double VolLow = h * h; /*set lower limit for volume reduction (as one element), or one square element the smallest an element can be*/
	/*double VolLow = 0.001 * Voltot; /*set lower limit for volume reduction as 0.1% of total volume */
	printf("\nHole volume reduction limit = %lf",Volin);
	/*printf("\nInital hole volume = %lf", vol1);*/
	
	if((vol1 - Volin) < 0.000001) /*if hole volume less than limit, return lambda, or zero if no holes are created*/
	{
		ArrayCopy(NumNodes, lsf_temp, lsfHole);
		free(lsf_temp);
		lam1 = (vol1 < VolLow) ? 0.0 : lam1; /*if volume reduction less than limit, then return zero lamT value*/
		printf("\nAfter 0 iterations Hole Lambda = %lf, vol=%lf",lam1,vol1);
		return(lam1);
	}
	
	/*if iniital volume less than minimum, then set target to minimum
	if(vol1 < VolLow)
	{
		Volin = VolLow;
		VolLow = h * h;
	}*/

	double grad, inter, ftemp;
	int i,count;
	
	/*set intial limits on lambda based on max nodal strain eneergy*/
	double lamlim, lammax, lammin;
	
	lamlim = 0.0;
	
	for(i=0;i<NumNodes;i++)
	{
		if((InSet(i,fixed,NumFix) == 0) && (InSet(i,active,NumAct) == 0) && (lsf[i] > -0.000001))
		{
			ftemp = Nenergy[i];
			lamlim = (ftemp > lamlim) ? ftemp : lamlim;
		}
	}
	
	lammin = 0.0;	/*set to zero to ensure volume reduction, i.e +ve velocity*/
	lammax = (lamlim > lambda) ? lamlim : lambda;	/*set to maximum internal strain energy, as this would mean the entire inside disappears*/
	/*printf("\nlamlim = %lf", lamlim);*/
	
	/*update min and max values if necessary*/
	if((vol1 - Volin) < 0.000001)
	{
		lammin = (lam1 > lammin) ? lam1 : lammin;
	}
	else if((vol1 - Volin) > - 0.000001)
	{
		lammax = (lam1 < lammax) ? lam1 : lammax;
	}
	
	lam2 = (vol1 > Volin) ? (0.5 * lam1) : (2.0 * lam1); /*choose intial value for lam2 based on vol1*/
	
	count = 0;
	do {
			free(lsf_temp);
			lsf_temp = malloc(NumNodes * sizeof(double)); /*temp array for lsfHole*/
			ArrayCopy(NumNodes, lsfHole, lsf_temp);

			vol2 = LamCalc4a(lam2, dt, NumNodes, Nenergy, lsf, lsf_temp, NumFix, fixed, NumAct, active, h, hbar, NodeX, NodeY, Nodes2, NodeStat, elemX, hxz);
			printf("\nlam1=%lf, lam2=%lf, vol1=%lf, vol2=%lf",lam1,lam2,vol1,vol2);
			
			/*work out percentage difference of current guess from target (Volin) */
			ftemp = vol2 - Volin;
			ftemp /= Volin;
			
			count++;
			/*if volume chnage close to target, then break the loop*/
			/*printf("\nftemp=%lf",ftemp);*/
			if(fabs(ftemp) < 0.01)
			{
				printf("\nAfter %i iterations Hole lambda = %lf", count, lam2);
				ArrayCopy(NumNodes, lsf_temp, lsfHole);
				free(lsf_temp);
				return(lam2); /*if converged then return lambda*/
			}
			
			/*update min and max values if necessary*/
			else if((vol2 - Volin) < 0.000001)
			{
				lammin = (lam2 > lammin) ? lam2 : lammin;
			}
			else if((vol2 - Volin) > -0.000001)
			{
				lammax = (lam2 < lammax) ? lam2 : lammax;
			}
			
			printf("\nlammax = %lf, lammin = %lf", lammax, lammin);
	
			/*work out difference between lammin and lammax*/
			ftemp = lammax - lammin;
			ftemp /= 0.5 * fabs(lammax + lammin);
			
			if(ftemp < 0.000001)
			{
				if((vol2 < VolLow) || (vol2 > Volin))
				{
					free(lsf_temp);
					printf("\nERROR! LamT not converged!");
					return(0.0); /*if not converged on an aceptable volume, then do not update lsfHole, and return lamT = 0.0*/
				}
				
				else /*otherwise volume is acceptable*/
				{
					ArrayCopy(NumNodes, lsf_temp, lsfHole);
					free(lsf_temp);
					return(lam2); /*if converged on a volume smaller than target, then return accepted valuey*/
				}
			}
			
			/*otherwise find the next guess by a linear search*/
			grad = vol2 - vol1;
			grad /= lam2 - lam1; /*calculate gradient*/
			
			inter = vol2 - (grad * lam2); /*calcualte intercept*/
			
			/*update previous lambda & volume values*/
			lam1 = lam2;
			vol1 = vol2;
			
			ftemp = grad / (h * h);
			ftemp /= h * h;
			
			if(fabs(ftemp) > 0.0001) /*if gradient isn't zero*/
			{
				lam2 = (Volin - inter) / grad; /*lambda value of next guess*/
			}
			else
			{
				lam2 = 0.5 * (lammax + lammin);
			}
			
			/*check for limits, if violated, use bi-section update*/
			if((lam2 < lammin) || (lam2 > lammax))
			{
				lam2 = 0.5 * (lammax + lammin);
			}
				
		} while (count < 25); /*limit 25 interations*/
		
		if(count == 25)
		{
			/*compute volume for lammin*/
			vol2 = LamCalc4a(lammin, dt, NumNodes, Nenergy, lsf, lsf_temp, NumFix, fixed, NumAct, active, h, hbar, NodeX, NodeY, Nodes2, NodeStat, elemX, hxz);

			if(vol2 < Volin) /*if volume for lammmin < target volume, then return lammin*/
			{
				ArrayCopy(NumNodes, lsf_temp, lsfHole);
				free(lsf_temp);
				return(lammin);
			}
			printf("\nERROR! LamT not converged, due to iteration limit, lsfHole not updated");
			free(lsf_temp);
			return(0.0);
		}
    return (0.0);
}

/*function to calcualte volume reduction due to hole insertion*/
double LamCalc4a(double lambda, double dt, int NumNodes, double *Nenergy, double *lsf, double *lsfHole, int NumFix, int *fixed,
					int NumAct, int *active, double h, double hbar, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], short *NodeStat, int elemX, double hxz[elemX])
{
	int i,j,k,temp,num,count;
	double ftemp, xtemp;
	
	double *lsf_temp;
	lsf_temp = malloc(NumNodes * sizeof(double));
	
	/*Run algorithm to calcualte volume change due to hole creation*/
	k = 0;
	for(i=0;i<NumNodes;i++)
	{
		/*if((InSet(i,fixed,NumFix) == 0) && (InSet(i,active,NumAct) == 0) && (lsfHole[i] > -h))*/
		if( (InSet(i,fixed,NumFix) == 0) && (NodeStat[i] == 1) )
		{
			ftemp = lambda - Nenergy[i]; /*calcualte velocity*/
			ftemp *= dt; /*needs to be realtive to boundary movement*/
			lsf_temp[i] = lsfHole[i] - ftemp; /*update lsf_Hole*/
			
			/*if hole created indicate this for later*/
			if( (lsf_temp[i] < -0.000001) && (InSet(i,active,NumAct) == 0) )
			{
				k = 1;
			}
			/*restrict maximum & minimum lsfHole value*/
			if(fabs(lsf_temp[i]) > hbar)
			{
				lsf_temp[i] = (lsf_temp[i] > 0.0) ? hbar : -hbar;
			}
		}
		else {
			lsf_temp[i] = lsfHole[i];
		}

	}
	
	if(k == 0) /*if no holes are created, return zero volume reduction*/
	{
		ArrayCopy(NumNodes,lsf_temp,lsfHole); /*copy back the updated lsf_Hole array*/
		free(lsf_temp);
		printf("\nk=0, return zero volume");
		return(0.0);
	}
	
	/*-------------smooth out the lsf_hole function-------------*/
	double *lsf_smooth;
	lsf_smooth = malloc(NumNodes * sizeof(double));
	double *lsfVals;
	lsfVals = malloc(4 * sizeof(double)); /*array to store lsfHole values for surrounding nodes*/
	
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++)
		{
			num = Nodes2[i][j]; /*read in node number*/
			/*lsf_smooth[num] = (InSet(num,active,NumAct) == 0) ? lsf_temp[num] : lsf[num];*/
			
			/*If node part of lsf_Hole*/
			if( (NodeStat[num] == 1) && (InSet(num,fixed,NumFix) == 0) 
					&& (InSet(num,active,NumAct) == 0) && (lsf_temp[num] < 0.000001)) /* */
			{			
				/*read in lsfHole value for neigbouring nodes (taking account of nodes inside narrow band) 
				temp = Nodes2[i-1][j];
				lsfVals[0] = (InSet(temp,active,NumAct) == 0) ? lsf_temp[temp] : lsf[temp];
				temp = Nodes2[i+1][j];
				lsfVals[1] = (InSet(temp,active,NumAct) == 0) ? lsf_temp[temp] : lsf[temp];
				temp = Nodes2[i][j-1];
				lsfVals[2] = (InSet(temp,active,NumAct) == 0) ? lsf_temp[temp] : lsf[temp];
				temp = Nodes2[i][j+1];
				lsfVals[3] = (InSet(temp,active,NumAct) == 0) ? lsf_temp[temp] : lsf[temp];

				lsf_smooth[num] = 0.25 * (lsfVals[0] + lsfVals[1] + lsfVals[2] + lsfVals[3]);*/
				
				lsf_smooth[num] = (InSet(num,active,NumAct) == 0) ? lsf_temp[num] : lsf[num];
			}
			else
			{
				lsf_smooth[num] = (InSet(num,active,NumAct) == 0) ? lsf_temp[num] : lsf[num];
			}
		}
	}
	
	double VolHole = 0.0;	/*initalise to zero*/
	/*now estimate volume reduction from updated lsfHole array*/
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++)
		{
			num = Nodes2[i][j]; /*read in node number*/
			
			/*If hole created at current node, then estimate volume reductiony*/
			if( (lsf_smooth[num] < -0.000001) && (NodeStat[num] == 1) 
					 && (InSet(num,active,NumAct) == 0) && (InSet(num,fixed,NumFix) == 0) ) /**/
			{		
				/*read in smoothed lsfHole values for neigbouring nodes*/
				lsfVals[0] = lsf_smooth[Nodes2[i-1][j]];
				lsfVals[1] = lsf_smooth[Nodes2[i+1][j]];
				lsfVals[2] = lsf_smooth[Nodes2[i][j-1]];
				lsfVals[3] = lsf_smooth[Nodes2[i][j+1]];
				
				count = 0;
				ftemp = 0.0;
				for(k=0;k<4;k++)
				{
					if(lsfVals[k] < -0.000001) {
						VolHole += 1.0; /*if neighbouring node is out*/
						}
					else if(lsfVals[k] > -0.000001) {
						xtemp = lsf_smooth[num] / (lsf_smooth[num] - lsfVals[k]); /*using relative distance*/
						count++;
						VolHole += xtemp; /*if neighbouring node is in, estimate volume by interpolation*/
						ftemp += xtemp;
						}
				}
				
				/*if((xtemp < 0.000001) || (count == 4)) /*if hole is just created at a single node, then don't create a hole here
				{
					VolHole -= xtemp;
					lsfHole[num] = 0.001 * h;
				}*/
				
				if(count == 4) /*if hole is just created at a single node, then don't create a hole here*/
				{
					lsf_temp[num] = 0.00001;
					VolHole -= ftemp;
					printf("\nHole is only at node %i!!!",num+1);
				}
				
				/*printf("\nVolHole = %lf",VolHole);*/
			}
		}
	}
	
	VolHole *= h * h * 0.25; /*scaling factor for volume calc*/
	ArrayCopy(NumNodes,lsf_temp,lsfHole); /*copy back the updated lsf_Hole array*/
	free(lsfVals);
	free(lsf_temp);
	free(lsf_smooth);
	return(VolHole);
}

/*Function to calcualte volume change for a given lambda value*/
double LamCalc5a(double lambda, double *Vnorm, Aux *Lbound, Aux *auxNode, double h, double alt, double *Nenergy,
					int NIOtotal, int NumNodes, Coord NodeCoord[NumNodes], double *VmaxA, int itt, double maxX, double maxY,
						double Blen_in, double Uave_in, int Ntot, int NumFix, int *fixed, Astn *EmaxA, int Hflag, int elemX, double hxz[elemX])
{
	int bsn = Lbound[0].n; /*read number of boundary nodes*/
	int i,j,num;
	int temp,temp2; /*continuation flags for loops*/
	double ftemp,xCrd,yCrd,a1,del1,del2,maxT,minT;
	double dt;	/*time step*/
	short *lamfix;	/*array to store additional fixed nodes for lambda calc*/
	lamfix = calloc(bsn,sizeof(short)); /*assign memory*/
	int itta = 0;
	
	double Vmax = -10000000.0;
	double Vmin = 10000000.0;
	
	/*calculate normal velocites on the boundary using Vn = L - U for all boundary points*/
	for(i=1;i<bsn;i++)
	{
		num = Lbound[i].n;
		if(InSet(num,fixed,NumFix) == 0) /*if node not fixed*/
		{
			Vnorm[num] = lambda - Nenergy[num]; /*normal velocity to boundary = lambda - strain energy*/
			/*find Vmax and Vmin*/
			Vmax = (Vnorm[num] > Vmax) ? Vnorm[num] : Vmax;
			Vmin = (Vnorm[num] < Vmin) ? Vnorm[num] : Vmin;
		}
	}
	
	/*printf("\tVmax = %f", Vmax);
	printf("\tVmin = %f", Vmin);*/
	/*find |V|max*/
	double Vabs = (fabs(Vmax) > fabs(Vmin)) ? fabs(Vmax) : fabs(Vmin);
	
	/*do loop to find time step*/
	do
	{
		temp = 0;
		dt = h / Vabs; /*Since h will be smaller then hxz[i] so I will use this as the limit*/
		dt *= alt;
		/*printf("\nVabs=%lf, dt=%lf",Vabs, dt);*/
		
		/*may need to re-set some velocities from previous iteration*/
		for(i=1;i<bsn;i++)
		{
			if(lamfix[i] == 1)
			{
				num = Lbound[i].n;
				Vnorm[num] = lambda - Nenergy[num]; 
			}
		}
		
		free(lamfix);
		lamfix = calloc(bsn,sizeof(short)); /*assign memory*/
		
		temp2 = 0;
		/*check that Vabs is unaffected by dt value (mainly for boundary points near domain boundary)*/
		for(i=1;i<bsn;i++)
		{
			num = Lbound[i].n;
			if(InSet(num,fixed,NumFix) == 0) /*if node not fixed*/
			{
			/*first find x & y co-ords of current node*/
				if(num >= NumNodes)
				{
					for(j=0;j<(2 * NIOtotal);j++)
					{
						if(auxNode[j].n == num+1)
						{
							xCrd = auxNode[j].xz;
							yCrd = auxNode[j].y;
							break;
						}
					}
				}
				else
				{
					xCrd = NodeCoord[num].xz;
					yCrd = NodeCoord[num].y;
				}
				
				/*printf("\nNode %i, x=%f, y=%f",num+1,xCrd,yCrd);
				printf("\tMaxX = %f", maxX);*/
				
				/*Now work out if that node lies within one grid spacing of the boundary of the domain*/
				if((xCrd < hxz[0]) || (xCrd > (maxX - hxz[elemX-1])) || (yCrd < h) || (yCrd > (maxY - h)))
				{
					if(Vnorm[num] < 0.0) /*if velocity -ve (i.e. boundary moves outward)*/
					{
						/*find shortest distance from node to domain boundary*/
						a1 = maxX-xCrd;
						del1 = (xCrd < a1) ? xCrd : a1;
						a1 = maxY - yCrd;
						del2 = (yCrd < a1) ? yCrd : a1;
						
						a1 = (del1 < del2) ? del1 : del2;

						/*printf("\n a1 = %f", a1);
						printf("\t Vnorm[%i] * dt = %f X %f = %f", num, Vnorm[num], dt, Vnorm[num] * dt);*/
					
						if((Vnorm[num] * dt) < (-a1-0.0001))	/*if boundary would move beyond domain, then set velocity to only reach boundary*/
						{
							/*printf("\tLamfix = 1");*/
							Vnorm[num] = (-a1 / dt); /*set velocity so that boundary only extends to domain edge*/
							lamfix[i] = 1;
							/*printf("\nremove node %i",num+1);*/
							temp2 = 1; /*signal a change has occured*/
							/*printf("\nNode %i on boundary, Changed Vnorm = %lf",num+1,Vnorm[num]);*/
						}
					}
				}
			}

		}
		/*printf("\ntemp2 = %i", temp2);*/
		
		/*check to see if Vabs has changed*/
		if(temp2 == 1) /*if a change has occured to vnorm*/
		{
			minT = 10000000.0;
			maxT = -10000000.0; /*initalise temp max.min values*/
			
			for(i=1;i<bsn;i++)
			{
				/*printf("\tA");*/
				num = Lbound[i].n;
				/*printf("\tnum = %i, %i, %i", num, lamfix[i], InSet(num,fixed,NumFix));*/
				/*check possible new Vmin vlaue*/
				if((InSet(num,fixed,NumFix) == 0) && (lamfix[i] == 0)) /*if node not fixed and not constrained on boundary*/
				{
					maxT = (Vnorm[num] > maxT) ? Vnorm[num] : maxT;
					minT = (Vnorm[num] < minT) ? Vnorm[num] : minT; /*update temp max/min values*/
					/*printf("\nVnorm[%i] = %f", num, Vnorm[num]);
					printf("\nMaxT = %f \t MinT = %f", maxT, minT);*/
				}
			}
			/*printf("\nMaxT = %f \t MinT = %f", maxT, minT);*/
			
			ftemp = (fabs(maxT) > fabs(minT)) ? fabs(maxT) : fabs(minT); /*calc current Vabs value*/
			a1 = (ftemp - Vabs) / Vabs; /*work out % change of Vabs*/
			/*printf("\nftemp = %f", ftemp);
			printf("\na1 = %f", a1);*/
			/*if Vabs has changed*/
			if(fabs(a1) > 0.01)
			{
				Vabs = ftemp; /*upadte Vabs*/
				temp = 1; /*flag to continue the loop*/
			}
			
		}
		/*printf("\ntemp = %i", temp);*/
		/*if(itta == 4)
		{	
			double X,Y;
			X = 7/0;
		}*/
		itta++;
	}while(temp == 1);
	
	/*printf("\nOUT OF THE LOOP!!!!");*/
	
	/*increase convergent by using non-linear map of velocity function
	if(Hflag == 0) /*only if not using hole creation method
	{
		NlMap(Ntot, Vnorm, Vabs);
	}
	/*Now calculate volume change due to boundary movement (including near domain edge movement)*/
	double VolB = 0.0;
	double VolA = 0.0;
	
	for(i=1;i<bsn;i++)
	{
		num = Lbound[i].n;
		if((InSet(num,fixed,NumFix) == 0))
		{				
			ftemp = Lbound[i].x * Vnorm[num];
			VolB += ftemp;
			VolA += fabs(ftemp);
		}
	}
	
	VolB *= dt; /*finally need to multiply by the time step*/	
	*VmaxA = Vabs; /*store absolute max boundary velocity*/
	
	free(lamfix);
	return(VolB);
}

/*Function to find lambda value by iteration*/
double LamCalc5b(double delVol, double *Vnorm, Aux *Lbound, int Ntot, Aux *auxNode, double h, Astn *EmaxA, double alt, int NumFix, int *fixed,
					double *Nenergy, int NIOtotal, int NumNodes, Coord NodeCoord[NumNodes], double *VmaxA, int itt, double maxX, double maxY, int Hflag, int elemX, double hxz[elemX])
{
	int temp,i,count; /*Incrementors*/
	double ftemp;
	double lam1, lam2; /*lambda values*/
	double vol0, vol1, vol2; /*volume reduction values*/
	double grad, inter; /*Gradient and intercept*/
	double *Vtemp;
	Vtemp = calloc(Ntot,sizeof(double));
	
	double Blen = 0.0; /*boundary length*/
	double Uave = 0.0; /*total boundary strain energy*/
	double lamlim = 0.0; /*initalise to zero*/
	
	/*calculate total boundary length and average boundary strain energy*/
	int bsn = Lbound[0].n; /*read number of boundary nodes*/
	for(i=1;i<bsn;i++)
	{
		/*first find out if node is fixed*/
		if(InSet(Lbound[i].n,fixed,NumFix) == 0) /*if not fixed*/
		{
			Blen += Lbound[i].x;
			Uave += Lbound[i].y;
			temp = Lbound[i].n;
			ftemp = Nenergy[temp];
			lamlim = (ftemp > lamlim) ? ftemp : lamlim;
			/*Emax = (Nenergy[num] > Emax) ? Nenergy[num] : Emax;
			Emin = (Nenergy[num] < Emin) ? Nenergy[num] : Emin;*/
		}
	}
	
	/*printf("\nBlen=%lf, Ub=%lf, Uave=%lf",Blen,Uave,(Uave/Blen));
	printf("\nlamlim = %lf",lamlim);*/
	
	ftemp = Blen * 0.2 * h * alt;	/*set a limit on total boundary movement, based on boundary length, grid spacing, & time step (h is smallest edge lenght so use that)*/
	ftemp *= (delVol > 0.0) ? 1.0 : -1.0; /*ensure correct sign*/
	/*printf("\nftemp = %f",ftemp);
	printf("\ndelVol = %f", delVol);*/
	double Volin = (fabs(delVol) < fabs(ftemp)) ? delVol : ftemp; /*change volume target if greater than limit*/
	
	printf("\nBoundary volume reduction target = %lf",Volin);
	
	/*calcualte volume change for inital guess*/
	double lam0; /*inital lambda value*/
	lam0 = (itt == 0) ? (Uave / Blen) : EmaxA[itt-1].e;  /*set inital guess as average boudary strain energy for 1st iteration*/
															/*Otherwise use previous value as inital guess*/
	/*printf("\nlam0 = %f", lam0);*/
	vol0 = LamCalc5a(lam0, Vtemp, Lbound, auxNode, h, alt, Nenergy, NIOtotal, NumNodes, NodeCoord, VmaxA, itt, maxX, maxY, Blen, Uave, Ntot, NumFix, fixed, EmaxA, Hflag, elemX, hxz);
	/*printf("\nVolin = %f", Volin);
	printf("\nVol0 = %f", vol0);*/

	/*work out percentage difference of inital guess from target (Volin) */
	ftemp = vol0 - Volin;
	ftemp /= Volin;
	/*printf("\nftemp = %f", ftemp);*/
	/*if volume chnage close to target, then break the loop*/
	if(fabs(ftemp) < 0.01)
	{
		/*copy the temp velocities to actual array after convergence*/
		ArrayCopy(Ntot, Vtemp, Vnorm);
		free(Vtemp);
		EmaxA[itt].e = lam0;		
		printf("\nAfter 0 iterations: Boundary lambda = %lf",lam0);
		
		return(lam0);
	}
	
	/*work out a limit on lambda
	for(i=0;i<Ntot;i++)
	{
		ftemp = Nenergy[i];
		lamlim = (ftemp > lamlim) ? ftemp : lamlim;
	}*/

	double lammin;;
	double lammax;

	do{
		lam1 = lam0; /*set first guess*/
		vol1 = vol0;
		lam2 = (vol1 < Volin) ? (lam0 * 2.0) : (lam0 * 0.5); /*set 2nd guess as average boudary strian energy modified depending on vol1*/
		
		lammin = (Volin > vol1) ? lam1 : -lamlim; /*if volume to be reduced then lambda must be +ve*/
		lammax = (Volin < vol1) ? lam1 : lamlim; /*if volume to be increased then lambda must be -ve*/
		/*printf("\nlammax = %f,  lammin = %f", lammax, lammin);*/
		
		count = 0;
		temp = 0;
		do{
			free(Vtemp);
			Vtemp = calloc(Ntot,sizeof(double)); /*re-set the Vtemp array*/
			/*printf("\nlam2 = %f", lam2);*/
			vol2 = LamCalc5a(lam2, Vtemp, Lbound, auxNode, h, alt, Nenergy, NIOtotal, NumNodes, NodeCoord, VmaxA, itt, maxX, maxY, Blen, Uave, Ntot, NumFix, fixed, EmaxA, Hflag, elemX, hxz);
			
			/*work out percentage difference of current guess from target (Volin) */
			ftemp = vol2 - Volin;
			ftemp /= Volin;
			/*printf("\nVolin = %f", Volin);
			printf("\nVol2 = %f", vol2);
			printf("\nftemp = %f", ftemp);*/
			count++;
			/*if volume chnage close to target, then break the loop*/
			if(fabs(ftemp) < 0.01)
			{
				temp = 1; /*signal the break*/
				break;
			}
			
			/*update min and max values if necessary*/
			else if((vol2 - Volin) < 0.000001)
			{
				lammin = (lam2 > lammin) ? lam2 : lammin;
			}
			else if((vol2 - Volin) > -0.000001)
			{
				lammax = (lam2 < lammax) ? lam2 : lammax;
			}
			
			/*printf("\nlam1=%lf, lam2=%lf, vol1=%lf, vol2=%lf",lam1,lam2,vol1,vol2);*/
			
			/*work out difference between lammin and lammax*/
			ftemp = lammax - lammin;
			ftemp /= 0.5 * fabs(lammax + lammin);
			
			if(ftemp < 0.000001) /*if lammax = lammin, then break the loop*/
			{
				printf("\nlammax(%lf) = lammin(%lf): break",lammax,lammin);
				break;
			}
			
			/*otherwise find the next guess by a linear search*/
			grad = vol2 - vol1;
			grad /= lam2 - lam1; /*calculate gradient*/
			
			inter = vol2 - (grad * lam2); /*calcualte intercept*/
			
			/*printf("\ngrad = %lf, inter = %lf",grad,inter);*/
			
			/*update previous lambda & volume values*/
			lam1 = lam2;
			vol1 = vol2;
			
			ftemp = grad / (h * h);
			ftemp /= h * h; /*modify gradient by grid spacing*/
			
			if(fabs(ftemp) > 0.0001) /*if gradient isn't zero*/
			{
				lam2 = (Volin - inter) / grad; /*lambda value of next guess*/
			}

			else /*otherwise use bi-section update*/
			{
				lam2 = 0.5 * (lammax + lammin);
			}
			
			/*check for limits, if violated, use bi-section method*/
			if((lam2 < lammin) || (lam2 > lammax))
			{
				/*printf("\nlam2=%lf, out of range, using bi-section");*/
				lam2 = 0.5 * (lammax + lammin);
			}
			/*printf("\nlammax = %f,  lammin = %f", lammax, lammin);*/
			
		} while(count < 20);
		
		if(temp == 0) /* if a natural break has not occured, then work out volume associated with lammin or lammax*/
		{
			/*free(Vtemp);
			Vtemp = calloc(Ntot,sizeof(double)); /*re-set the Vtemp array*/
			if(delVol > 0.0)
			{
				if(vol2 <= delVol) /*if volume reduces somewhere between 0 and delVol then accept lambda value*/
				{
					temp = 1;
				}
			}
			else
			{
				if(vol2 >= delVol) /*if volume increases somewhere between 0 and delVol then accept lambda value*/
				{
					temp = 1;
				}
			}
		}
		
		if(temp == 0) /*if solution still not found then reduce volume target for next loop*/
		{
			Volin *= 0.75;
			printf("\nReducing Volin to %lf:", Volin); /*NB: this shouldn't actually get used, but is here for stability reasons*/
		}
		
		if(fabs(Volin) < 0.000001)
		{
			temp = 1;
			printf("\nERROR! Lambda calculation aborted! Volin = %lf",Volin);
		}
		
	} while(temp == 0);
	
	/*copy the temp velocities to actual array after convergence*/
	ArrayCopy(Ntot, Vtemp, Vnorm);
	free(Vtemp);
	
	/*store lambda value*/
	EmaxA[itt].e = lam2;
	
	printf("\nAfter %i iterations: Boundary lambda = %lf",count,lam2);
	
	return(lam2);
}

/*Non-linear mapping funciton for velcoity to improve convergence*/
void NlMap(int Ntot, double *Vnorm, double Vmax)
{
	int i;
	double ftemp;
	
	for(i=0;i<Ntot;i++)
	{
		if(fabs(Vnorm[i]) > 0.000001)
		{
			/*ftemp = (Vmax - fabs(Vnorm[i])) / Vmax; /*normalized velocity difference*/
			/*Vnorm[i] *= exp(ftemp); /*increase intermediate velocites using the non-linear map*/
			ftemp = 2.0 - ( fabs(Vnorm[i]) / Vmax );
			Vnorm[i] *= ftemp;
		}
	}
}

/*Function to calculate extension velocities using the fast marching method*/
void Vext(int NodeX, int NodeY, int Nodes2[NodeX][NodeY], short *NodeStat, double *Vnorm, double *lsf, int *active, int NumAct, int NumNodes, int Ntot, 
			Coord NodeCoord[NumNodes], double h, Aux *AuxNodes, int sign, double dt, int elemX, double hxz[elemX])
{
int i,j,k,count,NS,temp; /*incrementors*/
int Nnum,Nnum2,Nnum3; /*node numbers*/
int Ci,Cj; /*variables to store x,y location of node with minimum lsf_temp values in trial set*/
double lsf_min; /*variable to track minimum lsf_temp value in trial set*/
double Af, Bf, Cf; /*variables for lsf_temp calculation by solution of a quadratic equn.*/
double vt1,vt2;	/*variables for tempary velocities*/
double ftemp, dtemp; /*tempory variables*/
double Lsum, Fsum; /*varibales to calcualte Vext by summing contributions from surrounding nodes*/

short *known;
known = calloc(NumNodes,sizeof(short));
short *trial;
trial = calloc(NumNodes,sizeof(short));
double *lsf_temp;
lsf_temp = calloc(NumNodes,sizeof(double));
double *Vtemp;
Vtemp = calloc(NumNodes,sizeof(double));

    
    
    
    
/*Initalise Known set*/
for(i=0;i<NumAct;i++)
{
    
   // Vtemp[i]=0;
	j = active[i];
    //printf("active %d", j);
	/*If velocity already determined or node is on the boundary, then add to known set*/
	if(  (fabs(Vnorm[j]) > 0.000001)||NodeStat[j] == 2)
	{
		known[j] = 1;
		lsf_temp[j] = fabs(lsf[j]); /*this should be boundary values i,e lsf = 0.0*/
		Vtemp[j] = Vnorm[j];
		//printf("\nNode %i = known",i+1);
       // printf("vtemp: %f \n", Vtemp[j]);
	}
}

short side = (sign > 0) ? 1 : 0;

/*Initalise Trial set*/
for(j=1;j<NodeY-1;j++)
{
	for(i=1;i<NodeX-1;i++)
	{
		Nnum = Nodes2[i][j]; /*read in node number*/
		ftemp = 1000.0 * lsf[Nnum];		
		/*printf("\nNode Num = %i",Nnum);*/
		/*If current node not in known set, but is in the narrow band active set and is the correct side of the boundary*/
		if((known[Nnum] != 1) && (NodeStat[Nnum] == side) && (InSet(Nnum,active,NumAct) == 1))	
		{
			if(NodeStat[Nodes2[i-1][j]] != side) {
				trial[Nnum] = 1;
				}
			else if(NodeStat[Nodes2[i+1][j]] != side) {
				trial[Nnum] = 1;
				}
			else if(NodeStat[Nodes2[i][j-1]] != side) {
				trial[Nnum] = 1;
				}
			else if(NodeStat[Nodes2[i][j+1]] != side) {
				trial[Nnum] = 1;
				}
		}
	}
}

printf("\nInital Trial Set");
    int kkk=0;
for(i=0;i<NumNodes;i++)
{
	if(trial[i] == 1)
	{
		//printf("\nNode %i",i+1);
        kkk++;
	}
    //printf("%d", active[i]);
}
//printf("kkk is %d",kkk);
    

int stop = 4*(Ntot-NumNodes); /*number of auxillary nodes*/

/*calculate lsf_temp and Vnorm for all current trial functions*/
LocalVext2(NodeX, NodeY, Nodes2, NodeStat, Vnorm, Vtemp, lsf, lsf_temp, known, trial, NumNodes, NodeCoord, stop, h, AuxNodes, sign, elemX, hxz);

    /*update trial set by considering neighbours of nodes in known set*/
    
   /*
    for(j=1;j<NodeY-1;j++)
    {
        for(i=1;i<NodeX-1;i++)
        {
            Nnum = Nodes2[i][j]; //read in node number
            //If current node not in known set, but is in the narrow band active set and is the correct side of the boundary
            if((known[Nnum] != 1) && (NodeStat[Nnum] == side) && (InSet(Nnum,active,NumAct) == 1))
            {
                //if any neighbouring node is in known set then curret node is a trial node
                if(known[Nodes2[i-1][j]] == 1) {
                    trial[Nnum] = 1;
                }
                else if(known[Nodes2[i+1][j]] == 1) {
                    trial[Nnum] = 1;
                }
                else if(known[Nodes2[i][j-1]] == 1) {
                    trial[Nnum] = 1;
                }
                else if(known[Nodes2[i][j+1]] == 1) {
                    trial[Nnum] = 1;
                }
            }
        }
    }
    */
    
    /*Do until trial set is empty*/
    short flag; /*used to detect when all extension velocities for active node set have been calcualted*/
    double *lsf_trial;
    lsf_trial = calloc(NumNodes,sizeof(double)); /*array to store trial values of the lsf_temp during the fast marching method*/
    
    count = 0;
    /*set a maximum lsf value to larger than maximum domain length*/
    double lsf_max = (NodeX > NodeY) ? (100.0 * h * (double)NodeX) : (100.0 * h * (double)NodeY);
    
    /*set two arrays to store nodes to be updated during each iteration of the fast marching method*/
    int *xind, *yind, ncount;
    xind = malloc(NumNodes * sizeof(int));
    yind = malloc(NumNodes * sizeof(int));
    
    do /*do unitl all velocities are calculated*/
    {
        /*printf("\nUpdated Sets");
         for(i=0;i<NumNodes;i++)
         {
         if(trial[i] == 1)
         {
         printf("\nNode %i trial",i+1);
         }
         if(known[i] == 1)
         {
         printf("\nNode %i known",i+1);
         printf(" - lsf_temp = %lf",lsf_temp[i]);
         }
         }*/
        
        flag = 0;
        lsf_min = lsf_max; /*initalize to a value much bigger than the element edge length*/
        
        /*For all trial nodes calculate lsf temp*/
        for(j=1;j<NodeY-1;j++)
        {
            for(i=1;i<NodeX-1;i++)
            {
                Nnum = Nodes2[i][j]; /*read in node number*/
                /*If current node not in known set, but is in the narrow band active set and is the correct side of the boundary*/
                if(trial[Nnum] == 1)
                {
                    flag = 1; /*indicates that more velocites required calculating*/
                    /*printf("\nNode %i is trial",Nnum+1);*/
                    
                    if(lsf_trial[Nnum] < 0.000001) /*if trial values needs to be calculated*/
                    {
                        Af = 0.0;
                        Bf = 0.0;
                        Cf = 0.0;
                        ftemp = 0.0;  /*Initalise values*/
                        
                        /*look at all neighbouring nodes*/
                        Nnum2 = Nodes2[i][j-1];
                        Nnum3 = Nodes2[i][j+1]; /*read in node number (above and below)*/
                        
                        /*If neigbouring node is in known set update co-efficents for lsf_temp calc*/
                        /*ensure only closest node to boundary is used for upwind scheme*/
                        if((known[Nnum2] == 1) && (known[Nnum3] == 1))
                        {
                            Af += 1.0;
                            ftemp =  (lsf_temp[Nnum2] < lsf_temp[Nnum3]) ? lsf_temp[Nnum2] : lsf_temp[Nnum3]; /*choose closest known node*/
                        }
                        else if(known[Nnum2] == 1)
                        {
                            Af += 1.0;
                            ftemp = lsf_temp[Nnum2];
                        }
                        else if(known[Nnum3] == 1)
                        {
                            Af += 1.0;
                            ftemp = lsf_temp[Nnum3];
                        }
                        
                        Bf += ftemp;
                        Cf += ftemp * ftemp;	/*Update quadratic co-efficients*/
                        
                        Nnum2 = Nodes2[i-1][j];
                        Nnum3 = Nodes2[i+1][j]; /*read in node number (left and right)*/
                        ftemp = 0.0; /*re-initalise for other direction*/
                        
                        /*If neigbouring node is in known set update co-efficents for lsf_temp calc*/
                        /*ensure only closest node to boundary is used for upwind scheme*/
                        if((known[Nnum2] == 1) && (known[Nnum3] == 1))
                        {
                            Af += 1.0;
                            ftemp =  (lsf_temp[Nnum2] < lsf_temp[Nnum3]) ? lsf_temp[Nnum2] : lsf_temp[Nnum3]; /*choose closest known node*/
                        }
                        else if(known[Nnum2] == 1)
                        {
                            Af += 1.0;
                            ftemp = lsf_temp[Nnum2];
                        }
                        else if(known[Nnum3] == 1)
                        {
                            Af += 1.0;
                            ftemp = lsf_temp[Nnum3];
                        }
                        
                        Bf += ftemp;
                        Cf += ftemp * ftemp;	/*Update quadratic co-efficients*/
                        
                        /*Need to do final update of the Bf and Cf co-efficients*/
                        Bf *= 2.0;
                        Cf -= h * h;
                        
                        if(Bf < 0.0)
                        {
                            printf("\nERROR! In Vext for node %i, Bf = %lf",Nnum+1,Bf);
                        }
                        
                        /*Calculate lsf_temp by solving the quadratic equation*/
                        ftemp = Bf * Bf;
                        ftemp -= 4 * Af * Cf;
                        dtemp = sqrt(ftemp);
                        ftemp = Bf + dtemp;
                        lsf_trial[Nnum] = ftemp / (2.0 * Af);
                        
                        if( (lsf_trial[Nnum] < 0.000001) ) /*|| ((lsf_trial[Nnum] - lsf_max) > -0.00001) )*/
                        {
                            printf("\nERROR! In Vext for node %i, lsf_trial = %lf",Nnum+1,lsf_trial[Nnum]);
                            printf("\nAf=%lf,Bf=%lf,Cf=%lf",Af,Bf,Cf);
                        }
                    }
                    
                    /*If current lsf_temp distance is less than current minimum, update minimum and store node number*/
                    if((lsf_trial[Nnum] - lsf_min) < 0.0)
                    {
                        lsf_min = lsf_trial[Nnum];
                    }
                }
            }
        }
        
        if(flag == 0) /*if there are no trial nodes, then end the loop*/
        {
            break;
        }
        
        /*printf("\nlsf_min = %lf", lsf_min);*/
        
        ncount = 0; /*re-initialize*/
        /*check to see which nodes need updating this iteration*/
        for(j=1;j<NodeY-1;j++)
        {
            for(i=1;i<NodeX-1;i++)
            {
                Nnum = Nodes2[i][j]; /*read in node number*/
                
                if(lsf_trial[Nnum] < 0.0)
                {
                    printf("\nERROR IN ReInt: lsf_trial = %lf, for node %i",lsf_trial[Nnum],Nnum+1);
                }
                
                /*if a trial node has a trail lsf value <= minimum then update*/
                if( (trial[Nnum] == 1) && (fabs(lsf_trial[Nnum] - lsf_min) < 0.000001) )
                {
                    /*printf("\nlsf_trail = %lf, lsf_min = %lf",lsf_trial[Nnum],fabs(lsf_min));
                     /*-----------For lowest node with lowest lsf temp value------------*/
                    Nnum = Nodes2[i][j]; /*read in node number*/
                    /*printf("\nTrial EX = %i",Nnum+1);*/
                    lsf_temp[Nnum] = lsf_trial[Nnum]; /*update lsf_temp for node with its trial value*/
                    /*look at all neighbouring nodes*/
                    Lsum = 0.0;
                    Fsum = 0.0; /*initaize sum variables to zero*/
                    
                    /*look at all neighbouring nodes*/
                    Nnum2 = Nodes2[i][j-1];
                    Nnum3 = Nodes2[i][j+1]; /*read in node number (above and below)*/
                    
                    /*If neigbouring node is in known set- potentially use for Vext calc*/
                    /*ensure only closest node to boundary is used for upwind scheme*/
                    if((known[Nnum2] == 1) && (known[Nnum3] == 1))
                    {
                        Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum2]); /*absolute diff to node below*/
                        Bf = fabs(lsf_temp[Nnum] - lsf_temp[Nnum3]); /*absolute diff to node above*/
                        Cf = lsf_temp[Nnum3] - lsf_temp[Nnum2];		 /*diff between above and below nodes*/
                        
                        /*choose values in shortest direction*/
                        if(fabs(Cf) < 0.000001)  /*if distance is same both sides*/
                        {
                            Lsum += Af; /*could be Bf, doesn't matter*/
                            vt1 = Vtemp[Nnum2];
                            vt2 = Vtemp[Nnum3]; /*find absolute values of neighbouring velocities*/
                            Fsum += ((sign * vt1) > (sign * vt2)) ? (Af * vt1) : (Bf * vt2); /*choose velocity that will move lsf towards zero, or least distance away from zero*/
                        }
                        
                        else if(Cf > -0.000001) /*else if Af < Bf, i.e. node below is closer*/
                        {
                            Lsum += Af;
                            Fsum += Af * Vtemp[Nnum2];
                        }
                        
                        else if(Cf < 0.000001)  /*else if Bf < Af, i.e. node above is closer*/
                        {
                            Lsum += Bf;
                            Fsum += Bf * Vtemp[Nnum3];
                        }
                        
                        else
                        {
                            printf("\nERROR! Can't calculate Vext for Node %i",Nnum);
                        }
                    }
                    else if(known[Nnum2] == 1)
                    {
                        Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum2]);
                        Lsum += Af;
                        Fsum += Af * Vtemp[Nnum2];
                    }
                    else if(known[Nnum3] == 1)
                    {
                        Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum3]);
                        Lsum += Af;
                        Fsum += Af * Vtemp[Nnum3];
                    }
                    
                    /*else
                     {
                     printf("\nERROR! Node %i identified as Trial has no Y neighbours",Nnum+1);
                     }*/
                    
                    Nnum2 = Nodes2[i-1][j];
                    Nnum3 = Nodes2[i+1][j]; /*read in node number (left and right)*/
                    
                    /*If neigbouring node is in known set- potentially use for Fext calc*/
                    /*ensure only closest node to boundary is used for upwind scheme*/
                    if((known[Nnum2] == 1) && (known[Nnum3] == 1))
                    {
                        Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum2]); /*absolute diff to left node*/
                        Bf = fabs(lsf_temp[Nnum] - lsf_temp[Nnum3]); /*absolute diff to right node*/
                        Cf = lsf_temp[Nnum3] - lsf_temp[Nnum2];		 /*diff between left and right nodes*/
                        
                        /*choose values in shortest direction*/
                        if(fabs(Cf) < 0.000001)  /*if distance is same both sides*/
                        {
                            Lsum += Af; /*could be Bf, doesn't matter*/
                            vt1 = Vtemp[Nnum2];
                            vt2 = Vtemp[Nnum3]; /*find absolute values of neighbouring velocities*/
                            Fsum += ((sign * vt1) > (sign * vt2)) ? (Af * vt1) : (Bf * vt2); /*choose velocity that will move lsf towards zero, or least distance away from zero*/
                        }
                        
                        else if(Cf > -0.000001) /*else if Af < Bf, i.e. left node is closer*/
                        {
                            Lsum += Af;
                            Fsum += Af * Vtemp[Nnum2];
                        }
                        
                        else if(Cf < 0.000001)  /*else if Bf < Af, i.e. right node is closer*/
                        {
                            Lsum += Bf;
                            Fsum += Bf * Vtemp[Nnum3];
                        }
                        
                        else
                        {
                            printf("\nERROR! Can't calculate Vext for Node %i",Nnum);
                        }
                    }
                    else if(known[Nnum2] == 1)
                    {
                        Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum2]);
                        Lsum += Af;
                        Fsum += Af * Vtemp[Nnum2];
                    }
                    else if(known[Nnum3] == 1)
                    {
                        Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum3]);
                        Lsum += Af;
                        Fsum += Af * Vtemp[Nnum3];
                    }
                    
                    /*else
                     {
                     printf("\nERROR! Node %i identified as Trial has no X neighbours",Nnum+1);
                     }*/
                    
                    /*Calculate Extension velocity and move node to known set*/
                    /*printf("\nNode %i lsf_temp=%lf, Lsum=%lf, Fsum=%lf",Nnum+1,lsf_min,Lsum,Fsum);*/
                    Vtemp[Nnum] = Fsum / Lsum;
                    
                    xind[ncount] = i;
                    yind[ncount++] = j; /*store update node for next bit*/
                    
                    /*printf("\nlsf_min = %lf",lsf_min);*/
                    count++;
                    
                }
            }
        }
        
        /*update known set for next iteration*/
        for(i=0;i<ncount;i++)
        {
            Ci = xind[i];
            Cj = yind[i]; /*read in updated node indicators*/
            
            /*read in node number*/				
            Nnum = Nodes2[Ci][Cj];
            
            trial[Nnum] = 0;
            known[Nnum] = 1; /*update known set*/
        }
        
        /*update trial set for next iteration*/
        for(i=0;i<ncount;i++)
        {
            Ci = xind[i];
            Cj = yind[i]; /*read in updated node indicators*/
            
            /*look at all neighbouring nodes*/				
            Nnum2 = Nodes2[Ci][Cj-1];
            Nnum3 = Nodes2[Ci][Cj+1]; /*read in node number (above and below)*/
            
            /*Check to see if node below needs to be added to trial set*/
            if((known[Nnum2] != 1) && (InSet(Nnum2,active,NumAct) == 1) && (NodeStat[Nnum2] == side))
            {
                trial[Nnum2] = 1;
                lsf_trial[Nnum2] = 0.0; /*reset trial value*/
            }
            
            /*Check to see if node above needs to be added to trial set*/
            if((known[Nnum3] != 1) && (InSet(Nnum3,active,NumAct) == 1) && (NodeStat[Nnum3] == side))
            {
                trial[Nnum3] = 1;
                lsf_trial[Nnum3] = 0.0; /*reset trial value*/
            }
            
            Nnum2 = Nodes2[Ci-1][Cj];
            Nnum3 = Nodes2[Ci+1][Cj]; /*read in node number (left and right)*/
            
            /*Check to see if node to left needs to be added to trial set*/
            if((known[Nnum2] != 1) && (InSet(Nnum2,active,NumAct) == 1) && (NodeStat[Nnum2] == side))
            {
                trial[Nnum2] = 1;
                lsf_trial[Nnum2] = 0.0; /*reset trial value*/
            }
            
            /*Check to see if node to right needs to be added to trial set*/
            if((known[Nnum3] != 1) && (InSet(Nnum3,active,NumAct) == 1) && (NodeStat[Nnum3] == side))
            {
                trial[Nnum3] = 1;
                lsf_trial[Nnum3] = 0.0; /*reset trial value*/
            }
        }
        
        if(ncount == 0)
        {
            printf("\nERROR! ncount=0 in Vext! Aborting");
        }
    }
    while(ncount > 0);
    
    /*FILE *outfile;	/*File varible for output files
     char plotname[20];	/*variable to change names of plotting output files
     sprintf(plotname,"Phi_temp%i.txt",sign); /*set name for current topology output file*/
    
    /* Write temp Node signed distance information file
     outfile = fopen(plotname, "w");
     if(outfile == NULL){
     printf("\nFailed to open Phi_temp Implicit writefile\n");
     }
     else{
     /*fprintf(outfile,"Node Num\tx\ty\tphi\n"); /*column headings
     for(i=0;i<NumNodes;i++)
     {
     fprintf(outfile,"%i\t",(i+1));
     fprintf(outfile,"%lf\t",NodeCoord[i].x);
     fprintf(outfile,"%lf\t",NodeCoord[i].y);
     fprintf(outfile,"%lf\t",lsf_temp[i]);
     
     fprintf(outfile,"\n");
     }
     }
     
     fclose(outfile);
     printf("\nlsf_temp Info File written\n");*/
    


/*finally copy calcualted Vtemp values to Vnorm array*/
for(i=0;i<NumNodes;i++)
{
	if(known[i] == 1)
	{
		Vnorm[i] = Vtemp[i];
        
	}
    //printf("%d", known[i]);
    //printf("Vnorm: %f", Vnorm[i]);
}

/*Free memory used*/
free(Vtemp);
free(known);
free(trial);
free(lsf_temp);

		
}

/*function that works out Velocity for nodes close to the boundary*/
void LocalVext2(int NodeX, int NodeY, int Nodes[NodeX][NodeY], short *NodeStat, double *Vnorm, double *Vtemp, double *lsf, double *lsf_temp, short *known, short *trial,
					int NumNodes, Coord NodeCoord[NumNodes], int stop, double h, Aux *AuxNodes, int sign, int elemX, double hxz[elemX])
{	
	double s,s1,s2,t,t1,t2; /*distance variables*/
	double vs,vs1,vs2,vt,vt1,vt2; /*velocity variables*/
	int i,j,k; /*incrementors*/
	double ftemp; /*temp varibale*/
	int Nnum,Nnum2; /*node number variables*/
	int check; /*variable to check aux node has been correctly identified*/
	double h2xz[elemX];
	double h3xz;
	double h2 = 0.5 * h;
	for(i=0;i<elemX;i++)
	{
		h2xz[i] = hxz[i]*0.5;
	}
	
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++)
		{
			Nnum = Nodes[i][j];
			if(trial[Nnum] == 1)
			{
				/*look at each neighbouring node in turn*/
				
				/*-----------look at node below-----------*/
				Nnum2 = Nodes[i][j-1];
				/*if( (fabs(lsf[Nnum2]) < 0.0000001) && (NodeStat[Nnum2] == 2) ) /*If neigbouring node is on the boundary*/
				if(NodeStat[Nnum2] == 2) /*If neigbouring node is on the boundary*/
				{
					t1 = h; /*distance to boundary = element edge length*/
					vt1 = Vnorm[Nnum2]; /*store velocity of the node*/
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) /*otherwise check to see is grid line is intersected*/
				{
					t1 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					t1 += 1.0;
					t1 *= h2;
					check = -1;
					
					if(NodeStat[Nnum] == 2)
					{
						vt1 = Vnorm[Nnum];
						check = 1;
					}
					
					else
					{
						for(k=0;k<stop;k++) /*For all auxillary nodes (contained in the auxNode array)*/
						{
							if(fabs(NodeCoord[Nnum].x - AuxNodes[k].x) < 0.000001)	/*If aux node has same x-coord*/
							{
								/*If aux node lies between the two grid nodes*/
								if( ((AuxNodes[k].y - NodeCoord[Nnum].y) < 0.000001) && ((AuxNodes[k].y - NodeCoord[Nnum2].y) > -0.000001) )
								{
									vt1 = Vnorm[AuxNodes[k].n -1]; /*store velocity of the auxillary node*/
									check = 1;
									break;
								}
							}
						}
					}
					
					if(check == -1)
					{
						printf("\nERROR! couldn't find aux node below in LocalVext near node %i",Nnum+1);
						printf("\t ..... %lf!",(lsf[Nnum] * lsf[Nnum2]));
					}
					
					/*check to see if t1 > h or = 0.0, i.e. something has gone wrong in update*/
					/*If so then approximate*/
					if( ((t1 - h) > -0.000001) || (t1 < 0.000001) )
					{
						if( ((sign * -lsf[Nnum2]) - h) > -0.000001 )
						{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							t1 = (lsf[Nnum] * h) / ftemp;
						}
						else
						{
							t1 = h + (sign * lsf[Nnum2]);
						}
						
						printf("\nERROR in LocalVext for node %i, t1 approx as %lf",Nnum+1,t1);
					}

				}
				else
				{
					t1 = h + 1.0; /*set distance to h+1.0 to make algorithm work later*/
					vt1 = 0.0; /*set velocity to 0.0, also for later*/
				}
				
				if(t1 < 0.000001)
				{
					printf("\nERROR! for node %i in LocalVext t1 = %lf",Nnum+1,t1);
				}
				
				/*----------look at node above----------*/
				Nnum2 = Nodes[i][j+1];
				/*if( (fabs(lsf[Nnum2]) < 0.0000001) && (NodeStat[Nnum2] == 2) ) /*If neigbouring node is on the boundary*/
				if(NodeStat[Nnum2] == 2) /*If neigbouring node is on the boundary*/
				{
					t2 = h; /*distance to boundary = element edge length*/
					vt2 = Vnorm[Nnum2]; /*store velocity of the node*/
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) /*otherwise check to see is grid line is intersected*/
				{
					t2 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					t2 += 1.0;
					t2 *= h2;
					check = -1;	
							
					if(NodeStat[Nnum] == 2)
					{
						vt2 = Vnorm[Nnum];
						check = 1;
					}
					
					else
					{
						for(k=0;k<stop;k++) /*For all auxillary nodes (contained in the auxNode array)*/
						{
							if( fabs(NodeCoord[Nnum].x - AuxNodes[k].x) < 0.000001 )	/*If aux node has same x-coord*/
							{	
								if( ((AuxNodes[k].y - NodeCoord[Nnum2].y) < 0.000001) && ((AuxNodes[k].y - NodeCoord[Nnum].y) > -0.000001) ) /*If aux node lies between the two grid nodes*/
								{
									vt2 = Vnorm[AuxNodes[k].n -1]; /*store velocity of the auxillary node*/
									check = 1;
									break;
								}
							}
						}
					}
					
					if(check == -1)
					{
						printf("\nERROR! couldn't find aux node above in LocalVext near node %i",Nnum+1);
						printf("\t ..... %lf!",(lsf[Nnum] * lsf[Nnum2]));
					}
					
					/*check to see if t2 > h or = 0.0, i.e. something has gone wrong in update*/
					/*If so then approximate*/	
					if( ((t2 - h) > -0.000001) || (t2 < 0.000001) )
					{		
						if( ((sign * -lsf[Nnum2]) - h) > -0.000001 )
						{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							t2 = (lsf[Nnum] * h) / ftemp;
						}
						else
						{
							t2 = h + (sign * lsf[Nnum2]);
						}
						
						printf("\nERROR in LocalVext for node %i, t2 approx as %lf",Nnum+1,t2);
					}

				}
				else
				{
					t2 = h + 1.0; /*set distance to h+1.0 to make algorithm work later*/
					vt2 = 0.0; /*set velocity to 0.0, also for later*/
				}
				
				if(t2 < 0.000001)
				{
					printf("\nERROR! for node %i in LocalVext t2 = %lf",Nnum+1,t2);
				}
				
				/*----------look at node to left----------*/
				Nnum2 = Nodes[i-1][j];
				/*if( (fabs(lsf[Nnum2]) < 0.0000001) && (NodeStat[Nnum2] == 2) ) /*If neigbouring node is on the boundary*/
				if(NodeStat[Nnum2] == 2) /*If neigbouring node is on the boundary*/
				{
					s1 = hxz[i-2]; /*distance to boundary = element edge length*/
					vs1 = Vnorm[Nnum2]; /*store velocity of the node*/
					/*printf("\t s1, bN");
					printf("\t hxz[%i] = %f", i-2, hxz[i-2]);*/
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) /*otherwise check to see is grid line is intersected*/
				{
					s1 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					s1 += 1.0;
					s1 *= h2xz[i-2];
					check = -1;
			
					if(NodeStat[Nnum] == 2)
					{
						vs1 = Vnorm[Nnum];
						check = 1;
					}
					
					else
					{
						for(k=0;k<stop;k++) /*For all auxillary nodes (contained in the auxNode array)*/
						{
							if( fabs(NodeCoord[Nnum].y - AuxNodes[k].y) < 0.000001 )/*If aux node has same y-coord*/
							{
								/*If aux node lies between the two grid nodes*/
								if( ((AuxNodes[k].x - NodeCoord[Nnum].x) < 0.000001) && ((AuxNodes[k].x - NodeCoord[Nnum2].x) > -0.000001) )
								{
									vs1 = Vnorm[AuxNodes[k].n -1]; /*store velocity of the auxillary node*/
									check = 1;
									break;
								}
							}
						}
					}
					
					if(check == -1)
					{
						printf("\nERROR! couldn't find aux node to left in LocalVext near node %i",Nnum+1);
						printf("\t ..... %lf!",(lsf[Nnum] * lsf[Nnum2]));
					}
					
					/*check to see if s1 > h or s1 = 0.0, i.e. something has gone wrong in update*/
					/*If so then approximate*/
					if( ((s1 - hxz[i-2]) > -0.000001) || (s1 < 0.000001) )
					{
						if( ((sign * -lsf[Nnum2]) - hxz[i-2]) > -0.000001 )
						{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							s1 = (lsf[Nnum] * hxz[i-2]) / ftemp;
						}
						else
						{
							s1 = hxz[i-2] + (sign * lsf[Nnum2]);
						}
						
						printf("\nERROR in LocalVext for node %i, s1 approx as %lf",Nnum+1,s1);
					}

				}
				else
				{
					s1 = hxz[i-2] + 1.0; /*set distance to h+1.0 to make algorithm work later*/
					vs1 = 0.0; /*set velocity to 0.0, also for later*/
				}
				
				if(s1 < 0.000001)
				{
					printf("\nERROR! for node %i in LocalVext s1 = %lf",Nnum+1,s1);
				}
				
				/*----------look at node to right----------*/
				Nnum2 = Nodes[i+1][j];
				/*if( (fabs(lsf[Nnum2]) < 0.0000001) && (NodeStat[Nnum2] == 2) ) /*If neigbouring node is on the boundary*/
				if(NodeStat[Nnum2] == 2) /*If neigbouring node is on the boundary*/
				{
					s2 = hxz[i-1]; /*distance to boundary = element edge length*/
					vs2 = Vnorm[Nnum2]; /*store velocity of the node*/
					/*printf("\t s2, bN");
					printf("\t hxz[%i] = %f", i-1, hxz[i-1]);*/
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) /*otherwise check to see is grid line is intersected*/
				{
					s2 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					s2 += 1.0;
					s2 *= h2xz[i-1];
					check = -1;
							
					if(NodeStat[Nnum] == 2)
					{
						vs2 = Vnorm[Nnum];
						check = 1;
					}
					
					else
					{
						for(k=0;k<stop;k++) /*For all auxillary nodes (contained in the auxNode array)*/
						{
							if( fabs(NodeCoord[Nnum].y - AuxNodes[k].y) < 0.000001 )/*If aux node has same y-coord*/
							{
								/*If aux node lies between the two grid nodes*/
								if( ((AuxNodes[k].x - NodeCoord[Nnum2].x) < 0.000001) && ((AuxNodes[k].x - NodeCoord[Nnum].x) > -0.000001) )
								{
									vs2 = Vnorm[AuxNodes[k].n -1]; /*store velocity of the auxillary node*/
									check = 1;
									break;
								}
							}
						}
					}
					
					if(check == -1)
					{
						printf("\nERROR! couldn't find aux node to right in LocalVext near node %i",Nnum+1);
						printf("\t ..... %lf!",(lsf[Nnum] * lsf[Nnum2]));
					}
					
					/*check to see if s2 > h or s2 = 0.0, i.e. something has gone wrong in update*/
					/*If so then approximate*/
					if( ((s2 - hxz[i-1]) > -0.000001) || (s2 < 0.000001) )
					{
						if( ((sign * -lsf[Nnum2]) - hxz[i-1]) > -0.000001 )
						{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							s2 = (lsf[Nnum] * hxz[i-1]) / ftemp;
						}
						else
						{
							s2 = hxz[i-1] + (sign * lsf[Nnum2]);
						}
						
						printf("\nERROR in LocalVext for node %i, s2 approx as %lf",Nnum+1,s2);
					}
				}
				else
				{
					s2 = hxz[i-1] + 1.0; /*set distance to h+1.0 to make algorithm work later*/
					vs2 = 0.0; /*set velocity to 0.0, also for later*/
				}
				
				if(s2 < 0.000001)
				{
					printf("\nERROR! for node %i in LocalVext s2 = %lf",Nnum+1,s2);
				}
				/*printf("\n s1 = %f  \ts2 = %f \tt1 = %f \tt2 = %f", s1, s2, t1, t2);*/
				
				/*----------------------Now calculate lsf temp and Vnorm for node---------------------------*/
				/*choose lowest s and t values, i.e. nearest to boundary*/
				if(fabs(s1 - s2) < 0.000001) /*if distances are the same then choose greatest velocity*/
				{
					s = s1; /*could be s2, doesn't make much difference*/
					h3xz = hxz[i-2];
					vs = ((sign * vs1) > (sign * vs2)) ? vs1 : vs2; /*choose velocity that will move lsf towards zero, or least distance away from zero*/
					/*printf("\t A  h3xz[%i] = %f",i-2, h3xz);*/
				}
				else if(s1 < s2)	{
					s = s1;
					vs = vs1;
					h3xz = hxz[i-2];
					/*printf("\t B   h3xz[%i] = %f",i-2, h3xz);*/
				}
				else if(s1 > s2) {
					s = s2;
					vs = vs2;
					h3xz = hxz[i-1];
					/*printf("\t C   h3xz[%i] = %f",i-1, h3xz);*/
				}
				/*printf("\t h3xz = %f", h3xz);*/
				/* t values */
				if(fabs(t1 - t2) < 0.000001) /*if distances are the same then choose greatest velocity*/
				{
					t = t1;  /*could be t2, doesn't make much difference*/
					vt = ((sign * vt1) > (sign * vt2)) ? vt1 : vt2; /*choose velocity that will move lsf towards zero, or least distance away from zero*/
				}
				else if(t1 < t2)	{
					t = t1;
					vt = vt1;
				}
				else if(t1 > t2) {
					t = t2;
					vt = vt2;
				}
				
				/*Calculate Extension velocity*/
				s2 = ((s - h3xz) > 0.1) ? 0.0 : (1.0 / (s * s));
				t2 = ((t - h) > 0.1) ? 0.0 : (1.0 / (t * t)); /*weighting for vs and vt velocities*/
				/*printf("\ns = %lf, t = %lf",s,t);*/
				ftemp = s2 * vs;
				ftemp += t2 * vt;
				ftemp /= (s2 + t2); /*Vnorm = Vext, by using inverse square weighting of distances*/
				Vtemp[Nnum] = ftemp;
				
				/*printf("\nVnorm = %lf",ftemp);*/
				
				/*Now calculate lsf temp*/
				if((s - h3xz) > 0.1) /*If horizontal distance > edge length only use vertical distance*/
				{
					lsf_temp[Nnum] = t;
				}
				
				else if((t - h) > 0.1)	/*If vertical distance > edge length only use horizontal distance*/
				{
					lsf_temp[Nnum] = s;
				}
				
				else
				{
					s2 = s * s;
					t2 = t * t; /*squared distances*/
					
					ftemp = s * t;
					ftemp /= sqrt(s2 + t2);
					lsf_temp[Nnum] = ftemp; /*calculate perpendicular distance to the boundary (by pythag)*/
				}
				if((s > h3xz) && (t > h)) /*print error if both horizontal and vertical distances greater than edge length*/
				{
					printf("\nERROR! For node %i, both s & t > h in LocalVext, s=%lf, t=%lf", Nnum+1,s,t);
				}
				
				/*Move node into known set*/
				trial[Nnum] = 0;
				known[Nnum] = 1;
			}
		}
	}
	/*printf("\nLocal Known Set");
	for(i=0;i<NumNodes;i++)
	{
		if(trial[i] == 1)
		{
			printf("\nNode %i",i+1);
		}
	}*/
}

/*Function to calcualte gradient using WENO scheme*/
double GradWENO(int Xi, int Yj, int num, double *lsf, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int sign, 
				double Vn, double dt, double h, int elemX, double hxz[elemX])
{
	double v1,v2,v3,v4,v5; /*varibales for finite differences used*/
	double grad;	/*variable for gradient approx*/
	double ftemp;	/*tempory variable*/
	double htemp;
	int flag = 0;
	double phi = lsf[num]; /*read in lsf value of current node*/
	
    
	if(Xi == 1) /*If the node is on the left edge*/
	{	
		if(Yj == 1) /*node is bottom left corner*/
		{
			/*printf("\nNode %i bottom left",Nodes2[Xi][Yj]+1);
			/*if lsf of node above and to right have the same value as the corner node then use node on diagonal*/
			if( (fabs(lsf[Nodes2[Xi][Yj+1]] - phi) < 0.000001) && (fabs(lsf[Nodes2[Xi+1][Yj]] - phi) < 0.000001) )
			{
				/*printf("\nUse Diagonal");*/
				grad = fabs(phi - lsf[Nodes2[Xi+1][Yj+1]]); /*distance to diagonal node*/
                htemp = 1.414213562;//(h*h)+(hxz[0]*hxz[0]);
				grad *= sqrt(htemp); /*multiply by sqrt(of distance from corner to oposite corner node)*/
				grad /= h;		/*divide by edge length*/
				flag = 1;	/*flag gradient already calcualted*/
			}
		}
		
		else if(Yj == NodeY -2) /*node is top left corner*/
		{
			/*printf("\nNode %i top left",Nodes2[Xi][Yj]+1);
			/*if lsf of node below and to right have the same value as the corner node then use node on diagonal*/
			if( (fabs(lsf[Nodes2[Xi][Yj-1]] - phi) < 0.000001) && (fabs(lsf[Nodes2[Xi+1][Yj]] - phi) < 0.000001) )
			{
				/*printf("\nUse Diagonal");*/
				grad = fabs(phi - lsf[Nodes2[Xi+1][Yj-1]]); /*distance to diagonal node*/
                htemp = 1.414213562;//(h*h)+(hxz[0]*hxz[0]);
				grad *= sqrt(htemp); /*multiply by sqrt(of distance from corner to oposite corner node)*/
				grad /= h;		/*divide by edge length*/
				flag = 1;	/*flag gradient already calcualted*/
			} 
		}
	}
	
	else if(Xi == NodeX -2) /*If the node is on the right edge*/
	{
		if(Yj == 1) /*node is bottom right corner*/
		{
			/*printf("\nNode %i bottom right",Nodes2[Xi][Yj]+1);
			/*if lsf of node above and to left have the same value as the corner node then use node on diagonal*/
			if( (fabs(lsf[Nodes2[Xi][Yj+1]] - phi) < 0.000001) && (fabs(lsf[Nodes2[Xi-1][Yj]] - phi) < 0.000001) )
			{
				/*printf("\nUse Diagonal");*/
				grad = fabs(phi - lsf[Nodes2[Xi-1][Yj+1]]); /*distance to diagonal node*/
                htemp = 1.414213562;//(h*h)+(hxz[elemX-1]*hxz[elemX-1]);
				grad *= sqrt(htemp); /*multiply by sqrt(of distance from corner to oposite corner node)*/
				grad /= h;		/*divide by edge length*/
				flag = 1;	/*flag gradient already calcualted*/
			} 
		}
		
		else if(Yj == NodeY -2) /*node is top right corner*/
		{
			/*printf("\nNode %i top right",Nodes2[Xi][Yj]+1);
			/*if lsf of node below and to left have the same value as the corner node then use node on diagonal*/
			if( (fabs(lsf[Nodes2[Xi][Yj-1]] - phi) < 0.000001) && (fabs(lsf[Nodes2[Xi-1][Yj]] - phi) < 0.000001) )
			{
				/*printf("\nUse Diagonal");*/
				grad = fabs(phi - lsf[Nodes2[Xi-1][Yj-1]]); /*distance to diagonal node*/
                htemp = 1.414213562;//(h*h)+(hxz[elemX-1]*hxz[elemX-1]);
				grad *= sqrt(htemp); /*multiply by sqrt(of distance from corner to oposite corner node)*/
				grad /= h;		/*divide by edge length*/
				flag = 1;	/*flag gradient already calcualted*/
			}  
			
		}
	}
	
	if(flag != 1)
	{
	
	/*------------first assume info travelling to the right----------*/
	if(Xi == 1) /*if node on left edge*/
	{
		v1 = (lsf[Nodes2[4][Yj]] - lsf[Nodes2[3][Yj]]) / hxz[2];
		v2 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / hxz[1];
		v3 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / hxz[0];
		v4 = v3; /*approx gradient outside of domain*/
		v5 = v3; /*approx gradient outside of domain*/
	}	
	else if(Xi == 2)
	{
		v1 = (lsf[Nodes2[5][Yj]] - lsf[Nodes2[4][Yj]]) / hxz[4];
		v2 = (lsf[Nodes2[4][Yj]] - lsf[Nodes2[3][Yj]]) / hxz[3];
		v3 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / hxz[2];
		v4 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / hxz[1];
		v5 = v4; /*approx gradient outside of domain*/
	}	
	else if(Xi == NodeX -2) /*if node on right edge*/
	{
		v5 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / hxz[elemX-2];
		v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / hxz[elemX-1];
		v3 = v4;
		v2 = v4;
		v1 = v4;	/*approx gradients outside of domain*/
	}	
	else if(Xi == NodeX -3)
	{
		v5 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / hxz[elemX-3];
		v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / hxz[elemX-2];
		v3 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / hxz[elemX-1];
		v2 = v3;
		v1 = v3;	/*approx gradients outside of domain*/
	}	
	else if(Xi == NodeX -4)
	{
		v5 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / hxz[elemX-4];
		v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / hxz[elemX-3];
		v3 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / hxz[elemX-2];
		v2 = (lsf[Nodes2[Xi+2][Yj]] - lsf[Nodes2[Xi+1][Yj]]) / hxz[elemX-1];
		v1 = v2;
	}
	else
	{
		v1 = (lsf[Nodes2[Xi+3][Yj]] - lsf[Nodes2[Xi+2][Yj]]) / hxz[Xi+1];
		v2 = (lsf[Nodes2[Xi+2][Yj]] - lsf[Nodes2[Xi+1][Yj]]) / hxz[Xi];
		v3 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / hxz[Xi-1];
		v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / hxz[Xi-2];
		v5 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / hxz[Xi-3];
	}
	
	double Xr = (double)sign * GWsub(v1, v2, v3, v4, v5); /*calcualte gradient, info going to the right*/
	
	/*------------Now assume info travelling to the left----------*/
	if(Xi == NodeX -2) /*if node on right edge*/
	{
		v1 = (lsf[Nodes2[Xi-2][Yj]] - lsf[Nodes2[Xi-3][Yj]]) / hxz[elemX-3];
		v2 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / hxz[elemX-2];
		v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / hxz[elemX-1];
		v4 = v3; /*approx gradient outside of domain*/
		v5 = v3; /*approx gradient outside of domain*/
	}
	else if(Xi == NodeX -3)
	{
		v1 = (lsf[Nodes2[Xi-2][Yj]] - lsf[Nodes2[Xi-3][Yj]]) / hxz[elemX-4];
		v2 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / hxz[elemX-3];
		v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / hxz[elemX-2];
		v4 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / hxz[elemX-1];
		v5 = v4; /*approx gradient outside of domain*/
	}	
	else if(Xi == 1) /*if node on left edge*/
	{
		v5 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / hxz[1];
		v4 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / hxz[0];
		v3 = v4;
		v2 = v4;
		v1= v4;	/*approx gradients outside of domain*/
	}	
	else if(Xi == 2)
	{
		v5 = (lsf[Nodes2[4][Yj]] - lsf[Nodes2[3][Yj]]) / hxz[2];
		v4 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / hxz[1];
		v3 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / hxz[0];
		v2 = v3;
		v1= v3;	/*approx gradients outside of domain*/
	}	
	else if(Xi == 3)
	{
		v5 = (lsf[Nodes2[5][Yj]] - lsf[Nodes2[4][Yj]]) / hxz[3];
		v4 = (lsf[Nodes2[4][Yj]] - lsf[Nodes2[3][Yj]]) / hxz[2];
		v3 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / hxz[1];
		v2 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / hxz[0];
		v1 = v2;	/*approx gradients outside of domain*/
	}
	else
	{
		v5 = (lsf[Nodes2[Xi+2][Yj]] - lsf[Nodes2[Xi+1][Yj]]) / hxz[Xi];
		v4 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / hxz[Xi-1];
		v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / hxz[Xi-2];
		v2 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / hxz[Xi-3];
		v1 = (lsf[Nodes2[Xi-2][Yj]] - lsf[Nodes2[Xi-3][Yj]]) / hxz[Xi-4];
	}
	
	double Xl = (double)sign * GWsub(v1, v2, v3, v4, v5); /*calcualte gradient info going to the left*/
	
	/*------------Now assume info travelling upward----------*/
	if(Yj == 1) /*if node on bottom edge*/
	{
		v1 = (lsf[Nodes2[Xi][4]] - lsf[Nodes2[Xi][3]]) / h;
		v2 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
		v3 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
		v4 = v3; /*approx gradient outside of domain*/
		v5 = v3; /*approx gradient outside of domain*/
	}	
	else if(Yj == 2)
	{
		v1 = (lsf[Nodes2[Xi][5]] - lsf[Nodes2[Xi][4]]) / h;
		v2 = (lsf[Nodes2[Xi][4]] - lsf[Nodes2[Xi][3]]) / h;
		v3 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
		v4 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
		v5 = v4; /*approx gradient outside of domain*/
	}	
	else if(Yj == NodeY -2) /*if node on top edge*/
	{
		v5 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
		v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
		v3 = v4;
		v2 = v4;
		v1 = v4;	/*approx gradients outside of domain*/
	}	
	else if(Yj == NodeY -3)
	{
		v5 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
		v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
		v3 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
		v2 = v3;
		v1 = v3;	/*approx gradients outside of domain*/
	}	
	else if(Yj == NodeY -4)
	{
		v5 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
		v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
		v3 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
		v2 = (lsf[Nodes2[Xi][Yj+2]] - lsf[Nodes2[Xi][Yj+1]]) / h;
		v1 = v2;
	}
	else
	{
		v5 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
		v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
		v3 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
		v2 = (lsf[Nodes2[Xi][Yj+2]] - lsf[Nodes2[Xi][Yj+1]]) / h;
		v1 = (lsf[Nodes2[Xi][Yj+3]] - lsf[Nodes2[Xi][Yj+2]]) / h;
	}
	
	double Yu = (double)sign * GWsub(v1, v2, v3, v4, v5); /*calcualte gradient info going up*/
	
	/*------------Now assume info travelling down----------*/
	if(Yj == NodeY -2) /*if node on right edge*/
	{
		v1 = (lsf[Nodes2[Xi][Yj-2]] - lsf[Nodes2[Xi][Yj-3]]) / h;
		v2 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
		v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
		v4 = v3; /*approx gradient outside of domain*/
		v5 = v3; /*approx gradient outside of domain*/
	}	
	else if(Yj == NodeY -3)
	{
		v1 = (lsf[Nodes2[Xi][Yj-2]] - lsf[Nodes2[Xi][Yj-3]]) / h;
		v2 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
		v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
		v4 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
		v5 = v4; /*approx gradient outside of domain*/
	}	
	else if(Yj == 1) /*if node on bottom edge*/
	{
		v5 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
		v4 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
		v3 = v4;
		v2 = v4;
		v1= v4;	/*approx gradients outside of domain*/
	}	
	else if(Yj == 2)
	{
		v5 = (lsf[Nodes2[Xi][4]] - lsf[Nodes2[Xi][3]]) / h;
		v4 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
		v3 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
		v2 = v3;
		v1= v3;	/*approx gradients outside of domain*/
	}	
	else if(Yj == 3)
	{
		v5 = (lsf[Nodes2[Xi][5]] - lsf[Nodes2[Xi][4]]) / h;
		v4 = (lsf[Nodes2[Xi][4]] - lsf[Nodes2[Xi][3]]) / h;
		v3 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
		v2 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
		v1 = v2;	/*approx gradients outside of domain*/
	}
	else
	{
		v1 = (lsf[Nodes2[Xi][Yj-2]] - lsf[Nodes2[Xi][Yj-3]]) / h;
		v2 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
		v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
		v4 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
		v5 = (lsf[Nodes2[Xi][Yj+2]] - lsf[Nodes2[Xi][Yj+1]]) / h;
	}
	
	double Yd = (double)sign * GWsub(v1, v2, v3, v4, v5); /*calcualte gradient info going down*/
	
	/*now calculate final gradinet using upwind scheme*/
	ftemp = (Xl > 0.0) ? Xl : 0.0;
	ftemp *= ftemp;
	grad = ftemp;
	
	ftemp = (Xr < 0.0) ? Xr : 0.0;
	ftemp *= ftemp;
	grad += ftemp;
	
	ftemp = (Yd > 0.0) ? Yd : 0.0;
	ftemp *= ftemp;
	grad += ftemp;
	
	ftemp = (Yu < 0.0) ? Yu : 0.0;
	ftemp *= ftemp;
	grad += ftemp;
	
	ftemp = grad;
	grad = sqrt(ftemp); /*need to square root answer*/
	
	}
	
	else
	{
		printf("\nError gradient approimxated in corner");
	}
	
	/*need to check that gradient won't take the boundary outside the domain*/
	/*if( (Yj == 1) || (Yj == NodeY -2) || (Xi == 1) || (Xi == NodeX -2) ) //if node on edge of domain
	{
		if( (phi - (grad * Vn * dt)) > 0.0 ) //if updated hpi value is positive
		{
			grad = phi;
			grad /= (Vn * dt);	//ensure boundary only moves to the edge of the domain
		}
	}
	*/
	return(grad);
}

/*sub-function for GradWENO*/
double GWsub(double v1,double v2,double v3,double v4,double v5)
{
	double ftemp;	/*tempory variable*/
	
	/*Calculate the smoothness of approx 1*/
	ftemp = (v1 - (2.0 * v2) + v3);
	ftemp *= ftemp;
	double s1 = 13.0 * ftemp;
	ftemp = (v1 - (4.0 * v2) + (3.0 * v3));
	ftemp *= ftemp;
	s1 += 3.0 * ftemp;
	
	/*Calculate the smoothness of approx 2*/
	ftemp = (v2 - (2.0 * v3) + v4);
	ftemp *= ftemp;
	double s2 = 13.0 * ftemp;
	ftemp = v2 - v4;
	ftemp *= ftemp;
	s2 += 3.0 * ftemp;
	
	/*Calculate the smoothness of approx 3*/
	ftemp = (v3 - (2.0 * v4) + v5);
	ftemp *= ftemp;
	double s3 = 13.0 * ftemp;
	ftemp = (v5 - (4.0 * v4) + (3.0 * v3));
	ftemp *= ftemp;
	s3 += 3.0 * ftemp;
	
	/*printf("\ns1=%lf, s2=%lf, s3=%lf",s1,s2,s3);*/
	
	/*find minimum s value
	double min = s1;
	min = (s2 < min) ? s2 : min;
	min = (s3 < min) ? s3 : min;*/
	
	/*find maximum v^2 value 
	double max = v1 * v1;
	ftemp = v2 * v2;
	max = (ftemp > max) ? ftemp : max;
	ftemp = v3 * v3;
	max = (ftemp > max) ? ftemp : max;
	ftemp = v4 * v4;
	max = (ftemp > max) ? ftemp : max;
	ftemp = v5 * v5;
	max = (ftemp > max) ? ftemp : max;*/
	
	/*double eps = (0.00001 * max);*/
	/*eps = (eps < 0.0000001) ? 0.000001 : eps;*/
	/*double eps = (min < 0.00001) ? 0.00001 : 0.0;*/
	double eps = 0.00001;
	
	/*calcualte alpha weighting values*/
	/*ftemp = (s1 + eps);
	ftemp *= ftemp;*/
	ftemp = (s1 < eps) ? (eps * eps) : (s1 * s1);
	double a1 = 1.0 / ftemp;
	
	/*ftemp = (s2 + eps);
	ftemp *= ftemp;*/
	ftemp = (s2 < eps) ? (eps * eps) : (s2 * s2);
	double a2 = 6.0 / ftemp;
	
	/*ftemp = (s3 + eps);
	ftemp *= ftemp;*/
	ftemp = (s3 < eps) ? (eps * eps) : (s3 * s3);
	double a3 = 3.0 / ftemp;
	
	double asum = a1 + a2 + a3;
	
	/*calculate normalized weightings*/
	a1 /= asum;
	a2 /= asum;
	a3 /= asum;
	
	/*Finally calculate the approximate gradient*/
	double grad = a1 * ( (2.0 * v1) - (7.0 * v2) + (11.0 * v3) );
	grad += a2 * ( (2.0 * v4) - (v2) + (5.0 * v3) );
	grad += a3 * ( (2.0 * v3) - (v5) + (5.0 * v4) );
	grad *= 0.166666666667;
	
	return(grad);
}

/*Function to compute avergae and max difference between two level set functions*/
void PhiComp(int NumNodes, double *lsf1, double *lsf2, int itt, Coord *Delta)
{
	int i; /*incrementor*/
	double ftemp;
	
	double diff = 0.0; /*initialize difference to zero*/
	double max = 0.0; /*initlaize maximum difference*/
	
	for(i=0;i<NumNodes;i++)
	{
		ftemp = fabs(lsf1[i] - lsf2[i]); /*compute difference at current node*/
		diff += ftemp; /*sum up difference values*/
		max = (ftemp > max) ? ftemp : max; /*update maximum value*/
	}
	
	diff /= (double)NumNodes; /*average difference by number of nodes*/
	
	Delta[itt].x = max;
	Delta[itt].y = diff; /*store difference values*/
	
	printf("Max Phi Difference = %lf\nAverage Phi Difference =%lf\n",max,diff);
}

/*Function to smooth level set function using simple linear approach*/
void SmoothPhi(int NumNodes, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], double *lsf, int *fixed, int NumFix)
{
	int i,j,temp; /*incrementors etc.*/
	double p1,p2,p3,p4,div; /*variables to perform smoothing*/
	
	double *lsf_temp;
	lsf_temp = malloc(NumNodes * sizeof(double));
	
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++) /*For all nodes*/
		{
			temp = Nodes2[i][j]; /*read in node number*/
			
			/*If lsf value not fixed and node is outside the boundary, then smooth*/
			if((InSet(temp,fixed,NumFix) == 0) && (lsf[temp] < 0.000001))
			{
				/*first initlaize values to zero*/
				div = 0.0;
				p1 = 0.0;
				p2 = 0.0;
				p3 = 0.0;
				p4 = 0.0;
				
				if(j > 1) /*if node not on bottom edge*/
				{
					p1 = lsf[Nodes2[i][j-1]];
					div += 1.0;
				}
				
				if(j < NodeY-2) /*if node not on top edge*/
				{
					p3 = lsf[Nodes2[i][j+1]];
					div += 1.0;
				}
				
				if(i > 1) /*if node not on left edge*/
				{
					p4 = lsf[Nodes2[i-1][j]];
					div += 1.0;
				}
				
				if(i < NodeX-2) /*if node not on right edge*/
				{
					p2 = lsf[Nodes2[i+1][j]];
					div += 1.0;
				}
				
				lsf_temp[temp] = (p1 + p2 + p3 + p4) / div; /*linearly smooth lsf value*/
			}
			
			/*Otherwise just copy value*/
			else
			{
				lsf_temp[temp] = lsf[temp];
			}
		}
	}
	
	/*finally copy smoothed lsf to original lsf and free memory*/
	ArrayCopy(NumNodes, lsf_temp, lsf);
	free(lsf_temp);
}
				
/*Function to re-initalise the lsf as a signed distance function - similar to Vext above*/
void ReInt(int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int NumNodes, double *lsf, double h, int sign, 
			double *lsf_temp, int *fixed, int NumFix, int elemX, double hxz[elemX])
{

int i,j,k,count; /*incrementors*/
int Nnum,Nnum2,Nnum3; /*node numbers*/
int Ci,Cj; /*variables to store x,y location of node with minimum lsf_temp values in trial set*/
double lsf_min,NS; /*variable to track minimum lsf_temp value in trial set*/
double Af, Bf, Cf; /*variables for lsf_temp calculation by solution of a quadratic equn.*/
double ftemp, dtemp; /*tempory variables*/
double Lsum, Fsum; /*varibales to calcualte Vext by summing contributions from surrounding nodes*/

short *known;
known = calloc(NumNodes,sizeof(short));
short *trial;
trial = calloc(NumNodes,sizeof(short));

/*first ensure all nodes on domain boundary have lsf <= 0.0*/
/*bottom edge*/
for(i=1;i<NodeX-1;i++)
{
	Nnum = Nodes2[i][1]; /*read in node number*/
	if(lsf[Nnum] > -0.000001)
	{
		lsf[Nnum] = 0.0; /*set lsf value to zero if greater than zero*/
	}
}

/*top edge*/
k = NodeY-2;
for(i=1;i<NodeX-1;i++)
{
	Nnum = Nodes2[i][k]; /*read in node number*/
	if(lsf[Nnum] > -0.000001)
	{
		lsf[Nnum] = 0.0; /*set lsf value to zero if greater than zero*/
	}
}

/*left edge*/
for(i=1;i<NodeY-1;i++)
{
	Nnum = Nodes2[1][i]; /*read in node number*/
	if(lsf[Nnum] > -0.000001)
	{
		lsf[Nnum] = 0.0; /*set lsf value to zero if greater than zero*/
	}
}

/*right edge*/
k = NodeX-2;
for(i=1;i<NodeY-1;i++)
{
	Nnum = Nodes2[k][i]; /*read in node number*/
	if(lsf[Nnum] > -0.000001)
	{
		lsf[Nnum] = 0.0; /*set lsf value to zero if greater than zero*/
	}
}

/*Initalise Known set*/
for(i=0;i<NumNodes;i++)
{	
	/*If node is fixed or on boundary, then retain the lsf*/
	if(fabs(lsf[i]) < 0.000001)
	{
		known[i] = 1;
		lsf_temp[i] = 0.0;
	}
	/*else if(InSet(i,fixed,NumFix) == 1)
	{
		known[i] = 1;
		lsf_temp[i] = lsf[i]; /*add node to known set*/
		/*printf("\nNode %i = known",i+1);
	}*/
}

/*Initalise Trial set*/
for(j=1;j<NodeY-1;j++)
{
	for(i=1;i<NodeX-1;i++)
	{
		Nnum = Nodes2[i][j]; /*read in node number*/
		
		/*If current node not in known set and is the correct side of the boundary*/
		if((known[Nnum] != 1) && ((sign * lsf[Nnum]) > 0.0))	
		{
			NS = 1000.0 * lsf[Nnum];
			/*printf("\nNode %i, lsf=%lf",Nnum+1,NS);*/
			/*If any neighbouring node is on opposite side of the boundary, or on boundary, then it is a trial node*/
			if((lsf[Nodes2[i-1][j]] * NS) < 0.0) {
				/*printf("\n left lsf* = %lf",(lsf[Nodes2[i-1][j]] * NS));*/
				trial[Nnum] = 1;
				}
			else if((lsf[Nodes2[i+1][j]] * NS) < 0.0) {
				/*printf("\n right lsf* = %lf",(lsf[Nodes2[i-1][j]] * NS));*/
				trial[Nnum] = 1;
				}
			else if((lsf[Nodes2[i][j-1]] * NS) < 0.0) {
				/*printf("\n below lsf* = %lf",(lsf[Nodes2[i-1][j]] * NS));*/
				trial[Nnum] = 1;
				}
			else if((lsf[Nodes2[i][j+1]] * NS) < 0.0) {
				/*printf("\n above lsf* = %lf",(lsf[Nodes2[i-1][j]] * NS));*/
				trial[Nnum] = 1;
				}
		}
	}
}

/*printf("\nInital Trial Set");
for(i=0;i<NumNodes;i++)
{
	if(trial[i] == 1)
	{
		printf("\nNode %i",i+1);
	}
}*/

/*calculate lsf_temp for all current trial nodes*/
LocalInt(NodeX, NodeY, Nodes2, lsf, lsf_temp, known, trial, h, sign, elemX, hxz);

/*update trial set by considering neighbours of nodes in known set*/
for(j=1;j<NodeY-1;j++)
{
	for(i=1;i<NodeX-1;i++)
	{
		Nnum = Nodes2[i][j]; /*read in node number*/
		/*If current node not in known set and is the correct side of the boundary*/
		if((known[Nnum] != 1)  /*&& (InSet(Nnum,fixed,NumFix) == 0)*/ && ((sign * lsf[Nnum]) > 0.0))	
		{
			/*If any neighbouring node is in known set then curret node is a trial node*/
			if(known[Nodes2[i-1][j]] == 1) {
				trial[Nnum] = 1;
				}
			else if(known[Nodes2[i+1][j]] == 1) {
				trial[Nnum] = 1;
				}
			else if(known[Nodes2[i][j-1]] == 1) {
				trial[Nnum] = 1;
				}
			else if(known[Nodes2[i][j+1]] == 1) {
				trial[Nnum] = 1;
				}
		}
	}
}

/*printf("\nUpdated Trial Set");
for(i=0;i<NumNodes;i++)
{
	if(trial[i] == 1)
	{
		printf("\nNode %i",i+1);
	}
}
	
/*Do until trial set is empty*/
short flag;
double *lsf_trial;
lsf_trial = calloc(NumNodes,sizeof(double)); /*array to store trial values of the lsf_temp during the fast marching method*/

/*set two arrays to store nodes to be updatedduring each iteration of the fast marching method*/
int *xind, *yind, ncount;
xind = malloc(NumNodes * sizeof(int));
yind = malloc(NumNodes * sizeof(int));

count = 0;
double lsf_max = (NodeX > NodeY) ? (100.0 * h * NodeX) : (100.0 * h * NodeY);
do
{
	/*printf("\nUpdated Sets");
	for(i=0;i<NumNodes;i++)
	{
		if(trial[i] == 1)
		{
			printf("\nNode %i trial",i+1);
		}
		if(known[i] == 1)
		{
			printf("\nNode %i known",i+1);
			printf(" - lsf_temp = %lf",lsf_temp[i]);
		}
	}*/

	flag = 0;
	lsf_min = lsf_max; /*initalize to a value much bigger than the element edge length*/
	/*For all trial nodes calculate lsf temp*/
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++)
		{
			Nnum = Nodes2[i][j]; /*read in node number*/
			/*If current node not in known set, but is in the narrow band active set and is the correct side of the boundary*/
			if(trial[Nnum] == 1)
			{
				flag = 1;
				/*printf("\nNode %i is trial",Nnum+1);*/
				
				if(fabs(lsf_trial[Nnum]) < 0.000001)
				{	
					Af = 0.0;
					Bf = 0.0;
					Cf = 0.0; /*Initalise values*/
					ftemp = 0.0;
					
					/*look at all neighbouring nodes*/				
					Nnum2 = Nodes2[i][j-1];
					Nnum3 = Nodes2[i][j+1]; /*read in node number (above and below)*/
					
					/*If neigbouring node is in known set update co-efficents for lsf_temp calc*/
					/*ensure only closest node to boundary is used for upwind scheme*/
					if((known[Nnum2] == 1) && (known[Nnum3] == 1))
					{
						Af += 1.0;
						ftemp =  ( (fabs(lsf_temp[Nnum2])) < (fabs(lsf_temp[Nnum3])) ) ? lsf_temp[Nnum2] : lsf_temp[Nnum3]; /*choose closest node*/
					}
					else if(known[Nnum2] == 1)
					{
						Af += 1.0;
						ftemp = lsf_temp[Nnum2];
					}
					else if(known[Nnum3] == 1)
					{
						Af += 1.0;
						ftemp = lsf_temp[Nnum3];
					}

					Bf += sign * ftemp;
					Cf += ftemp * ftemp;	/*Update quadratic co-efficients*/
					
					Nnum2 = Nodes2[i-1][j];
					Nnum3 = Nodes2[i+1][j]; /*read in node number (left and right)*/
					ftemp = 0.0; /*re-initalise for other direction*/
					/*If neigbouring node is in known set update co-efficents for lsf_temp calc*/
					/*ensure only closest node to boundary is used for upwind scheme*/
					if((known[Nnum2] == 1) && (known[Nnum3] == 1))
					{
						Af += 1.0;
						ftemp =  ( (fabs(lsf_temp[Nnum2])) < (fabs(lsf_temp[Nnum3])) ) ? lsf_temp[Nnum2] : lsf_temp[Nnum3];  /*choose closest node*/
					}
					else if(known[Nnum2] == 1)
					{
						Af += 1.0;
						ftemp = lsf_temp[Nnum2];
					}
					else if(known[Nnum3] == 1)
					{
						Af += 1.0;
						ftemp = lsf_temp[Nnum3];
					}

					Bf += sign * ftemp;
					Cf += ftemp * ftemp;	/*Update quadratic co-efficients*/
					
					/*Need to do final update of the Bf and Cf co-efficients*/
					Bf *= 2.0;
					Cf -= h * h; 
					
					/*Calculate lsf_temp by solving the quadratic equation*/
					ftemp = Bf * Bf;
					ftemp -= 4.0 * Af * Cf;
					dtemp = sqrt(ftemp);
					ftemp = Bf + dtemp;
					lsf_trial[Nnum] = ftemp / (2.0 * Af);
					
					if( (lsf_trial[Nnum] < 0.000001) || (lsf_trial[Nnum] > lsf_max)) /*|| ((lsf_trial[Nnum] - lsf_max) > -0.00001) )*/
					{
						printf("\nERROR! In ReInt lsf_trial = %lf",lsf_trial[Nnum]);
						printf("\nAf=%lf,Bf=%lf,Cf=%lf",Af,Bf,Cf);
					}
				}
				
				/*If current lsf_temp distance is less than current minimum, update minimum and store node number*/
				if(lsf_trial[Nnum] < lsf_min)
				{
					lsf_min = lsf_trial[Nnum];
				}
			}
		}
	}
	
	if(flag == 0) /*if there are no trial nodes, then end the loop*/
	{
		break;
	}
	
	ncount = 0; /*re-initialize*/
	/*check to see which nodes need updating this iteration*/
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++)
		{
			Nnum = Nodes2[i][j]; /*read in node number*/
			
			if(lsf_trial[Nnum] < 0.0)
			{
				printf("\nERROR IN ReInt: lsf_trial = %lf, for node %i",lsf_trial[Nnum],Nnum+1);
			}
			
			/*if a trial node has a trail lsf value <= minimum then update*/
			if( (trial[Nnum] == 1) && (fabs(lsf_trial[Nnum] - lsf_min) < 0.000001) )
			{
				/*printf("\nTrial EX = %i, lsf=%lf",Nnum+1,lsf_min);*/
				lsf_temp[Nnum] = sign * lsf_trial[Nnum]; /*update lsf_temp for node with lsf_min*/
				
				/*Move node to known set*/
				/*printf("\nNode %i lsf_temp=%lf, Lsum=%lf, Fsum=%lf",Nnum+1,lsf_min,Lsum,Fsum);*/
				xind[ncount] = i;
				yind[ncount++] = j; /*store update node for next bit*/	
				count++; /*update number of re-initalized nodes*/
			}
		}
	}
	
	/*printf("\nCount = %i, lsf_min=%f",count,lsf_min);*/
	
	/*update known set for next iteration*/
	for(i=0;i<ncount;i++)
	{
		Ci = xind[i];
		Cj = yind[i]; /*read in updated node indicators*/
		
		/*read in node number*/				
		Nnum = Nodes2[Ci][Cj];
		
		trial[Nnum] = 0;
		known[Nnum] = 1; /*update known set*/
	}
	
	/*update trial set for next iteration*/
	for(i=0;i<ncount;i++)
	{
		Ci = xind[i];
		Cj = yind[i]; /*read in updated node indicators*/
		
		/*look at all neighbouring nodes*/				
		Nnum2 = Nodes2[Ci][Cj-1];
		Nnum3 = Nodes2[Ci][Cj+1]; /*read in node number (above and below)*/
				
		/*Check to see if node below needs to be added to trial set*/
		if((known[Nnum2] != 1) && /*(InSet(Nnum2,fixed,NumFix) == 0) &&*/ ((sign * lsf[Nnum2]) > 0.0))
		{
			trial[Nnum2] = 1;
			lsf_trial[Nnum2] = 0.0; /*reset trial value*/
		}

		/*Check to see if node above needs to be added to trial set*/
		if((known[Nnum3] != 1) &&  ((sign * lsf[Nnum3]) > 0.0))
		{
			trial[Nnum3] = 1;
			lsf_trial[Nnum3] = 0.0; /*reset trial value*/
		}
		
		Nnum2 = Nodes2[Ci-1][Cj];
		Nnum3 = Nodes2[Ci+1][Cj]; /*read in node number (left and right)*/

		/*Check to see if node to left needs to be added to trial set*/
		if((known[Nnum2] != 1) && ((sign * lsf[Nnum2]) > 0.0))
		{
			trial[Nnum2] = 1;
			lsf_trial[Nnum2] = 0.0; /*reset trial value*/
		}

		/*Check to see if node to right needs to be added to trial set*/
		if((known[Nnum3] != 1) && ((sign * lsf[Nnum3]) > 0.0))
		{
			trial[Nnum3] = 1;
			lsf_trial[Nnum3] = 0.0; /*reset trial value*/
		}
	}

	if(ncount == 0)
	{
		printf("\nERROR! ncount=0 in ReInt! Aborting");
	}
}
while(ncount > 0);

free(known);
free(trial);
free(lsf_trial);
free(xind);
free(yind);
		
}

/*Function to re-initalise the lsf as a signed distance function - similar to Vext above*/
void ReIntlsf2(int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int NumNodes, double *lsf, double h, int sign, 
			double *lsf_temp, int *fixed, int NumFix, int elemX, double hxz[elemX])
{

int i,j,k,count; /*incrementors*/
int Nnum,Nnum2,Nnum3; /*node numbers*/
int Ci,Cj; /*variables to store x,y location of node with minimum lsf_temp values in trial set*/
double lsf_min,NS; /*variable to track minimum lsf_temp value in trial set*/
double Af, Bf, Cf; /*variables for lsf_temp calculation by solution of a quadratic equn.*/
double ftemp, dtemp; /*tempory variables*/
double Lsum, Fsum; /*varibales to calcualte Vext by summing contributions from surrounding nodes*/
double h3;

short *known;
known = calloc(NumNodes,sizeof(short));
short *trial;
trial = calloc(NumNodes,sizeof(short));

/*first ensure all nodes on domain boundary have lsf <= 0.0*/
/*bottom edge*
for(i=1;i<NodeX-1;i++)
{
	Nnum = Nodes2[i][1]; /*read in node number*
	if(lsf[Nnum] > -0.000001)
	{
		lsf[Nnum] = 0.0; /*set lsf value to zero if greater than zero*
	}
}

/*top edge*
k = NodeY-2;
for(i=1;i<NodeX-1;i++)
{
	Nnum = Nodes2[i][k]; /*read in node number*
	if(lsf[Nnum] > -0.000001)
	{
		lsf[Nnum] = 0.0; /*set lsf value to zero if greater than zero*
	}
}

/*left edge*
for(i=1;i<NodeY-1;i++)
{
	Nnum = Nodes2[1][i]; /*read in node number*
	if(lsf[Nnum] > -0.000001)
	{
		lsf[Nnum] = 0.0; /*set lsf value to zero if greater than zero*
	}
}

/*right edge*
k = NodeX-2;
for(i=1;i<NodeY-1;i++)
{
	Nnum = Nodes2[k][i]; /*read in node number*
	if(lsf[Nnum] > -0.000001)
	{
		lsf[Nnum] = 0.0; /*set lsf value to zero if greater than zero*
	}
}*/

/*Initalise Known set*/
for(i=0;i<NumNodes;i++)
{	
	/*If node is fixed or on boundary, then retain the lsf*/
	if(fabs(lsf[i]) < 0.000001)
	{
		known[i] = 1;
		lsf_temp[i] = 0.0;
	}
	/*else if(InSet(i,fixed,NumFix) == 1)
	{
		known[i] = 1;
		lsf_temp[i] = lsf[i]; /*add node to known set*/
		/*printf("\nNode %i = known",i+1);
	}*/
}

/*Initalise Trial set*/
for(j=1;j<NodeY-1;j++)
{
	for(i=1;i<NodeX-1;i++)
	{
		Nnum = Nodes2[i][j]; /*read in node number*/
		
		/*If current node not in known set and is the correct side of the boundary*/
		if((known[Nnum] != 1) && ((sign * lsf[Nnum]) > 0.0))	
		{
			NS = 1000.0 * lsf[Nnum];
			/*printf("\nNode %i, lsf=%lf",Nnum+1,NS);*/
			/*If any neighbouring node is on opposite side of the boundary, or on boundary, then it is a trial node*/
			if((lsf[Nodes2[i-1][j]] * NS) < 0.0) {
				/*printf("\n left lsf* = %lf",(lsf[Nodes2[i-1][j]] * NS));*/
				trial[Nnum] = 1;
				}
			else if((lsf[Nodes2[i+1][j]] * NS) < 0.0) {
				/*printf("\n right lsf* = %lf",(lsf[Nodes2[i-1][j]] * NS));*/
				trial[Nnum] = 1;
				}
			else if((lsf[Nodes2[i][j-1]] * NS) < 0.0) {
				/*printf("\n below lsf* = %lf",(lsf[Nodes2[i-1][j]] * NS));*/
				trial[Nnum] = 1;
				}
			else if((lsf[Nodes2[i][j+1]] * NS) < 0.0) {
				/*printf("\n above lsf* = %lf",(lsf[Nodes2[i-1][j]] * NS));*/
				trial[Nnum] = 1;
				}
		}
	}
}

/*printf("\nInital Trial Set");

for(i=0;i<NumNodes;i++)
{
	if(trial[i] == 1)
	{

		printf("\nNode %i",i+1);
	}
}*/

/*calculate lsf_temp for all current trial nodes*/
LocalInt(NodeX, NodeY, Nodes2, lsf, lsf_temp, known, trial, h, sign, elemX, hxz);

/*update trial set by considering neighbours of nodes in known set*/
for(j=1;j<NodeY-1;j++)
{
	for(i=1;i<NodeX-1;i++)
	{
		Nnum = Nodes2[i][j]; /*read in node number*/
		/*If current node not in known set and is the correct side of the boundary*/
		if((known[Nnum] != 1)  /*&& (InSet(Nnum,fixed,NumFix) == 0)*/ && ((sign * lsf[Nnum]) > 0.0))	
		{
			/*If any neighbouring node is in known set then curret node is a trial node*/
			if(known[Nodes2[i-1][j]] == 1) {
				trial[Nnum] = 1;
				}
			else if(known[Nodes2[i+1][j]] == 1) {
				trial[Nnum] = 1;
				}
			else if(known[Nodes2[i][j-1]] == 1) {
				trial[Nnum] = 1;
				}
			else if(known[Nodes2[i][j+1]] == 1) {
				trial[Nnum] = 1;
				}
		}
	}
}

/*printf("\nUpdated Trial Set");
for(i=0;i<NumNodes;i++)
{
	if(trial[i] == 1)

	{
		printf("\nNode %i",i+1);
	}
}
	

/*Do until trial set is empty*/
short flag;
double *lsf_trial;
lsf_trial = calloc(NumNodes,sizeof(double)); /*array to store trial values of the lsf_temp during the fast marching method*/

/*set two arrays to store nodes to be updatedduring each iteration of the fast marching method*/
int *xind, *yind, ncount;
xind = malloc(NumNodes * sizeof(int));
yind = malloc(NumNodes * sizeof(int));

count = 0;
double lsf_max = (NodeX > NodeY) ? (100.0 * h * NodeX) : (100.0 * h * NodeY);
do
{
	/*printf("\nUpdated Sets");
	for(i=0;i<NumNodes;i++)

	{
		if(trial[i] == 1)
		{
			printf("\nNode %i trial",i+1);

		}
		if(known[i] == 1)
		{
			printf("\nNode %i known",i+1);
			printf(" - lsf_temp = %lf",lsf_temp[i]);

		}
	}*/

	flag = 0;
	lsf_min = lsf_max; /*initalize to a value much bigger than the element edge length*/
	/*For all trial nodes calculate lsf temp*/
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++)
		{
			Nnum = Nodes2[i][j]; /*read in node number*/
			/*If current node not in known set, but is in the narrow band active set and is the correct side of the boundary*/
			if(trial[Nnum] == 1)
			{
				flag = 1;
				/*printf("\nNode %i is trial",Nnum+1);*/
				
				if(fabs(lsf_trial[Nnum]) < 0.000001)
				{	
					Af = 0.0;
					Bf = 0.0;
					Cf = 0.0; /*Initalise values*/
					ftemp = 0.0;
					
					/*look at all neighbouring nodes*/				
					Nnum2 = Nodes2[i][j-1];
					Nnum3 = Nodes2[i][j+1]; /*read in node number (above and below)*/
					
					/*If neigbouring node is in known set update co-efficents for lsf_temp calc*/
					/*ensure only closest node to boundary is used for upwind scheme*/
					if((known[Nnum2] == 1) && (known[Nnum3] == 1))
					{
						Af += 1.0;
						ftemp =  ( (fabs(lsf_temp[Nnum2])) < (fabs(lsf_temp[Nnum3])) ) ? lsf_temp[Nnum2] : lsf_temp[Nnum3]; /*choose closest node*/
						h3 = h;
					}
					else if(known[Nnum2] == 1)
					{
						Af += 1.0;
						ftemp = lsf_temp[Nnum2];
						h3 = h;
					}
					else if(known[Nnum3] == 1)
					{
						Af += 1.0;
						ftemp = lsf_temp[Nnum3];
						h3 = h;
					}

					Bf += sign * ftemp;
					Cf += ftemp * ftemp;	/*Update quadratic co-efficients*/
					
					Nnum2 = Nodes2[i-1][j];
					Nnum3 = Nodes2[i+1][j]; /*read in node number (left and right)*/
					ftemp = 0.0; /*re-initalise for other direction*/
					/*If neigbouring node is in known set update co-efficents for lsf_temp calc*/
					/*ensure only closest node to boundary is used for upwind scheme*/
					if((known[Nnum2] == 1) && (known[Nnum3] == 1))
					{
						Af += 1.0;
						ftemp =  ( (fabs(lsf_temp[Nnum2])) < (fabs(lsf_temp[Nnum3])) ) ? lsf_temp[Nnum2] : lsf_temp[Nnum3];  /*choose closest node*/
						h3 = ( (fabs(lsf_temp[Nnum2])) < (fabs(lsf_temp[Nnum3])) ) ? hxz[i-2]:hxz[i-1];
					}
					else if(known[Nnum2] == 1)
					{
						Af += 1.0;
						ftemp = lsf_temp[Nnum2];
						h3 = hxz[i-2];
					}
					else if(known[Nnum3] == 1)
					{
						Af += 1.0;
						ftemp = lsf_temp[Nnum3];
						h3 = hxz[i-1];
					}

					Bf += sign * ftemp;
					Cf += ftemp * ftemp;	/*Update quadratic co-efficients*/
					
					/*Need to do final update of the Bf and Cf co-efficients*/
					Bf *= 2.0;
					Cf -= h3 * h3; 
					
					/*Calculate lsf_temp by solving the quadratic equation*/
					ftemp = Bf * Bf;
					ftemp -= 4.0 * Af * Cf;
					dtemp = sqrt(ftemp);
					ftemp = Bf + dtemp;
					lsf_trial[Nnum] = ftemp / (2.0 * Af);
					
					if( (lsf_trial[Nnum] < 0.000001) || (lsf_trial[Nnum] > lsf_max)) /*|| ((lsf_trial[Nnum] - lsf_max) > -0.00001) )*/
					{
						printf("\nERROR! In ReInt lsf_trial = %lf",lsf_trial[Nnum]);
						printf("\nAf=%lf,Bf=%lf,Cf=%lf",Af,Bf,Cf);
					}
				}
				
				/*If current lsf_temp distance is less than current minimum, update minimum and store node number*/
				if(lsf_trial[Nnum] < lsf_min)
				{
					lsf_min = lsf_trial[Nnum];
				}
			}
		}
	}
	
	if(flag == 0) /*if there are no trial nodes, then end the loop*/
	{
		break;
	}
	
	ncount = 0; /*re-initialize*/
	/*check to see which nodes need updating this iteration*/
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++)
		{
			Nnum = Nodes2[i][j]; /*read in node number*/
			
			if(lsf_trial[Nnum] < 0.0)
			{
				printf("\nERROR IN ReInt: lsf_trial = %lf, for node %i",lsf_trial[Nnum],Nnum+1);
			}
			
			/*if a trial node has a trail lsf value <= minimum then update*/
			if( (trial[Nnum] == 1) && (fabs(lsf_trial[Nnum] - lsf_min) < 0.000001) )
			{
				/*printf("\nTrial EX = %i, lsf=%lf",Nnum+1,lsf_min);*/
				lsf_temp[Nnum] = sign * lsf_trial[Nnum]; /*update lsf_temp for node with lsf_min*/
				
				/*Move node to known set*/
				/*printf("\nNode %i lsf_temp=%lf, Lsum=%lf, Fsum=%lf",Nnum+1,lsf_min,Lsum,Fsum);*/
				xind[ncount] = i;
				yind[ncount++] = j; /*store update node for next bit*/	
				count++; /*update number of re-initalized nodes*/
			}
		}
	}
	
	/*printf("\nCount = %i, lsf_min=%f",count,lsf_min);*/
	
	/*update known set for next iteration*/
	for(i=0;i<ncount;i++)
	{
		Ci = xind[i];
		Cj = yind[i]; /*read in updated node indicators*/
		
		/*read in node number*/				
		Nnum = Nodes2[Ci][Cj];
		
		trial[Nnum] = 0;
		known[Nnum] = 1; /*update known set*/
	}
	
	/*update trial set for next iteration*/
	for(i=0;i<ncount;i++)
	{
		Ci = xind[i];
		Cj = yind[i]; /*read in updated node indicators*/
		
		/*look at all neighbouring nodes*/				
		Nnum2 = Nodes2[Ci][Cj-1];
		Nnum3 = Nodes2[Ci][Cj+1]; /*read in node number (above and below)*/
				
		/*Check to see if node below needs to be added to trial set*/
		if((known[Nnum2] != 1) && /*(InSet(Nnum2,fixed,NumFix) == 0) &&*/ ((sign * lsf[Nnum2]) > 0.0))
		{
			trial[Nnum2] = 1;
			lsf_trial[Nnum2] = 0.0; /*reset trial value*/
		}

		/*Check to see if node above needs to be added to trial set*/
		if((known[Nnum3] != 1) &&  ((sign * lsf[Nnum3]) > 0.0))
		{
			trial[Nnum3] = 1;
			lsf_trial[Nnum3] = 0.0; /*reset trial value*/
		}
		
		Nnum2 = Nodes2[Ci-1][Cj];
		Nnum3 = Nodes2[Ci+1][Cj]; /*read in node number (left and right)*/

		/*Check to see if node to left needs to be added to trial set*/
		if((known[Nnum2] != 1) && ((sign * lsf[Nnum2]) > 0.0))
		{
			trial[Nnum2] = 1;
			lsf_trial[Nnum2] = 0.0; /*reset trial value*/
		}

		/*Check to see if node to right needs to be added to trial set*/
		if((known[Nnum3] != 1) && ((sign * lsf[Nnum3]) > 0.0))
		{
			trial[Nnum3] = 1;
			lsf_trial[Nnum3] = 0.0; /*reset trial value*/
		}
	}

	if(ncount == 0)
	{
		printf("\nERROR! ncount=0 in ReInt! Aborting");
	}
}
while(ncount > 0);

free(known);
free(trial);
free(lsf_trial);
free(xind);
free(yind);
		
}

/*Function to reinitalise lsf for inital set of trial nodes - similar to LocalVext above*/
void LocalInt(int NodeX, int NodeY, int Nodes[NodeX][NodeY], double *lsf, double *lsf_temp, short *known, short *trial, double h, int sign, int elemX, double hxz[elemX])
{
	
	double s,s1,s2,t,t1,t2; /*distance variables*/
	int i,j,k; /*incrementors*/
	double ftemp; /*temp varibale*/
	int Nnum,Nnum2; /*node number variables*/
	double h2 = h * 0.5;
	double h2xz[elemX];
	double h3xz;
	for(i=0;i<elemX;i++)
	{
		h2xz[i] = hxz[i]*0.5;
	}
	
	for(j=1;j<NodeY-1;j++)
	{
		for(i=1;i<NodeX-1;i++)
		{
			Nnum = Nodes[i][j];
			if(trial[Nnum] == 1)
			{
				/*look at each neighbouring node in turn*/
				
				/*-----------look at node below-----------*/
				Nnum2 = Nodes[i][j-1];
				if(fabs(lsf[Nnum2]) < 0.000001) /*If neigbouring node is on the boundary*/
				{
					t1 = h; /*distance to boundary = element edge length*/
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) /*otherwise check to see if grid line is intersected*/
				{
					/*t1 = sign * lsf[Nnum];*/
					t1 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					t1 += 1.0;
					t1 *= h2;
					/*check to see if t1 > h, i.e. something has gone wrong in update*/
					/*If so then approximate*/
					if((t1 - h) > -0.000001)
					{
						/*printf("\nERROR! t1(%lf) > h",t1);*/
						if( ((sign * -lsf[Nnum2]) - h) > -0.000001 )
						{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							t1 = (lsf[Nnum] * h) / ftemp;
							/*printf("\tApproxed 1 as %lf",t1);*/
						}
						else
						{
							t1 = (sign * lsf[Nnum2]) + h;
							/*printf("\tApproxed 2 as %lf",t1);*/
						}
						
						printf("\nERROR in LocalInt for node %i, t1 approx as %lf",Nnum+1,t1);
					}
				}
				else
				{
					t1 = h + 1.0; /*set distance to h+1.0 to make algorithm work later*/
				}
				
				/*----------look at node above----------*/
				Nnum2 = Nodes[i][j+1];
				if(fabs(lsf[Nnum2]) < 0.000001) /*If neigbouring node is on the boundary*/
				{
					t2 = h; /*distance to boundary = element edge length*/
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) /*otherwise check to see is grid line is intersected*/
				{
					/*t2 = sign * lsf[Nnum];*/
					t2 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					t2 += 1.0;
					t2 *= h2;
					/*check to see if t2 > h, i.e. something has gone wrong in update*/
					/*If so then approximate*/
					if((t2 - h) > -0.000001)
					{
						/*printf("\nERROR! t2(%lf) > h",t2);*/
						if( ((sign * -lsf[Nnum2]) - h) > -0.000001 )
						{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							t2 = (lsf[Nnum] * h) / ftemp;
							/*printf("\tApproxed 1 as %lf",t2);*/
						}
						else
						{
							t2 = (sign * lsf[Nnum2]) + h;
							/*printf("\tApproxed 2 as %lf",t2);*/
						}
						
						printf("\nERROR in LocalInt for node %i, t2 approx as %lf",Nnum+1,t2);
					}
				}
				else
				{
					t2 = h + 1.0; /*set distance to h+1.0 to make algorithm work later*/
				}
				
				/*----------look at node to left----------*/
				Nnum2 = Nodes[i-1][j];
				if(fabs(lsf[Nnum2]) < 0.000001) /*If neigbouring node is on the boundary*/
				{
					s1 = hxz[i-2]; /*distance to boundary = element edge length*/
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) /*otherwise check to see is grid line is intersected*/
				{
					/*s1 = sign * lsf[Nnum];*/
					s1 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					s1 += 1.0;
					s1 *= h2xz[i-2];
					/*check to see if s1 > h, i.e. something has gone wrong in update*/
					/*If so then approximate*/
					if((s1 - hxz[i-2]) > -0.000001)
					{
						/*printf("\nERROR! s1(%lf) > h",s1);*/
						if( ((sign * -lsf[Nnum2]) - hxz[i-2]) > -0.000001 )
						{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							s1 = (lsf[Nnum] * hxz[i-2]) / ftemp;
							/*printf("\tApproxed 1 as %lf",s1);*/
						}
						else
						{
							s1 = (sign * lsf[Nnum2]) + hxz[i-2];
							/*printf("\tApproxed 2 as %lf",s1);*/
						}
						
						printf("\nERROR in LocalInt for node %i, s1 approx as %lf",Nnum+1,s1);
					}
				}
				else
				{
					s1 = hxz[i-2] + 2.0; /*set distance to h+1.0 to make algorithm work later*/
				}
				
				/*----------look at node to right----------*/
				Nnum2 = Nodes[i+1][j];
				if(fabs(lsf[Nnum2]) < 0.000001) /*If neigbouring node is on the boundary*/
				{
					s2 = hxz[i-1]; /*distance to boundary = element edge length*/
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) /*otherwise check to see is grid line is intersected*/
				{
					/*s2 = sign * lsf[Nnum];*/
					s2 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					s2 += 1.0;
					s2 *= h2xz[i-1];
					/*check to see if s2 > h, i.e. something has gone wrong in update*/
					/*If so then approximate*/
					if((s2 - hxz[i-1]) > -0.000001)
					{
						/*printf("\nERROR! s2(%lf) > h",s2);*/
						if( ((sign * -lsf[Nnum2]) - hxz[i-1]) > -0.000001 )
						{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							s2 = (lsf[Nnum] * hxz[i-1]) / ftemp;
							/*printf("\tApproxed 1 as %lf",s2);*/
						}
						else
						{
							s2 = (sign * lsf[Nnum2]) + hxz[i-1];
							/*printf("\tApproxed 2 as %lf",s2);*/
						}
						
						printf("\nERROR in LocalInt for node %i, s2 approx as %lf",Nnum+1,s2);
					}

				}
				else
				{
					s2 = hxz[i-1] + 2.0; /*set distance to h+1.0 to make algorithm work later*/
				}
				
				/*---------Now calculate lsf temp for node---------*/
				/*choose lowest s and t values, i.e. nearest to boundary*/				
				s = (s1 < s2) ? s1 : s2;
				h3xz = (s1 < s2) ? hxz[i-2]:hxz[i-1];
				t = (t1 < t2) ? t1 : t2;

				/*Now calculate lsf temp*/
				if((s - h3xz) > 0.9) /*If horizontal distance > edge length only use vertical distance*/
				{
					lsf_temp[Nnum] = (double)sign * t;
				}
				
				else if((t - h) > 0.9)	/*If vertical distance > edge length only use horizontal distance*/
				{
					lsf_temp[Nnum] = (double)sign * s;
				}
				
				else
				{		
					s2 = s * s;
					t2 = t * t; /*squared distances*/
					
					ftemp = s * t;
					ftemp /= sqrt(s2 + t2);
					lsf_temp[Nnum] = (double)sign * ftemp; /*calculate perpendicular distance to the boundary (by pythag)*/
				}
				
				if((s > h3xz) && (t > h)) /*print error if both horizontal and vertical distances greater than edge length*/
				{
					printf("\nError! For node %i, both s & t > h in LocalInt, s=%lf, t=%lf", Nnum+1,s,t);
				}
				
				/*Move node into known set*/
				trial[Nnum] = 0;
				known[Nnum] = 1;
			}
		}
	}
	/*printf("\nLocal Known Set");
	/*for(i=0;i<NumNodes;i++)
	{
		if(trial[i] == 1)
		{
			printf("\nNode %i",i+1);
		}
	}*/
}

/*function to determine if a number is in a set*/
int InSet(int num, int *array, int length)
{
	int i;
	
	for(i=0;i<length;i++)
	{
		if(array[i] == num) /*if number in set return 1*/
		{
			return(1);
		}
	}
	
	return(0); /*ifnumber not in set return 0*/
}


