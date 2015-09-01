/*
 *  FixedGrid.c
 *  
 *  Created by Peter Dunning on 26/02/2007.
 *	Modified for BLES program on 13/02/2009
 *
 */

#include "FixedGrid.h"
#include "Strain.h"
#include "ls_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Function to calcualte lsf around a rectangular hole*/
void RectHole(int NumNodes, Coord NodeCoord[NumNodes], int NumRect, Coord *Rect, double *lsf, int elemX, int elemY, double hxz[elemX], double maxXZ)
{
	int i,j,k,o;
	double xtemp,ytemp,ftemp,dist;
	double minX, minY, maxX, maxY;
	double Thxz[elemX+1];
	Thxz[0] = 0;
	/*printf("\nThxz[0] = %f", Thxz[0]);*/
	for(o=1;o<elemX+1;o++)
	{	
		Thxz[o] = Thxz[o-1] + hxz[o-1];
		/*printf("\nThxz[%i] = %f", o, Thxz[o]);*/
	}
	
	/*---------For all rectangular holes change signed distance fucntion if less that it current value---------*/
	for(i=0;i<NumRect;i++)
	{
		/*read in rectangular hole dimensions*/
		k = 2 * i;

		for(o=0;o<elemX+1;o++)
		{
			/*Get minx in XZ corrds (so around the plate)*/
			if(Rect[k].x<0.0)
			{
				minX = Rect[k].x;
			}
			else if(Rect[k].x>maxXZ)
			{
				minX = (Rect[k].x - elemX) + maxXZ;
			}
			else if(Rect[k].x == o)
			{
				minX = Thxz[o];
			}
			else if((Rect[k].x > o)&&(Rect[k].x < o+1))
			{
				minX = Thxz[o]+(Rect[k].x-o)*hxz[o];
			}
			/*Get maxx in XZ corrds (so around the plate)*/
			if(Rect[k+1].x<0.0)
			{
				maxX = Rect[k+1].x;
			}
			else if(Rect[k+1].x>maxXZ)
			{
				maxX = (Rect[k+1].x - elemX) + maxXZ;
			}
			else if(Rect[k+1].x == o)
			{
				maxX = Thxz[o];
			}
			else if((Rect[k+1].x > o)&&(Rect[k+1].x < o+1))
			{
				maxX = Thxz[o]+(Rect[k+1].x-o)*hxz[o];
			}
		}
		
		/*printf("\nminX = %f", minX);
		printf("\nmaxX = %f", maxX);*/


		minY = Rect[k].y;
		maxY = Rect[k+1].y;
		
		for(j=0;j<NumNodes;j++)
		{
			/*read in node coords*/
			xtemp = NodeCoord[j].xz;
			ytemp = NodeCoord[j].y;
			
			/*if node to left or right of the rectangle*/
			if(((xtemp > maxX) || (xtemp < minX)) && (ytemp <= maxY) && (ytemp >= minY))
			{
				dist = minX - xtemp;
				ftemp = xtemp - maxX;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
			}
			
			/*else if node above or below the rectangle*/
			else if(((ytemp > maxY) || (ytemp < minY)) && (xtemp <= maxX) && (xtemp >= minX))
			{
				dist = minY - ytemp;
				ftemp = ytemp - maxY;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
			}
			
			/*If node inside the rectagle, find closest edge*/
			else if((xtemp <= maxX) && (xtemp >= minX) && (ytemp <= maxY) && (ytemp >= minY))
			{
				dist = minX - xtemp;
				ftemp = minY - ytemp;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
				ftemp = minY - ytemp;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
				ftemp = xtemp - maxX;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
				ftemp = ytemp - maxY;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
			}
			
			/*if node to bottom left of rectangle*/
			else if ((xtemp < minX) && (ytemp < minY))
			{
				ftemp = (minX - xtemp) * (minX - xtemp);
				ftemp += (minY - ytemp) * (minY - ytemp);
				dist = sqrt(ftemp);
			}
			
			/*if node to bottom right of rectangle*/
			else if ((xtemp > maxX) && (ytemp < minY))
			{
				ftemp = (maxX - xtemp) * (maxX - xtemp);
				ftemp += (minY - ytemp) * (minY - ytemp);
				dist = sqrt(ftemp);
			}
			
			/*if node to top left of rectangle*/
			else if ((xtemp < minX) && (ytemp > maxY))
			{
				ftemp = (minX - xtemp) * (minX - xtemp);
				ftemp += (maxY - ytemp) * (maxY - ytemp);
				dist = sqrt(ftemp);
			}
			
			/*if node to top right of rectangle*/
			else if ((xtemp > maxX) && (ytemp > maxY))
			{
				ftemp = (maxX - xtemp) * (maxX - xtemp);
				ftemp += (maxY - ytemp) * (maxY - ytemp);
				dist = sqrt(ftemp);
			}
			
			else
			{
				printf("\nERROR! cant locate node %i in RectHole function!",j+1);
			}
			
			/*If signed distance to hole edge is less than current value, then update*/
			if(dist < lsf[j])
			{
				lsf[j] = dist;
			}
		}
	}
}

/*Function that caculates the area of any Polygon
/ NB: vertices have to be numbered anti-clockwise*/
double PolyArea(int N,Coord *point)	
{
	int i,j,temp;
	double area = 0.0;
	double pix,pjx;
	/*printf("\nPoints recieved by POLYGON func:  \n");*/
	/*for(i=0;i<N;i++)
	{
		printf("point%i xz = %lf, y = %lf, x = %lf, z = %lf\n",i, point[i].xz, point[i].x, point[i].y, point[i].z);
	}*/
		
	for (i=0;i<N;i++) /*For all verticies of the polygon*/
	{
		j = (i == N-1) ? 0 : (i+1);	/*Need to loop back to first point to complete algorithm*/
		/*temp = ((point[i].x * point[i].x) + (point[i].z * point[i].z));
		pix = sqrt(temp);
		temp = ((point[j].x * point[j].x) + (point[j].z * point[j].z));
		pjx = sqrt(temp);*/
		area += point[i].xz * point[j].y;
		area -= point[j].xz * point[i].y;
		/*printf("\ni=%i, j=%i Current AREA = %lf",i,j,area);*/
	}
	
	area *= 0.5;
	/*printf("POLYGON func called: AREA = %lf  N = %i\n\n", area, N);*/
	return(fabs(area)); /*return absolute value just in case*/
}

double PolyArea2(double *lsf, int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], short *NodeStat, short *ElemStat, double *alpha0, double *hxz, double h, double *alpha)
{
	int i,j,temp,N,num,p,q,count,n,m;
	double dtemp;
	double area;
	double pix,pjx;
	Coord *point;
	point = malloc(5*sizeof(Coord));
	int *tnode;
	double *Tlsf;
	tnode = malloc(4*sizeof(int));
	Tlsf = malloc(4*sizeof(double));
	double aTotal = 0;

	for(n=0;n<elemX;n++)
	{
		for(m=0;m<elemY;m++)
		{
			num = Number[n][m].n -1;
			if(ElemStat[num]==0){area = 0.0;}
			else if(ElemStat[num]==4){area = 1.0;}
			else
			{
				area = 0.0;
				N = ElemStat[num]+2;
				area = 0;
				tnode[0] = Number[n][m].a -1;
				tnode[1] = Number[n][m].b -1;
				tnode[2] = Number[n][m].c -1;
				tnode[3] = Number[n][m].d -1;
	
				Tlsf[0] = lsf[tnode[0]];
				Tlsf[1] = lsf[tnode[1]];
				Tlsf[2] = lsf[tnode[2]];
				Tlsf[3] = lsf[tnode[3]];
				
				
				count = 0;
				for(p=0;p<4;p++)	/*Check the status of each edge of the elements*/
				{
					q = (p==3) ? 0:p+1;

					if(((NodeStat[tnode[p]]==1)&&(NodeStat[tnode[q]]==1))||((NodeStat[tnode[p]]+NodeStat[tnode[q]])>2))	/*Edge is "in"*/
					{
						point[count].xz = NodeCoord[tnode[p]].xz;
						point[count].y = NodeCoord[tnode[p]].y;
						count++;
					}
					else if((NodeStat[tnode[p]]+NodeStat[tnode[q]])==1)	/*Edge is "cut*/
					{
						if(NodeStat[tnode[p]] == 1)
						{
							point[count].xz = NodeCoord[tnode[p]].xz;
							point[count].y = NodeCoord[tnode[p]].y;
							count++;
						}

						if(p==0)
						{
							point[count].y = NodeCoord[tnode[p]].y;
							/*printf("Tlsf[%i] - Tlsf[%i] = %f - %f = %f", p, q, Tlsf[p], Tlsf[q], Tlsf[p] - Tlsf[q]);
							printf("\nXz");*/
							/*if(num==382){printf("\nCut p = %i", p);}*/
							dtemp = (Tlsf[p])/(Tlsf[p] - Tlsf[q]);
							point[count].xz = NodeCoord[tnode[p]].xz + hxz[n]*dtemp;
							count++;
						}	
						if(p==1)
						{
							point[count].xz = NodeCoord[tnode[p]].xz;
							/*printf("Tlsf[%i] - Tlsf[%i] = %f - %f = %f", p, q, Tlsf[p], Tlsf[q], Tlsf[p] - Tlsf[q]);
							printf("\nXz");*/
							/*if(num==382){printf("\nCut p = %i", p);}*/
							dtemp = (Tlsf[p])/(Tlsf[p] - Tlsf[q]);
							point[count].y = NodeCoord[tnode[p]].y + h*dtemp;
							count++;
						}
						if(p==2)
						{
							point[count].y = NodeCoord[tnode[p]].y;
							/*printf("Tlsf[%i] - Tlsf[%i] = %f - %f = %f", p, q, Tlsf[p], Tlsf[q], Tlsf[p] - Tlsf[q]);
							printf("\nXz");*/
							/*if(num==382){printf("\nCut p = %i", p);}*/
							dtemp = (Tlsf[q])/(Tlsf[q] - Tlsf[p]);
							/*printf("\tcount = %i  dtemp = %f   Tlsf[%i] = %f   Tlsf[%i] = %f", count, dtemp, p, Tlsf[p], q, Tlsf[q]);*/
							point[count].xz = NodeCoord[tnode[q]].xz + hxz[n]*dtemp;
							count++;
						}
						if(p==3)
						{
							point[count].xz = NodeCoord[tnode[p]].xz;
							/*printf("Tlsf[%i] - Tlsf[%i] = %f - %f = %f", p, q, Tlsf[p], Tlsf[q], Tlsf[p] - Tlsf[q]);
							printf("\nXz");*/
							/*if(num==382){printf("\nCut p = %i", p);}*/
							dtemp = (Tlsf[q])/(Tlsf[q] - Tlsf[p]);
							/*if(num==382){printf("\tcount = %i  dtemp = %f   Tlsf[%i] = %f   Tlsf[%i] = %f", count, dtemp, p, Tlsf[p], q, Tlsf[q]);}*/
							point[count].y = NodeCoord[tnode[q]].y + h*dtemp;
							count++;
						}				
					}
					else if((NodeStat[tnode[p]]==2)&&((NodeStat[tnode[p]]+NodeStat[tnode[q]])==2))	/*Edge is touching at one point*/
					{
						point[count].xz = NodeCoord[tnode[p]].xz;
						point[count].y = NodeCoord[tnode[p]].y;
						count++;
					}
				}
				if(count!=N){printf("\nERROR %i points selected when %i points needed!!!!", count, N);}
	
				/*printf("\nPoints recieved by POLYGON func:  \n");*/
				/*if(num==382)
				{
				for(i=0;i<N;i++)
				{
					printf("\npoint%i xz = %lf, y = %lf",i, point[i].xz, point[i].y);
				}
				}*/
					
				for (i=0;i<N;i++) /*For all verticies of the polygon*/
				{
					j = (i == N-1) ? 0 : (i+1);	/*Need to loop back to first point to complete algorithm*/
					/*temp = ((point[i].x * point[i].x) + (point[i].z * point[i].z));
			
					pix = sqrt(temp);
					temp = ((point[j].x * point[j].x) + (point[j].z * point[j].z));
					pjx = sqrt(temp);*/
					area += point[i].xz * point[j].y;
					area -= point[j].xz * point[i].y;
					/*printf("\ni=%i, j=%i Current AREA = %lf",i,j,area);*/
				}
		
				area *= 0.5;
			}
			aTotal += (area==0.0) ? 0.0:(area/(alpha0[num]));
			alpha[num] = (area==0.0) ? 0.0:(area/(alpha0[num]));
			if(alpha[num]>1.0){printf("\nError, element %i area ratio > 1.0", num);}
			/*printf("\nAREA[%i] = %f \t Atotal = %f", num, area, aTotal);*/
		}
	}
	/*printf("POLYGON func called: AREA = %lf  N = %i\n\n", area, N);*/
	return(fabs(aTotal)); /*return absolute value just in case*/
}

/*Function to determine if an auxillary node is within a fixed boundary*/
void AuxBound2(int bsn, Aux *Afix, int NumNodes, int Ntot, Aux *auxNode, int NIOtotal, int *TotFix2, int *FixDofs2)
{
	int i,j,k,temp,dtemp,num;
	int Asize = 2 * NIOtotal; /*size of auxNode array*/
	int max = NumNodes; /*maximum node number for original mesh*/
	double xa,ya; /*co-ordinate variables*/
	int count = *TotFix2; /*initalized fixed dof count to original value*/
	
	for(j=0;j<Asize;j++) /*for all entries in the auxNode array*/
	{
		num = auxNode[j].n; /*read in auxillary node number*/
		
		if(num > max)
		{
			xa = auxNode[j].x;
			ya = auxNode[j].y; /*read in auxillary node co-ordinates*/
			
			/*printf("\txa=%lf, ya=%lf",xa,ya);*/
			
			for(i=0;i<bsn;i++) /*for all fixed areas*/
			{
				temp = 2 * i; /*to access Afix array properly*/
				if( ((xa - Afix[temp+1].x) < 0.000001) && ((xa - Afix[temp].x) > -0.000001) ) /*If within x bounds*/
				{
					if( ((ya - Afix[temp+1].y) < 0.000001) && ((ya - Afix[temp].y) > -0.000001) ) /*If also within y bounds*/
					{
						for(k=0;k<2;k++) /*for each direction (x then y)*/
						{
							if(Afix[temp + k].n == 1) /*If direction is to be fixed*/
							{
								dtemp = (num - 1) * 2;
								dtemp += k + 1; /*this is the dof to be fixed*/
								FixDofs2[count++] = dtemp; /*add to the new array*/
								/*printf("\n%i aux dof fixed",dtemp);*/
							}
						}
						/*break; /*just in case auxillary node is in two areas!*/
					}
				}
			}
			max = num; /*update maximum auxnode number looked at so far*/
		}
	}

	FixDofs2[count] = (100 * Ntot); /*for IFG_Solve*/
	*TotFix2 = ++count; /*update count of fixed dofs*/
}



/*Function to add the auxillary nodes around the boundary*/
void auxAdd(double h, int NIOn, int Xi, int Yj, int elemX, int elemY, short *ElemStat, Elem Number[elemX][elemY], short *NodeStat,
			double *lsf, int *Acount, Aux *auxNode, int NumNodes, Coord NodeCoord[NumNodes], int *NIOnums, double hxz[elemX],double hz[elemX])
{
	/*printf("\nrecieved Acount = %i",*Acount);*/
	int i,k,eNeigh,temp,sumB;
	int n1,n2; /*node status variables*/
	double Xa,Za,XZa,node1,node2,ftemp;
	int *nodes;
	nodes = malloc(4 * sizeof(int));
	int I2 = NIOn * 2; /*Constant to access auxNode properly*/
	
	/*Read in node numbers*/
	nodes[0] = Number[Xi][Yj].a - 1;
	nodes[1] = Number[Xi][Yj].b - 1;
	nodes[2] = Number[Xi][Yj].c - 1;
	nodes[3] = Number[Xi][Yj].d - 1;
	
	int AN = 0; /*used to count number of aux nodes for the element*/
	
	/*Look at edge 1*/
	n1 = NodeStat[nodes[0]];
	n2 = NodeStat[nodes[1]];	
	if((n1 + n2) == 1) /*Edge is cut if 1 node IN & 1 OUT */
	{		
		if((Yj == 0) || (ElemStat[Number[Xi][Yj-1].n-1] == 0)) /*If element is on the bottom or element below is OUT, create new Aux node*/
		{
			temp = *Acount;
			*Acount = ++temp; /*update auxillary node count*/
			/*printf("\nAcount_1=%i",*Acount);*/
			auxNode[I2+AN].n = *Acount;		/*Create a new auxillary node*/
			auxNode[I2+AN].y = NodeCoord[nodes[0]].y;	/*Y coord is same as grid y*/
			node1 = lsf[nodes[0]];
			node2 = lsf[nodes[1]]; /*read in lsf vlaues for the two nodes*/
			Xa = 0.5 * h * (node1 + node2);
			Xa /= (node1 - node2);	/*Local co-ordinate*/
			Za = 0.5 * hz[Xi] * (node1 + node2);
			Za /= (node1 - node2);	/*Local co-ordinate*/
			XZa = 0.5 * hxz[Xi] * (node1 + node2);
			XZa /= (node1 - node2);	/*Local co-ordinate*/
			/*printf("\n1 XZA = %f, node1: %f, node2: %f, Hxz: %f",XZa, node1, node2, hxz[Xi]);*/
			if(fabs(Xa) >= (0.5 * h))
			{
				printf("\nERROR!, local Aux node co-ord in element %i = %lf", Number[Xi][Yj].n, Xa);
			}
			Xa += NodeCoord[nodes[0]].x + (0.5 * h);	/*Global co-ordinate*/
			auxNode[I2+AN].x = Xa;
			Za += NodeCoord[nodes[0]].z + (0.5 * hz[Xi]);	/*Global co-ordinate*/
			auxNode[I2+AN].z = Za;
			XZa += NodeCoord[nodes[0]].xz + (0.5 * hxz[Xi]);	/*Global co-ordinate*/
			auxNode[I2+AN].xz = XZa;
			AN++;
		}
		
		else /*Otherwise copy the aux node from below*/
		{
			eNeigh = Number[Xi][Yj-1].n; /*neigbouring element number*/
			for(i=0;i<NIOn;i++)
			{
				if(NIOnums[i] == eNeigh) /*If auxillary nodes already calculated for neigbouing element*/
				{
					i *= 2;
					k = ( fabs(auxNode[i].y - NodeCoord[nodes[0]].y) < 0.000001 ) ? 0 : 1;	/*work out if aux node is 1st or 2nd of one below*/
					/*Copy node data from aux node of element below*/
					auxNode[I2+AN].n = auxNode[i+k].n;
					auxNode[I2+AN].x = auxNode[i+k].x;
					auxNode[I2+AN].y = auxNode[i+k].y;
					auxNode[I2+AN].z = auxNode[i+k].z;
					auxNode[I2+AN].xz = auxNode[i+k].xz;
					AN++;
					break;
				}
			}
		}

	}
	
	/*Look at edge 2*/
	n1 = NodeStat[nodes[1]];
	n2 = NodeStat[nodes[2]];	
	if((n1 + n2) == 1) /*Edge is cut if 1 node IN & 1 OUT */
	{
		temp = *Acount;
		*Acount = ++temp; /*update auxillary node count*/
		/*printf("\nAcount_2=%i",*Acount);*/
		auxNode[I2+AN].n = *Acount;		/*Create a new auxillary node*/
		auxNode[I2+AN].x = NodeCoord[nodes[1]].x;	/*X coord is same as grid x*/
		auxNode[I2+AN].z = NodeCoord[nodes[1]].z;	/*Z coord is same as grid z*/
		auxNode[I2+AN].xz = NodeCoord[nodes[1]].xz;	/*XZ coord is same as grid xz*/
		node1 = lsf[nodes[1]];
		node2 = lsf[nodes[2]]; /*read in lsf vlaues for the two nodes*/
		Xa = 0.5 * h * (node1 + node2);
		Xa /= (node1 - node2);	/*Local co-ordinate*/
		if(fabs(Xa) >= (0.5 * h))
		{
			printf("\nERROR!, local Aux node co-ord in element %i = %lf", Number[Xi][Yj].n, Xa);
		}
		Xa += NodeCoord[nodes[1]].y + (0.5 * h);	/*Global co-ordinate*/
		auxNode[I2+AN].y = Xa;
		AN++;
	}
	
	/*Look at edge 3*/
	n1 = NodeStat[nodes[2]];
	n2 = NodeStat[nodes[3]];	
	if((n1 + n2) == 1) /*Edge is cut if 1 node IN & 1 OUT */
	{
		temp = *Acount;
		*Acount = ++temp; /*update auxillary node count*/
		/*printf("\nAcount_3=%i",*Acount);*/
		auxNode[I2+AN].n = *Acount;		/*Create a new auxillary node*/
		auxNode[I2+AN].y = NodeCoord[nodes[2]].y;	/*Y coord is same as grid y*/
		node1 = lsf[nodes[2]];
		node2 = lsf[nodes[3]]; /*read in lsf vlaues for the two nodes*/
		Xa = 0.5 * h * (node1 + node2);
		Xa /= (node2 - node1);	/*Local co-ordinate*/
		Za = 0.5 * hz[Xi] * (node1 + node2);
		Za /= (node2 - node1);	/*Local co-ordinate*/
		XZa = 0.5 * hxz[Xi] * (node1 + node2);
		XZa /= (node2 - node1);	/*Local co-ordinate*/
		/*printf("\n2 XZA = %f, node1: %f, node2: %f, Hxz: %f",XZa, node1, node2, hxz[Xi]);
		printf("\n n2+n1 = %f, n2-n1 = %f", node1 + node2, node2 - node1);*/
		if(fabs(Xa) >= (0.5 * h))
		{
			printf("\nERROR!, local Aux node co-ord in element %i = %lf", Number[Xi][Yj].n, Xa);
		}
		Xa += NodeCoord[nodes[0]].x + (0.5 * h);	/*Global co-ordinate*/
		auxNode[I2+AN].x = Xa;
		Za += NodeCoord[nodes[0]].z + (0.5 * hz[Xi]);	/*Global co-ordinate*/
		auxNode[I2+AN].z = Za;
		XZa += NodeCoord[nodes[0]].xz + (0.5 * hxz[Xi]);	/*Global co-ordinate*/
		auxNode[I2+AN].xz = XZa;
		AN++;
	}
	
	/*Look at edge 4*/
	n1 = NodeStat[nodes[0]];
	n2 = NodeStat[nodes[3]];	
	if((n1 + n2) == 1) /*Edge is cut if 1 node IN & 1 OUT */
	{
		if((Xi == 0) || (ElemStat[Number[Xi-1][Yj].n-1] == 0)) /*If element is on the left or element to left is OUT, create new Aux node*/
		{
			temp = *Acount;
			*Acount = ++temp; /*update auxillary node count*/
			/*printf("\nAcount_4=%i",*Acount);*/
			auxNode[I2+AN].n = *Acount;		/*Create a new auxillary node*/
			auxNode[I2+AN].x = NodeCoord[nodes[0]].x;	/*X coord is same as grid x*/
			auxNode[I2+AN].z = NodeCoord[nodes[0]].z;	/*Z coord is same as grid z*/
			auxNode[I2+AN].xz = NodeCoord[nodes[0]].xz;	/*XZ coord is same as grid xz*/
			node1 = lsf[nodes[0]];
			node2 = lsf[nodes[3]]; /*read in lsf vlaues for the two nodes*/
			Xa = 0.5 * h * (node1 + node2);
			Xa /= (node1 - node2);	/*Local co-ordinate*/
			if(fabs(Xa) >= (0.5 * h))
			{
				printf("\nERROR!, local Aux node co-ord in element %i = %lf", Number[Xi][Yj].n, Xa);
			}
			Xa += NodeCoord[nodes[0]].y + (0.5 * h);	/*Global co-ordinate*/
			auxNode[I2+AN].y = Xa;
			AN++;
		}
		
		else /*Otherwise copy the aux node from the left*/
		{
			eNeigh = Number[Xi-1][Yj].n; /*neigbouring element number*/
			for(i=0;i<NIOn;i++)
			{
				if(NIOnums[i] == eNeigh) /*If auxillary nodes already calculated for neigbouing element*/
				{
					i *= 2;
					k = ( fabs(auxNode[i].x - NodeCoord[nodes[0]].x) < 0.000001 ) ? 0 : 1;	/*work out if aux node is 1st or 2nd of one left*/
					/*Copy node data from aux node of element to left*/
					auxNode[I2+AN].n = auxNode[i+k].n;
					auxNode[I2+AN].x = auxNode[i+k].x;
					auxNode[I2+AN].y = auxNode[i+k].y;
					auxNode[I2+AN].z = auxNode[i+k].z;
					auxNode[I2+AN].xz = auxNode[i+k].xz;
					AN++;
					break;
				}
			}
		}
	}
	
	if(AN != 2)
	{
		if(AN > 2)
		{
			printf("\nERROR! In AuxAdd. Element %i has %i aux nodes", Number[Xi][Yj].n, AN);
		}
		
		else
		{
			for(i=AN;i<2;i++) /*if less than required aux nodes were created then, set some dummy values*/
			{
				auxNode[I2+i].n = -1;
				auxNode[I2+i].x = -1.0;
				auxNode[I2+i].y = -1.0;
				auxNode[I2+i].z = -1.0;
				auxNode[I2+i].xz = -1.0;
			}
		}
	}
	
	free(nodes);
	
	/*printf("\nAcount Now= %i",*Acount);*/
	/*for(i=0;i<AN;i++)
	{
		printf("\nauxNode[%i][%i].n = %i, x=%lf, y=%lf",NIOn,i+1,auxNode[I2+i].n,auxNode[I2+i].x,auxNode[I2+i].y);
	}*/
}

/*Function to count the number of various element types*/
void elemType(int NumElem, short *ElemStat, int *Itotal, int *Ototal, int *NIOtotal, int *tCount, int *qCount, int *pCount)
{
	int i, temp; /*incermentors*/
	
	/*first re-set all values to zero*/
	*tCount = 0;
	*qCount = 0;
	*pCount = 0;
	*Itotal = 0;
	*NIOtotal = 0;
	*Ototal = 0;
	
	/*Update values*/
	for(i=0;i<NumElem;i++)
	{
		if(ElemStat[i] == 4) {
			temp = *Itotal;
			*Itotal = ++temp;
			}
		else if(ElemStat[i] == 0) {
			temp = *Ototal;
			*Ototal = ++temp;
			}
		else {
			temp = *NIOtotal;
			*NIOtotal = ++temp;
			if(ElemStat[i] == 1) {
				temp = *tCount;
				*tCount = ++temp;
				}
			else if(ElemStat[i] == 2) {
				temp = *qCount;
				*qCount = ++temp;
				}
			else if(ElemStat[i] == 3) {
				temp = *pCount;
				*pCount = ++temp;
				}
			else {
				printf("\nERROR! Element Status %i not 0->4",i+1);
				}
			}
	}
	
	/*store data for current iteration
	Etypes[t] = *Itotal;
	Etypes[t + 1] = *pCount;
	Etypes[t + 2] = *qCount;
	Etypes[t + 3] = *tCount;
	Etypes[t + 4] = *Ototal;*/
	
	printf("\nI Elems = %i, NIO Elems = %i, O Elems = %i",*Itotal,*NIOtotal,*Ototal);
	printf("\nNo. Tri = %i No. Quad = %i No. of Pen = %i",*tCount,*qCount,*pCount);
}

/*Function to read IN and auxillary node numbers for an NIO element*/
void readNodes(int nNod, int NIOnum, int i, int j, int elemX, int elemY, short *NodeStat, Elem Number[elemX][elemY], int *tnodes, Aux *auxNode)
{

	int node1,node2,node3,n,m,s;
	int *nodes;
	nodes = malloc(4 * sizeof(int));
	int N2 = NIOnum * 2; /*Constant to access auxNode properly*/
	/*printf("\nRn:nNod=%i,NIOnum=%i",nNod,NIOnum);

	/*Read in node numbers*/
	nodes[0] = Number[i][j].a - 1;
	nodes[1] = Number[i][j].b - 1;
	nodes[2] = Number[i][j].c - 1;
	nodes[3] = Number[i][j].d - 1;
	
	n = 0; /*used to count number of nodes found*/
	m = 0; /*used to count number of aux nodes for the element*/
	
	for(s=0;s<4;s++)
	{
		node2 = (s == 3) ? 0 : (s+1); /*Next node round*/
		node3 = (s == 0) ? 3 : (s-1); /*previous node*/
		/*Look at current node*/
		if(NodeStat[nodes[s]] == 1)
		{
			/*printf("\nNode %i is IN",s+1);*/
			tnodes[n++] = nodes[s];
		}
		
		else if(NodeStat[nodes[s]] == 2) /*If node is on boundary*/
		{
			/*Check to see if boundary node attahced to some part of the current element*/
			if((NodeStat[nodes[node2]] != 0) || (NodeStat[nodes[node3]] != 0))
			{
				tnodes[n++] = nodes[s];
			}
		}
		
		/*Look at next edge*/
		if((NodeStat[nodes[s]] + NodeStat[nodes[node2]]) == 1)
		{
			/*printf("\nEdge %i is cut",s+1);*/
			tnodes[n++] = auxNode[N2+m].n -1;
			m++;
		}
	}
	
	if(n != nNod)
	{
		printf("\nERROR! in read nodes. Element %i, node read=%i, nNod=%i", Number[i][j].n, n, nNod);
	}
	
	/*Fix for some quadrilaterals*/
	if(nNod == 4)
	{
		/*If node 1 is OUT and node 4 is IN, then re-order nodes properly*/
		if((NodeStat[nodes[0]] == 0) && (NodeStat[nodes[3]] == 1))
		{
			node1 = tnodes[3];
			tnodes[3] = tnodes[2];
			tnodes[2] = tnodes[1];
			tnodes[1] = tnodes[0];
			tnodes[0] = node1;
		}
	}
	
	/*Fix for some triangles*/
	else if(nNod == 3)
	{
		/*If node 1 is OUT and node 4 is IN, then re-order nodes properly*/
		if((NodeStat[nodes[3]] == 1) && (NodeStat[nodes[0]] == 0))
		{
			node1 = tnodes[2];
			tnodes[2] = tnodes[1];
			tnodes[1] = tnodes[0];
			tnodes[0] = node1;
		}
	}
	
	free(nodes);
	
	for(s=0;s<nNod;s++)
	{
		if(tnodes[s] < 0)
		{
			printf("\nERROR! IN read nodes, node num %i for element %i = %i", s+1, Number[i][j].n, tnodes[s]);
		}
	}
}

/*Function to read IN and auxillary node co-ordinates for an NIO element*/
void readCoord(int nNod, int N2, int i, int j, int elemX, int elemY, short *NodeStat, Elem Number[elemX][elemY], Aux *auxNode, 
				Coord *tcords, int NumNodes, Coord NodeCoord[NumNodes], int *tnodes)
{
	int node2,s,k;
	int *nodes;
	nodes = malloc(4 * sizeof(int));
	/*printf("\nRn:nNod=%i,NIOnum=%i",nNod,NIOnum);*/
	
	/*Read in node numbers*/
	nodes[0] = Number[i][j].a - 1;
	nodes[1] = Number[i][j].b - 1;
	nodes[2] = Number[i][j].c - 1;
	nodes[3] = Number[i][j].d - 1;
	
	for(s=0;s<nNod;s++)
	{
		/*If node isn't an auxillary node then read co-ordinates straight in*/
		if(tnodes[s] < 0)
		{
			printf("\nERROR! In read coord. auxillary node doesn't exist for Element %i", Number[i][j].n);
		}
		
		else if(tnodes[s] < NumNodes)
		{
			/*printf("\nNode %i is IN",s+1);*/
			tcords[s].x = NodeCoord[tnodes[s]].x;
			tcords[s].y = NodeCoord[tnodes[s]].y;
			tcords[s].z = NodeCoord[tnodes[s]].z;
			tcords[s].xz = NodeCoord[tnodes[s]].xz;
		}
		/*Otherwise find the co-ordinates of the auxillary node*/
		else
		{
			for(k=0;k<2;k++)
			{
				if((auxNode[N2+k].n -1) == tnodes[s])
				{
					tcords[s].x = auxNode[N2+k].x;
					tcords[s].y = auxNode[N2+k].y;
					tcords[s].z = auxNode[N2+k].z;
					tcords[s].xz = auxNode[N2+k].xz;
					break;
				}
			}
			
		}
		
	}
	
	free(nodes);
}

	
/*Edge calculation function for auxillary node formulation*/	
void EdgeCalc3(int Enum, int NIOnum, short *NodeStat, int NumNodes, int i, int j, int elemX, int elemY, Elem Number[elemX][elemY], 
				double h, Coord NodeCoord[NumNodes], Aux *auxNode, double *Edges)
{
	int *nodes;
	nodes = malloc(4 * sizeof(int));
	int n = 0; /*used to count element aux nodes*/
	int N2 = NIOnum * 2; /*constant to access auxNode properly*/
	
	/*read grid element node numbers*/
	nodes[0] = Number[i][j].a-1;
	nodes[1] = Number[i][j].b-1;
	nodes[2] = Number[i][j].c-1;
	nodes[3] = Number[i][j].d-1;
	
	/*look at edge 1*/	
	if((NodeStat[nodes[0]] * NodeStat[nodes[1]]) == 0) /*If one node is OUT */
	{
		if((NodeStat[nodes[0]] + NodeStat[nodes[1]]) != 1) /*If edge not cut then edge length = 0.0 */
		{
			Edges[0] = 0.0;
		}
		else /*otherwise work out edge length from aux node and grid node co-ordinates*/
		{
			Edges[0] = (NodeStat[nodes[0]] == 0) ? (NodeCoord[nodes[1]].x - auxNode[N2+n].x) : (auxNode[N2+n].x - NodeCoord[nodes[0]].x);
			n++; /*update aux node count*/
		}
	}
	else /*If neither node OUT, then entire edge is IN*/
	{
		Edges[0] = h;
	}
	
	
	/*look at edge 2*/
	if((NodeStat[nodes[1]] * NodeStat[nodes[2]]) == 0) /*If one node is OUT */
	{
		if((NodeStat[nodes[1]] + NodeStat[nodes[2]]) != 1) /*If edge not cut then edge length = 0.0 */
		{
			Edges[1] = 0.0;
		}
		else /*otherwise work out edge length from aux node and grid node co-ordinates*/
		{
			Edges[1] = (NodeStat[nodes[1]] == 0) ? (NodeCoord[nodes[2]].y - auxNode[N2+n].y) : (auxNode[N2+n].y - NodeCoord[nodes[1]].y);
			n++; /*update aux node count*/
		}
	}
	else /*If neither node OUT, then entire edge is IN*/
	{
		Edges[1] = h;
	}
	
	/*look at edge 3*/
	if((NodeStat[nodes[2]] * NodeStat[nodes[3]]) == 0) /*If one node is OUT */
	{
		if((NodeStat[nodes[2]] + NodeStat[nodes[3]]) != 1) /*If edge not cut then edge length = 0.0 */
		{
			Edges[2] = 0.0;
		}
		else /*otherwise work out edge length from aux node and grid node co-ordinates*/
		{
			Edges[2] = (NodeStat[nodes[2]] == 0) ? (auxNode[N2+n].x - NodeCoord[nodes[3]].x) : (NodeCoord[nodes[2]].x - auxNode[N2+n].x);
			n++; /*update aux node count*/
		}
	}
	else /*If neither node OUT, then entire edge is IN*/
	{
		Edges[2] = h;
	}
	
	/*look at edge 4*/
	if((NodeStat[nodes[3]] * NodeStat[nodes[0]]) == 0) /*If one node is OUT */
	{
		if((NodeStat[nodes[3]] + NodeStat[nodes[0]]) != 1) /*If edge not cut then edge length = 0.0 */
		{
			Edges[3] = 0.0;
		}
		else
		{
			Edges[3] = (NodeStat[nodes[3]] == 0) ? (auxNode[N2+n].y - NodeCoord[nodes[0]].y) : (NodeCoord[nodes[3]].y - auxNode[N2+n].y);
		}
	}
	else /*If neither node OUT, then entire edge is IN*/
	{
		Edges[3] = h;
	}
	
	free(nodes);
}
	
/*Function to calculate geometric and center co-ordinates for a pentagonal element*/
void PentCoord(int Enum, double h, double a1, double a2, double b1, double b2, short *NodeStat, int i, int j, 
				int elemX, int elemY, Elem Number[elemX][elemY], int pnum, Coord *Pcrds)
{
	int n; /*Incermentor*/
	int N4 = Enum * 4;
	int P6 = pnum * 6; /*to access Pcrds correctly*/
	double h2 = 0.5 * h;
	double *xCrd;
	double *yCrd; /*Arrays for nodes co-ords realtive to IN element center*/
	xCrd = malloc(5 * sizeof(double));
	yCrd = malloc(5 * sizeof(double));
	
	/*Indetify the shorter edge lengths*/
	double a = (a1 < a2) ? a1 : a2;
	double b = (b1 < b2) ? b1 : b2;	
	/*calculate element area*/
	double Area = 2 * h2 * h2;
	Area += (a + b) * h2;
	Area -= 0.5 * a * b;
	
	/*Define some co-ordinates at the start*/
	xCrd[0] = -h2;
	xCrd[1] = h2;
	xCrd[2] = h2;
	xCrd[4] = -h2;
	
	yCrd[0] = -h2;
	yCrd[1] = -h2;
	yCrd[2] = h2;
	yCrd[3] = h2;
	yCrd[4] = h2;
	
	/*Now work out missing co-ords and center co-ords*/
	if(NodeStat[Number[i][j].a - 1] == 0)
	{
		xCrd[0] = h2 - a1;
		xCrd[3] = -h2;
		yCrd[4] = h2 - b1;
		Pcrds[P6].x = Cent(a,b,h2,Area,1,-1);
		Pcrds[P6].y = Cent(b,a,h2,Area,1,-1); /*First entry is center co-ordinate*/
	}
	
	else if(NodeStat[Number[i][j].b - 1] == 0)
	{
		xCrd[1] = a1 - h2;
		yCrd[2] = h2 - b2;
		xCrd[3] = h2;
		Pcrds[P6].x = Cent(a,b,h2,Area,-1,1);
		Pcrds[P6].y = Cent(b,a,h2,Area,1,-1);
	}
	
	else if(NodeStat[Number[i][j].c - 1] == 0)
	{
		yCrd[2] = b2 - h2;
		xCrd[3] = a2 - h2;
		Pcrds[P6].x = Cent(a,b,h2,Area,-1,1);
		Pcrds[P6].y = Cent(b,a,h2,Area,-1,1);
	}
	
	else if(NodeStat[Number[i][j].d - 1] == 0)
	{
		xCrd[3] = h2 - a2;
		yCrd[4] = b1 - h2;
		Pcrds[P6].x = Cent(a,b,h2,Area,1,-1);
		Pcrds[P6].y = Cent(b,a,h2,Area,-1,1);
	}
	
	else
	{
		printf("\nERROR! could not find OUT node for pentagonal Element %i",Enum+1);
	}
	
	/*Calculate node co-ords relative to pentagonal element geometric center*/
	for(n=1;n<6;n++)
	{
		Pcrds[P6+n].x = xCrd[n-1] - Pcrds[P6].x;
		Pcrds[P6+n].y = yCrd[n-1] - Pcrds[P6].y;
	}
	
	free(xCrd);
	free(yCrd);
	
	/*printf("\nPent X coords:\n");
	for(n=0;n<6;n++)
	{
		printf("\t%lf",Pcrds[P6+n].x);
	}
	printf("\nPent Y coords:\n");
	for(n=0;n<6;n++)
	{
		printf("\t%lf",Pcrds[P6+n].y);
	}*/
}

/*Function for working out center co-ordinates of a pentagon relative to IN element center*/
double Cent(double a, double b, double h, double A, short p, short q)
{
	/*printf("\nCent called");*/
	double C; /*variable for center co-ordinate*/
	double D; /*Other variable*/
	
	D = 2.0 * h * h* b;
	D += a * h * b;
	D += (a -h) * 2.0 * h * a;
	D *= (double)q;
	
	C = 4.0 * h * h * h;
	C += a * a * b;
	C *= (double)p;
	
	C += D;
	C /= 6.0 * A;
	
	/*printf("\nC=%lf",C);*/
	return(C);
}

/*Function to find the nearest grid node number to a set of co-ordinates*/
int closeNode(double h, double xp, double yp, int NumNodes, Coord NodeCoord[NumNodes])
{
	/*printf("\nclose Node called");*/
	int n,m; /*incrementors*/
	int nodenum;
	double dtemp, d2, dx, dy;	/*variables to keep track of distnace closest node*/
	double xn, yn;				/*Variables for current node co-ordinates*/
	/*double lim = 0.00001;		/*small limit for rounding errors*/
	dtemp =  (2.0 * h);
	dtemp *= dtemp;			/*set initial distance to twice an elements edge length squared*/
	
	for(n=0;n<NumNodes;n++) /*For all nodes*/
	{
		xn = NodeCoord[n].x;
		yn = NodeCoord[n].y; /*read in node co-ords*/
		/*printf("\nxn,yn,zn = %lf, %lf, %lf",xn,yn,zn);*/
		dx = xp - xn;
		dy = yp - yn; /*find difference between node and target co-ords*/
		
		d2 = (dx * dx) + (dy * dy);	/*Squared distance between node and point*/
		/*printf("\n d2 = %lf, num = %i",d2,n);*/
		if(fabs(d2) < 0.000001)	/*If co-ordinates are same as current node, that that is closest*/
		{
			return(n);
		}
		
		if(d2 < dtemp)	/*If distance of current node is less than previous shortest distance then update*/
		{
			nodenum = n;
			dtemp = d2;
		}
	}
	
	/*printf("\n end of close node, node = %i", nodenum+1);*/
	return(nodenum);
}
