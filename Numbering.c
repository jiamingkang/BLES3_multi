/*
 *  Numbering.c
 *  
 *  Created by Peter Dunning on 23/02/2007.
 *
 */

#include "Numbering.h"
#include "ls_types.h"
#include <stdio.h>
#include <math.h>

/*--------Function that numbers all elements and nodes in the FG domain--------*/
void Numbering(int elemX,int elemY,Elem Number[elemX][elemY])
{

int i,j;
int num = 0; /*Incrementors*/

/* Element Numbering loop */
for(j=0;j<elemY;j++)
{
	for(i=0;i<elemX;i++)
	{
		Number[i][j].n = ++num;
	}
}

/*printf("\nElements Numbered\n");*/

/* State first element node numbers */
num = 0;

Number[0][0].a = ++num;
Number[0][0].b = ++num;
Number[0][0].c = ++num;
Number[0][0].d = ++num;

/*printf("Nodes of Element 1 Numbered\n");*/

/*Number first X row (Y = 0)*/

if(elemX > 1) /*if there is more than one element in the x direction number first row*/
{
	for(i=1;i<elemX;i++)
		{
			Number[i][0].a = Number[i-1][0].b;
			Number[i][0].b = ++num;
			Number[i][0].c = ++num;
			Number[i][0].d = Number[i-1][0].c;
		}
	/*printf("First Row Numbered\n");*/
}

/* Number rest of FG mesh */
if(elemY > 1)
{
	for(j=1;j<elemY;j++)
	{
		for(i=0;i<elemX;i++)
		{
			Number[i][j].a = Number[i][j-1].d;
			Number[i][j].b = Number[i][j-1].c;
			
			if(i==0)
			{
				Number[i][j].c = ++num;
				Number[i][j].d = ++num;
			}
			else
			{
				Number[i][j].c = ++num;
				Number[i][j].d = Number[i-1][j].c;
			}	
		}
	}
}

/*printf("Rest of Mesh Numbered\n");*/
}


/*-----------------Node Co-ordinate calculation function------------------*/
void Coordinates(int elemX, int elemY, int NumNodes, double hx, double hy, Elem Number[elemX][elemY], Coord NodeCoord[NumNodes])
{
	int i,j; /*incrementors*/
	int num;
	
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].a - 1;
				if((NodeCoord[num].x + NodeCoord[num].y) < 0.000001){   /*If node co-ordinate hasn't been calculated already*/
				NodeCoord[num].x = i * hx;
				NodeCoord[num].y = j * hy;
				}
				num = Number[i][j].b - 1;
				if((NodeCoord[num].x + NodeCoord[num].y) < 0.000001){
				NodeCoord[num].x = (i * hx) + hx;
				NodeCoord[num].y = j * hy;
				}
				num = Number[i][j].c - 1;
				/*printf("\nnum=%i", num);
				printf("\ntotal=%lf x=%lf y=%lf z=%lf",(NodeCoord[num].x + NodeCoord[num].y + NodeCoord[num].z),NodeCoord[num].x,NodeCoord[num].y,NodeCoord[num].z);*/
				if((NodeCoord[num].x + NodeCoord[num].y) < 0.000001){
				NodeCoord[num].x = (i * hx) + hx;
				NodeCoord[num].y = (j * hy) + hy;
				}
				num = Number[i][j].d - 1;
				if((NodeCoord[num].x + NodeCoord[num].y) < 0.000001){
				NodeCoord[num].x = i * hx;
				NodeCoord[num].y = (j * hy) + hy;
				}
		
			}
		}
}

/*----------------Function that orders node numbers into a 2D based on their relative positions--------------*/
void NodeNums2(double h, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int NumNodes, Coord NodeCoord[NumNodes])
{
	int i; /*incrementor*/
	int X,Y; /*node position variables*/
	
	/*for all grid node determine position in ordered 2D array by their co-ordinates and element edge length*/
	for(i=0;i<NumNodes;i++)
	{
		X = floor((NodeCoord[i].x / h)+0.0001) + 1;
		Y = floor((NodeCoord[i].y / h)+0.0001) + 1; /*+1 for layer of ghost nodes*/
		
		Nodes2[X][Y] = i;
	}
	
	/*now copy node numbers to create a layer of ghost nodes around the grid*/
	
	/*fill in the four corners*/
	Nodes2[0][0] = Nodes2[1][1];
	Nodes2[0][NodeY-1] = Nodes2[1][NodeY-2];
	Nodes2[NodeX-1][0] = Nodes2[NodeX-2][1];
	Nodes2[NodeX-1][NodeY-1] = Nodes2[NodeX-2][NodeY-2];
	
	/*now fill in top and bottom rows*/
	for(i=1;i<NodeX-1;i++)
	{
		Nodes2[i][0] = Nodes2[i][1];
		Nodes2[i][NodeY-1] = Nodes2[i][NodeY-2];
	}
	
	/*Finally fill in left and right columns*/
	for(i=1;i<NodeY-1;i++)
	{
		Nodes2[0][i] = Nodes2[1][i];
		Nodes2[NodeX-1][i] = Nodes2[NodeX-2][i];
	}
	
}
