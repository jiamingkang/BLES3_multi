/*
 *  EMatrix.c
 *  
 *  Created by Peter Dunning on 20/02/2007.
 *
 */

#include "EMatrix.h"
#include "ls_types.h"
#include "Solve.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*Array storing data for location of the 4 nodes in the square element*/
static Pos Q4[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};

/*look up tables for quadrilateral BFG element stiffness matrix calculation*/
static Three alOdd[4][4] = {{{1,3,-1},{0,1,2},{1,2,0},{0,2,1}},
					{{1,3,-1},{0,1,2},{1,2,0},{0,2,1}},
					{{1,2,0},{0,2,1},{2,1,0},{-1,3,1}},
					{{1,2,0},{0,2,1},{2,1,0},{-1,3,1}}};

static Three beOdd[4][4] = {{{1,3,-1},{1,3,-1},{1,2,0},{1,2,0}},
					{{0,1,2},{0,1,2},{0,2,1},{0,2,1}},
					{{1,2,0},{1,2,0},{2,1,0},{2,1,0}},
					{{0,2,1},{0,2,1},{-1,3,1},{-1,3,1}}};
				
static Four beEven[4][4] = {{{1,3,-1,1},{0,1,2,-1},{1,2,-1,0},{0,2,2,0}},
					{{0,1,2,-1},{0,0,2,2},{0,2,2,0},{0,-1,2,1}},
					{{1,2,-1,0},{0,2,2,0},{2,2,0,0},{-1,2,1,0}},
					{{0,2,2,0},{0,-1,2,1},{-1,2,1,0},{1,-1,3,1}}};
					
/*row and column indices for a 4x4 symmetric matric*/
static int jcnp[10] = {1,2,3,4,2,3,4,3,4,4};
static int irnp[10] = {1,1,1,1,2,2,2,3,3,4};

/*The following 2 functions calculate the values in each 2 x 2 block in the Element Stiffness Matrix 
	These functions integrate the matrix exactly. Trust me it works
	--------------------------------------------------------------------------------------------------*/

/*For odd (i + j = odd) Matrix enrites*/	
double Aodd(int i,int j,double e,double g)
{
	short a = (Q4[i].y * Q4[j].x);	
	short b = (Q4[i].x * Q4[j].y);
	
	double value = (0.25 * ((b * e) + (a * g)));
	
	return (value);
}

/*For even (i + j = even) Matrix enrites*/	
double Aeven(int i,int j,double e,double g)
{
	double a = 1.5 + (0.5 * (Q4[i].y * Q4[j].y));
	a *= (Q4[i].x * Q4[j].x);
	
	double b = 1.5 + (0.5 * (Q4[i].x * Q4[j].x));
	b *= (Q4[i].y * Q4[j].y);
	
	double value = (0.166666666667 * ((e * a) + (g * b)));
	printf("\n A: %f \t B: %f \t Value: %f", a, b, value);
		
	return (value);
}

/*---------------------Iso Matrix entry calculation funcions--------------*/

/*alpha value calculation*/
double Ievenv1(int i, int j, double a1, double a2, int flag)
{
	double value;
	
	if(flag == 1)
	{
		value = a1 * (double)(2 + Q4[i].y);
		value += a2 * (double)(2 - Q4[j].y);
	}
	else if(flag == 2)
	{
		if((Q4[i].x == 1) && (Q4[j].x == -1)) {
			value = a1 + a2;
		}
		
		else {
			value = a1 * (double)(2 + Q4[i].x);
			value += a2 * (double)(2 - Q4[j].x);
		}
	}
	
	return(value);
}

/*---------------------------------Time to calculate Element Stiffness Matirx------------------------------*/

void KEMatrix(double *KE, double c, double d, double g)
{
	int n,m; /*Incrementors*/

	/* Populate entire Matrix with the following*/
	int row,col,temp;

	/*For each 2 x 2 block*/
	for(n=0;n<4;n++)
	{
		for(m=0;m<4;m++)
		{
			row = 2 * n;
			col = 2 * m;
			
			temp = (8 * row) + col;
			KE[temp] = Aeven(n,m,c,g);
			KE[temp+1] = Aodd(n,m,d,g);
			KE[temp+8] = Aodd(n,m,g,d);
			KE[temp+9] = Aeven(n,m,g,c);
		}
	}
}

/*-----------------------------------Iso-Parametric Matrix--------------------------------*/
/*New function that calculates isoparametric matrices correctly*/
void IsoMatrix(double *Iso, double c, double d, double g, double a1, double a2, double b1, double b2, double h)
{	
	double bt = h * 0.5;
	double as = bt;
	
	int n,m,row,col; /*Incrementors*/
	int temp;
	short *Count; /*Array to ensure values are determined correctly*/
	Count = malloc(4 * sizeof(short));
	
	/*variables for algrebric Isoparametric calcualtion*/
	double alphaE, alphaO, betaE, betaO, denom;
	
	/*If both a1 & a2 do not equal zero*/
	if((a1 > 0.000001) && (a2 > 0.000001))
	{
		/*Setup count array for gamma value determination*/
		if( ((b1 - h) > -0.000001) && ((b1 - h) < 0.000001) )	/*if left edge is completely IN, b1=h*/
		{
			Count[0] = 0;
			Count[1] = 1;
			Count[2] = 2;
			Count[3] = 3;
		}
		
		else if( ((b2 - h) > -0.000001) && ((b2 - h) < 0.000001) )	 /*if right edge is completely IN, b2=h*/
		{
			Count[0] = 1;
			Count[1] = 0;
			Count[2] = 3;
			Count[3] = 2;
		}
		
		else
		{
			printf("\nERROR! In IsoMatrix: can't find orientation! -> both b1 & b2 != h");
		}

		denom = a1 * a1;
		denom += a2 * a2;
		denom += 4.0 * a1 * a2;
		
		for(n=0;n<4;n++)
			{
				for(m=n;m<4;m++)
				{
					/*Calculate alpha value for even entry*/
					alphaE = Ievenv1(n, m, a1, a2, 1);
					alphaE *= bt * bt;
					alphaE *= Q4[n].x * Q4[m].x;
					
					/*Calculate beta value for even entry*/
					betaE = (double)beEven[Count[n]][Count[m]].a * a1 * a1 * a1;
					betaE += (double)beEven[Count[n]][Count[m]].b * a1 * a1 * a2;
					betaE += (double)beEven[Count[n]][Count[m]].c * a1 * a2 * a2;
					betaE += (double)beEven[Count[n]][Count[m]].d * a2 * a2 * a2;
					betaE *= 0.25;
					betaE *= Q4[n].y * Q4[m].y;
					
					/*Calculate alpha value for odd entry*/
					alphaO = (double)alOdd[Count[n]][Count[m]].a * a1 * a1;
					alphaO += (double)alOdd[Count[n]][Count[m]].b * a1 * a2;
					alphaO += (double)alOdd[Count[n]][Count[m]].c * a2 * a2;
					alphaO *= Q4[n].x * Q4[m].y;
					
					/*Calculate beta value for odd entry*/
					betaO = (double)beOdd[Count[n]][Count[m]].a * a1 * a1;
					betaO += (double)beOdd[Count[n]][Count[m]].b * a1 * a2;
					betaO += (double)beOdd[Count[n]][Count[m]].c * a2 * a2;
					betaO *= Q4[n].y * Q4[m].x;
					
					/*printf("\nalphaE = %lf, betaE = %lf, alphaO = %lf, betaO = %lf\ndenom = %lf", alphaE, betaE, alphaO, betaO, denom);*/
					
					row = n * 2;
					col = m * 2;
					temp = (8 * row) + col;
					
					Iso[temp] = (c * alphaE) + (g * betaE);
					Iso[temp] /= (bt * denom);
					Iso[temp + 9] = (g * alphaE) + (c * betaE);
					Iso[temp + 9] /= (bt * denom);
					
					Iso[temp + 1] = (d * alphaO) + (g * betaO);
					Iso[temp + 1] /= 2.0 * denom;
					Iso[temp + 8] = (g * alphaO) + (d * betaO);
					Iso[temp + 8] /= 2.0 * denom;
				}
			}
		}
		
	/*If both b1 & b2 do not equal zero*/
	else if((b1 > 0.000001) && (b2 > 0.000001))
	{
		/*Setup count array for gamma value determination*/
		if( ((a1 - h) > -0.000001) && ((a1 - h) < 0.000001) )	 /*if bottom edge is completely IN, a1=h*/
		{
			Count[0] = 0;
			Count[1] = 3;
			Count[2] = 2;
			Count[3] = 1;
		}
		
		else if( ((a2 - h) > -0.000001) && ((a2 - h) < 0.000001) )	 /*if top edge is completely IN, a2=h*/
		{
			Count[0] = 1;
			Count[1] = 2;
			Count[2] = 3;
			Count[3] = 0;
		}
		
		else
		{
			printf("\nERROR! In IsoMatrix: can't find orientation! -> both a1 & a2 != h");
		}
		
		denom = b1 * b1;
		denom += b2 * b2;
		denom += 4.0 * b1 * b2;
		
		for(n=0;n<4;n++)
			{
				for(m=n;m<4;m++)
				{
					/*Calculate alpha value for even entry*/
					alphaE = Ievenv1(n, m, b1, b2, 2);
					alphaE *= as * as;
					alphaE *= Q4[n].y * Q4[m].y;
					
					/*Calculate beta value for even entry*/
					betaE = beEven[Count[n]][Count[m]].a * b1 * b1 * b1;
					betaE += beEven[Count[n]][Count[m]].b * b1 * b1 * b2;
					betaE += beEven[Count[n]][Count[m]].c * b1 * b2 * b2;
					betaE += beEven[Count[n]][Count[m]].d * b2 * b2 * b2;
					betaE *= 0.25;
					betaE *= Q4[n].x * Q4[m].x;
					
					/*Calculate alpha value for odd entry*/
					alphaO = alOdd[Count[n]][Count[m]].a * b1 * b1;
					alphaO += alOdd[Count[n]][Count[m]].b * b1 * b2;
					alphaO += alOdd[Count[n]][Count[m]].c * b2 * b2;
					alphaO *= Q4[n].y * Q4[m].x;
					
					/*Calculate beta value for odd entry*/
					betaO = beOdd[Count[n]][Count[m]].a * b1 * b1;
					betaO += beOdd[Count[n]][Count[m]].b * b1 * b2;
					betaO += beOdd[Count[n]][Count[m]].c * b2 * b2;
					betaO *= Q4[n].x * Q4[m].y;
					
					row = n * 2;
					col = m * 2;
					temp = (8 * row) + col;

					Iso[temp] = (g * alphaE) + (c * betaE);
					Iso[temp] /= (as * denom);
					Iso[temp + 9] = (c * alphaE) + (g * betaE);
					Iso[temp + 9] /= (as * denom);
					
					Iso[temp + 1] = (g * alphaO) + (d * betaO);
					Iso[temp + 1] /= 2.0 * denom;
					Iso[temp + 8] = (d * alphaO) + (g * betaO);
					Iso[temp + 8] /= 2.0 * denom;
				}
			}
		}
		
	else
	{
		printf("\nERROR! In IsoMatrix: can't find orientation!");
	}
	
	/*Now copy Upper triangle to lower triangle*/
	for(m=1;m<8;m++)
	{
		for(n=0;n<m;n++)
		{
			Iso[(8 * m)+n] = Iso[(8 * n)+m];
		}
	}
	
	free(Count);
	
	/*Print matrix to screen
	printf("\nQuad Matrix\n");
	for(n=0;n<8;n++)
	{
		for(m=0;m<8;m++)
		{
			printf("%lf\t",Iso[(n * 8) + m]);
		}
		printf("\n");
	}
	
	printf("\nEnter No to cont: ");
	scanf("%i",&row);*/
}

/*------------------Constant Strain Triangle Stiffness martrix-------------------*/
void Triangle(int Enum, short *NodeStat, int Xi, int Yj, int elemX, int elemY, Elem Number[elemX][elemY], double *Tri, 
				double a1, double b1, double e1, double v1, double g1, double h)
{
	int i,j,row,col; /*Incrementors*/
	int temp;
	short *nodes;
	nodes = malloc(4 * sizeof(short));
	/*Read in node's status*/
	nodes[0] = NodeStat[Number[Xi][Yj].a - 1];
	nodes[1] = NodeStat[Number[Xi][Yj].b - 1];
	nodes[2] = NodeStat[Number[Xi][Yj].c - 1];
	nodes[3] = NodeStat[Number[Xi][Yj].d - 1];
	
	/*for(i=0;i<4;i++)
	{
		printf("\nnodes[%i]=%i",i,nodes[i]);
	}*/
	
	/*Array for varibles that calculate stiffness matrix entries*/
	Coord *nodals;
	nodals = calloc(3, sizeof(Coord));
	
	double area = a1 * b1 * 2.0; /* actually 4 x area of the triangle*/
	/*printf("\nNodals Allocated");
	/*Material Constants - Including area factor! - NB: plane Stress*/
	double c = e1 / area;
	double g = g1 / area;
	double d = v1 / area;
	/*printf("\narea=%lf, E=%lf, V=%lf, G=%lf",area,c,d,g);*/
	
	/*Identify the In node and assign appropriate edge lengths*/
	if(nodes[0] == 1)
	{
		nodals[0].x = -a1;
		nodals[0].y = -b1;
		nodals[2].x = a1;
		nodals[1].y = b1;
	}
	
	else if(nodes[1] == 1)
	{
		nodals[1].x = -a1;
		nodals[0].y = -b1;
		nodals[2].x = a1;
		nodals[1].y = b1;
	}
	
	else if(nodes[2] == 1)
	{
		nodals[0].x = -a1;
		nodals[2].y = -b1;
		nodals[1].x = a1;
		nodals[1].y = b1;
	}
	
	else if(nodes[3] == 1)
	{
		nodals[0].x = -a1;
		nodals[2].y = -b1;
		nodals[2].x = a1;
		nodals[1].y = b1;
	}
	
	else
	{
		printf("\nERROR! Triangle func could not find an IN node for Element %i",Enum+1);
	}
	
	/*for(i=0;i<3;i++)
	{
		printf("\nx%i = %lf, y%i = %lf",i+1,nodals[i].x,i+1,nodals[i].y);
	}*/
	
	/*Calculate upper triangle of matrix*/
	for(i=0;i<3;i++)
	{
		for(j=i;j<3;j++)
		{
			row = 2 * i;
			col = 2 * j;
			temp = (6 * row) + col;
			
			Tri[temp] = nodals[i].y * nodals[j].y * c;
			Tri[temp] += nodals[i].x * nodals[j].x * g;
			Tri[temp + 1] = nodals[i].y * nodals[j].x * d;
			Tri[temp + 1] += nodals[i].x * nodals[j].y * g;	
			Tri[temp + 7] = nodals[i].x * nodals[j].x * c;
			Tri[temp + 7] += nodals[i].y * nodals[j].y * g;
			
			if(row < col) /*Only need to calculate upper triangle*/
			{
				Tri[temp + 6] = nodals[i].x * nodals[j].y * d;
				Tri[temp + 6] += nodals[i].y * nodals[j].x * g;
			}

		}
	}
	
	/*Now copy Upper triangle to lower triangle*/
	for(j=1;j<6;j++)
	{
		for(i=0;i<j;i++)
		{
			Tri[(6 * j)+i] = Tri[(6 * i)+j];
		}
	}
	
	free(nodes);
	free(nodals);
	
	/*Print matrix to screen
	printf("\nTri Matrix\n");
	for(i=0;i<6;i++)
	{
		for(j=0;j<6;j++)
		{
			printf("%lf\t",Tri[(6 * i)+j]);
		}
		printf("\n");
	}
	
	printf("\nEnter No to cont: ");
	scanf("%i",&row);*/
}

/*Function to calculate a pentagonal element matrix using the least squares method*/
void Pentagon(int pnum, Coord *Pcrds, double *Pent, double e1, double v1, double g1, Shape *Pshp)
{	
	int P6 = 6 * pnum; /*to access Pcrds correctly*/

	int i,i4,j,k,k2,row,col;
	int temp;
	double ftemp;
	
	/*double *Brhs;	/*Array to store rhs' of matrix equation*/
	double *Brhs;
	double *Entries;
	Brhs = malloc(20 * sizeof(double));
	Entries = calloc(10, sizeof(double));
	
	/*Now Use Least Squares Linear Regression to calculate Strain Field constants for element*/
	for(i=1;i<6;i++)
	{
		/*Update RHS by adding all entries in data thus..*/
		j = (i-1) * 4;
		Brhs[j] = 1.0;
		Brhs[j+1] = Pcrds[P6+i].x;
		Brhs[j+2] = Pcrds[P6+i].y;
		Brhs[j+3] = Pcrds[P6+i].x * Pcrds[P6+i].y;
		
		/*Update Entries by adding all entries in data thus..*/
		Entries[1] += Pcrds[P6+i].x;
		Entries[2] += Pcrds[P6+i].y;
		Entries[3] += Pcrds[P6+i].x * Pcrds[P6+i].y;
		Entries[4] += Pcrds[P6+i].x * Pcrds[P6+i].x;
		Entries[6] += Pcrds[P6+i].x * Pcrds[P6+i].x * Pcrds[P6+i].y;	
		Entries[7] += Pcrds[P6+i].y * Pcrds[P6+i].y;
		Entries[8] += Pcrds[P6+i].x * Pcrds[P6+i].y * Pcrds[P6+i].y;
		Entries[9] += Pcrds[P6+i].x * Pcrds[P6+i].x * Pcrds[P6+i].y * Pcrds[P6+i].y;
	}
	
	/*Complete Entries Arrray*/
	Entries[0] = 5.0;
	Entries[5] = Entries[3];

	/*printf("\n(RHS):\n");
	for(i=0;i<5;i++)
	{
		for(j=0;j<4;j++)
		{
			printf("\t%lf",Brhs[(4 * i) + j]);
		}
		printf("\n");
	}*/
	
	/*Solve Matrix Equation*/
	int MatOrd = 4;
	int NumEnt = 10;
	int Nrhs = 5;
	int pinfo = 0;
	
	/*solve(MatOrd,NumEnt,irnp,jcnp,Entries,Nrhs,Brhs,pinfo)*/

	/*printf("\n(RHS-Solved):\n");
	for(i=0;i<5;i++)
	{
		for(j=0;j<4;j++)
		{
			printf("\t%lf",Brhs[(4 * i) + j]);
		}
		printf("\n");
	}
	
	/*-----------Store Shape funciton co-efficients for later----------*/
	int p5 = pnum * 5; /*constant to store co-efficient properly*/
	for(i=0;i<5;i++)
	{
		i4 = 4 * i;
		Pshp[p5 + i].a = Brhs[i4];
		Pshp[p5 + i].b = Brhs[i4+1];
		Pshp[p5 + i].c = Brhs[i4+2];
		Pshp[p5 + i].d = Brhs[i4+3];
	}
	
	/*-------------Calculate area of each sub-domain-----------*/
	double *subArea; /*Array to store areas of subdomains*/
	subArea = malloc(5 * sizeof(double));
	
	for(i=1;i<6;i++)
	{
		j = (i==5) ? 1 : i+1;
		ftemp = Pcrds[P6+i].x * (Pcrds[P6+j].y - Pcrds[P6+i].y);
		ftemp -= Pcrds[P6+i].y * (Pcrds[P6+j].x - Pcrds[P6+i].x);
		ftemp *= 0.5;
		subArea[i-1] = fabs(ftemp);
		/*printf("\nArea %i = %lf",i,subArea[i-1]);*/
	}
	
	/*----------Calculate mid edge points & associated areas----------*/
	/*outside edges*/
	Coord *Edges; /*array for mid edge point co-ords*/
	double *iArea; /*array for intergration area weightings*/
	Edges = malloc(10 * sizeof(Coord));
	iArea = malloc(10 * sizeof(double));

	for(i=1;i<6;i++)
	{
		j = i - 1;
		k = (i == 5) ? 1 : i+1; /*loop back to start*/
		k2 = (j == 0) ? 4 : j-1; /*loop to end*/
		/*For outside mid edge point co-ords*/
		Edges[j].x = 0.5 * (Pcrds[P6+i].x + Pcrds[P6+k].x);
		Edges[j].y = 0.5 * (Pcrds[P6+i].y + Pcrds[P6+k].y);
		iArea[j] = subArea[j] * 0.33333333333; /*multiply by 1/3 now to save time later*/
		/*For inside mid edge point co-ords*/
		Edges[j+5].x = 0.5 * Pcrds[P6+i].x;
		Edges[j+5].y = 0.5 * Pcrds[P6+i].y;
		iArea[j+5] = (subArea[j] + subArea[k2]) * 0.33333333333; /*multiply by 1/3 now to save time later*/
	}
	
	/*for(i=0;i<10;i++)
	{
		printf("\nEdge%i x=%lf y=%lf area = %lf",i+1,Edges[i].x,Edges[i].y,iArea[i]);
	}*/
	
	/*-----------------------Calculate matrix------------------------*/
	double valXX,valXY,valYX,valYY; /*variables for matrix entry evaluation*/
	Coord *Bmat; /*array to store B matrix values during intergration*/
	Bmat = malloc(5 * sizeof(Coord));
	
	/*Need to re-initalize upper triangle Pent matrix to zero*/
	for(i=0;i<10;i++)
	{	
		for(j=i;j<10;j++)
		{
			Pent[(10 * i) + j] = 0.0;
		}
	}
	
	for(i=0;i<10;i++) /*for all integration points*/
	{		
		/*calculate strain displacement (B) matrix entries*/
		for(j=0;j<5;j++)
		{
			Bmat[j].x = Pshp[p5 + j].b + (Pshp[p5 + j].d * Edges[i].y);
			Bmat[j].y = Pshp[p5 + j].c + (Pshp[p5 + j].d * Edges[i].x);
		}
		
		/*evaluate stiffness matrix for current integration point - upper triangle only*/
		for(k=0;k<5;k++)
		{
			for(k2=k;k2<5;k2++)
			{
				row = 2 * k;
				col = 2 * k2;
				temp = (10 * row) + col;
				
				/*calculate entry values*/
				valXX = Bmat[k].x * Bmat[k2].x;
				valYY = Bmat[k].y * Bmat[k2].y;
				valXY = Bmat[k].x * Bmat[k2].y;
				valYX = Bmat[k].y * Bmat[k2].x;
				
				/*evaluate stiffness for integration point, weighted by area*/
				Pent[temp] += ((valXX * e1) + (valYY * g1)) * iArea[i];
				Pent[temp + 11] += ((valYY * e1) + (valXX * g1)) * iArea[i];
				Pent[temp + 1] += ((valXY * v1) + (valYX * g1)) * iArea[i];
				
				if(row < col) /*only need to evaluate upper triangle*/
				{
					Pent[temp + 10] += ((valXY * g1) + (valYX * v1)) * iArea[i];
				}
			}
		}
	}
	
	/*Now copy Upper triangle to lower triangle*/
	for(i=1;i<10;i++)
	{
		for(j=0;j<i;j++)
		{
			Pent[(10 * i)+j] = Pent[(10 * j)+i];
		}
	}
	
	free(Brhs);
	free(Entries);
	free(Bmat);
	free(Edges);
	free(iArea);
	free(subArea);
	
	/*Print matrix to screen
	printf("\nPent Matrix\n");
	for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			printf("%lf\t",Pent[(i * 10) + j]);
		}
		printf("\n");
	}
	
	printf("\nEnter No to cont: ");
	scanf("%i",&row);

	
	/*Print Pentagon Matrix to file
	outfile = fopen("Pentagon.txt", "w");
	if(outfile == NULL){
		printf("\nFailed to open Pentagon writefile\n");
		}
	else{
		for(i=0;i<10;i++)
		{
			for(j=0;j<10;j++)
			{
				fprintf(outfile,"%lf\t",Pent[i][j]);
			}
			fprintf(outfile,"\n");
		}
	}

	fclose(outfile);*/

}		
