/*
 *  AFG.c
 *  
 *
 *  Created by Peter Dunning on 05/05/2009.
 *
 */

#include "AFG.h"
#include "KMatrix.h"
#include "FixedGrid.h"
#include "Levels.h"
#include "ls_types.h"
#include "Strain.h"
#include "FibreSensitivity.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void AFG_Matrix(double e1, double e2, double v12, double v21, double g12, double g23, double g13, double *KE, double h, short *ElemStat, short *NodeStat, Aux *auxNode, int *irn, int *jcn, double *A, int elemX, int elemY, Elem Number[elemX][elemY], int NumElem, int NumNodes, Coord NodeCoord[NumNodes], double *alpha, short *redun, double hxz[elemX], double *theta, int numply, double plyt, double *Trans)
{

	/*Note during the asembly process the KE value used for each column (whih changes with x) is correct for an element of that size. The key thing is to
	get the relative area for that element and multiply that by the stiffness matrix. IF the element is fully in the KE slection can be used strait.*/
	FILE *outfile;	/*File varible for output files*/
	int n,m,o,i,j,num,temp;
	int nNod; /*Varaibale for number of nodes for a NIO element*/
	short NIO_ElemStat; /*Varibale for the status of an NIO Element*/
	int *tnodes;	/*tempery array for node numbers (upto 5 for pentagonal areas)*/
	tnodes = malloc(5 * sizeof(int));
	Coord *tcords;	/*array to store node (and interception point) co-ords*/
	tcords = malloc(5 * sizeof(Coord));
	int IsCount = 0;	/*Variable to track current NIO number*/
	int begin = 0; /*Begin variable used to track how many entries have been indexed so far*/
	double hxztemp;
	double *KEQuad; /*2D Array for an Quadrilateral Stiffness Matrix*/
	KEQuad = malloc(576 * sizeof(double));
	double *thetatemp;
	thetatemp = malloc(numply*sizeof(double));
	double atemp; /*area ratio variable*/
	double AreaElem;	/*element area*/
	double X;
	double C, S, C4, S4, C2S2, CS, C3S, CS3, C2, S2, C3, S3, C2S, CS2; 	/*Cosine and Sine varibles*/
	double Q11, Q22, Q12, Q16, Q26, Q66, Q44, Q45, Q55;			/*Qmatrix varibles*/
	double E11, E22, E12, E66, E44, E55;					/*Material Property varibles for strain*/

	/*For All Elements*/
	for(m=0;m<elemY;m++)
	{	
		for(n=0;n<elemX;n++)
		{
			num = Number[n][m].n - 1;	/*Number of current element, -1 for way C stores data in Arrays*/
			AreaElem = h * hxz[n];
			
			/*Perform Area calcualtion so update of zero lsf line will be stable*/
			
			tnodes[0] = Number[n][m].a - 1;
			tnodes[1] = Number[n][m].b - 1;
			tnodes[2] = Number[n][m].c - 1;
			tnodes[3] = Number[n][m].d - 1;

			alpha[num] = 1.0;

			
					
			for(o=0;o<4;o++)
			{
				temp = tnodes[o] * 6;
				redun[temp] = 1;
				redun[temp+1] = 1; /*indicate the dof's are not redundant*/
				redun[temp+2] = 1;
				redun[temp+3] = 1;
				redun[temp+4] = 1;
				redun[temp+5] = 1;
			}
		
			hxztemp = hxz[n];
			for(o=0;o<numply;o++)
			{
				thetatemp[o] = theta[numply*num + o];
				/*printf("\nTheta = %f", theta[numply*num + o]);
				printf("\nTheta = %f", thetatemp[o]);*/
			}
			CompositeKE(num,e1,e2,v12,v21,g12,g23,g13, KE, thetatemp, numply, plyt, h, hxztemp);

			Transform(KE, Trans, KEQuad, n);

			/*multiply In stiffness matrix by area ratio*/
			/*for(i=0;i<576;i++)
			{
				KEQuad[i] = KE[i]; /*0 atemp;*
			}*/
			/*Assemble area ratio weighted stiffness matrix into index arrays*/
			Assemble2(begin, 4, tnodes, KEQuad, irn, jcn, A);
			begin += 300;		
		}
	}
outfile = fopen("QuadMatrix.txt", "w");
if(outfile == NULL){
	printf("\nFailed to open Numbering writefile\n");
	}
else{
for(i=0;i<24;i++)
{
	for(j=0;j<24;j++)
	{
		fprintf(outfile,"%0.2f  ", KEQuad[j+i*24]);

	}
	fprintf(outfile,"\n");
}
}

fclose(outfile);

	free(tnodes);
	free(tcords);
	free(KEQuad);
	free(thetatemp);


}

void AFG_Strain_LS(int elemX, int elemY, Elem Number[elemX][elemY], short *ElemStat, short *NodeStat, double *alpha, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double h, Aux *auxNode, double *Nenergy, double *rhs1, double *rhs2, int NumNodes, Coord NodeCoord[NumNodes], int gpoints, int Acount, int NIOtotal, double *Eenergy, short Hflag, double wgt, double hxz[elemX], double hz[elemX], double *theta, double *Trans, int numply, double plyt)
{
	FILE *outfile;	/*File varible for output files*/
	double rad = 2.0 * h; /*hard coded for two elements around a node*/
	double AreaElem,ftemp,hxztemp,hztemp;

	/*printf("\nrad = %f",rad);*/
	int n,m,o,num,temp;
	int Gcount = 0;		/*varible to track number of points evaluated*/
	Gstn *Genergy;	/*array to store strain energy data for least squares filtering*/
	/*int gpoints;	/*variable for number of gauss & intergration points in structure*/
	double cx,cy;	/*variables for grid element center*/
	int *tnodes;	/*tempery array for node numbers (upto 5 for pentagonals)*/
	tnodes = malloc(4 * sizeof(int));
	double *thetatemp;
	thetatemp = malloc(numply*sizeof(double));
	
	/*gpoints = Itotal + qCount + tCount + pCount;	/*number of point = no. of non OUT elements*/
	Genergy = calloc(gpoints,sizeof(Gstn));	/*define memory for least squares method*/
	/*Genergy = calloc(gpoints,sizeof(double));	/*define memory for superconvergent method*/
	/*For All Elements*/
	for(m=0;m<elemY;m++)
	{	
		for(n=0;n<elemX;n++)
		{
			num = Number[n][m].n - 1;	/*Number of current element, -1 for way C stores data in Arrays*/
			/*printf("\nnum=%i n=%i, m=%i",num,n,m);*/
			AreaElem = h * hxz[n];
			hxztemp = hxz[n];
			hztemp = hz[n];
			if(ElemStat[num] != -1) /* For all elements*/
			{	
				/*printf("\nElememt %i is IN!!",num+1);*/
				tnodes[0] = Number[n][m].a - 1;
				tnodes[1] = Number[n][m].b - 1;
				tnodes[2] = Number[n][m].c - 1;
				tnodes[3] = Number[n][m].d - 1;
				for(o=0;o<numply;o++)
				{
					thetatemp[o] = theta[numply*num + o];
				}
				/*calcualte strain energy at Gauss points*/
				GaINstrainV2(tnodes, rhs1, rhs2, 1.0, e1, e2, v12, v21, g12, g23, g13, h, Gcount, Genergy, NumNodes, NodeCoord, num, Eenergy, AreaElem, wgt, hxztemp, hztemp, thetatemp, Trans, numply, plyt, n);

				Gcount += 4; /*update point count*/
				
				/*calcualte strain energy at element center
				CeINstrain(tnodes, rhs, alpha[num], e11, v12, g33, h, Gcount, Genergy,NumNodes, NodeCoord);
				Gcount ++; /*update point count*/
                
                
			}
            
            
            
		}
	}
	
	/* Write Gauss point Strain Energy information file 
	outfile = fopen("GaussStrain.txt", "w");
	if(outfile == NULL){
		printf("\nFailed to open Gauss Point Strain Energy writefile\n");
		}
	else{
	fprintf(outfile,"node\tX\tY\tUn\tarea\n"); /*column headings*/
	/*	for(n=0;n<gpoints;n++)
		{
			/*fprintf(outfile,"%i\t",(n+1));
			fprintf(outfile,"%f\t",Genergy[n].x);
			fprintf(outfile,"%f\t",Genergy[n].y);
			fprintf(outfile,"%f\t",Genergy[n].u);
			fprintf(outfile,"%f",Genergy[n].a);

			fprintf(outfile,"\n");
		}
	}

	fclose(outfile);*/

	int stat = (Hflag == 1) ? 0 : 1; /*variable to indicate calcualtion of boundary or all Inside & boundary node strain energies*/
	
	/*now calculate filtered strain energy for all non-OUT nodes*/
	for(n=0;n<Acount;n++)
	{
		if((n >= NumNodes) || (NodeStat[n] > -1/*stat*/))  /*if node auxillary or on boundary (or within structure) */
		{						  /*We want every node to be included in the node senstivties so we can use them as 
								  the implicit fuction in this tow optimization model*/
			cx = -1.0; /*set to monitor errors*/
			/*first need to find co-ordinates of current node*/
			if(n < NumNodes) /*If node is a grid node*/
			{
				cx = NodeCoord[n].xz;
				cy = NodeCoord[n].y; /*read in co-ordiantes directly*/
			}
			else /*Otherwise node is an auxillary node*/
			{
				for(m=0;m<(2 * NIOtotal);m++) /*for all entries in auxNode*/
				{
					if(n == auxNode[m].n-1) /*find an auxNode with matching node number*/
					{
						cx = auxNode[m].xz;
						cy = auxNode[m].y; /*read in auxillary node co-ordinates*/
						break;
					}
				}
			}
			
			if(cx == -1.0)
			{
				printf("\nERROR! could not find co-ordinates of node %i in Strain calc",n+1);
			}
			/*use least squares to calculate nodal strain energy*/
			ftemp = LstrainV22(cx, cy, rad, gpoints, Genergy,4); /*weighted by alpha / distance*/
			Nenergy[n] += wgt * ftemp; /* multiply smothed strain energy by weight */
		}
	}
	free(tnodes);
	free(Genergy);
}

void CompositeKE(int num, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double *KE, double *theta,int numply, double plyt, int h, double hxz)
{
/*NOTE For the sake of the dgemm functions set up here is column strong (arrays go down thecolumns first, then along the rows). The fina out put however will be the other way around so that we can use it in assmeble*/

int n,m,i,j,p,q;
double dn1dx,dn1dy,dn2dx,dn2dy,dn3dx,dn3dy,dn4dx,dn4dy,nx,ny,n1j,n2j,n3j,n4j;
/*Intialize all the matrices we'll be using*/
double *A, *D, *B, *C0matrix, *Cmatrix, *F, *Q0matrix, *Qmatrix;
/*Material property matrices*/
A = malloc(9*sizeof(double));
D = malloc(9*sizeof(double));
B = malloc(9*sizeof(double));
F = malloc(4*sizeof(double));
Cmatrix = malloc(4*sizeof(double));
C0matrix = malloc(4*sizeof(double));
Q0matrix = malloc(9*sizeof(double));
Qmatrix = malloc(9*sizeof(double));
/*Shape function Matrices*/
double *BB,*BM,*BS;
BM = malloc(24*sizeof(double));
BB = malloc(36*sizeof(double));
BS = malloc(24*sizeof(double));

/*StiffnessMatrixes*/
double *Km, *Kp, *tempKm, *tempKp, *tempKps;
Km = malloc(64*sizeof(double));
Kp = malloc(144*sizeof(double));
tempKm = malloc(24*sizeof(double));
tempKp = malloc(36*sizeof(double));
tempKps = malloc(36*sizeof(double));

char transA;
char transB;
int Arows;
int Bcol;
int Acol;
double alpha = 1.0;
double beta =0.0;
int LDA;
int LDB;
int LDC;

Q0matrix[0] = e1/(1-v12*v21);
Q0matrix[4] = e2/(1-v12*v21);
Q0matrix[1] = (v21*e1)/(1-v12*v21);
Q0matrix[8] = g12;

C0matrix[0] = g13;
C0matrix[3] = g23;

double c,s,c2,s2,c3,s3,c4,s4;

	for(n=0;n<3;n++)
	{
		for(m=0;m<3;m++)
		{
			A[n+3*m] = 0;
			D[n+3*m] = 0;
			B[n+3*m] = 0;
		}
	}
	for(p=0;p<8;p++)
	{
		for(q=0;q<8;q++)
		{
			Km[p+8*q] = 0;
		}
	}
	for(p=0;p<12;p++)
	{
		for(q=0;q<12;q++)
		{
			Kp[p+12*q] = 0;
		}
	}				
		
	for(q=0;q<numply;q++)
	{
		/*if((num == 361)||(num==378)){printf("\ntheta[%i] = %f : %f",q, theta[q], theta[q]*(180/3.14));}*/
		c = cos(theta[q]);
		s = sin(theta[q]);
		/*if((num == 361)||(num==378)){printf("\nc = %f, s = %f", c, s);}*/
		c2 = c*c;
		s2 = s*s;
		c3 = c2*c;
		s3 = s2*s;
		c4 = c2*c2;
		s4 = s2*s2;
		Qmatrix[0] = c4*Q0matrix[0] + 2*c2*s2*Q0matrix[1] + s4*Q0matrix[4] + 4*c2*s2*Q0matrix[8];
		Qmatrix[3] = c2*s2*Q0matrix[0] + (c4+s4)*Q0matrix[1] + c2*s2*Q0matrix[4] - 4*c2*s2*Q0matrix[8];
		Qmatrix[1] = c2*s2*Q0matrix[0] + (c4+s4)*Q0matrix[1] + c2*s2*Q0matrix[4] - 4*c2*s2*Q0matrix[8];
		Qmatrix[4] = s4*Q0matrix[0] + 2*c2*s2*Q0matrix[1] + c4*Q0matrix[4] + 4*c2*s2*Q0matrix[8];
		Qmatrix[2] = c3*s*Q0matrix[0] - c*s*(c2-s2)*Q0matrix[1] - c*s3*Q0matrix[4] - 2*c*s*(c2-s2)*Q0matrix[8];
		Qmatrix[6] = c3*s*Q0matrix[0] - c*s*(c2-s2)*Q0matrix[1] - c*s3*Q0matrix[4] - 2*c*s*(c2-s2)*Q0matrix[8];
		Qmatrix[7] = c*s3*Q0matrix[0] + c*s*(c2-s2)*Q0matrix[1] - c3*s*Q0matrix[4] + 2*c*s*(c2-s2)*Q0matrix[8];
		Qmatrix[5] = c*s3*Q0matrix[0] + c*s*(c2-s2)*Q0matrix[1] - c3*s*Q0matrix[4] + 2*c*s*(c2-s2)*Q0matrix[8];
		Qmatrix[8] = c2*s2*Q0matrix[0] - 2*c2*s2*Q0matrix[1] + c2*s2*Q0matrix[4] + (c2-s2)*(c2-s2)*Q0matrix[8];

		Cmatrix[0] = s2*C0matrix[3] + C0matrix[0]*c2;
		Cmatrix[3] = c2*C0matrix[3] + C0matrix[0]*s2;
		Cmatrix[1] = (C0matrix[3] - C0matrix[0])*c*s;
		Cmatrix[2] = (C0matrix[3] - C0matrix[0])*c*s;

		for(n=0;n<3;n++)
		{
			for(m=0;m<3;m++)
			{
				A[n+3*m] += plyt*Qmatrix[n+3*m];
				B[n+3*m] += -1*0.5*plyt*plyt*Qmatrix[n+3*m];
				D[n+3*m] += (plyt*plyt*plyt*Qmatrix[n+3*m])/12.0;
				if((n<2)&&(m<2)){F[n+2*m] = (5.0/6.0)*plyt*Cmatrix[n+2*m];}
			}
		}				
	}

	/*printf("\nC0matrix\n");
	printf("%f\t %f\n", C0matrix[0], C0matrix[2]);
	printf("%f\t %f\n", C0matrix[1], C0matrix[3]);

	printf("\nCmatrix\n");
	printf("%f\t %f\n", Cmatrix[0], Cmatrix[2]);
	printf("%f\t %f\n", Cmatrix[1], Cmatrix[3]);

	printf("\nQ0matrix\n");
	printf("%f\t %f\t %f\n", Q0matrix[0], Q0matrix[3], Q0matrix[6]);
	printf("%f\t %f\t %f\n", Q0matrix[1], Q0matrix[4], Q0matrix[7]);
	printf("%f\t %f\t %f\n", Q0matrix[2], Q0matrix[5], Q0matrix[8]);


	printf("\nQmatrix\n");
	printf("%f\t %f\t %f\n", Qmatrix[0], Qmatrix[3], Qmatrix[6]);
	printf("%f\t %f\t %f\n", Qmatrix[1], Qmatrix[4], Qmatrix[7]);
	printf("%f\t %f\t %f\n", Qmatrix[2], Qmatrix[5], Qmatrix[8]);*/

	/*printf("\nAmatrix\n");
	printf("%f\t %f\t %f\n", A[0], A[3], A[6]);
	printf("%f\t %f\t %f\n", A[1], A[4], A[7]);
	printf("%f\t %f\t %f\n", A[2], A[5], A[8]);

	/*printf("\nBmatrix\n");
	printf("%f\t %f\t %f\n", Bmatrix[0], Bmatrix[3], Bmatrix[6]);
	printf("%f\t %f\t %f\n", Bmatrix[1], Bmatrix[4], Bmatrix[7]);
	printf("%f\t %f\t %f\n", Bmatrix[2], Bmatrix[5], Bmatrix[8]);*

	printf("\nDmatrix\n");
	printf("%f\t %f\t %f\n", D[0], D[3], D[6]);
	printf("%f\t %f\t %f\n", D[1], D[4], D[7]);
	printf("%f\t %f\t %f\n", D[2], D[5], D[8]);

	printf("\nFmatrix\n");
	printf("%f\t %f\n", F[0], F[2]);
	printf("%f\t %f\n", F[1], F[3]);*/



for(i=0;i<4;i++)
{
		if((i==0)||(i==3)){nx = -0.5773502692;}
		else{nx = 0.5773502692;}

		if((i==0)||(i==1)){ny = -0.5773502692;}
		else{ny = 0.5773502692;}

		n1j = 0.125*(1-nx)*(1-ny);
		n2j = 0.125*(1+nx)*(1-ny);
		n3j = 0.125*(1+nx)*(1+ny);
		n4j = 0.125*(1-nx)*(1+ny);

		dn1dx = -0.25*(1-ny)/hxz;
		dn2dx = 0.25*(1-ny)/hxz;
		dn3dx = 0.25*(1+ny)/hxz;
		dn4dx = -0.25*(1+ny)/hxz;

		dn1dy = -0.25*(1-nx)/h;
		dn2dy = -0.25*(1+nx)/h;
		dn3dy = 0.25*(1+nx)/h;
		dn4dy = 0.25*(1-nx)/h;

		/*printf("\n dn1dx = %f", dn1dx);
		printf("\n dn1dx = %f", dn2dx);
		printf("\n dn1dx = %f", dn3dx);
		printf("\n dn1dx = %f", dn4dx);

		printf("\n dn1dx = %f", dn1dy);
		printf("\n dn1dx = %f", dn2dy);
		printf("\n dn1dx = %f", dn3dy);
		printf("\n dn1dx = %f", dn4dy);

		printf("\n n1j = %f", n1j);
		printf("\n n2j = %f", n2j);
		printf("\n n3j = %f", n3j);
		printf("\n n4j = %f\n", n4j);*/

		BM[0] = dn1dx;
		BM[3] = 0;
		BM[6] = dn2dx;
		BM[9] = 0;
		BM[12] = dn3dx;
		BM[15] = 0;
		BM[18] = dn4dx;
		BM[21] = 0;
		
		BM[1] = 0;
		BM[4] = dn1dy;
		BM[7] = 0;
		BM[10] = dn2dy;
		BM[13] = 0;
		BM[16] = dn3dy;
		BM[19] = 0;
		BM[22] = dn4dy;

		BM[2] = dn1dy;
		BM[5] = dn1dx;
		BM[8] = dn2dy;
		BM[11] = dn2dx;
		BM[14] = dn3dy;
		BM[17] = dn3dx;
		BM[20] = dn4dy;
		BM[23] = dn4dx;

		BB[0] = 0;
		BB[3] = 0;
		BB[6] = -dn1dx;
		BB[9] = 0;
		BB[12] = 0;
		BB[15] = -dn2dx;
		BB[18] = 0;
		BB[21] = 0;
		BB[24] = -dn3dx;
		BB[27] = 0;
		BB[30] = 0;
		BB[33] = -dn4dx;

		BB[1] = 0;
		BB[4] = dn1dy;
		BB[7] = 0;
		BB[10] = 0;
		BB[13] = dn2dy;
		BB[16] = 0;
		BB[19] = 0;
		BB[22] = dn3dy;
		BB[25] = 0;
		BB[28] = 0;
		BB[31] = dn4dy;
		BB[34] = 0;

		BB[2] = 0;
		BB[5] = dn1dx;
		BB[8] = -dn1dy;
		BB[11] = 0;
		BB[14] = dn2dx;
		BB[17] = -dn2dy;
		BB[20] = 0;
		BB[23] = dn3dx;
		BB[26] = -dn3dy;
		BB[29] = 0;
		BB[32] = dn4dx;
		BB[35] = -dn4dy;

/*printf("BMshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[0], BM[3], BM[6], BM[9], BM[12], BM[15], BM[18], BM[21]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[1], BM[4], BM[7], BM[10], BM[13], BM[16], BM[19], BM[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[2], BM[5], BM[8], BM[11], BM[14], BM[17], BM[20], BM[23]);


printf("BBshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[0], BB[3], BB[6], BB[9], BB[12], BB[15], BB[18], BB[21], BB[24], BB[27], BB[30], BB[33]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[1], BB[4], BB[7], BB[10], BB[13], BB[16], BB[19], BB[22], BB[25], BB[28], BB[31], BB[34]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[2], BB[5], BB[8], BB[11], BB[14], BB[17], BB[20], BB[23], BB[26], BB[29], BB[32], BB[35]);*/



/*printf("\nA");
printf("\nA");*/

/*Membrain terms*/		
Km[0] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[0] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[1] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[2]);
Km[1] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[0] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[1] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[2]);
Km[2] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[0] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[1] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[2]);
Km[3] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[0] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[1] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[2]);
Km[4] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[0] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[1] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[2]);
Km[5] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[0] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[1] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[2]);
Km[6] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[0] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[1] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[2]);
Km[7] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[0] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[1] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[2]);

Km[8] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[3] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[4] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[5]);
Km[9] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[3] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[4] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[5]);
Km[10] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[3] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[4] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[5]);
Km[11] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[3] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[4] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[5]);
Km[12] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[3] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[4] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[5]);
Km[13] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[3] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[4] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[5]);
Km[14] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[3] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[4] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[5]);
Km[15] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[3] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[4] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[5]);

Km[16] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[6] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[7] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[8]);
Km[17] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[6] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[7] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[8]);
Km[18] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[6] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[7] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[8]);
Km[19] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[6] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[7] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[8]);
Km[20] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[6] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[7] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[8]);
Km[21] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[6] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[7] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[8]);
Km[22] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[6] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[7] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[8]);
Km[23] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[6] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[7] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[8]);

Km[24] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[9] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[10] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[11]);
Km[25] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[9] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[10] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[11]);
Km[26] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[9] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[10] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[11]);
Km[27] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[9] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[10] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[11]);
Km[28] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[9] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[10] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[11]);
Km[29] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[9] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[10] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[11]);
Km[30] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[9] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[10] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[11]);
Km[31] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[9] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[10] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[11]);
		
Km[32] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[12] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[13] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[14]);
Km[33] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[12] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[13] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[14]);
Km[34] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[12] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[13] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[14]);
Km[35] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[12] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[13] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[14]);
Km[36] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[12] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[13] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[14]);
Km[37] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[12] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[13] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[14]);
Km[38] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[12] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[13] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[14]);
Km[39] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[12] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[13] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[14]);

Km[40] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[15] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[16] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[17]);
Km[41] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[15] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[16] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[17]);
Km[42] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[15] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[16] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[17]);
Km[43] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[15] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[16] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[17]);
Km[44] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[15] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[16] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[17]);
Km[45] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[15] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[16] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[17]);
Km[46] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[15] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[16] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[17]);
Km[47] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[15] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[16] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[17]);

Km[48] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[18] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[19] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[20]);
Km[49] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[18] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[19] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[20]);
Km[50] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[18] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[19] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[20]);
Km[51] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[18] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[19] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[20]);
Km[52] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[18] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[19] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[20]);
Km[53] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[18] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[19] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[20]);
Km[54] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[18] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[19] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[20]);
Km[55] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[18] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[19] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[20]);

Km[56] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[21] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[22] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[23]);
Km[57] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[21] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[22] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[23]);
Km[58] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[21] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[22] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[23]);
Km[59] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[21] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[22] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[23]);
Km[60] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[21] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[22] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[23]);
Km[61] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[21] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[22] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[23]);
Km[62] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[21] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[22] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[23]);
Km[63] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[21] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[22] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[23]);


/*Bending terms*/
Kp[0] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[0] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[1] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[2]);
Kp[1] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[0] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[1] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[2]);
Kp[2] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[0] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[1] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[2]);
Kp[3] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[0] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[1] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[2]);
Kp[4] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[0] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[1] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[2]);
Kp[5] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[0] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[1] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[2]);
Kp[6] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[0] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[1] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[2]);
Kp[7] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[0] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[1] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[2]);
Kp[8] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[0] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[1] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[2]);
Kp[9] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[0] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[1] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[2]);
Kp[10] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[0] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[1] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[2]);
Kp[11] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[0] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[1] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[2]);

Kp[12] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[3] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[4] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[5]);
Kp[13] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[3] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[4] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[5]);
Kp[14] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[3] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[4] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[5]);
Kp[15] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[3] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[4] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[5]);
Kp[16] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[3] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[4] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[5]);
Kp[17] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[3] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[4] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[5]);
Kp[18] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[3] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[4] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[5]);
Kp[19] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[3] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[4] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[5]);
Kp[20] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[3] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[4] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[5]);
Kp[21] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[3] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[4] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[5]);
Kp[22] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[3] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[4] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[5]);
Kp[23] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[3] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[4] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[5]);

Kp[24] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[6] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[7] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[8]);
Kp[25] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[6] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[7] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[8]);
Kp[26] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[6] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[7] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[8]);
Kp[27] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[6] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[7] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[8]);
Kp[28] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[6] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[7] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[8]);
Kp[29] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[6] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[7] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[8]);
Kp[30] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[6] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[7] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[8]);
Kp[31] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[6] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[7] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[8]);
Kp[32] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[6] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[7] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[8]);
Kp[33] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[6] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[7] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[8]);
Kp[34] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[6] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[7] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[8]);
Kp[35] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[6] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[7] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[8]);

Kp[36] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[9] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[10] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[11]);
Kp[37] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[9] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[10] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[11]);
Kp[38] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[9] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[10] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[11]);
Kp[39] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[9] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[10] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[11]);
Kp[40] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[9] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[10] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[11]);
Kp[41] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[9] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[10] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[11]);
Kp[42] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[9] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[10] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[11]);
Kp[43] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[9] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[10] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[11]);
Kp[44] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[9] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[10] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[11]);
Kp[45] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[9] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[10] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[11]);
Kp[46] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[9] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[10] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[11]);
Kp[47] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[9] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[10] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[11]);

Kp[48] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[12] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[13] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[14]);
Kp[49] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[12] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[13] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[14]);
Kp[50] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[12] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[13] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[14]);
Kp[51] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[12] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[13] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[14]);
Kp[52] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[12] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[13] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[14]);
Kp[53] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[12] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[13] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[14]);
Kp[54] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[12] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[13] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[14]);
Kp[55] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[12] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[13] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[14]);
Kp[56] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[12] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[13] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[14]);
Kp[57] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[12] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[13] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[14]);
Kp[58] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[12] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[13] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[14]);
Kp[59] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[12] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[13] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[14]);

Kp[60] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[15] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[16] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[17]);
Kp[61] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[15] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[16] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[17]);
Kp[62] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[15] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[16] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[17]);
Kp[63] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[15] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[16] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[17]);
Kp[64] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[15] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[16] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[17]);
Kp[65] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[15] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[16] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[17]);
Kp[66] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[15] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[16] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[17]);
Kp[67] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[15] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[16] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[17]);
Kp[68] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[15] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[16] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[17]);
Kp[69] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[15] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[16] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[17]);
Kp[70] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[15] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[16] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[17]);
Kp[71] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[15] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[16] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[17]);

Kp[72] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[18] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[19] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[20]);
Kp[73] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[18] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[19] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[20]);
Kp[74] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[18] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[19] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[20]);
Kp[75] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[18] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[19] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[20]);
Kp[76] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[18] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[19] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[20]);
Kp[77] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[18] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[19] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[20]);
Kp[78] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[18] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[19] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[20]);
Kp[79] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[18] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[19] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[20]);
Kp[80] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[18] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[19] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[20]);
Kp[81] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[18] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[19] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[20]);
Kp[82] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[18] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[19] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[20]);
Kp[83] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[18] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[19] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[20]);

Kp[84] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[21] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[22] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[23]);
Kp[85] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[21] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[22] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[23]);
Kp[86] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[21] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[22] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[23]);
Kp[87] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[21] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[22] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[23]);
Kp[88] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[21] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[22] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[23]);
Kp[89] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[21] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[22] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[23]);
Kp[90] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[21] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[22] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[23]);
Kp[91] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[21] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[22] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[23]);
Kp[92] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[21] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[22] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[23]);
Kp[93] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[21] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[22] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[23]);
Kp[94] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[21] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[22] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[23]);
Kp[95] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[21] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[22] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[23]);

Kp[96] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[24] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[25] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[26]);
Kp[97] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[24] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[25] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[26]);
Kp[98] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[24] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[25] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[26]);
Kp[99] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[24] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[25] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[26]);
Kp[100] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[24] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[25] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[26]);
Kp[101] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[24] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[25] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[26]);
Kp[102] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[24] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[25] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[26]);
Kp[103] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[24] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[25] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[26]);
Kp[104] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[24] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[25] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[26]);
Kp[105] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[24] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[25] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[26]);
Kp[106] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[24] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[25] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[26]);
Kp[107] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[24] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[25] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[26]);

Kp[108] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[27] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[28] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[29]);
Kp[109] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[27] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[28] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[29]);
Kp[110] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[27] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[28] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[29]);
Kp[111] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[27] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[28] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[29]);
Kp[112] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[27] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[28] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[29]);
Kp[113] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[27] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[28] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[29]);
Kp[114] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[27] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[28] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[29]);
Kp[115] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[27] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[28] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[29]);
Kp[116] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[27] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[28] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[29]);
Kp[117] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[27] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[28] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[29]);
Kp[118] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[27] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[28] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[29]);
Kp[119] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[27] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[28] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[29]);

Kp[120] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[30] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[31] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[32]);
Kp[121] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[30] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[31] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[32]);
Kp[122] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[30] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[31] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[32]);
Kp[123] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[30] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[31] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[32]);
Kp[124] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[30] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[31] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[32]);
Kp[125] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[30] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[31] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[32]);
Kp[126] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[30] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[31] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[32]);
Kp[127] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[30] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[31] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[32]);
Kp[128] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[30] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[31] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[32]);
Kp[129] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[30] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[31] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[32]);
Kp[130] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[30] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[31] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[32]);
Kp[131] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[30] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[31] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[32]);

Kp[132] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[33] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[34] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[35]);
Kp[133] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[33] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[34] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[35]);
Kp[134] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[33] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[34] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[35]);
Kp[135] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[33] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[34] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[35]);
Kp[136] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[33] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[34] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[35]);
Kp[137] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[33] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[34] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[35]);
Kp[138] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[33] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[34] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[35]);
Kp[139] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[33] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[34] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[35]);
Kp[140] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[33] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[34] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[35]);
Kp[141] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[33] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[34] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[35]);
Kp[142] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[33] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[34] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[35]);
Kp[143] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[33] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[34] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[35]);
}
/*Shear terms go through single guass intergration only*/
nx = 0;
ny = 0;

n1j = 0.25*(1-nx)*(1-ny);
n2j = 0.25*(1+nx)*(1-ny);
n3j = 0.25*(1+nx)*(1+ny);
n4j = 0.25*(1-nx)*(1+ny);

dn1dx = -0.5*(1-ny)/hxz;
dn2dx = 0.5*(1-ny)/hxz;
dn3dx = 0.5*(1+ny)/hxz;
dn4dx = -0.5*(1+ny)/hxz;

dn1dy = -0.5*(1-nx)/h;
dn2dy = -0.5*(1+nx)/h;
dn3dy = 0.5*(1+nx)/h;
dn4dy = 0.5*(1-nx)/h;


BS[0] = dn1dx;
BS[2] = 0;
BS[4] = n1j;
BS[6] = dn2dx;
BS[8] = 0;
BS[10] = n2j;
BS[12] = dn3dx;
BS[14] = 0;
BS[16] = n3j;
BS[18] = dn4dx;
BS[20] = 0;
BS[22] = n4j;

BS[1] = dn1dy;
BS[3] = -n1j;
BS[5] = 0;
BS[7] = dn2dy;
BS[9] = -n2j;
BS[11] = 0;
BS[13] = dn3dy;
BS[15] = -n3j;
BS[17] = 0;
BS[19] = dn4dy;
BS[21] = -n4j;
BS[23] = 0;

/*printf("BSshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[0], BS[2], BS[4], BS[6], BS[8], BS[10], BS[12], BS[14], BS[16], BS[18], BS[20], BS[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[1], BS[3], BS[5], BS[7], BS[9], BS[11], BS[13], BS[15], BS[17], BS[19], BS[21], BS[23]);*/

/*Shear terms*/
Kp[0] += ((BS[0]*F[0]+BS[1]*F[1])*BS[0] + (BS[0]*F[2]+BS[1]*F[3])*BS[1]);
Kp[1] += ((BS[2]*F[0]+BS[3]*F[1])*BS[0] + (BS[2]*F[2]+BS[3]*F[3])*BS[1]);
Kp[2] += ((BS[4]*F[0]+BS[5]*F[1])*BS[0] + (BS[4]*F[2]+BS[5]*F[3])*BS[1]);
Kp[3] += ((BS[6]*F[0]+BS[7]*F[1])*BS[0] + (BS[6]*F[2]+BS[7]*F[3])*BS[1]);
Kp[4] += ((BS[8]*F[0]+BS[9]*F[1])*BS[0] + (BS[8]*F[2]+BS[9]*F[3])*BS[1]);
Kp[5] += ((BS[10]*F[0]+BS[11]*F[1])*BS[0] + (BS[10]*F[2]+BS[11]*F[3])*BS[1]);
Kp[6] += ((BS[12]*F[0]+BS[13]*F[1])*BS[0] + (BS[12]*F[2]+BS[13]*F[3])*BS[1]);
Kp[7] += ((BS[14]*F[0]+BS[15]*F[1])*BS[0] + (BS[14]*F[2]+BS[15]*F[3])*BS[1]);
Kp[8] += ((BS[16]*F[0]+BS[17]*F[1])*BS[0] + (BS[16]*F[2]+BS[17]*F[3])*BS[1]);
Kp[9] += ((BS[18]*F[0]+BS[19]*F[1])*BS[0] + (BS[18]*F[2]+BS[19]*F[3])*BS[1]);
Kp[10] += ((BS[20]*F[0]+BS[21]*F[1])*BS[0] + (BS[20]*F[2]+BS[21]*F[3])*BS[1]);
Kp[11] += ((BS[22]*F[0]+BS[23]*F[1])*BS[0] + (BS[22]*F[2]+BS[23]*F[3])*BS[1]);

Kp[12] += ((BS[0]*F[0]+BS[1]*F[1])*BS[2] + (BS[0]*F[2]+BS[1]*F[3])*BS[3]);
Kp[13] += ((BS[2]*F[0]+BS[3]*F[1])*BS[2] + (BS[2]*F[2]+BS[3]*F[3])*BS[3]);
Kp[14] += ((BS[4]*F[0]+BS[5]*F[1])*BS[2] + (BS[4]*F[2]+BS[5]*F[3])*BS[3]);
Kp[15] += ((BS[6]*F[0]+BS[7]*F[1])*BS[2] + (BS[6]*F[2]+BS[7]*F[3])*BS[3]);
Kp[16] += ((BS[8]*F[0]+BS[9]*F[1])*BS[2] + (BS[8]*F[2]+BS[9]*F[3])*BS[3]);
Kp[17] += ((BS[10]*F[0]+BS[11]*F[1])*BS[2] + (BS[10]*F[2]+BS[11]*F[3])*BS[3]);
Kp[18] += ((BS[12]*F[0]+BS[13]*F[1])*BS[2] + (BS[12]*F[2]+BS[13]*F[3])*BS[3]);
Kp[19] += ((BS[14]*F[0]+BS[15]*F[1])*BS[2] + (BS[14]*F[2]+BS[15]*F[3])*BS[3]);
Kp[20] += ((BS[16]*F[0]+BS[17]*F[1])*BS[2] + (BS[16]*F[2]+BS[17]*F[3])*BS[3]);
Kp[21] += ((BS[18]*F[0]+BS[19]*F[1])*BS[2] + (BS[18]*F[2]+BS[19]*F[3])*BS[3]);
Kp[22] += ((BS[20]*F[0]+BS[21]*F[1])*BS[2] + (BS[20]*F[2]+BS[21]*F[3])*BS[3]);
Kp[23] += ((BS[22]*F[0]+BS[23]*F[1])*BS[2] + (BS[22]*F[2]+BS[23]*F[3])*BS[3]);

Kp[24] += ((BS[0]*F[0]+BS[1]*F[1])*BS[4] + (BS[0]*F[2]+BS[1]*F[3])*BS[5]);
Kp[25] += ((BS[2]*F[0]+BS[3]*F[1])*BS[4] + (BS[2]*F[2]+BS[3]*F[3])*BS[5]);
Kp[26] += ((BS[4]*F[0]+BS[5]*F[1])*BS[4] + (BS[4]*F[2]+BS[5]*F[3])*BS[5]);
Kp[27] += ((BS[6]*F[0]+BS[7]*F[1])*BS[4] + (BS[6]*F[2]+BS[7]*F[3])*BS[5]);
Kp[28] += ((BS[8]*F[0]+BS[9]*F[1])*BS[4] + (BS[8]*F[2]+BS[9]*F[3])*BS[5]);
Kp[29] += ((BS[10]*F[0]+BS[11]*F[1])*BS[4] + (BS[10]*F[2]+BS[11]*F[3])*BS[5]);
Kp[30] += ((BS[12]*F[0]+BS[13]*F[1])*BS[4] + (BS[12]*F[2]+BS[13]*F[3])*BS[5]);
Kp[31] += ((BS[14]*F[0]+BS[15]*F[1])*BS[4] + (BS[14]*F[2]+BS[15]*F[3])*BS[5]);
Kp[32] += ((BS[16]*F[0]+BS[17]*F[1])*BS[4] + (BS[16]*F[2]+BS[17]*F[3])*BS[5]);
Kp[33] += ((BS[18]*F[0]+BS[19]*F[1])*BS[4] + (BS[18]*F[2]+BS[19]*F[3])*BS[5]);
Kp[34] += ((BS[20]*F[0]+BS[21]*F[1])*BS[4] + (BS[20]*F[2]+BS[21]*F[3])*BS[5]);
Kp[35] += ((BS[22]*F[0]+BS[23]*F[1])*BS[4] + (BS[22]*F[2]+BS[23]*F[3])*BS[5]);

Kp[36] += ((BS[0]*F[0]+BS[1]*F[1])*BS[6] + (BS[0]*F[2]+BS[1]*F[3])*BS[7]);
Kp[37] += ((BS[2]*F[0]+BS[3]*F[1])*BS[6] + (BS[2]*F[2]+BS[3]*F[3])*BS[7]);
Kp[38] += ((BS[4]*F[0]+BS[5]*F[1])*BS[6] + (BS[4]*F[2]+BS[5]*F[3])*BS[7]);
Kp[39] += ((BS[6]*F[0]+BS[7]*F[1])*BS[6] + (BS[6]*F[2]+BS[7]*F[3])*BS[7]);
Kp[40] += ((BS[8]*F[0]+BS[9]*F[1])*BS[6] + (BS[8]*F[2]+BS[9]*F[3])*BS[7]);
Kp[41] += ((BS[10]*F[0]+BS[11]*F[1])*BS[6] + (BS[10]*F[2]+BS[11]*F[3])*BS[7]);
Kp[42] += ((BS[12]*F[0]+BS[13]*F[1])*BS[6] + (BS[12]*F[2]+BS[13]*F[3])*BS[7]);
Kp[43] += ((BS[14]*F[0]+BS[15]*F[1])*BS[6] + (BS[14]*F[2]+BS[15]*F[3])*BS[7]);
Kp[44] += ((BS[16]*F[0]+BS[17]*F[1])*BS[6] + (BS[16]*F[2]+BS[17]*F[3])*BS[7]);
Kp[45] += ((BS[18]*F[0]+BS[19]*F[1])*BS[6] + (BS[18]*F[2]+BS[19]*F[3])*BS[7]);
Kp[46] += ((BS[20]*F[0]+BS[21]*F[1])*BS[6] + (BS[20]*F[2]+BS[21]*F[3])*BS[7]);
Kp[47] += ((BS[22]*F[0]+BS[23]*F[1])*BS[6] + (BS[22]*F[2]+BS[23]*F[3])*BS[7]);

Kp[48] += ((BS[0]*F[0]+BS[1]*F[1])*BS[8] + (BS[0]*F[2]+BS[1]*F[3])*BS[9]);
Kp[49] += ((BS[2]*F[0]+BS[3]*F[1])*BS[8] + (BS[2]*F[2]+BS[3]*F[3])*BS[9]);
Kp[50] += ((BS[4]*F[0]+BS[5]*F[1])*BS[8] + (BS[4]*F[2]+BS[5]*F[3])*BS[9]);
Kp[51] += ((BS[6]*F[0]+BS[7]*F[1])*BS[8] + (BS[6]*F[2]+BS[7]*F[3])*BS[9]);
Kp[52] += ((BS[8]*F[0]+BS[9]*F[1])*BS[8] + (BS[8]*F[2]+BS[9]*F[3])*BS[9]);
Kp[53] += ((BS[10]*F[0]+BS[11]*F[1])*BS[8] + (BS[10]*F[2]+BS[11]*F[3])*BS[9]);
Kp[54] += ((BS[12]*F[0]+BS[13]*F[1])*BS[8] + (BS[12]*F[2]+BS[13]*F[3])*BS[9]);
Kp[55] += ((BS[14]*F[0]+BS[15]*F[1])*BS[8] + (BS[14]*F[2]+BS[15]*F[3])*BS[9]);
Kp[56] += ((BS[16]*F[0]+BS[17]*F[1])*BS[8] + (BS[16]*F[2]+BS[17]*F[3])*BS[9]);
Kp[57] += ((BS[18]*F[0]+BS[19]*F[1])*BS[8] + (BS[18]*F[2]+BS[19]*F[3])*BS[9]);
Kp[58] += ((BS[20]*F[0]+BS[21]*F[1])*BS[8] + (BS[20]*F[2]+BS[21]*F[3])*BS[9]);
Kp[59] += ((BS[22]*F[0]+BS[23]*F[1])*BS[8] + (BS[22]*F[2]+BS[23]*F[3])*BS[9]);

Kp[60] += ((BS[0]*F[0]+BS[1]*F[1])*BS[10] + (BS[0]*F[2]+BS[1]*F[3])*BS[11]);
Kp[61] += ((BS[2]*F[0]+BS[3]*F[1])*BS[10] + (BS[2]*F[2]+BS[3]*F[3])*BS[11]);
Kp[62] += ((BS[4]*F[0]+BS[5]*F[1])*BS[10] + (BS[4]*F[2]+BS[5]*F[3])*BS[11]);
Kp[63] += ((BS[6]*F[0]+BS[7]*F[1])*BS[10] + (BS[6]*F[2]+BS[7]*F[3])*BS[11]);
Kp[64] += ((BS[8]*F[0]+BS[9]*F[1])*BS[10] + (BS[8]*F[2]+BS[9]*F[3])*BS[11]);
Kp[65] += ((BS[10]*F[0]+BS[11]*F[1])*BS[10] + (BS[10]*F[2]+BS[11]*F[3])*BS[11]);
Kp[66] += ((BS[12]*F[0]+BS[13]*F[1])*BS[10] + (BS[12]*F[2]+BS[13]*F[3])*BS[11]);
Kp[67] += ((BS[14]*F[0]+BS[15]*F[1])*BS[10] + (BS[14]*F[2]+BS[15]*F[3])*BS[11]);
Kp[68] += ((BS[16]*F[0]+BS[17]*F[1])*BS[10] + (BS[16]*F[2]+BS[17]*F[3])*BS[11]);
Kp[69] += ((BS[18]*F[0]+BS[19]*F[1])*BS[10] + (BS[18]*F[2]+BS[19]*F[3])*BS[11]);
Kp[70] += ((BS[20]*F[0]+BS[21]*F[1])*BS[10] + (BS[20]*F[2]+BS[21]*F[3])*BS[11]);

Kp[71] += ((BS[22]*F[0]+BS[23]*F[1])*BS[10] + (BS[22]*F[2]+BS[23]*F[3])*BS[11]);
Kp[72] += ((BS[0]*F[0]+BS[1]*F[1])*BS[12] + (BS[0]*F[2]+BS[1]*F[3])*BS[13]);
Kp[73] += ((BS[2]*F[0]+BS[3]*F[1])*BS[12] + (BS[2]*F[2]+BS[3]*F[3])*BS[13]);
Kp[74] += ((BS[4]*F[0]+BS[5]*F[1])*BS[12] + (BS[4]*F[2]+BS[5]*F[3])*BS[13]);
Kp[75] += ((BS[6]*F[0]+BS[7]*F[1])*BS[12] + (BS[6]*F[2]+BS[7]*F[3])*BS[13]);
Kp[76] += ((BS[8]*F[0]+BS[9]*F[1])*BS[12] + (BS[8]*F[2]+BS[9]*F[3])*BS[13]);
Kp[77] += ((BS[10]*F[0]+BS[11]*F[1])*BS[12] + (BS[10]*F[2]+BS[11]*F[3])*BS[13]);
Kp[78] += ((BS[12]*F[0]+BS[13]*F[1])*BS[12] + (BS[12]*F[2]+BS[13]*F[3])*BS[13]);
Kp[79] += ((BS[14]*F[0]+BS[15]*F[1])*BS[12] + (BS[14]*F[2]+BS[15]*F[3])*BS[13]);
Kp[80] += ((BS[16]*F[0]+BS[17]*F[1])*BS[12] + (BS[16]*F[2]+BS[17]*F[3])*BS[13]);
Kp[81] += ((BS[18]*F[0]+BS[19]*F[1])*BS[12] + (BS[18]*F[2]+BS[19]*F[3])*BS[13]);
Kp[82] += ((BS[20]*F[0]+BS[21]*F[1])*BS[12] + (BS[20]*F[2]+BS[21]*F[3])*BS[13]);
Kp[83] += ((BS[22]*F[0]+BS[23]*F[1])*BS[12] + (BS[22]*F[2]+BS[23]*F[3])*BS[13]);

Kp[84] += ((BS[0]*F[0]+BS[1]*F[1])*BS[14] + (BS[0]*F[2]+BS[1]*F[3])*BS[15]);
Kp[85] += ((BS[2]*F[0]+BS[3]*F[1])*BS[14] + (BS[2]*F[2]+BS[3]*F[3])*BS[15]);
Kp[86] += ((BS[4]*F[0]+BS[5]*F[1])*BS[14] + (BS[4]*F[2]+BS[5]*F[3])*BS[15]);
Kp[87] += ((BS[6]*F[0]+BS[7]*F[1])*BS[14] + (BS[6]*F[2]+BS[7]*F[3])*BS[15]);
Kp[88] += ((BS[8]*F[0]+BS[9]*F[1])*BS[14] + (BS[8]*F[2]+BS[9]*F[3])*BS[15]);
Kp[89] += ((BS[10]*F[0]+BS[11]*F[1])*BS[14] + (BS[10]*F[2]+BS[11]*F[3])*BS[15]);
Kp[90] += ((BS[12]*F[0]+BS[13]*F[1])*BS[14] + (BS[12]*F[2]+BS[13]*F[3])*BS[15]);
Kp[91] += ((BS[14]*F[0]+BS[15]*F[1])*BS[14] + (BS[14]*F[2]+BS[15]*F[3])*BS[15]);
Kp[92] += ((BS[16]*F[0]+BS[17]*F[1])*BS[14] + (BS[16]*F[2]+BS[17]*F[3])*BS[15]);
Kp[93] += ((BS[18]*F[0]+BS[19]*F[1])*BS[14] + (BS[18]*F[2]+BS[19]*F[3])*BS[15]);
Kp[94] += ((BS[20]*F[0]+BS[21]*F[1])*BS[14] + (BS[20]*F[2]+BS[21]*F[3])*BS[15]);
Kp[95] += ((BS[22]*F[0]+BS[23]*F[1])*BS[14] + (BS[22]*F[2]+BS[23]*F[3])*BS[15]);

Kp[96] += ((BS[0]*F[0]+BS[1]*F[1])*BS[16] + (BS[0]*F[2]+BS[1]*F[3])*BS[17]);
Kp[97] += ((BS[2]*F[0]+BS[3]*F[1])*BS[16] + (BS[2]*F[2]+BS[3]*F[3])*BS[17]);
Kp[98] += ((BS[4]*F[0]+BS[5]*F[1])*BS[16] + (BS[4]*F[2]+BS[5]*F[3])*BS[17]);
Kp[99] += ((BS[6]*F[0]+BS[7]*F[1])*BS[16] + (BS[6]*F[2]+BS[7]*F[3])*BS[17]);
Kp[100] += ((BS[8]*F[0]+BS[9]*F[1])*BS[16] + (BS[8]*F[2]+BS[9]*F[3])*BS[17]);
Kp[101] += ((BS[10]*F[0]+BS[11]*F[1])*BS[16] + (BS[10]*F[2]+BS[11]*F[3])*BS[17]);
Kp[102] += ((BS[12]*F[0]+BS[13]*F[1])*BS[16] + (BS[12]*F[2]+BS[13]*F[3])*BS[17]);
Kp[103] += ((BS[14]*F[0]+BS[15]*F[1])*BS[16] + (BS[14]*F[2]+BS[15]*F[3])*BS[17]);
Kp[104] += ((BS[16]*F[0]+BS[17]*F[1])*BS[16] + (BS[16]*F[2]+BS[17]*F[3])*BS[17]);
Kp[105] += ((BS[18]*F[0]+BS[19]*F[1])*BS[16] + (BS[18]*F[2]+BS[19]*F[3])*BS[17]);
Kp[106] += ((BS[20]*F[0]+BS[21]*F[1])*BS[16] + (BS[20]*F[2]+BS[21]*F[3])*BS[17]);
Kp[107] += ((BS[22]*F[0]+BS[23]*F[1])*BS[16] + (BS[22]*F[2]+BS[23]*F[3])*BS[17]);

Kp[108] += ((BS[0]*F[0]+BS[1]*F[1])*BS[18] + (BS[0]*F[2]+BS[1]*F[3])*BS[19]);
Kp[109] += ((BS[2]*F[0]+BS[3]*F[1])*BS[18] + (BS[2]*F[2]+BS[3]*F[3])*BS[19]);
Kp[110] += ((BS[4]*F[0]+BS[5]*F[1])*BS[18] + (BS[4]*F[2]+BS[5]*F[3])*BS[19]);
Kp[111] += ((BS[6]*F[0]+BS[7]*F[1])*BS[18] + (BS[6]*F[2]+BS[7]*F[3])*BS[19]);
Kp[112] += ((BS[8]*F[0]+BS[9]*F[1])*BS[18] + (BS[8]*F[2]+BS[9]*F[3])*BS[19]);
Kp[113] += ((BS[10]*F[0]+BS[11]*F[1])*BS[18] + (BS[10]*F[2]+BS[11]*F[3])*BS[19]);
Kp[114] += ((BS[12]*F[0]+BS[13]*F[1])*BS[18] + (BS[12]*F[2]+BS[13]*F[3])*BS[19]);
Kp[115] += ((BS[14]*F[0]+BS[15]*F[1])*BS[18] + (BS[14]*F[2]+BS[15]*F[3])*BS[19]);
Kp[116] += ((BS[16]*F[0]+BS[17]*F[1])*BS[18] + (BS[16]*F[2]+BS[17]*F[3])*BS[19]);
Kp[117] += ((BS[18]*F[0]+BS[19]*F[1])*BS[18] + (BS[18]*F[2]+BS[19]*F[3])*BS[19]);
Kp[118] += ((BS[20]*F[0]+BS[21]*F[1])*BS[18] + (BS[20]*F[2]+BS[21]*F[3])*BS[19]);
Kp[119] += ((BS[22]*F[0]+BS[23]*F[1])*BS[18] + (BS[22]*F[2]+BS[23]*F[3])*BS[19]);

Kp[120] += ((BS[0]*F[0]+BS[1]*F[1])*BS[20] + (BS[0]*F[2]+BS[1]*F[3])*BS[21]);
Kp[121] += ((BS[2]*F[0]+BS[3]*F[1])*BS[20] + (BS[2]*F[2]+BS[3]*F[3])*BS[21]);
Kp[122] += ((BS[4]*F[0]+BS[5]*F[1])*BS[20] + (BS[4]*F[2]+BS[5]*F[3])*BS[21]);
Kp[123] += ((BS[6]*F[0]+BS[7]*F[1])*BS[20] + (BS[6]*F[2]+BS[7]*F[3])*BS[21]);
Kp[124] += ((BS[8]*F[0]+BS[9]*F[1])*BS[20] + (BS[8]*F[2]+BS[9]*F[3])*BS[21]);
Kp[125] += ((BS[10]*F[0]+BS[11]*F[1])*BS[20] + (BS[10]*F[2]+BS[11]*F[3])*BS[21]);
Kp[126] += ((BS[12]*F[0]+BS[13]*F[1])*BS[20] + (BS[12]*F[2]+BS[13]*F[3])*BS[21]);
Kp[127] += ((BS[14]*F[0]+BS[15]*F[1])*BS[20] + (BS[14]*F[2]+BS[15]*F[3])*BS[21]);
Kp[128] += ((BS[16]*F[0]+BS[17]*F[1])*BS[20] + (BS[16]*F[2]+BS[17]*F[3])*BS[21]);
Kp[129] += ((BS[18]*F[0]+BS[19]*F[1])*BS[20] + (BS[18]*F[2]+BS[19]*F[3])*BS[21]);
Kp[130] += ((BS[20]*F[0]+BS[21]*F[1])*BS[20] + (BS[20]*F[2]+BS[21]*F[3])*BS[21]);
Kp[131] += ((BS[22]*F[0]+BS[23]*F[1])*BS[20] + (BS[22]*F[2]+BS[23]*F[3])*BS[21]);

Kp[132] += ((BS[0]*F[0]+BS[1]*F[1])*BS[22] + (BS[0]*F[2]+BS[1]*F[3])*BS[23]);
Kp[133] += ((BS[2]*F[0]+BS[3]*F[1])*BS[22] + (BS[2]*F[2]+BS[3]*F[3])*BS[23]);
Kp[134] += ((BS[4]*F[0]+BS[5]*F[1])*BS[22] + (BS[4]*F[2]+BS[5]*F[3])*BS[23]);
Kp[135] += ((BS[6]*F[0]+BS[7]*F[1])*BS[22] + (BS[6]*F[2]+BS[7]*F[3])*BS[23]);
Kp[136] += ((BS[8]*F[0]+BS[9]*F[1])*BS[22] + (BS[8]*F[2]+BS[9]*F[3])*BS[23]);
Kp[137] += ((BS[10]*F[0]+BS[11]*F[1])*BS[22] + (BS[10]*F[2]+BS[11]*F[3])*BS[23]);
Kp[138] += ((BS[12]*F[0]+BS[13]*F[1])*BS[22] + (BS[12]*F[2]+BS[13]*F[3])*BS[23]);
Kp[139] += ((BS[14]*F[0]+BS[15]*F[1])*BS[22] + (BS[14]*F[2]+BS[15]*F[3])*BS[23]);
Kp[140] += ((BS[16]*F[0]+BS[17]*F[1])*BS[22] + (BS[16]*F[2]+BS[17]*F[3])*BS[23]);
Kp[141] += ((BS[18]*F[0]+BS[19]*F[1])*BS[22] + (BS[18]*F[2]+BS[19]*F[3])*BS[23]);
Kp[142] += ((BS[20]*F[0]+BS[21]*F[1])*BS[22] + (BS[20]*F[2]+BS[21]*F[3])*BS[23]);
Kp[143] += ((BS[22]*F[0]+BS[23]*F[1])*BS[22] + (BS[22]*F[2]+BS[23]*F[3])*BS[23]);

/*if((num == 361)||(num==378))
{
printf("\n Membrain Matrix");
for(p=0;p<8;p++)
{
	printf("\n");
	for(j=0;j<8;j++)
	{
		printf("%0.3f\t", Km[j*8+p]);
	}
}
printf("\n");
}*/

/*Check the reuslts so far*/
/*printf("\n Membrain Matrix");
for(p=0;p<8;p++)
{
	printf("\n");
	for(j=0;j<8;j++)
	{
		printf("%0.3f\t", Km[j*8+p]);
	}
}
printf("\n");
printf("\n Plate Matrix");
for(i=0;i<12;i++)
{
	printf("\n");
	for(j=0;j<12;j++)
	{
		printf("%0.3f\t", Kp[j*12+i]);
	}
}*/

/*Assemble Full Matirx (Try to transpose it so that it is correct in C matrix form by not hugely important since it should be symetrical*/
for(i=0;i<576;i++)
{
	KE[i] = 0.00000000000000001;	/*Intialize every value to zero*/
}

/*Fill in the values*/
for(j=0;j<4;j++)
{
	p =24*6*j;
	KE[p] = Km[0+2*j];
	KE[1+p] = Km[8+2*j];
	KE[6+p] = Km[16+2*j];
	KE[7+p] = Km[24+2*j];
	KE[12+p] = Km[32+2*j];
	KE[13+p] = Km[40+2*j];
	KE[18+p] = Km[48+2*j];
	KE[19+p] = Km[56+2*j];
	KE[24+p] = Km[1+2*j];
	KE[25+p] = Km[9+2*j];
	KE[30+p] = Km[17+2*j];
	KE[31+p] = Km[25+2*j];
	KE[36+p] = Km[33+2*j];
	KE[37+p] = Km[41+2*j];
	KE[42+p] = Km[49+2*j];
	KE[43+p] = Km[57+2*j];

	q =24*6*j+24*2;
	KE[2+q] = Kp[0+3*j];
	KE[3+q] = Kp[12+3*j];
	KE[4+q] = Kp[24+3*j];
	KE[8+q] = Kp[36+3*j];
	KE[9+q] = Kp[48+3*j];
	KE[10+q] = Kp[60+3*j];
	KE[14+q] = Kp[72+3*j];
	KE[15+q] = Kp[84+3*j];
	KE[16+q] = Kp[96+3*j];
	KE[20+q] = Kp[108+3*j];

	KE[21+q] = Kp[120+3*j];
	KE[22+q] = Kp[132+3*j];
	q =24*6*j+24*3;
	KE[2+q] = Kp[1+3*j];
	KE[3+q] = Kp[13+3*j];
	KE[4+q] = Kp[25+3*j];
	KE[8+q] = Kp[37+3*j];
	KE[9+q] = Kp[49+3*j];
	KE[10+q] = Kp[61+3*j];
	KE[14+q] = Kp[73+3*j];
	KE[15+q] = Kp[85+3*j];
	KE[16+q] = Kp[97+3*j];
	KE[20+q] = Kp[109+3*j];
	KE[21+q] = Kp[121+3*j];
	KE[22+q] = Kp[133+3*j];
	q =24*6*j+24*4;
	KE[2+q] = Kp[2+3*j];
	KE[3+q] = Kp[14+3*j];
	KE[4+q] = Kp[26+3*j];
	KE[8+q] = Kp[38+3*j];
	KE[9+q] = Kp[50+3*j];
	KE[10+q] = Kp[62+3*j];
	KE[14+q] = Kp[74+3*j];
	KE[15+q] = Kp[86+3*j];
	KE[16+q] = Kp[98+3*j];
	KE[20+q] = Kp[110+3*j];
	KE[21+q] = Kp[122+3*j];
	KE[22+q] = Kp[134+3*j];

}

/*printf("\n Element Matrix");
for(i=0;i<24;i++)
{
	printf("\n");
	for(j=0;j<24;j++)
	{
		printf("%0.2f  ", KE[j+i*24]);

	}
} */

free(A);
free(D);
free(B);
free(F);
free(Cmatrix);
free(C0matrix);
free(Q0matrix);
free(Qmatrix);
free(BB);
free(BM);
free(BS);
free(tempKm);
free(tempKp);
free(tempKps);
free(Km);
free(Kp);
/*printf("\n Temp Membrain Matrix");*/
/*for(p=0;p<12;p++)
{
	printf("\n");
	for(j=0;j<3;j++)
	{
		printf("%0.3f\t", tempKp[j*12+p]);
	}
}*/
/*printf("\nB");
printf("\nB");	*/	

}

void CompositeKEGauss(int GP, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double *KE, double *theta, int numply, double plyt, double h, double hxz)
{
/*NOTE For the sake of the dgemm functions set up here is column strong (arrays go down thecolumns first, then along the rows). The fina out put however will be the other way around so that we can use it in assmeble*/


int n,m,i,j,p,q;
double dn1dx,dn1dy,dn2dx,dn2dy,dn3dx,dn3dy,dn4dx,dn4dy,nx,ny,n1j,n2j,n3j,n4j;
/*Intialize all the matrices we'll be using*/
double *A, *D, *B, *C0matrix, *Cmatrix, *F, *Q0matrix, *Qmatrix;
/*Material property matrices*/
A = malloc(9*sizeof(double));
D = malloc(9*sizeof(double));
B = malloc(9*sizeof(double));
F = malloc(4*sizeof(double));
Cmatrix = malloc(4*sizeof(double));
C0matrix = malloc(4*sizeof(double));
Q0matrix = malloc(9*sizeof(double));
Qmatrix = malloc(9*sizeof(double));
/*Shape function Matrices*/
double *BB,*BM,*BS;
BM = malloc(24*sizeof(double));
BB = malloc(36*sizeof(double));
BS = malloc(24*sizeof(double));

/*StiffnessMatrixes*/
double *Km, *Kp, *tempKm, *tempKp, *tempKps;
Km = malloc(64*sizeof(double));
Kp = malloc(144*sizeof(double));
tempKm = malloc(24*sizeof(double));
tempKp = malloc(36*sizeof(double));
tempKps = malloc(36*sizeof(double));

char transA;
char transB;
int Arows;
int Bcol;
int Acol;
double alpha = 1.0;
double beta =0.0;
int LDA;
int LDB;
int LDC;

Q0matrix[0] = e1/(1-v12*v21);
Q0matrix[4] = e2/(1-v12*v21);
Q0matrix[1] = (v21*e1)/(1-v12*v21);
Q0matrix[8] = g12;

C0matrix[0] = g23;
C0matrix[3] = g13;

double c,s,c2,s2,c3,s3,c4,s4;

	for(n=0;n<3;n++)
	{
		for(m=0;m<3;m++)
		{
			A[n+3*m] = 0;
			D[n+3*m] = 0;
			B[n+3*m] = 0;
		}
	}
	for(p=0;p<8;p++)
	{
		for(q=0;q<8;q++)
		{
			Km[p+8*q] = 0;
		}
	}
	for(p=0;p<12;p++)
	{
		for(q=0;q<12;q++)
		{
			Kp[p+12*q] = 0;
		}
	}				
		
	for(q=0;q<numply;q++)
	{
		/*if((num == 361)||(num==378)){printf("\ntheta[%i] = %f : %f",q, theta[q], theta[q]*(180/3.14));}*/
		c = cos(theta[q]);
		s = sin(theta[q]);
		/*if((num == 361)||(num==378)){printf("\nc = %f, s = %f", c, s);}*/
		c2 = c*c;
		s2 = s*s;
		c3 = c2*c;
		s3 = s2*s;
		c4 = c2*c2;
		s4 = s2*s2;
		Qmatrix[0] = c4*Q0matrix[0] + 2*c2*s2*Q0matrix[1] + s4*Q0matrix[4] + 4*c2*s2*Q0matrix[8];
		Qmatrix[3] = c2*s2*Q0matrix[0] + (c4+s4)*Q0matrix[1] + c2*s2*Q0matrix[4] - 4*c2*s2*Q0matrix[8];
		Qmatrix[1] = c2*s2*Q0matrix[0] + (c4+s4)*Q0matrix[1] + c2*s2*Q0matrix[4] - 4*c2*s2*Q0matrix[8];
		Qmatrix[4] = s4*Q0matrix[0] + 2*c2*s2*Q0matrix[1] + c4*Q0matrix[4] + 4*c2*s2*Q0matrix[8];
		Qmatrix[2] = c3*s*Q0matrix[0] - c*s*(c2-s2)*Q0matrix[1] - c*s3*Q0matrix[4] - 2*c*s*(c2-s2)*Q0matrix[8];
		Qmatrix[6] = c3*s*Q0matrix[0] - c*s*(c2-s2)*Q0matrix[1] - c*s3*Q0matrix[4] - 2*c*s*(c2-s2)*Q0matrix[8];
		Qmatrix[7] = c*s3*Q0matrix[0] + c*s*(c2-s2)*Q0matrix[1] - c3*s*Q0matrix[4] + 2*c*s*(c2-s2)*Q0matrix[8];
		Qmatrix[5] = c*s3*Q0matrix[0] + c*s*(c2-s2)*Q0matrix[1] - c3*s*Q0matrix[4] + 2*c*s*(c2-s2)*Q0matrix[8];
		Qmatrix[8] = c2*s2*Q0matrix[0] - 2*c2*s2*Q0matrix[1] + c2*s2*Q0matrix[4] + (c2-s2)*(c2-s2)*Q0matrix[8];

		Cmatrix[0] = s2*C0matrix[3] + C0matrix[0]*c2;
		Cmatrix[3] = c2*C0matrix[3] + C0matrix[0]*s2;
		Cmatrix[1] = (C0matrix[0] - C0matrix[3])*c*s;
		Cmatrix[2] = (C0matrix[0] - C0matrix[3])*c*s;

		for(n=0;n<3;n++)
		{
			for(m=0;m<3;m++)
			{
				A[n+3*m] += plyt*Qmatrix[n+3*m];
				B[n+3*m] += -1*0.5*plyt*plyt*Qmatrix[n+3*m];
				D[n+3*m] += (plyt*plyt*plyt*Qmatrix[n+3*m])/12.0;
				if((n<2)&&(m<2)){F[n+2*m] = (5.0/6.0)*plyt*Cmatrix[n+2*m];}
			}
		}				
	}
	/*if((num == 361)||(num==378))

	{
	printf("\nC0matrix\n");

	printf("%f\t %f\n", C0matrix[0], C0matrix[2]);

	printf("%f\t %f\n", C0matrix[1], C0matrix[3]);


	printf("\nCmatrix\n");

	printf("%f\t %f\n", Cmatrix[0], Cmatrix[2]);
	printf("%f\t %f\n", Cmatrix[1], Cmatrix[3]);


	printf("\nQ0matrix\n");

	printf("%f\t %f\t %f\n", Q0matrix[0], Q0matrix[3], Q0matrix[6]);
	printf("%f\t %f\t %f\n", Q0matrix[1], Q0matrix[4], Q0matrix[7]);
	printf("%f\t %f\t %f\n", Q0matrix[2], Q0matrix[5], Q0matrix[8]);




	printf("\nQmatrix\n");
	printf("%f\t %f\t %f\n", Qmatrix[0], Qmatrix[3], Qmatrix[6]);
	printf("%f\t %f\t %f\n", Qmatrix[1], Qmatrix[4], Qmatrix[7]);

	printf("%f\t %f\t %f\n", Qmatrix[2], Qmatrix[5], Qmatrix[8]);


	printf("\nAmatrix\n");
	printf("%f\t %f\t %f\n", Amatrix[0], Amatrix[3], Amatrix[6]);
	printf("%f\t %f\t %f\n", Amatrix[1], Amatrix[4], Amatrix[7]);
	printf("%f\t %f\t %f\n", Amatrix[2], Amatrix[5], Amatrix[8]);

	printf("\nBmatrix\n");
	printf("%f\t %f\t %f\n", Bmatrix[0], Bmatrix[3], Bmatrix[6]);
	printf("%f\t %f\t %f\n", Bmatrix[1], Bmatrix[4], Bmatrix[7]);
	printf("%f\t %f\t %f\n", Bmatrix[2], Bmatrix[5], Bmatrix[8]);


	printf("\nDmatrix\n");
	printf("%f\t %f\t %f\n", Dmatrix[0], Dmatrix[3], Dmatrix[6]);
	printf("%f\t %f\t %f\n", Dmatrix[1], Dmatrix[4], Dmatrix[7]);
	printf("%f\t %f\t %f\n", Dmatrix[2], Dmatrix[5], Dmatrix[8]);

	printf("\nFmatrix\n");
	printf("%f\t %f\n", Fmatrix[0], Fmatrix[2]);
	printf("%f\t %f\n", Fmatrix[1], Fmatrix[3]);
	}*/

		if((GP==0)||(GP==3)){nx = -0.5773502692;}
		else{nx = 0.5773502692;}

		if((GP==0)||(GP==1)){ny = -0.5773502692;}
		else{ny = 0.5773502692;}

		/*printf("\n GP = %i, nx = %f, ny = %f", GP, nx, ny);*/


		n1j = 0.125*(1-nx)*(1-ny);
		n2j = 0.125*(1+nx)*(1-ny);
		n3j = 0.125*(1+nx)*(1+ny);
		n4j = 0.125*(1-nx)*(1+ny);

		dn1dx = -0.25*(1-ny)/hxz;
		dn2dx = 0.25*(1-ny)/hxz;
		dn3dx = 0.25*(1+ny)/hxz;
		dn4dx = -0.25*(1+ny)/hxz;

		dn1dy = -0.25*(1-nx)/h;
		dn2dy = -0.25*(1+nx)/h;
		dn3dy = 0.25*(1+nx)/h;
		dn4dy = 0.25*(1-nx)/h;

		/*printf("\n dn1dx = %f", dn1dx);
		printf("\n dn1dx = %f", dn2dx);
		printf("\n dn1dx = %f", dn3dx);
		printf("\n dn1dx = %f", dn4dx);

		printf("\n dn1dx = %f", dn1dy);
		printf("\n dn1dx = %f", dn2dy);
		printf("\n dn1dx = %f", dn3dy);

		printf("\n dn1dx = %f", dn4dy);


		printf("\n n1j = %f", n1j);
		printf("\n n2j = %f", n2j);
		printf("\n n3j = %f", n3j);

		printf("\n n4j = %f\n", n4j);*/

		BM[0] = dn1dx;
		BM[3] = 0;
		BM[6] = dn2dx;
		BM[9] = 0;
		BM[12] = dn3dx;
		BM[15] = 0;
		BM[18] = dn4dx;
		BM[21] = 0;
		
		BM[1] = 0;
		BM[4] = dn1dy;
		BM[7] = 0;
		BM[10] = dn2dy;
		BM[13] = 0;
		BM[16] = dn3dy;
		BM[19] = 0;
		BM[22] = dn4dy;

		BM[2] = dn1dy;
		BM[5] = dn1dx;
		BM[8] = dn2dy;
		BM[11] = dn2dx;
		BM[14] = dn3dy;
		BM[17] = dn3dx;
		BM[20] = dn4dy;
		BM[23] = dn4dx;

		BB[0] = 0;
		BB[3] = 0;
		BB[6] = -dn1dx;
		BB[9] = 0;
		BB[12] = 0;
		BB[15] = -dn2dx;
		BB[18] = 0;
		BB[21] = 0;
		BB[24] = -dn3dx;
		BB[27] = 0;
		BB[30] = 0;
		BB[33] = -dn4dx;

		BB[1] = 0;
		BB[4] = dn1dy;
		BB[7] = 0;
		BB[10] = 0;
		BB[13] = dn2dy;
		BB[16] = 0;
		BB[19] = 0;
		BB[22] = dn3dy;
		BB[25] = 0;
		BB[28] = 0;
		BB[31] = dn4dy;
		BB[34] = 0;

		BB[2] = 0;
		BB[5] = dn1dx;
		BB[8] = -dn1dy;
		BB[11] = 0;
		BB[14] = dn2dx;
		BB[17] = -dn2dy;
		BB[20] = 0;
		BB[23] = dn3dx;
		BB[26] = -dn3dy;
		BB[29] = 0;
		BB[32] = dn4dx;
		BB[35] = -dn4dy;

		BS[0] = dn1dx;
		BS[2] = 0;
		BS[4] = n1j;
		BS[6] = dn2dx;
		BS[8] = 0;
		BS[10] = n2j;
		BS[12] = dn3dx;
		BS[14] = 0;
		BS[16] = n3j;
		BS[18] = dn4dx;
		BS[20] = 0;
		BS[22] = n4j;

		BS[1] = dn1dy;
		BS[3] = -n1j;
		BS[5] = 0;
		BS[7] = dn2dy;
		BS[9] = -n2j;
		BS[11] = 0;
		BS[13] = dn3dy;
		BS[15] = -n3j;
		BS[17] = 0;
		BS[19] = dn4dy;
		BS[21] = -n4j;
		BS[23] = 0;

/*if((num == 361)||(num==378))

{
printf("BMshape\n");

printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[0], BM[3], BM[6], BM[9], BM[12], BM[15], BM[18], BM[21]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[1], BM[4], BM[7], BM[10], BM[13], BM[16], BM[19], BM[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[2], BM[5], BM[8], BM[11], BM[14], BM[17], BM[20], BM[23]);
}

/*printf("BBshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[0], BB[3], BB[6], BB[9], BB[12], BB[15], BB[18], BB[21], BB[24], BB[27], BB[30], BB[33]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[1], BB[4], BB[7], BB[10], BB[13], BB[16], BB[19], BB[22], BB[25], BB[28], BB[31], BB[34]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[2], BB[5], BB[8], BB[11], BB[14], BB[17], BB[20], BB[23], BB[26], BB[29], BB[32], BB[35]);

printf("BSshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[0], BS[2], BS[4], BS[6], BS[8], BS[10], BS[12], BS[14], BS[16], BS[18], BS[20], BS[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[1], BS[3], BS[5], BS[7], BS[9], BS[11], BS[13], BS[15], BS[17], BS[19], BS[21], BS[23]);*/



		/*Multiply the shape functions and the material matrices together using dgemm*/
/*printf("\nA");

printf("\nA");*/

/*Membrain terms*/		
Km[0] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[0] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[1] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[2]);
Km[1] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[0] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[1] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[2]);
Km[2] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[0] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[1] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[2]);
Km[3] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[0] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[1] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[2]);
Km[4] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[0] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[1] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[2]);
Km[5] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[0] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[1] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[2]);
Km[6] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[0] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[1] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[2]);
Km[7] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[0] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[1] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[2]);

Km[8] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[3] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[4] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[5]);
Km[9] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[3] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[4] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[5]);
Km[10] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[3] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[4] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[5]);
Km[11] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[3] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[4] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[5]);
Km[12] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[3] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[4] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[5]);
Km[13] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[3] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[4] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[5]);
Km[14] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[3] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[4] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[5]);
Km[15] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[3] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[4] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[5]);

Km[16] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[6] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[7] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[8]);
Km[17] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[6] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[7] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[8]);
Km[18] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[6] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[7] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[8]);
Km[19] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[6] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[7] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[8]);
Km[20] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[6] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[7] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[8]);
Km[21] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[6] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[7] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[8]);
Km[22] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[6] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[7] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[8]);
Km[23] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[6] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[7] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[8]);

Km[24] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[9] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[10] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[11]);
Km[25] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[9] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[10] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[11]);
Km[26] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[9] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[10] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[11]);
Km[27] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[9] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[10] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[11]);
Km[28] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[9] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[10] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[11]);
Km[29] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[9] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[10] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[11]);
Km[30] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[9] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[10] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[11]);
Km[31] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[9] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[10] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[11]);
		
Km[32] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[12] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[13] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[14]);
Km[33] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[12] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[13] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[14]);
Km[34] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[12] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[13] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[14]);
Km[35] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[12] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[13] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[14]);
Km[36] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[12] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[13] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[14]);
Km[37] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[12] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[13] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[14]);
Km[38] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[12] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[13] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[14]);
Km[39] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[12] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[13] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[14]);

Km[40] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[15] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[16] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[17]);
Km[41] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[15] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[16] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[17]);
Km[42] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[15] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[16] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[17]);
Km[43] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[15] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[16] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[17]);
Km[44] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[15] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[16] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[17]);
Km[45] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[15] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[16] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[17]);
Km[46] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[15] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[16] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[17]);
Km[47] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[15] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[16] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[17]);

Km[48] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[18] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[19] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[20]);
Km[49] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[18] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[19] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[20]);
Km[50] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[18] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[19] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[20]);
Km[51] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[18] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[19] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[20]);
Km[52] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[18] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[19] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[20]);
Km[53] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[18] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[19] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[20]);
Km[54] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[18] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[19] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[20]);
Km[55] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[18] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[19] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[20]);

Km[56] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[21] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[22] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[23]);
Km[57] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[21] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[22] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[23]);
Km[58] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[21] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[22] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[23]);
Km[59] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[21] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[22] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[23]);
Km[60] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[21] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[22] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[23]);
Km[61] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[21] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[22] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[23]);
Km[62] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[21] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[22] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[23]);
Km[63] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[21] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[22] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[23]);


/*Bending terms*/
Kp[0] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[0] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[1] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[2]);
Kp[1] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[0] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[1] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[2]);
Kp[2] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[0] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[1] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[2]);
Kp[3] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[0] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[1] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[2]);
Kp[4] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[0] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[1] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[2]);
Kp[5] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[0] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[1] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[2]);
Kp[6] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[0] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[1] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[2]);
Kp[7] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[0] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[1] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[2]);
Kp[8] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[0] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[1] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[2]);
Kp[9] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[0] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[1] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[2]);
Kp[10] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[0] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[1] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[2]);
Kp[11] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[0] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[1] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[2]);

Kp[12] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[3] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[4] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[5]);
Kp[13] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[3] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[4] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[5]);
Kp[14] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[3] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[4] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[5]);
Kp[15] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[3] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[4] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[5]);
Kp[16] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[3] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[4] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[5]);
Kp[17] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[3] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[4] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[5]);
Kp[18] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[3] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[4] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[5]);
Kp[19] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[3] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[4] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[5]);
Kp[20] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[3] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[4] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[5]);
Kp[21] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[3] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[4] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[5]);
Kp[22] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[3] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[4] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[5]);
Kp[23] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[3] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[4] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[5]);

Kp[24] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[6] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[7] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[8]);
Kp[25] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[6] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[7] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[8]);
Kp[26] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[6] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[7] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[8]);
Kp[27] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[6] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[7] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[8]);
Kp[28] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[6] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[7] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[8]);
Kp[29] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[6] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[7] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[8]);
Kp[30] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[6] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[7] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[8]);
Kp[31] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[6] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[7] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[8]);
Kp[32] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[6] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[7] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[8]);
Kp[33] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[6] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[7] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[8]);
Kp[34] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[6] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[7] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[8]);
Kp[35] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[6] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[7] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[8]);

Kp[36] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[9] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[10] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[11]);
Kp[37] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[9] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[10] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[11]);
Kp[38] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[9] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[10] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[11]);
Kp[39] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[9] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[10] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[11]);
Kp[40] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[9] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[10] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[11]);
Kp[41] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[9] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[10] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[11]);
Kp[42] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[9] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[10] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[11]);
Kp[43] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[9] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[10] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[11]);
Kp[44] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[9] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[10] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[11]);
Kp[45] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[9] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[10] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[11]);
Kp[46] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[9] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[10] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[11]);
Kp[47] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[9] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[10] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[11]);

Kp[48] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[12] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[13] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[14]);
Kp[49] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[12] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[13] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[14]);
Kp[50] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[12] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[13] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[14]);
Kp[51] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[12] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[13] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[14]);
Kp[52] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[12] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[13] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[14]);
Kp[53] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[12] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[13] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[14]);
Kp[54] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[12] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[13] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[14]);
Kp[55] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[12] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[13] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[14]);
Kp[56] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[12] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[13] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[14]);
Kp[57] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[12] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[13] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[14]);
Kp[58] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[12] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[13] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[14]);
Kp[59] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[12] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[13] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[14]);

Kp[60] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[15] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[16] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[17]);
Kp[61] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[15] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[16] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[17]);
Kp[62] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[15] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[16] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[17]);
Kp[63] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[15] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[16] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[17]);
Kp[64] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[15] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[16] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[17]);
Kp[65] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[15] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[16] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[17]);
Kp[66] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[15] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[16] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[17]);
Kp[67] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[15] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[16] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[17]);
Kp[68] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[15] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[16] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[17]);
Kp[69] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[15] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[16] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[17]);
Kp[70] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[15] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[16] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[17]);
Kp[71] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[15] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[16] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[17]);

Kp[72] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[18] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[19] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[20]);
Kp[73] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[18] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[19] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[20]);
Kp[74] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[18] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[19] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[20]);
Kp[75] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[18] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[19] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[20]);
Kp[76] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[18] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[19] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[20]);
Kp[77] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[18] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[19] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[20]);
Kp[78] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[18] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[19] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[20]);
Kp[79] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[18] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[19] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[20]);
Kp[80] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[18] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[19] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[20]);
Kp[81] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[18] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[19] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[20]);
Kp[82] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[18] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[19] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[20]);
Kp[83] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[18] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[19] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[20]);

Kp[84] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[21] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[22] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[23]);
Kp[85] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[21] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[22] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[23]);
Kp[86] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[21] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[22] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[23]);
Kp[87] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[21] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[22] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[23]);
Kp[88] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[21] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[22] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[23]);
Kp[89] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[21] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[22] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[23]);
Kp[90] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[21] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[22] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[23]);
Kp[91] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[21] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[22] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[23]);
Kp[92] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[21] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[22] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[23]);
Kp[93] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[21] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[22] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[23]);
Kp[94] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[21] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[22] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[23]);
Kp[95] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[21] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[22] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[23]);

Kp[96] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[24] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[25] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[26]);
Kp[97] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[24] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[25] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[26]);
Kp[98] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[24] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[25] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[26]);
Kp[99] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[24] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[25] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[26]);
Kp[100] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[24] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[25] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[26]);
Kp[101] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[24] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[25] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[26]);
Kp[102] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[24] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[25] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[26]);
Kp[103] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[24] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[25] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[26]);
Kp[104] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[24] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[25] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[26]);
Kp[105] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[24] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[25] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[26]);
Kp[106] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[24] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[25] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[26]);
Kp[107] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[24] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[25] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[26]);

Kp[108] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[27] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[28] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[29]);
Kp[109] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[27] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[28] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[29]);
Kp[110] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[27] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[28] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[29]);
Kp[111] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[27] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[28] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[29]);
Kp[112] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[27] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[28] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[29]);
Kp[113] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[27] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[28] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[29]);
Kp[114] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[27] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[28] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[29]);
Kp[115] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[27] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[28] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[29]);
Kp[116] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[27] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[28] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[29]);
Kp[117] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[27] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[28] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[29]);
Kp[118] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[27] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[28] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[29]);
Kp[119] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[27] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[28] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[29]);

Kp[120] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[30] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[31] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[32]);
Kp[121] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[30] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[31] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[32]);
Kp[122] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[30] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[31] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[32]);
Kp[123] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[30] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[31] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[32]);
Kp[124] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[30] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[31] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[32]);
Kp[125] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[30] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[31] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[32]);
Kp[126] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[30] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[31] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[32]);
Kp[127] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[30] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[31] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[32]);
Kp[128] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[30] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[31] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[32]);
Kp[129] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[30] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[31] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[32]);
Kp[130] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[30] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[31] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[32]);
Kp[131] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[30] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[31] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[32]);

Kp[132] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[33] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[34] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[35]);
Kp[133] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[33] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[34] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[35]);
Kp[134] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[33] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[34] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[35]);
Kp[135] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[33] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[34] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[35]);
Kp[136] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[33] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[34] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[35]);
Kp[137] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[33] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[34] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[35]);
Kp[138] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[33] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[34] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[35]);
Kp[139] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[33] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[34] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[35]);
Kp[140] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[33] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[34] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[35]);
Kp[141] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[33] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[34] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[35]);
Kp[142] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[33] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[34] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[35]);
Kp[143] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[33] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[34] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[35]);

/*Shear terms go through single guass intergration only*/
nx = 0;
ny = 0;

n1j = 0.25*(1-nx)*(1-ny);
n2j = 0.25*(1+nx)*(1-ny);
n3j = 0.25*(1+nx)*(1+ny);
n4j = 0.25*(1-nx)*(1+ny);

dn1dx = -0.5*(1-ny)/hxz;
dn2dx = 0.5*(1-ny)/hxz;
dn3dx = 0.5*(1+ny)/hxz;
dn4dx = -0.5*(1+ny)/hxz;

dn1dy = -0.5*(1-nx)/h;
dn2dy = -0.5*(1+nx)/h;
dn3dy = 0.5*(1+nx)/h;
dn4dy = 0.5*(1-nx)/h;


BS[0] = dn1dx;
BS[2] = 0;
BS[4] = n1j;
BS[6] = dn2dx;
BS[8] = 0;
BS[10] = n2j;
BS[12] = dn3dx;
BS[14] = 0;
BS[16] = n3j;
BS[18] = dn4dx;
BS[20] = 0;
BS[22] = n4j;

BS[1] = dn1dy;
BS[3] = -n1j;
BS[5] = 0;
BS[7] = dn2dy;
BS[9] = -n2j;
BS[11] = 0;
BS[13] = dn3dy;
BS[15] = -n3j;
BS[17] = 0;
BS[19] = dn4dy;
BS[21] = -n4j;
BS[23] = 0;

/*printf("BSshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[0], BS[2], BS[4], BS[6], BS[8], BS[10], BS[12], BS[14], BS[16], BS[18], BS[20], BS[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[1], BS[3], BS[5], BS[7], BS[9], BS[11], BS[13], BS[15], BS[17], BS[19], BS[21], BS[23]);*/

/*Shear terms*/
Kp[0] += ((BS[0]*F[0]+BS[1]*F[1])*BS[0] + (BS[0]*F[2]+BS[1]*F[3])*BS[1]);
Kp[1] += ((BS[2]*F[0]+BS[3]*F[1])*BS[0] + (BS[2]*F[2]+BS[3]*F[3])*BS[1]);
Kp[2] += ((BS[4]*F[0]+BS[5]*F[1])*BS[0] + (BS[4]*F[2]+BS[5]*F[3])*BS[1]);
Kp[3] += ((BS[6]*F[0]+BS[7]*F[1])*BS[0] + (BS[6]*F[2]+BS[7]*F[3])*BS[1]);
Kp[4] += ((BS[8]*F[0]+BS[9]*F[1])*BS[0] + (BS[8]*F[2]+BS[9]*F[3])*BS[1]);
Kp[5] += ((BS[10]*F[0]+BS[11]*F[1])*BS[0] + (BS[10]*F[2]+BS[11]*F[3])*BS[1]);
Kp[6] += ((BS[12]*F[0]+BS[13]*F[1])*BS[0] + (BS[12]*F[2]+BS[13]*F[3])*BS[1]);
Kp[7] += ((BS[14]*F[0]+BS[15]*F[1])*BS[0] + (BS[14]*F[2]+BS[15]*F[3])*BS[1]);
Kp[8] += ((BS[16]*F[0]+BS[17]*F[1])*BS[0] + (BS[16]*F[2]+BS[17]*F[3])*BS[1]);
Kp[9] += ((BS[18]*F[0]+BS[19]*F[1])*BS[0] + (BS[18]*F[2]+BS[19]*F[3])*BS[1]);
Kp[10] += ((BS[20]*F[0]+BS[21]*F[1])*BS[0] + (BS[20]*F[2]+BS[21]*F[3])*BS[1]);
Kp[11] += ((BS[22]*F[0]+BS[23]*F[1])*BS[0] + (BS[22]*F[2]+BS[23]*F[3])*BS[1]);

Kp[12] += ((BS[0]*F[0]+BS[1]*F[1])*BS[2] + (BS[0]*F[2]+BS[1]*F[3])*BS[3]);
Kp[13] += ((BS[2]*F[0]+BS[3]*F[1])*BS[2] + (BS[2]*F[2]+BS[3]*F[3])*BS[3]);
Kp[14] += ((BS[4]*F[0]+BS[5]*F[1])*BS[2] + (BS[4]*F[2]+BS[5]*F[3])*BS[3]);
Kp[15] += ((BS[6]*F[0]+BS[7]*F[1])*BS[2] + (BS[6]*F[2]+BS[7]*F[3])*BS[3]);
Kp[16] += ((BS[8]*F[0]+BS[9]*F[1])*BS[2] + (BS[8]*F[2]+BS[9]*F[3])*BS[3]);
Kp[17] += ((BS[10]*F[0]+BS[11]*F[1])*BS[2] + (BS[10]*F[2]+BS[11]*F[3])*BS[3]);
Kp[18] += ((BS[12]*F[0]+BS[13]*F[1])*BS[2] + (BS[12]*F[2]+BS[13]*F[3])*BS[3]);
Kp[19] += ((BS[14]*F[0]+BS[15]*F[1])*BS[2] + (BS[14]*F[2]+BS[15]*F[3])*BS[3]);
Kp[20] += ((BS[16]*F[0]+BS[17]*F[1])*BS[2] + (BS[16]*F[2]+BS[17]*F[3])*BS[3]);
Kp[21] += ((BS[18]*F[0]+BS[19]*F[1])*BS[2] + (BS[18]*F[2]+BS[19]*F[3])*BS[3]);
Kp[22] += ((BS[20]*F[0]+BS[21]*F[1])*BS[2] + (BS[20]*F[2]+BS[21]*F[3])*BS[3]);
Kp[23] += ((BS[22]*F[0]+BS[23]*F[1])*BS[2] + (BS[22]*F[2]+BS[23]*F[3])*BS[3]);

Kp[24] += ((BS[0]*F[0]+BS[1]*F[1])*BS[4] + (BS[0]*F[2]+BS[1]*F[3])*BS[5]);
Kp[25] += ((BS[2]*F[0]+BS[3]*F[1])*BS[4] + (BS[2]*F[2]+BS[3]*F[3])*BS[5]);
Kp[26] += ((BS[4]*F[0]+BS[5]*F[1])*BS[4] + (BS[4]*F[2]+BS[5]*F[3])*BS[5]);
Kp[27] += ((BS[6]*F[0]+BS[7]*F[1])*BS[4] + (BS[6]*F[2]+BS[7]*F[3])*BS[5]);
Kp[28] += ((BS[8]*F[0]+BS[9]*F[1])*BS[4] + (BS[8]*F[2]+BS[9]*F[3])*BS[5]);
Kp[29] += ((BS[10]*F[0]+BS[11]*F[1])*BS[4] + (BS[10]*F[2]+BS[11]*F[3])*BS[5]);
Kp[30] += ((BS[12]*F[0]+BS[13]*F[1])*BS[4] + (BS[12]*F[2]+BS[13]*F[3])*BS[5]);
Kp[31] += ((BS[14]*F[0]+BS[15]*F[1])*BS[4] + (BS[14]*F[2]+BS[15]*F[3])*BS[5]);

Kp[32] += ((BS[16]*F[0]+BS[17]*F[1])*BS[4] + (BS[16]*F[2]+BS[17]*F[3])*BS[5]);
Kp[33] += ((BS[18]*F[0]+BS[19]*F[1])*BS[4] + (BS[18]*F[2]+BS[19]*F[3])*BS[5]);
Kp[34] += ((BS[20]*F[0]+BS[21]*F[1])*BS[4] + (BS[20]*F[2]+BS[21]*F[3])*BS[5]);
Kp[35] += ((BS[22]*F[0]+BS[23]*F[1])*BS[4] + (BS[22]*F[2]+BS[23]*F[3])*BS[5]);

Kp[36] += ((BS[0]*F[0]+BS[1]*F[1])*BS[6] + (BS[0]*F[2]+BS[1]*F[3])*BS[7]);
Kp[37] += ((BS[2]*F[0]+BS[3]*F[1])*BS[6] + (BS[2]*F[2]+BS[3]*F[3])*BS[7]);
Kp[38] += ((BS[4]*F[0]+BS[5]*F[1])*BS[6] + (BS[4]*F[2]+BS[5]*F[3])*BS[7]);
Kp[39] += ((BS[6]*F[0]+BS[7]*F[1])*BS[6] + (BS[6]*F[2]+BS[7]*F[3])*BS[7]);
Kp[40] += ((BS[8]*F[0]+BS[9]*F[1])*BS[6] + (BS[8]*F[2]+BS[9]*F[3])*BS[7]);
Kp[41] += ((BS[10]*F[0]+BS[11]*F[1])*BS[6] + (BS[10]*F[2]+BS[11]*F[3])*BS[7]);
Kp[42] += ((BS[12]*F[0]+BS[13]*F[1])*BS[6] + (BS[12]*F[2]+BS[13]*F[3])*BS[7]);
Kp[43] += ((BS[14]*F[0]+BS[15]*F[1])*BS[6] + (BS[14]*F[2]+BS[15]*F[3])*BS[7]);
Kp[44] += ((BS[16]*F[0]+BS[17]*F[1])*BS[6] + (BS[16]*F[2]+BS[17]*F[3])*BS[7]);
Kp[45] += ((BS[18]*F[0]+BS[19]*F[1])*BS[6] + (BS[18]*F[2]+BS[19]*F[3])*BS[7]);
Kp[46] += ((BS[20]*F[0]+BS[21]*F[1])*BS[6] + (BS[20]*F[2]+BS[21]*F[3])*BS[7]);
Kp[47] += ((BS[22]*F[0]+BS[23]*F[1])*BS[6] + (BS[22]*F[2]+BS[23]*F[3])*BS[7]);

Kp[48] += ((BS[0]*F[0]+BS[1]*F[1])*BS[8] + (BS[0]*F[2]+BS[1]*F[3])*BS[9]);
Kp[49] += ((BS[2]*F[0]+BS[3]*F[1])*BS[8] + (BS[2]*F[2]+BS[3]*F[3])*BS[9]);
Kp[50] += ((BS[4]*F[0]+BS[5]*F[1])*BS[8] + (BS[4]*F[2]+BS[5]*F[3])*BS[9]);
Kp[51] += ((BS[6]*F[0]+BS[7]*F[1])*BS[8] + (BS[6]*F[2]+BS[7]*F[3])*BS[9]);
Kp[52] += ((BS[8]*F[0]+BS[9]*F[1])*BS[8] + (BS[8]*F[2]+BS[9]*F[3])*BS[9]);
Kp[53] += ((BS[10]*F[0]+BS[11]*F[1])*BS[8] + (BS[10]*F[2]+BS[11]*F[3])*BS[9]);
Kp[54] += ((BS[12]*F[0]+BS[13]*F[1])*BS[8] + (BS[12]*F[2]+BS[13]*F[3])*BS[9]);
Kp[55] += ((BS[14]*F[0]+BS[15]*F[1])*BS[8] + (BS[14]*F[2]+BS[15]*F[3])*BS[9]);
Kp[56] += ((BS[16]*F[0]+BS[17]*F[1])*BS[8] + (BS[16]*F[2]+BS[17]*F[3])*BS[9]);
Kp[57] += ((BS[18]*F[0]+BS[19]*F[1])*BS[8] + (BS[18]*F[2]+BS[19]*F[3])*BS[9]);
Kp[58] += ((BS[20]*F[0]+BS[21]*F[1])*BS[8] + (BS[20]*F[2]+BS[21]*F[3])*BS[9]);
Kp[59] += ((BS[22]*F[0]+BS[23]*F[1])*BS[8] + (BS[22]*F[2]+BS[23]*F[3])*BS[9]);

Kp[60] += ((BS[0]*F[0]+BS[1]*F[1])*BS[10] + (BS[0]*F[2]+BS[1]*F[3])*BS[11]);
Kp[61] += ((BS[2]*F[0]+BS[3]*F[1])*BS[10] + (BS[2]*F[2]+BS[3]*F[3])*BS[11]);
Kp[62] += ((BS[4]*F[0]+BS[5]*F[1])*BS[10] + (BS[4]*F[2]+BS[5]*F[3])*BS[11]);
Kp[63] += ((BS[6]*F[0]+BS[7]*F[1])*BS[10] + (BS[6]*F[2]+BS[7]*F[3])*BS[11]);
Kp[64] += ((BS[8]*F[0]+BS[9]*F[1])*BS[10] + (BS[8]*F[2]+BS[9]*F[3])*BS[11]);
Kp[65] += ((BS[10]*F[0]+BS[11]*F[1])*BS[10] + (BS[10]*F[2]+BS[11]*F[3])*BS[11]);
Kp[66] += ((BS[12]*F[0]+BS[13]*F[1])*BS[10] + (BS[12]*F[2]+BS[13]*F[3])*BS[11]);
Kp[67] += ((BS[14]*F[0]+BS[15]*F[1])*BS[10] + (BS[14]*F[2]+BS[15]*F[3])*BS[11]);
Kp[68] += ((BS[16]*F[0]+BS[17]*F[1])*BS[10] + (BS[16]*F[2]+BS[17]*F[3])*BS[11]);
Kp[69] += ((BS[18]*F[0]+BS[19]*F[1])*BS[10] + (BS[18]*F[2]+BS[19]*F[3])*BS[11]);
Kp[70] += ((BS[20]*F[0]+BS[21]*F[1])*BS[10] + (BS[20]*F[2]+BS[21]*F[3])*BS[11]);

Kp[71] += ((BS[22]*F[0]+BS[23]*F[1])*BS[10] + (BS[22]*F[2]+BS[23]*F[3])*BS[11]);

Kp[72] += ((BS[0]*F[0]+BS[1]*F[1])*BS[12] + (BS[0]*F[2]+BS[1]*F[3])*BS[13]);
Kp[73] += ((BS[2]*F[0]+BS[3]*F[1])*BS[12] + (BS[2]*F[2]+BS[3]*F[3])*BS[13]);
Kp[74] += ((BS[4]*F[0]+BS[5]*F[1])*BS[12] + (BS[4]*F[2]+BS[5]*F[3])*BS[13]);
Kp[75] += ((BS[6]*F[0]+BS[7]*F[1])*BS[12] + (BS[6]*F[2]+BS[7]*F[3])*BS[13]);
Kp[76] += ((BS[8]*F[0]+BS[9]*F[1])*BS[12] + (BS[8]*F[2]+BS[9]*F[3])*BS[13]);
Kp[77] += ((BS[10]*F[0]+BS[11]*F[1])*BS[12] + (BS[10]*F[2]+BS[11]*F[3])*BS[13]);
Kp[78] += ((BS[12]*F[0]+BS[13]*F[1])*BS[12] + (BS[12]*F[2]+BS[13]*F[3])*BS[13]);
Kp[79] += ((BS[14]*F[0]+BS[15]*F[1])*BS[12] + (BS[14]*F[2]+BS[15]*F[3])*BS[13]);
Kp[80] += ((BS[16]*F[0]+BS[17]*F[1])*BS[12] + (BS[16]*F[2]+BS[17]*F[3])*BS[13]);
Kp[81] += ((BS[18]*F[0]+BS[19]*F[1])*BS[12] + (BS[18]*F[2]+BS[19]*F[3])*BS[13]);
Kp[82] += ((BS[20]*F[0]+BS[21]*F[1])*BS[12] + (BS[20]*F[2]+BS[21]*F[3])*BS[13]);
Kp[83] += ((BS[22]*F[0]+BS[23]*F[1])*BS[12] + (BS[22]*F[2]+BS[23]*F[3])*BS[13]);

Kp[84] += ((BS[0]*F[0]+BS[1]*F[1])*BS[14] + (BS[0]*F[2]+BS[1]*F[3])*BS[15]);
Kp[85] += ((BS[2]*F[0]+BS[3]*F[1])*BS[14] + (BS[2]*F[2]+BS[3]*F[3])*BS[15]);
Kp[86] += ((BS[4]*F[0]+BS[5]*F[1])*BS[14] + (BS[4]*F[2]+BS[5]*F[3])*BS[15]);
Kp[87] += ((BS[6]*F[0]+BS[7]*F[1])*BS[14] + (BS[6]*F[2]+BS[7]*F[3])*BS[15]);
Kp[88] += ((BS[8]*F[0]+BS[9]*F[1])*BS[14] + (BS[8]*F[2]+BS[9]*F[3])*BS[15]);
Kp[89] += ((BS[10]*F[0]+BS[11]*F[1])*BS[14] + (BS[10]*F[2]+BS[11]*F[3])*BS[15]);
Kp[90] += ((BS[12]*F[0]+BS[13]*F[1])*BS[14] + (BS[12]*F[2]+BS[13]*F[3])*BS[15]);
Kp[91] += ((BS[14]*F[0]+BS[15]*F[1])*BS[14] + (BS[14]*F[2]+BS[15]*F[3])*BS[15]);
Kp[92] += ((BS[16]*F[0]+BS[17]*F[1])*BS[14] + (BS[16]*F[2]+BS[17]*F[3])*BS[15]);
Kp[93] += ((BS[18]*F[0]+BS[19]*F[1])*BS[14] + (BS[18]*F[2]+BS[19]*F[3])*BS[15]);
Kp[94] += ((BS[20]*F[0]+BS[21]*F[1])*BS[14] + (BS[20]*F[2]+BS[21]*F[3])*BS[15]);
Kp[95] += ((BS[22]*F[0]+BS[23]*F[1])*BS[14] + (BS[22]*F[2]+BS[23]*F[3])*BS[15]);

Kp[96] += ((BS[0]*F[0]+BS[1]*F[1])*BS[16] + (BS[0]*F[2]+BS[1]*F[3])*BS[17]);
Kp[97] += ((BS[2]*F[0]+BS[3]*F[1])*BS[16] + (BS[2]*F[2]+BS[3]*F[3])*BS[17]);
Kp[98] += ((BS[4]*F[0]+BS[5]*F[1])*BS[16] + (BS[4]*F[2]+BS[5]*F[3])*BS[17]);
Kp[99] += ((BS[6]*F[0]+BS[7]*F[1])*BS[16] + (BS[6]*F[2]+BS[7]*F[3])*BS[17]);
Kp[100] += ((BS[8]*F[0]+BS[9]*F[1])*BS[16] + (BS[8]*F[2]+BS[9]*F[3])*BS[17]);
Kp[101] += ((BS[10]*F[0]+BS[11]*F[1])*BS[16] + (BS[10]*F[2]+BS[11]*F[3])*BS[17]);
Kp[102] += ((BS[12]*F[0]+BS[13]*F[1])*BS[16] + (BS[12]*F[2]+BS[13]*F[3])*BS[17]);
Kp[103] += ((BS[14]*F[0]+BS[15]*F[1])*BS[16] + (BS[14]*F[2]+BS[15]*F[3])*BS[17]);
Kp[104] += ((BS[16]*F[0]+BS[17]*F[1])*BS[16] + (BS[16]*F[2]+BS[17]*F[3])*BS[17]);
Kp[105] += ((BS[18]*F[0]+BS[19]*F[1])*BS[16] + (BS[18]*F[2]+BS[19]*F[3])*BS[17]);
Kp[106] += ((BS[20]*F[0]+BS[21]*F[1])*BS[16] + (BS[20]*F[2]+BS[21]*F[3])*BS[17]);
Kp[107] += ((BS[22]*F[0]+BS[23]*F[1])*BS[16] + (BS[22]*F[2]+BS[23]*F[3])*BS[17]);

Kp[108] += ((BS[0]*F[0]+BS[1]*F[1])*BS[18] + (BS[0]*F[2]+BS[1]*F[3])*BS[19]);
Kp[109] += ((BS[2]*F[0]+BS[3]*F[1])*BS[18] + (BS[2]*F[2]+BS[3]*F[3])*BS[19]);
Kp[110] += ((BS[4]*F[0]+BS[5]*F[1])*BS[18] + (BS[4]*F[2]+BS[5]*F[3])*BS[19]);
Kp[111] += ((BS[6]*F[0]+BS[7]*F[1])*BS[18] + (BS[6]*F[2]+BS[7]*F[3])*BS[19]);
Kp[112] += ((BS[8]*F[0]+BS[9]*F[1])*BS[18] + (BS[8]*F[2]+BS[9]*F[3])*BS[19]);
Kp[113] += ((BS[10]*F[0]+BS[11]*F[1])*BS[18] + (BS[10]*F[2]+BS[11]*F[3])*BS[19]);
Kp[114] += ((BS[12]*F[0]+BS[13]*F[1])*BS[18] + (BS[12]*F[2]+BS[13]*F[3])*BS[19]);
Kp[115] += ((BS[14]*F[0]+BS[15]*F[1])*BS[18] + (BS[14]*F[2]+BS[15]*F[3])*BS[19]);
Kp[116] += ((BS[16]*F[0]+BS[17]*F[1])*BS[18] + (BS[16]*F[2]+BS[17]*F[3])*BS[19]);
Kp[117] += ((BS[18]*F[0]+BS[19]*F[1])*BS[18] + (BS[18]*F[2]+BS[19]*F[3])*BS[19]);
Kp[118] += ((BS[20]*F[0]+BS[21]*F[1])*BS[18] + (BS[20]*F[2]+BS[21]*F[3])*BS[19]);
Kp[119] += ((BS[22]*F[0]+BS[23]*F[1])*BS[18] + (BS[22]*F[2]+BS[23]*F[3])*BS[19]);

Kp[120] += ((BS[0]*F[0]+BS[1]*F[1])*BS[20] + (BS[0]*F[2]+BS[1]*F[3])*BS[21]);
Kp[121] += ((BS[2]*F[0]+BS[3]*F[1])*BS[20] + (BS[2]*F[2]+BS[3]*F[3])*BS[21]);
Kp[122] += ((BS[4]*F[0]+BS[5]*F[1])*BS[20] + (BS[4]*F[2]+BS[5]*F[3])*BS[21]);
Kp[123] += ((BS[6]*F[0]+BS[7]*F[1])*BS[20] + (BS[6]*F[2]+BS[7]*F[3])*BS[21]);
Kp[124] += ((BS[8]*F[0]+BS[9]*F[1])*BS[20] + (BS[8]*F[2]+BS[9]*F[3])*BS[21]);
Kp[125] += ((BS[10]*F[0]+BS[11]*F[1])*BS[20] + (BS[10]*F[2]+BS[11]*F[3])*BS[21]);
Kp[126] += ((BS[12]*F[0]+BS[13]*F[1])*BS[20] + (BS[12]*F[2]+BS[13]*F[3])*BS[21]);
Kp[127] += ((BS[14]*F[0]+BS[15]*F[1])*BS[20] + (BS[14]*F[2]+BS[15]*F[3])*BS[21]);
Kp[128] += ((BS[16]*F[0]+BS[17]*F[1])*BS[20] + (BS[16]*F[2]+BS[17]*F[3])*BS[21]);
Kp[129] += ((BS[18]*F[0]+BS[19]*F[1])*BS[20] + (BS[18]*F[2]+BS[19]*F[3])*BS[21]);
Kp[130] += ((BS[20]*F[0]+BS[21]*F[1])*BS[20] + (BS[20]*F[2]+BS[21]*F[3])*BS[21]);
Kp[131] += ((BS[22]*F[0]+BS[23]*F[1])*BS[20] + (BS[22]*F[2]+BS[23]*F[3])*BS[21]);

Kp[132] += ((BS[0]*F[0]+BS[1]*F[1])*BS[22] + (BS[0]*F[2]+BS[1]*F[3])*BS[23]);
Kp[133] += ((BS[2]*F[0]+BS[3]*F[1])*BS[22] + (BS[2]*F[2]+BS[3]*F[3])*BS[23]);
Kp[134] += ((BS[4]*F[0]+BS[5]*F[1])*BS[22] + (BS[4]*F[2]+BS[5]*F[3])*BS[23]);
Kp[135] += ((BS[6]*F[0]+BS[7]*F[1])*BS[22] + (BS[6]*F[2]+BS[7]*F[3])*BS[23]);
Kp[136] += ((BS[8]*F[0]+BS[9]*F[1])*BS[22] + (BS[8]*F[2]+BS[9]*F[3])*BS[23]);
Kp[137] += ((BS[10]*F[0]+BS[11]*F[1])*BS[22] + (BS[10]*F[2]+BS[11]*F[3])*BS[23]);
Kp[138] += ((BS[12]*F[0]+BS[13]*F[1])*BS[22] + (BS[12]*F[2]+BS[13]*F[3])*BS[23]);
Kp[139] += ((BS[14]*F[0]+BS[15]*F[1])*BS[22] + (BS[14]*F[2]+BS[15]*F[3])*BS[23]);
Kp[140] += ((BS[16]*F[0]+BS[17]*F[1])*BS[22] + (BS[16]*F[2]+BS[17]*F[3])*BS[23]);
Kp[141] += ((BS[18]*F[0]+BS[19]*F[1])*BS[22] + (BS[18]*F[2]+BS[19]*F[3])*BS[23]);
Kp[142] += ((BS[20]*F[0]+BS[21]*F[1])*BS[22] + (BS[20]*F[2]+BS[21]*F[3])*BS[23]);
Kp[143] += ((BS[22]*F[0]+BS[23]*F[1])*BS[22] + (BS[22]*F[2]+BS[23]*F[3])*BS[23]);

/*if((num == 361)||(num==378))
{
printf("\n Membrain Matrix");
for(p=0;p<8;p++)
{
	printf("\n");
	for(j=0;j<8;j++)
	{
		printf("%0.3f\t", Km[j*8+p]);
	}
}
printf("\n");
}*/

/*Check the reuslts so far*/
/*printf("\n Membrain Matrix");
for(p=0;p<8;p++)
{
	printf("\n");
	for(j=0;j<8;j++)
	{
		printf("%0.3f\t", Km[j*8+p]);
	}
}
printf("\n");
printf("\n Plate Matrix");
for(i=0;i<12;i++)
{
	printf("\n");
	for(j=0;j<12;j++)
	{
		printf("%0.3f\t", Kp[j*12+i]);
	}
}*/

/*Assemble Full Matirx (Try to transpose it so that it is correct in C matrix form by not hugely important since it should be symetrical*/
for(i=0;i<576;i++)
{
	KE[i] = 0.0000000000000000;	/*Intialize every value to zero*/
}

/*Fill in the values*/
for(j=0;j<4;j++)
{
	p =24*6*j;
	KE[p] = Km[0+2*j];
	KE[1+p] = Km[8+2*j];
	KE[6+p] = Km[16+2*j];
	KE[7+p] = Km[24+2*j];
	KE[12+p] = Km[32+2*j];
	KE[13+p] = Km[40+2*j];
	KE[18+p] = Km[48+2*j];
	KE[19+p] = Km[56+2*j];
	KE[24+p] = Km[1+2*j];

	KE[25+p] = Km[9+2*j];
	KE[30+p] = Km[17+2*j];
	KE[31+p] = Km[25+2*j];
	KE[36+p] = Km[33+2*j];
	KE[37+p] = Km[41+2*j];
	KE[42+p] = Km[49+2*j];
	KE[43+p] = Km[57+2*j];

	q =24*6*j+24*2;
	KE[2+q] = Kp[0+3*j];
	KE[3+q] = Kp[12+3*j];
	KE[4+q] = Kp[24+3*j];
	KE[8+q] = Kp[36+3*j];
	KE[9+q] = Kp[48+3*j];
	KE[10+q] = Kp[60+3*j];
	KE[14+q] = Kp[72+3*j];
	KE[15+q] = Kp[84+3*j];
	KE[16+q] = Kp[96+3*j];
	KE[20+q] = Kp[108+3*j];

	KE[21+q] = Kp[120+3*j];
	KE[22+q] = Kp[132+3*j];
	q =24*6*j+24*3;
	KE[2+q] = Kp[1+3*j];
	KE[3+q] = Kp[13+3*j];
	KE[4+q] = Kp[25+3*j];
	KE[8+q] = Kp[37+3*j];
	KE[9+q] = Kp[49+3*j];
	KE[10+q] = Kp[61+3*j];
	KE[14+q] = Kp[73+3*j];
	KE[15+q] = Kp[85+3*j];
	KE[16+q] = Kp[97+3*j];
	KE[20+q] = Kp[109+3*j];
	KE[21+q] = Kp[121+3*j];
	KE[22+q] = Kp[133+3*j];
	q =24*6*j+24*4;
	KE[2+q] = Kp[2+3*j];
	KE[3+q] = Kp[14+3*j];
	KE[4+q] = Kp[26+3*j];
	KE[8+q] = Kp[38+3*j];
	KE[9+q] = Kp[50+3*j];
	KE[10+q] = Kp[62+3*j];
	KE[14+q] = Kp[74+3*j];
	KE[15+q] = Kp[86+3*j];
	KE[16+q] = Kp[98+3*j];
	KE[20+q] = Kp[110+3*j];
	KE[21+q] = Kp[122+3*j];
	KE[22+q] = Kp[134+3*j];

}

/*printf("\n Element Matrix");
for(i=0;i<24;i++)
{
	printf("\n");
	for(j=0;j<24;j++)
	{
		printf("%0.2f  ", KE[j+i*24]);

	}
}*/


free(A);
free(D);
free(B);
free(F);
free(Cmatrix);
free(C0matrix);
free(Q0matrix);
free(Qmatrix);
free(BB);
free(BM);
free(BS);
free(tempKm);
free(tempKp);
free(tempKps);
free(Km);
free(Kp);
/*printf("\n Temp Membrain Matrix");*/
/*for(p=0;p<12;p++)
{
	printf("\n");
	for(j=0;j<3;j++)
	{
		printf("%0.3f\t", tempKp[j*12+p]);
	}
}*/
/*printf("\nB");
printf("\nB");	*/	

}
/*Calculates the tow path by setting it parrell to the part of the interger values of the implicit function*/
void ThetaUpdate(int elemX, int elemY,int NumNodes, Elem Number[elemX][elemY],Coord NodeCoord[NumNodes], double *theta, double *lsf, int tow, double *hxz, double h, double maxXZ)
{

	int Numtows, k, n, m, i, j, Tsum, minT, num, count;
	double mintow, lsfmin, dtemp, Ttow,ftemp,lsf1,lsf2;

	Numtows = 4*(elemY+elemX); /*Make sure it is large enough to get every element's tow, safey margine for the effects of a curve*/
	int *TNodeStat, *Tnodes, *done;
	double *Tlsf, *localX, *localY;
	TNodeStat = malloc(4*sizeof(int));
	Tnodes = malloc(4*sizeof(int));
	Tlsf = malloc(4*sizeof(double));
	done = malloc(elemX*elemY*sizeof(int));
	localX = malloc(4*sizeof(double));
	localY = malloc(4*sizeof(double));
	for(i=0;i<(elemX*elemY);i++)
	{
		done[i] = 0;
	}
	
	mintow = -1*(Numtows/2);


	for(k=0;k<Numtows;k++)
	{
		Ttow = mintow + (k * tow);
		for(n=0;n<elemX;n++)
		{
			for(m=0;m<elemY;m++)
			{

				/*See if this tow crosses this element*/
				num = Number[n][m].n-1;
				if(done[num]!=1)
				{
				Tnodes[0] = Number[n][m].a-1;
				Tnodes[1] = Number[n][m].b-1;
				Tnodes[2] = Number[n][m].c-1;
				Tnodes[3] = Number[n][m].d-1;
				Tsum = 0;				

				for(i=0;i<4;i++)
				{
					dtemp = lsf[Tnodes[i]]-Ttow;

					if(dtemp<0.00000000001){TNodeStat[i] = 0;}
					else{TNodeStat[i] = 1;}
					/*printf("\n lsf[%i] - Tow = %f - %f = %f", Tnodes[i], lsf[Tnodes[i]], Ttow, dtemp);*/

					Tsum += TNodeStat[i];
				}

				if((Tsum!=0)&&(Tsum!=4))	/*Then the element is cut by this Tow*/
				{
					Tlsf[0] = lsf[Tnodes[0]] - Ttow;
					Tlsf[1] = lsf[Tnodes[1]] - Ttow;
					Tlsf[2] = lsf[Tnodes[2]] - Ttow;
					Tlsf[3] = lsf[Tnodes[3]] - Ttow;

					count = 0;

					for(i=0;i<4;i++) /*Check for double cuts*/
					{
						j = (i == 3) ? 0 : (i + 1); 
						if(((TNodeStat[i]*TNodeStat[j]) == 0)&&((TNodeStat[i] + TNodeStat[j]) !=0)){count++;}
					}
					/*if((num==400)||(num==4)){printf("\ncount = %i", count);}*/

					if(count>2)
					{
						lsfmin = fabs(Tlsf[0]);
						minT = 0;
						if(fabs(Tlsf[1])<lsfmin){lsfmin = fabs(Tlsf[1]); minT = 1;}
						if(fabs(Tlsf[2])<lsfmin){lsfmin = fabs(Tlsf[2]); minT = 2;}
						if(fabs(Tlsf[3])<lsfmin){lsfmin = fabs(Tlsf[3]); minT = 3;}

						TNodeStat[minT] = 2;
					}


					/*Get the location of all the cut edge points, if element is cut twice by this tow deploy the corection proccedure*/
					count = 0;
					for(i=0;i<4;i++)
					{
						j = (i == 3) ? 0 : (i + 1); 
						/*if((num==1)||(num==389)){
							printf("\n I = %i, J = %i", i, j);
							printf("\n TnodeStat[0] = %i,   Tlsf[0] = %f", TNodeStat[0], Tlsf[0]);
							printf("\n TnodeStat[1] = %i,   Tlsf[1] = %f", TNodeStat[1], Tlsf[1]);
							printf("\n TnodeStat[2] = %i,   Tlsf[2] = %f", TNodeStat[2], Tlsf[2]);
							printf("\n TnodeStat[3] = %i,   Tlsf[3] = %f", TNodeStat[3], Tlsf[3]);
						}*/

						if(((TNodeStat[i]*TNodeStat[j]) == 0)&&((TNodeStat[i] + TNodeStat[j]) !=2)&&((TNodeStat[i] + TNodeStat[j]) !=0))	/*edge is cut*/
						{
							lsf1 = Tlsf[i];
							lsf2 = Tlsf[j];
							if(fabs(lsf1) < 0.0001){lsf1 = 0.0001;}
							if(fabs(lsf2) < 0.0001){lsf2 = 0.0001;}	/*Prevents both local points from being placed exactly at one node*/

							ftemp = fabs(lsf1 / (lsf1 - lsf2));
							/*if((num==1)||(num==389)){printf("\nlsf1 = %f, lsf2 = %f", lsf1, lsf2);}*/
							/*if((num==1)||(num==389)){printf("\ni = %i, ftemp = %f", i, ftemp);}
							if((num==1)||(num==389)){printf("\n lsf1-lsf2 = %f", lsf1-lsf2);}*/

							if((fabs(lsf1-lsf2))<0.00001)/*System to deal with elments whoes edges line up exactly with the */
							{
		/*printf("\nIn check");*/
		if(i==0){localX[0] = NodeCoord[Tnodes[0]].x;localX[1] = NodeCoord[Tnodes[1]].x;localY[0] = NodeCoord[Tnodes[0]].y;localY[1] = NodeCoord[Tnodes[1]].y;}
		if(i==1){localX[0] = NodeCoord[Tnodes[1]].x;localX[1] = NodeCoord[Tnodes[2]].x;localY[0] = NodeCoord[Tnodes[1]].y;localY[1] = NodeCoord[Tnodes[2]].y;}
		if(i==2){localX[0] = NodeCoord[Tnodes[2]].x;localX[1] = NodeCoord[Tnodes[3]].x;localY[0] = NodeCoord[Tnodes[2]].y;localY[1] = NodeCoord[Tnodes[3]].y;}
		if(i==3){localX[0] = NodeCoord[Tnodes[3]].x;localX[1] = NodeCoord[Tnodes[0]].x;localY[0] = NodeCoord[Tnodes[3]].y;localY[1] = NodeCoord[Tnodes[0]].y;}
		count = 2;
							}
							else
							{
								if(i==0) {localX[count] = NodeCoord[Tnodes[0]].x + ftemp*hxz[n]; localY[count] = NodeCoord[Tnodes[0]].y;}
								if(i==1) {localX[count] = NodeCoord[Tnodes[1]].x; localY[count] = NodeCoord[Tnodes[1]].y + ftemp*h;;}
								if(i==2) {localX[count] = NodeCoord[Tnodes[2]].x - ftemp*hxz[n]; localY[count] = NodeCoord[Tnodes[2]].y;}
								if(i==3) {localX[count] = NodeCoord[Tnodes[3]].x; localY[count] = NodeCoord[Tnodes[3]].y - ftemp*h;}
								count++;
							}
						}					
					}
					/*if((num==1)||(num==389)){printf("\nnum:%i Ttow: %f localX: %f, %f \t localY: %f, %f",num, Ttow, localX[1],localX[0],localY[1],localY[0]);}*/

					/*Update Theta (one ply only for now)*/
					dtemp =  (localY[1] - localY[0])/(localX[1] - localX[0]);
					/*if((num==1)||(num==389)){printf("\ndtemp: %f",dtemp);}*/
					theta[num] = atan(dtemp);
					if(theta[num] > 3.14159265){theta[num] -= 3.14159265;}
					else if(theta[num] < 0){theta[num] += 3.14159265;}
					/*if((num==1)||(num==389)){printf("\ttheta: %f",theta[num]);}*/
					done[num] = 1;
				}
				}
			}
		}
	}

	free(TNodeStat);
	free(Tnodes);
	free(Tlsf);
	free(done);
	free(localX);
	free(localY);

}
/*Updates the Theta values by intergrating the implicit funtion across the element to find the direction of the maximum implict function gradient. The tow is laided normal to this 
maximum gradient. Uses a Gradweno scheme to calcualte the graident*/
void ThetaUpdate2(int *Lnodes, Coord *LCoord, int num, double *lsf, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], double h, int elemX, double *hxz, double *theta)
{
	double v1,v2,v3,v4,v5; /*varibales for finite differences used*/
	double grad;	/*variable for gradient approx*/
	double *gradX, *gradY; /*Pointer to store the X and Y gradients at each node.*/
	gradX = malloc(4*sizeof(double));
	gradY = malloc(4*sizeof(double));
	double XGE, YGE;	/*mean gradient of the element in both directions*/
	double ftemp, dtemp;	/*tempory variable*/
	double htemp;
	int flag = 0;
	double phi;		/*Current Node lsf value*/
	int Xi, Yj,i;		/*Current Node Cords*/

	for(i=0;i<4;i++)
	{
		phi = lsf[Lnodes[i]];
		Xi =  LCoord[i].x;
		Yj =  LCoord[i].y;

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

		double Xr = GWsub(v1, v2, v3, v4, v5); /*calcualte gradient, info going to the right*/

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

		double Xl = GWsub(v1, v2, v3, v4, v5); /*calcualte gradient info going to the left*/
		/*if((num==19)||(num==59))
		{
			printf("\n\n %i, %i", num, Xi);
			printf("\n Xr = %f\tXl = %f", Xr, Xl);
		}*/

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

		double Yu = GWsub(v1, v2, v3, v4, v5); /*calcualte gradient info going up*/

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

		double Yd = GWsub(v1, v2, v3, v4, v5); /*calcualte gradient info going down*/

		/*now calculate gradinet in each direction using upwind scheme!!!Old Wrong Method Kept around just in case*/
		/*Get X Gradient Magnitude*
		ftemp = Xl*Xl;
		gradX[i] = ftemp;
		ftemp = Xr*Xr;
		gradX[i] += ftemp;
		ftemp = gradX[i]/2;
		gradX[i] = sqrt(ftemp); /*need to square root answer*
		/*Now get direction correct*
		if(((Xl==0)&&(Xr==0))||((Xl-Xr) == 0)){  gradX[i] = 0;}
		else if((Xr<-0.000000001)&&(Xl<-0.000000001)){ gradX[i] = -1*gradX[i];}
		else if((Xr<-0.000000001)||(Xl<-0.000000001)){ gradX[i] = (Xr+Xl<0.0000000001) ? (-1*gradX[i]):(gradX[i]);}
		else {gradX[i] = gradX[i];}*/

		/*Get X Gradient from average of both estimates*/
		ftemp = Xr+Xl; 
		gradX[i] = ftemp/2;
		
		/*Get Y Gradient from average of both estimates*/
		ftemp = Yu+Yd; 
		gradY[i] = ftemp/2;


	}

	/*Now we have got the gradient over each of the nodes we get the elment gradient as the mean value*/
	XGE = (gradX[0] + gradX[1] + gradX[2] + gradX[3])/4;
	YGE = (gradY[0] + gradY[1] + gradY[2] + gradY[3])/4;
	/*if((num==19)||(num==59))
	{
		printf("\n\n %i", num);
		printf("\n GRADX:(%f, %f, %f, %f)", gradX[0], gradX[1], gradX[2], gradX[3]);
		printf("\n GRADY:(%f, %f, %f, %f)", gradY[0], gradY[1], gradY[2], gradY[3]);
		printf("\n XGE: %f", XGE);
		printf("\n YGE: %f", YGE);
	}*/

	/*Now we use these to find the angel of the maximum gradient direction*/
	if(XGE==0) {ftemp = 1.570796327;} /*Set theta to 90 since tan90 is infinate*/
	else if (YGE==0) {ftemp = 0;}
	else
	{
		dtemp = YGE/XGE;
		ftemp = atan(dtemp);
	}

	/*Up date theta to be the normal to this angle*/
	theta[num] = ftemp + 1.570796327;

	/*Tidy up theta to be between 0 and 180 degrees*/
	if(theta[num] > 3.14158265){theta[num] -= 3.14158265;}
	else if(theta[num] < 0){theta[num] += 3.14158265;}

	/*if((num==19)||(num==59))
	{
		printf("\ntheta[%i] = %f", num, theta[num]);
		printf("\t = %f", theta[num]*(180/3.14158265));
	}*/

	/*Free the tempoary memory functions*/
	free(gradX);
	free(gradY);
}
/*Function that uses shape functions to calculate the fibre orientation*/
void ThetaUpdate3(double *Theta, int elemX, int elemY, Elem Number[elemX][elemY], double *lsf, double h, double *hxz)
{
	int i,j,k,num;							/*interger varribles*/
	int *Lnodes;
	double *Tlsf;				/*Local node data*/
	double *Nx, *Ny;					/*Differnetial of shape functions*/
	double s, n, a, b, X, Y;			/*local shape funtion co-ordinate values*/
	double dx, dy;						/*Differnetial of lsf wrt to x and y*/
	double ftemp, dtemp;			/*tempoary varibles*/

	b = h;
	
	Lnodes = malloc(4*sizeof(int));
	Tlsf = malloc(4*sizeof(double));
	Nx = malloc(4*sizeof(double));
	Ny = malloc(4*sizeof(double));

	for(i=0;i<elemX;i++)
	{
		a = hxz[i];
		for(j=0;j<elemY;j++)
		{
			/*Get the local Node data*/
			num = Number[i][j].n - 1;
			Lnodes[0] = Number[i][j].a - 1;
			Lnodes[1] = Number[i][j].b - 1;
			Lnodes[2] = Number[i][j].c - 1;
			Lnodes[3] = Number[i][j].d - 1;

			Tlsf[0] = lsf[Lnodes[0]];
			Tlsf[1] = lsf[Lnodes[1]];
			Tlsf[2] = lsf[Lnodes[2]];
			Tlsf[3] = lsf[Lnodes[3]];

			/*Shape function local cordinates, for this problem a single gauss point at the centre of the element will be adequate*/
			s = 0.0;
			n = 0.0;

			X = ((1)/(4*a*a*b));
			Y = ((1)/(4*a*b*b));
			/*Shape function derivates*/
			Nx[0] = X*(-1+n);
			Nx[1] = X*(1-n);
			Nx[2] = X*(1+n);
			Nx[3] = X*(-1-n);

			Ny[0] = Y*(s-1);
			Ny[1] = Y*(-s-1);
			Ny[2] = Y*(s+1);
			Ny[3] = Y*(-s+1);
			
			dx = 0.0;
			dy = 0.0;

			/*Now use the shape function derivaties to calcualte x and y lsf gradients over the element*/
			for(k=0;k<4;k++)
			{
				dx += Nx[k]*Tlsf[k];
				dy += Ny[k]*Tlsf[k];
			}
			/*Now we have the gradients we need to get the angle theta from it*/

			if(fabs(dx)<0.00000001) {ftemp = 1.570796327;} /*Set theta to 90 since tan90 is infinate*/
			else if (fabs(dy)<0.00000001) {ftemp = 0.0;}
			else
			{
				dtemp = dy/dx;
				ftemp = atan(dtemp);
			}

			/*Up date theta to be the normal to this angle*/
			Theta[num] = ftemp + 1.570796327;
            
            
            
           
            //  This part is to verify the sensitivity
            double dT;
            dT=0.0000001;
            if (num==80)
            {
               // Theta[num]+=dT;
                
                // C = cos(Theta[num]+dT);
                //printf("C= %f", C);
                //S = sin(Theta[num]+dT);
            }
            
			/*Tidy up theta to be between 0 and 180 degrees*/
			if(Theta[num] > 1.570796327){Theta[num] -= 3.14159265;}
			else if(Theta[num] < -1.570796327){Theta[num] += 3.14159265;}
            printf("\n element %d theta is equal to %f", num, Theta[num] );
			/*if((num==19)||(num==400))
			{
				printf("\n\n %i", num);
				printf("\n lsfX:(%f, %f, %f, %f)", Tlsf[0], Tlsf[1], Tlsf[2], Tlsf[3]);
				printf("\n dx: %f", dx);
				printf("\n dy: %f", dy);
				printf("\n ftemp: %f", ftemp);
				printf("\n Theta: %f", Theta[num]);
			}*/

		}
	}
	free(Lnodes);
	free(Tlsf);
	free(Nx);
	free(Ny);
}

/*Function that uses shape functions to calculate the fibre orientation*/
void ThetaUpdate4(double *Theta, int elemX, int elemY, Elem Number[elemX][elemY], double *lsf, double *lsfb, int Numlsf, double h, double *hxz)
{
	int i,j,k,o,p1,p2,num,p,q;		/*interger varribles*/
	int *Lnodes;				/*Local Nodes*/
	double *Tlsf,*Tlsfb;				/*Local node lsf data*/
	double lsf_min, lsf_max;		/*max and min lsf in an elemenet varibles*/
	int minN, maxN;				/*max and min lsf number*/
	double ftemp, ftemp1, ftemp2, dtemp;	/*tempoary varibles*/
	int Scount,Gcount;			/*Count varibles*/
	double x,x1,x2,y,y1,y2,xa,xb,ya,yb,a,b;	/*Point varibles*/
	double Clsf,Clsfb;				/*level set function value at element centre*/
	Coord *Cut;				/*Point element edge is cut*/
	int *Lstat;				/*Local Node Status*/
	int *Point;				/*edge the cut is on*/
	int NumNodes = (elemX+1)*(elemY+1);
	Lnodes = malloc(4*sizeof(int));
	Point = malloc(4*sizeof(int));
	Lstat = malloc(4*sizeof(int));
	Tlsf = malloc(4*sizeof(double));
	Tlsfb = malloc(4*sizeof(double));
	Cut = malloc(4*sizeof(Coord));

	for(i=0;i<elemX;i++)
	{
		for(j=0;j<elemY;j++)
		{
			num = Number[i][j].n-1;
			Lnodes[0] = Number[i][j].a-1;
			Lnodes[1] = Number[i][j].b-1;
			Lnodes[2] = Number[i][j].c-1;
			Lnodes[3] = Number[i][j].d-1;

			/*First get the lsf value at the element centre from the closet lsf (the mean according to linear shape functions)*/
			/*If there are two level sets we need to choose the one closer to element to define the fibre angle*/
			Clsf = 10000.0;
			for(p=0;p<Numlsf;p++)
			{
				q = NumNodes*p;
				for(k=0;k<4;k++)
				{
					Tlsfb[k] = lsf[q+Lnodes[k]];
				}
				Clsfb = 0.25*(Tlsfb[0]+Tlsfb[1]+Tlsfb[2]+Tlsfb[3]);

				/*if(num==0){printf("\np = %i, q = %i \nClsfb = %f, Clsf = %f",p, q, Clsfb, Clsf);
					   printf("\n Tlsfb(%f, %f, %f, %f)", Tlsfb[0], Tlsfb[1], Tlsfb[2],Tlsfb[3]);}*/

				if(fabs(Clsfb)<fabs(Clsf))
				{
					Clsf = Clsfb;
					for(k=0;k<4;k++)
					{
						Tlsf[k] = Tlsfb[k];
					}
				}
			}
			

			/*if(num==0){printf("\nClsf = %f\n", Clsf); printf("\n Tlsf(%f, %f, %f, %f)", Tlsf[0], Tlsf[1], Tlsf[2],Tlsf[3]);}*/
			/*In order to actually go through the centre Clsf should be greater then two of the nodal lsf values and smaller the the other two*/
			lsf_min = 1000;
			lsf_max = -1000;
			Scount = 0;
			Gcount = 0;
			for(k=0;k<4;k++)
			{
				if((Tlsf[k]-Clsf) < 0.0){ ++Scount;}
				else{++Gcount;}
				
				/*Store minmum and maximum lsf values*/
				if(Tlsf[k]<lsf_min)
				{
					lsf_min = Tlsf[k];
					minN = k;
				}
				if(Tlsf[k]>lsf_max)
				{
					lsf_max = Tlsf[k];
					maxN = k;
				}

			}
			/*printf("\n Clsf = %f", Clsf);*/
			/*If the condition of two lsf nodal values either side of the centre is not reached set Clsf to the mean of the two middle nodal lsf values*/
			if(Scount!=2)
			{
				Clsf = 0.0;
				for(k=0;k<4;k++)
				{
					if((k!=minN)&&(k!=maxN)){Clsf += Tlsf[k];}
				}
				Clsf /= 2.0;
			}

			/*Now we have our Clsf value we need to find the points where it crosses the edges of the elements.*/
			/*Get the local Node Status*/
			for(k=0;k<4;k++)
			{
				if(fabs(Tlsf[k]-Clsf)< 0.00000001){Lstat[k] = 2;}
				else if((Tlsf[k]-Clsf)>0.00000001){Lstat[k] = 1;}
				else if((Tlsf[k]-Clsf)<0.00000001){Lstat[k] = 0;}
				else{printf("*******************ERROR No Node Stat*********************");}
			}

			/*Select the two nodes on each edge in turn and see if it is cut*/
			Scount = 0;
			for(k=0;k<4;k++){ Point[k] = -5; Cut[k].x =0.0; Cut[k].y =0.0; }
			/*Edge0*/
			if(Lstat[0]!=Lstat[1])
			{
				/*printf("\nAA");*/
				Cut[0].y = 0.0;
				a = fabs(Clsf-Tlsf[0]);
				b = fabs(Clsf-Tlsf[1]);
				Cut[0].x = a/(a+b);
				Point[Scount] = 0;
				++Scount;
			}
			/*Edge1*/
			if(Lstat[1]!=Lstat[2])
			{
				/*printf("\nBB");*/
				Cut[1].x = 1.0;
				a = fabs(Clsf-Tlsf[1]);
				b = fabs(Clsf-Tlsf[2]);
				Cut[1].y = a/(a+b);
				Point[Scount] = 1;
				++Scount;
				
			}
			/*Edge2*/
			if(Lstat[2]!=Lstat[3])
			{
				/*printf("\nCC");*/
				Cut[2].y = 1.0;
				a = fabs(Clsf-Tlsf[3]);
				b = fabs(Clsf-Tlsf[2]);
				Cut[2].x = a/(a+b);
				Point[Scount] = 2;
				++Scount;
			}
			/*Edge3*/
			if(Lstat[0]!=Lstat[3])
			{
				/*printf("\nDD");*/
				Cut[3].x = 0.0;
				a = fabs(Clsf-Tlsf[0]);
				b = fabs(Clsf-Tlsf[3]);
				Cut[3].y = a/(a+b);
				Point[Scount] = 3;
				++Scount;
			}

			/*Now we have all the cuts we calculate the fibre angle, note only 2 o four cuts should be possible, else something is wrong*/

			if(Scount == 2)
			{
				/*printf("\n\nNumber of Cuts = %i,  num: %i, Clsf: %f\n lsf(%f, %f, %f, %f)", Scount, num, Clsf, Tlsf[0], Tlsf[1], Tlsf[2], Tlsf[3]);*/

				/*Simple to get the fibre angel strait from the line path*/
				p1 = Point[0];
				p2 = Point[1];
				/*printf("\n p1: %i, p2: %i", p1, p2);*/
				x1 = Cut[p1].x;
				x2 = Cut[p2].x;
				y1 = Cut[p1].y;
				y2 = Cut[p2].y;
				/*printf("\n x1: %f, y1: %f, x2: %f, y2: %f,", x1, y1, x2, y2);*/
				x = x2 - x1;
				y = y2 - y1;
				/*printf("\nX: %f, Y: %f", x, y);*/
				dtemp = y/x;
				ftemp = atan(dtemp);
				/*printf("\n dtemp: %f", dtemp);*/
			}
			else if(Scount ==4)
			{
				/*In this case we need to pair up the right egdes, this is done by knowing that for four cuts to occur there must be a
				diagonal of lsf nodal values that is larger then Clsf and the Clsf line cannot cross.*/				
				
				if(Lstat[0]*Lstat[2]!=0) /*(((Tlsf[0]>=Clsf)&&(Tlsf[2]>=Clsf))||((Lstat[0]==2)&&(Lstat[2]==2))) /*Try diagonals 0 and 2*/
				{
					/*Line one*/
					x1 = Cut[0].x;
					x2 = Cut[1].x;
					y1 = Cut[0].y;
					y2 = Cut[1].y;
					x = x2 - x1;
					y = y2 - y1;
					dtemp = y/x;
					ftemp1 = atan(dtemp);
					if(ftemp1 > 1.570796327){ftemp1 -= 3.14159265;}
					else if(ftemp1 < -1.570796327){ftemp1 += 3.14159265;}
					/*if(num==0){printf("\n\na: (%f, %f)", xa, ya);}*/

/*if((num==0)||(num==899)){printf("\n\n ftemp 1 = %f\tCut (%f, %f)  %f", ftemp1*(180/3.141592654), x, y, dtemp);}*/

					/*Line two*/
					x1 = Cut[3].x;
					x2 = Cut[2].x;
					y1 = Cut[3].y;
					y2 = Cut[2].y;
					x = x2 - x1;
					y = y2 - y1;
					dtemp = y/x;
					ftemp2 = atan(dtemp);
					if(ftemp2 > 1.570796327){ftemp2 -= 3.14159265;}
					else if(ftemp2 < -1.570796327){ftemp2 += 3.14159265;}

/*if((num==0)||(num==899)){printf("\n ftemp 2 = %f\t Cut (%f, %f)  %f", ftemp2*(180/3.141592654), x, y, dtemp);}*/

					ftemp = 0.5*(ftemp1+ftemp2);
					/*if(num==0){printf("\nb: (%f, %f)", xb, yb);}*/

				}
				else if(Lstat[1]*Lstat[3]!=0) /*(((Tlsf[1]>=Clsf)&&(Tlsf[3]>=Clsf))||((Lstat[1]==2)&&(Lstat[3]==2))) /*Try diagonals 3 and 1*/
				{

					/*Line one*/
					x1 = Cut[0].x;
					x2 = Cut[3].x;
					y1 = Cut[0].y;
					y2 = Cut[3].y;
					x = x2 - x1;
					y = y2 - y1;
					dtemp = y/x;
					ftemp1 = atan(dtemp);
					if(ftemp1 > 1.570796327){ftemp1 -= 3.14159265;}
					else if(ftemp1 < -1.570796327){ftemp1 += 3.14159265;}
					
/*if((num==0)||(num==899)){printf("\n\n ftemp 1 = %f\tCut (%f, %f)  %f", ftemp1*(180/3.141592654), x, y, dtemp);}*/

					/*Line two*/
					x1 = Cut[1].x;
					x2 = Cut[2].x;
					y1 = Cut[1].y;
					y2 = Cut[2].y;
					x = x2 - x1;
					y = y2 - y1;
					dtemp = y/x;
					ftemp2 = atan(dtemp);
					if(ftemp2 > 1.570796327){ftemp2 -= 3.14159265;}
					else if(ftemp2 < -1.570796327){ftemp2 += 3.14159265;}

/*if((num==0)||(num==899)){printf("\n ftemp 2 = %f\t Cut (%f, %f)  %f", ftemp2*(180/3.141592654), x, y, dtemp);}*/

					ftemp = 0.5*(ftemp1+ftemp2);


				}
				else
				{
				printf("*******************************/tError/t*******************************");
				printf("\nFour Cuts impossibly aranged \nNumber of Cuts = %i,  num: %i, Clsf: %f\n lsf(%f, %f, %f, %f)", Scount, num, Clsf, Tlsf[0], Tlsf[1], Tlsf[2], Tlsf[3]);
				printf("\n Point:(%i, %i, %i, %i)", Point[0], Point[1], Point[2], Point[3]);
				printf("\nDist = (%0.15f, %0.15f, %0.15f, %0.15f)", Tlsf[0]-Clsf, Tlsf[1]-Clsf, Tlsf[2]-Clsf, Tlsf[3]-Clsf);
				printf("\n Lstat:(%i, %i, %i, %i)", Lstat[0], Lstat[1], Lstat[2], Lstat[3]);
				}

			}
			else
			{
				printf("\n*******************************/tError/t*******************************");
				printf("\n Too few cuts \nNumber of Cuts = %i,  num: %i, Clsf: %0.9f\n lsf(%0.9f, %0.9f, %0.9f, %0.9f)\n", Scount, num, Clsf, Tlsf[0], Tlsf[1], Tlsf[2], Tlsf[3]);
				printf("\n Point:(%i, %i, %i, %i)", Point[0], Point[1], Point[2], Point[3]);
				printf("\nDist = (%0.15f, %0.15f, %0.15f, %0.15f)", Tlsf[0]-Clsf, Tlsf[1]-Clsf, Tlsf[2]-Clsf, Tlsf[3]-Clsf);
				printf("\n Lstat:(%i, %i, %i, %i)", Lstat[0], Lstat[1], Lstat[2], Lstat[3]);
			}

			/*Up date theta to be this angle*/
             Theta[num] = ftemp;
         

			/*Tidy up theta to be between 0 and 180 degrees*/
			if(Theta[num] > 1.570796327){Theta[num] -= 3.14159265;}
			else if(Theta[num] < -1.570796327){Theta[num] += 3.14159265;}	
             printf("\n element %d theta is equal to %f", num, Theta[num] );
            
            if (num==80||num==81||num==40||num==41) {
               // printf("\n theta value of element lol %d: %.16lf", num, Theta[num] );
            }
			/*if((num==0)||(num==899))
			{
				printf("\n\n Number of Cuts = %i,  num: %i, Clsf: %0.9f\n lsf(%0.9f, %0.9f, %0.9f, %0.9f)", Scount, num, Clsf, Tlsf[0], Tlsf[1], Tlsf[2], Tlsf[3]);
				printf("\n Point:(%i, %i, %i, %i)", Point[0], Point[1], Point[2], Point[3]);
				printf("\n Dist = (%0.15f, %0.15f, %0.15f, %0.15f)", Tlsf[0]-Clsf, Tlsf[1]-Clsf, Tlsf[2]-Clsf, Tlsf[3]-Clsf);
				printf("\n Cuts: (%f,%f)\t(%f,%f)\t(%f,%f)\t(%f,%f)\t", Cut[0].x, Cut[0].y, Cut[1].x, Cut[1].y, Cut[2].x, Cut[2].y, Cut[3].x, Cut[3].y);
				printf("\n Lstat:(%i, %i, %i, %i)", Lstat[0], Lstat[1], Lstat[2], Lstat[3]);
				printf("\n Theta = %f", Theta[num]*(180/3.141592654));
			}	*/

		}
	}

free(Lnodes);
free(Cut);
free(Tlsf);
free(Point);
}

void ThetaUpdate5(double *Theta, int elemX, int elemY, Elem Number[elemX][elemY], double *lsf, double *lsfb, int Numlsf, double h, double *hxz)
{
    double s, n, a, b, X, Y;			/*local shape funtion co-ordinate values*/
    int i,j,k,o,p1,p2,num,p,q;		/*interger varribles*/
    int *Lnodes;				/*Local Nodes*/
    double *Tlsf,*Tlsfb;				/*Local node lsf data*/
    double lsf_min, lsf_max;		/*max and min lsf in an elemenet varibles*/
    int minN, maxN;				/*max and min lsf number*/
    double ftemp, ftemp1, ftemp2, dtemp;	/*tempoary varibles*/
    int Scount,Gcount;			/*Count varibles*/
    double x,x1,x2,y,y1,y2,xa,xb,ya,yb;	/*Point varibles*/
    double Clsf,Clsfb;				/*level set function value at element centre*/
    Coord *Cut;				/*Point element edge is cut*/
    int *Lstat;				/*Local Node Status*/
    int *Point;				/*edge the cut is on*/
    int NumNodes = (elemX+1)*(elemY+1);
   	double *Nx, *Ny;					/*Differnetial of shape functions*/

    double dx, dy;						/*Differnetial of lsf wrt to x and y*/

    
    b = h;
    
    Lnodes = malloc(4*sizeof(int));
    Point = malloc(4*sizeof(int));
    Lstat = malloc(4*sizeof(int));
    Tlsf = malloc(4*sizeof(double));
    Tlsfb = malloc(4*sizeof(double));
    Cut = malloc(4*sizeof(Coord));
    Nx = malloc(4*sizeof(double));
    Ny = malloc(4*sizeof(double));
   
    
    
    for(i=0;i<elemX;i++)
    {   a = hxz[i];
        for(j=0;j<elemY;j++)
        {
            num = Number[i][j].n-1;
            Lnodes[0] = Number[i][j].a-1;
            Lnodes[1] = Number[i][j].b-1;
            Lnodes[2] = Number[i][j].c-1;
            Lnodes[3] = Number[i][j].d-1;
            
            /*First get the lsf value at the element centre from the closet lsf (the mean according to linear shape functions)*/
            /*If there are two level sets we need to choose the one closer to element to define the fibre angle*/
            Clsf = 10000.0;
            for(p=0;p<Numlsf;p++)
            {
                q = NumNodes*p;
                for(k=0;k<4;k++)
                {
                    Tlsfb[k] = lsf[q+Lnodes[k]];
                }
                Clsfb = 0.25*(Tlsfb[0]+Tlsfb[1]+Tlsfb[2]+Tlsfb[3]);
                
                /*if(num==0){printf("\np = %i, q = %i \nClsfb = %f, Clsf = %f",p, q, Clsfb, Clsf);
                 printf("\n Tlsfb(%f, %f, %f, %f)", Tlsfb[0], Tlsfb[1], Tlsfb[2],Tlsfb[3]);}*/
                
                if(fabs(Clsfb)<fabs(Clsf))
                {
                    Clsf = Clsfb;
                    for(k=0;k<4;k++)
                    {
                        Tlsf[k] = Tlsfb[k];
                    }
                }
            }
            
            /*Shape function local cordinates, for this problem a single gauss point at the centre of the element will be adequate*/
            s = 0.0;
            n = 0.0;
            
            X = ((1)/(4*a*a*b));
            Y = ((1)/(4*a*b*b));
            /*Shape function derivates*/
            Nx[0] = X*(-1+n);
            Nx[1] = X*(1-n);
            Nx[2] = X*(1+n);
            Nx[3] = X*(-1-n);
            
            Ny[0] = Y*(s-1);
            Ny[1] = Y*(-s-1);
            Ny[2] = Y*(s+1);
            Ny[3] = Y*(-s+1);
            
            dx = 0.0;
            dy = 0.0;
            
            /*Now use the shape function derivaties to calcualte x and y lsf gradients over the element*/
            for(k=0;k<4;k++)
            {
                dx += Nx[k]*Tlsf[k];
                dy += Ny[k]*Tlsf[k];
            }
            /*Now we have the gradients we need to get the angle theta from it*/
            
            if(fabs(dx)<0.00000001) {ftemp = 1.570796327;} /*Set theta to 90 since tan90 is infinate*/
            else if (fabs(dy)<0.00000001) {ftemp = 0.0;}
            else
            {
                dtemp = dy/dx;
                ftemp = atan(dtemp);
            }
            
            /*Up date theta to be the normal to this angle*/
            Theta[num] = ftemp + 1.570796327;
            
            
            
            
            //  This part is to verify the sensitivity
            double dT;
            dT=0.0000001;
            if (num==80)
            {
                // Theta[num]+=dT;
                
                // C = cos(Theta[num]+dT);
                //printf("C= %f", C);
                //S = sin(Theta[num]+dT);
            }
            
            /*Tidy up theta to be between 0 and 180 degrees*/
            if(Theta[num] > 1.570796327){Theta[num] -= 3.14159265;}
            else if(Theta[num] < -1.570796327){Theta[num] += 3.14159265;}
            printf("\n element %d theta is equal to %f", num, Theta[num] );
            /*if((num==19)||(num==400))
             {
             printf("\n\n %i", num);
             printf("\n lsfX:(%f, %f, %f, %f)", Tlsf[0], Tlsf[1], Tlsf[2], Tlsf[3]);
             printf("\n dx: %f", dx);
             printf("\n dy: %f", dy);
             printf("\n ftemp: %f", ftemp);
             printf("\n Theta: %f", Theta[num]);
             }*/
           
            
            
            
            
        }
        
    }
    
    free(Lnodes);
    free(Cut);
    free(Tlsf);
    free(Point);
}



void Transform(double *KE,double *T, double *K, int n)
{
	/*printf("\n");
	printf("\n n = %i", n);
	printf("\n Trans[%i] = %f", 9*n+0, T[9*n+0]);
	printf("\n Trans[%i] = %f", 9*n+1, T[9*n+1]);
	printf("\n Trans[%i] = %f", 9*n+2, T[9*n+2]);
	printf("\n Trans[%i] = %f", 9*n+3, T[9*n+3]);
	printf("\n Trans[%i] = %f", 9*n+4, T[9*n+4]);
	printf("\n Trans[%i] = %f", 9*n+5, T[9*n+5]);
	printf("\n Trans[%i] = %f", 9*n+6, T[9*n+6]);
	printf("\n Trans[%i] = %f", 9*n+7, T[9*n+7]);
	printf("\n Trans[%i] = %f", 9*n+8, T[9*n+8]);*/
	int m = n*9;
	int p;
	int j=0;
	int k=0;
	int i;
	for(i=0;i<24;i++)
	{
	p = m+j;
	K[0+i] = (T[0+m]*KE[0+k]+T[3+m]*KE[24+k]+T[6+m]*KE[48+k])*T[0+p] + (T[0+m]*KE[1+k]+T[3+m]*KE[25+k]+T[6+m]*KE[49+k])*T[3+p] + (T[0+m]*KE[2+k]+T[3+m]*KE[26+k]+T[6+m]*KE[50+k])*T[6+p];
	K[24+i] = (T[1+m]*KE[0+k]+T[4+m]*KE[24+k]+T[7+m]*KE[48+k])*T[0+p] + (T[1+m]*KE[1+k]+T[4+m]*KE[25+k]+T[7+m]*KE[49+k])*T[3+p] + (T[1+m]*KE[2+k]+T[4+m]*KE[26+k]+T[7+m]*KE[50+k])*T[6+p];
	K[48+i] = (T[2+m]*KE[0+k]+T[5+m]*KE[24+k]+T[8+m]*KE[48+k])*T[0+p] + (T[2+m]*KE[1+k]+T[5+m]*KE[25+k]+T[8+m]*KE[49+k])*T[3+p] + (T[2+m]*KE[2+k]+T[5+m]*KE[26+k]+T[8+m]*KE[50+k])*T[6+p];

	K[72+i] = (T[0+m]*KE[72+k]+T[3+m]*KE[96+k]+T[6+m]*KE[120+k])*T[0+p] + (T[0+m]*KE[73+k]+T[3+m]*KE[97+k]+T[6+m]*KE[121+k])*T[3+p] + (T[0+m]*KE[74+k]+T[3+m]*KE[98+k]+T[6+m]*KE[122+k])*T[6+p];
	K[96+i] = (T[1+m]*KE[72+k]+T[4+m]*KE[96+k]+T[7+m]*KE[120+k])*T[0+p] + (T[1+m]*KE[73+k]+T[4+m]*KE[97+k]+T[7+m]*KE[121+k])*T[3+p] + (T[1+m]*KE[74+k]+T[4+m]*KE[98+k]+T[7+m]*KE[122+k])*T[6+p];
	K[120+i] = (T[2+m]*KE[72+k]+T[5+m]*KE[96+k]+T[8+m]*KE[120+k])*T[0+p] + (T[2+m]*KE[73+k]+T[5+m]*KE[97+k]+T[8+m]*KE[121+k])*T[3+p] + (T[2+m]*KE[74+k]+T[5+m]*KE[98+k]+T[8+m]*KE[122+k])*T[6+p];

	K[144+i] = (T[0+m]*KE[144+k]+T[3+m]*KE[168+k]+T[6+m]*KE[192+k])*T[0+p] + (T[0+m]*KE[145+k]+T[3+m]*KE[169+k]+T[6+m]*KE[193+k])*T[3+p] + (T[0+m]*KE[146+k]+T[3+m]*KE[170+k]+T[6+m]*KE[194+k])*T[6+p];
	K[168+i] = (T[1+m]*KE[144+k]+T[4+m]*KE[168+k]+T[7+m]*KE[192+k])*T[0+p] + (T[1+m]*KE[145+k]+T[4+m]*KE[169+k]+T[7+m]*KE[193+k])*T[3+p] + (T[1+m]*KE[146+k]+T[4+m]*KE[170+k]+T[7+m]*KE[194+k])*T[6+p];
	K[192+i] = (T[2+m]*KE[144+k]+T[5+m]*KE[168+k]+T[8+m]*KE[192+k])*T[0+p] + (T[2+m]*KE[145+k]+T[5+m]*KE[169+k]+T[8+m]*KE[193+k])*T[3+p] + (T[2+m]*KE[146+k]+T[5+m]*KE[170+k]+T[8+m]*KE[194+k])*T[6+p];

	K[216+i] = (T[0+m]*KE[216+k]+T[3+m]*KE[240+k]+T[6+m]*KE[264+k])*T[0+p] + (T[0+m]*KE[217+k]+T[3+m]*KE[241+k]+T[6+m]*KE[265+k])*T[3+p] + (T[0+m]*KE[218+k]+T[3+m]*KE[242+k]+T[6+m]*KE[266+k])*T[6+p];
	K[240+i] = (T[1+m]*KE[216+k]+T[4+m]*KE[240+k]+T[7+m]*KE[264+k])*T[0+p] + (T[1+m]*KE[217+k]+T[4+m]*KE[241+k]+T[7+m]*KE[265+k])*T[3+p] + (T[1+m]*KE[218+k]+T[4+m]*KE[242+k]+T[7+m]*KE[266+k])*T[6+p];
	K[264+i] = (T[2+m]*KE[216+k]+T[5+m]*KE[240+k]+T[8+m]*KE[264+k])*T[0+p] + (T[2+m]*KE[217+k]+T[5+m]*KE[241+k]+T[8+m]*KE[265+k])*T[3+p] + (T[2+m]*KE[218+k]+T[5+m]*KE[242+k]+T[8+m]*KE[266+k])*T[6+p];

	K[288+i] = (T[0+m]*KE[288+k]+T[3+m]*KE[312+k]+T[6+m]*KE[336+k])*T[0+p] + (T[0+m]*KE[289+k]+T[3+m]*KE[313+k]+T[6+m]*KE[337+k])*T[3+p] + (T[0+m]*KE[290+k]+T[3+m]*KE[314+k]+T[6+m]*KE[338+k])*T[6+p];
	K[312+i] = (T[1+m]*KE[288+k]+T[4+m]*KE[312+k]+T[7+m]*KE[336+k])*T[0+p] + (T[1+m]*KE[289+k]+T[4+m]*KE[313+k]+T[7+m]*KE[337+k])*T[3+p] + (T[1+m]*KE[290+k]+T[4+m]*KE[314+k]+T[7+m]*KE[338+k])*T[6+p];
	K[336+i] = (T[2+m]*KE[288+k]+T[5+m]*KE[312+k]+T[8+m]*KE[336+k])*T[0+p] + (T[2+m]*KE[289+k]+T[5+m]*KE[313+k]+T[8+m]*KE[337+k])*T[3+p] + (T[2+m]*KE[290+k]+T[5+m]*KE[314+k]+T[8+m]*KE[338+k])*T[6+p];

	K[360+i] = (T[0+m]*KE[360+k]+T[3+m]*KE[384+k]+T[6+m]*KE[408+k])*T[0+p] + (T[0+m]*KE[361+k]+T[3+m]*KE[385+k]+T[6+m]*KE[409+k])*T[3+p] + (T[0+m]*KE[362+k]+T[3+m]*KE[386+k]+T[6+m]*KE[410+k])*T[6+p];
	K[384+i] = (T[1+m]*KE[360+k]+T[4+m]*KE[384+k]+T[7+m]*KE[408+k])*T[0+p] + (T[1+m]*KE[361+k]+T[4+m]*KE[385+k]+T[7+m]*KE[409+k])*T[3+p] + (T[1+m]*KE[362+k]+T[4+m]*KE[386+k]+T[7+m]*KE[410+k])*T[6+p];
	K[408+i] = (T[2+m]*KE[360+k]+T[5+m]*KE[384+k]+T[8+m]*KE[408+k])*T[0+p] + (T[2+m]*KE[361+k]+T[5+m]*KE[385+k]+T[8+m]*KE[409+k])*T[3+p] + (T[2+m]*KE[362+k]+T[5+m]*KE[386+k]+T[8+m]*KE[410+k])*T[6+p];

	K[432+i] = (T[0+m]*KE[432+k]+T[3+m]*KE[456+k]+T[6+m]*KE[480+k])*T[0+p] + (T[0+m]*KE[433+k]+T[3+m]*KE[457+k]+T[6+m]*KE[481+k])*T[3+p] + (T[0+m]*KE[434+k]+T[3+m]*KE[458+k]+T[6+m]*KE[482+k])*T[6+p];
	K[456+i] = (T[1+m]*KE[432+k]+T[4+m]*KE[456+k]+T[7+m]*KE[480+k])*T[0+p] + (T[1+m]*KE[433+k]+T[4+m]*KE[457+k]+T[7+m]*KE[481+k])*T[3+p] + (T[1+m]*KE[434+k]+T[4+m]*KE[458+k]+T[7+m]*KE[482+k])*T[6+p];
	K[480+i] = (T[2+m]*KE[432+k]+T[5+m]*KE[456+k]+T[8+m]*KE[480+k])*T[0+p] + (T[2+m]*KE[433+k]+T[5+m]*KE[457+k]+T[8+m]*KE[481+k])*T[3+p] + (T[2+m]*KE[434+k]+T[5+m]*KE[458+k]+T[8+m]*KE[482+k])*T[6+p];

	K[504+i] = (T[0+m]*KE[504+k]+T[3+m]*KE[528+k]+T[6+m]*KE[552+k])*T[0+p] + (T[0+m]*KE[505+k]+T[3+m]*KE[529+k]+T[6+m]*KE[553+k])*T[3+p] + (T[0+m]*KE[506+k]+T[3+m]*KE[530+k]+T[6+m]*KE[554+k])*T[6+p];
	K[528+i] = (T[1+m]*KE[504+k]+T[4+m]*KE[528+k]+T[7+m]*KE[552+k])*T[0+p] + (T[1+m]*KE[505+k]+T[4+m]*KE[529+k]+T[7+m]*KE[553+k])*T[3+p] + (T[1+m]*KE[506+k]+T[4+m]*KE[530+k]+T[7+m]*KE[554+k])*T[6+p];
	K[552+i] = (T[2+m]*KE[504+k]+T[5+m]*KE[528+k]+T[8+m]*KE[552+k])*T[0+p] + (T[2+m]*KE[505+k]+T[5+m]*KE[529+k]+T[8+m]*KE[553+k])*T[3+p] + (T[2+m]*KE[506+k]+T[5+m]*KE[530+k]+T[8+m]*KE[554+k])*T[6+p];
	j++;
	if(j==3){j=0;k+=3;}
	}

}
