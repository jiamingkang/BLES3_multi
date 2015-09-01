/*
 *  Strain.c
 *  
 *  Created by Peter Dunning on 20/05/2008.
 *	Modified last 15/05/2009 - Changed the way strain is calculated for IN and quadrilateral elements
 */

#include "Strain.h"
#include "ls_types.h"
#include "FixedGrid.h"
#include "Solve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Array storing data for location of the 4 nodes in the square element*/
static Pos Po[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};

static int rows[10] = {1,1,1,1,2,2,2,3,3,4};
static int cols[10] = {1,2,3,4,2,3,4,3,4,4};

const double ga = 0.57735026919;	/*Constant for 2-point gauss rule*/
const double rs = 0.408248290464;	/*Square root of 1/6*/
const double rt = 1.414213562373;	/*Square root of 2*/

/*Function that calculates the strain energy of a node by a least squares (bi-linear order) filter of near-by gauss points, using MA44*/
double LstrainV2(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag)
{
	int i,j; /*Incrementor*/
	double r2 = rad * rad; /*radius squared*/
	double ftemp,x2,y2; /*variables for squared co-ordinates*/
	double a2;	/*variable for area squared*/
	double wgt;	/*varible for weight function*/
	int *pts;
	pts = malloc(gpoints * sizeof(int)); /*array to store gpoint numbers that lie within radius*/
	double *dist;
	dist = malloc(gpoints * sizeof(double)); /*array to store distance of points from point of interest*/
	
	double min = 0.0;
	int count = 0;	/*initialize count of number of points used to zero*/
	
	for(i=0;i<gpoints;i++)
	{
		ftemp = nx - Genergy[i].x;
		x2 = ftemp * ftemp;	/*squared x-distance from current gauss point to node of interest*/
		ftemp = ny - Genergy[i].y;
		y2 = ftemp * ftemp;	/*squared y-distance from current gauss point to node of interest*/
		
		ftemp = x2 + y2; /*squared distance between gauss point and node*/
		
		if((ftemp - r2) < 0.000001) /*if distance less than the radius squared, then add data to arrays*/
		{
			dist[count] = sqrt(ftemp);
			pts[count++] = i;
		}
	}
	
	if(count < 4)  /*if not enough points within radius then point must be on an island!*/
	{
		free(pts);
		free(dist);
		/*printf("\nNode at [%fl, %fl] considered part of an island", nx,ny);*/
		return(0.0); /*hence return zero strain energy*/
	}
		
	/*printf("\n%i found points",count);*/
	
	pts = realloc(pts,count * sizeof(int)); /*resize pts array to number of points inside radius*/
	dist = realloc(dist,count * sizeof(double)); /*resize pts array to number of points inside radius*/
	
	double *A; /*array to store gpoint values*/
	A = malloc( 4 * count * sizeof(double));
	double *B;	/*array to store rhs of equation*/
	B = malloc(count * sizeof(double));
	double *Coeff;
	Coeff = calloc(4,sizeof(double));
	
	double xb, yb; /*variables to calc x & y co-ords relative to point of interest*/

	for(i=0;i<count;i++)
	{
		j = pts[i];
		
		/*calculate the weight function*/
		switch (wFlag)
		{
			case 1:
					wgt = 1.0; /*no weighting*/
					break;
			case 2:
					wgt = 1.0 / sqrt(dist[i]); /*weighted by inverse distance*/
					break;
			case 3:
					wgt = sqrt(Genergy[j].a); /*weighted by area*/
					break;
			case 4:
					wgt = sqrt(Genergy[j].a / dist[i]); /*weighted by area & inverse distance*/
					break;
			default:
					printf("\nERROR!! wFlag out of range!");
		}
		
		xb = Genergy[j].x - nx; /*relative x-coord*/
		yb = Genergy[j].y - ny; /*relative y-coord*/
		
		A[i] = wgt;
		A[count + i] = xb * wgt;
		A[(2 * count) + i] = yb * wgt;
		A[(3 * count) + i] = xb * yb * wgt;
		
		B[i] = Genergy[j].u * wgt;
	}
	
	/*printf("\nindexing complete");
	/*solve least squares problem using MA44*/
	leastsq(count, 4, 0, A, B, Coeff, -1);
	
	/*Finally evaluate strain energy at the node using the co-efficients*/
	ftemp = Coeff[0]; /*only require first as relative co-ords were used*/
	
	free(A);
	free(B);
	free(Coeff);
	free(pts);
	free(dist);
	
	return(ftemp);
}

/*Function that calculates the strain energy of a node by a least squares (2nd order) filter of near-by gauss points*/
double LstrainV22(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag)
{
	int i,j; /*Incrementor*/
	double r2 = rad * rad; /*radius squared*/
	double ftemp,x2,y2,xb,yb; /*variables for squared co-ordinates*/
	double a2;	/*variable for area squared*/
	double wgt;	/*varible for weight function*/
	int *pts;
	pts = malloc(gpoints * sizeof(int)); /*array to store gpoint numbers that lie within radius*/
	double *dist;
	dist = malloc(gpoints * sizeof(double)); /*array to store distance of points from point of interest*/
	
	double min = 0.0;
	int count = 0;	/*initialize count of number of points used to zero*/
	
	for(i=0;i<gpoints;i++)
	{
		ftemp = nx - Genergy[i].x;
		x2 = ftemp * ftemp;	/*squared x-distance from current gauss point to node of interest*/
		ftemp = ny - Genergy[i].y;
		y2 = ftemp * ftemp;	/*squared y-distance from current gauss point to node of interest*/
		
		ftemp = x2 + y2; /*squared distance between gauss point and node*/
		
		if((ftemp < r2) && (r2 >= min))/*if distance less than the radius squared, then add data to arrays*/
		{
			dist[count] = sqrt(ftemp);
			pts[count++] = i;
		}
	}
	
	if(count < 8)  /*if not enough points within radius then point must be on an island!*/
	{
		free(pts);
		free(dist);
		printf("\nNode at [%fl, %fl] considered part of an island in LstrainV22", nx,ny);
		return(0.0); /*hence return zero strain energy*/
	}
	
	pts = realloc(pts,count * sizeof(int)); /*resize pts array to number of points inside radius*/
	dist = realloc(dist,count * sizeof(double)); /*resize pts array to number of points inside radius*/
	
	double *A; /*array to store gpoint values*/
	A = malloc( 6 * count * sizeof(double));
	double *B;	/*array to store rhs of equation*/
	B = malloc(count * sizeof(double));
	double *Coeff;
	Coeff = calloc(6,sizeof(double));
	
	for(i=0;i<count;i++)
	{
		j = pts[i];
		
		/*calculate the weight function*/
		switch (wFlag)
		{
			case 1:
					wgt = 1.0; /*no weighting*/
					break;
			case 2:
					wgt = 1.0 / sqrt(dist[i]); /*weighted by inverse distance*/
					break;
			case 3:
					wgt = sqrt(Genergy[j].a); /*weighted by area*/
					break;
			case 4:
					wgt = sqrt(Genergy[j].a / dist[i]); /*weighted by area & inverse distance*/
					break;
			default:
					printf("\nERROR!! wFlag out of range!");
		}
		
		xb = Genergy[j].x - nx; /*relative x-coord*/
		yb = Genergy[j].y - ny; /*relative y-coord*/
		
		A[i] = wgt;
		A[count + i] = xb * wgt;
		A[(2 * count) + i] = yb * wgt;
		A[(3 * count) + i] = xb * yb * wgt;
		A[(4 * count) + i] = xb * xb * wgt;
		A[(5 * count) + i] = yb * yb * wgt;
		
		B[i] = Genergy[j].u * wgt;
	}
	
	/*printf("\nindexing complete");
	/*solve least squares problem using MA44*/
	leastsq(count, 6, 0, A, B, Coeff, 0);
	
	/*Finally evaluate strain energy at the node using the co-efficients*/
	ftemp = Coeff[0]; /*only require first as relative co-ords were used*/
	
	free(A);
	free(B);
	free(Coeff);
	free(pts);
	free(dist);
	
	return(ftemp);
}

/*Function that calculates the strain energy of a node by a least squares (2nd order) filter of near-by gauss points*/
double LstrainV23(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag)
{
	int i,j; /*Incrementor*/
	double r2 = rad * rad; /*radius squared*/

	double ftemp,x2,y2,xb,yb; /*variables for squared co-ordinates*/
	double a2;	/*variable for area squared*/
	double wgt;	/*varible for weight function*/
	int *pts;
	pts = malloc(gpoints * sizeof(int)); /*array to store gpoint numbers that lie within radius*/
	double *dist;
	dist = malloc(gpoints * sizeof(double)); /*array to store distance of points from point of interest*/
	
	double min = 0.0;
	int count = 0;	/*initialize count of number of points used to zero*/
	
	for(i=0;i<gpoints;i++)
	{
		ftemp = nx - Genergy[i].x;
		x2 = ftemp * ftemp;	/*squared x-distance from current gauss point to node of interest*/
		ftemp = ny - Genergy[i].y;
		y2 = ftemp * ftemp;	/*squared y-distance from current gauss point to node of interest*/
		
		ftemp = x2 + y2; /*squared distance between gauss point and node*/
		
		if((ftemp < r2) && (r2 >= min))/*if distance less than the radius squared, then add data to arrays*/
		{
			dist[count] = sqrt(ftemp);
			pts[count++] = i;
		}
	}
	
	if(count < 8)  /*if not enough points within radius then point must be on an island!*/
	{
		free(pts);
		free(dist);
		printf("\nNode at [%fl, %fl] considered part of an island in LstrainV22", nx,ny);
		return(0.0); /*hence return zero strain energy*/
	}
	
	pts = realloc(pts,count * sizeof(int)); /*resize pts array to number of points inside radius*/
	dist = realloc(dist,count * sizeof(double)); /*resize pts array to number of points inside radius*/
	
	double *A; /*array to store gpoint values*/
	A = malloc( 6 * count * sizeof(double));
	double *B;	/*array to store rhs of equation*/
	B = malloc(count * sizeof(double));
	double *Coeff;
	Coeff = calloc(6,sizeof(double));
	
	for(i=0;i<count;i++)
	{
		j = pts[i];
		
		/*calculate the weight function*/
		switch (wFlag)
		{
			case 1:
					wgt = 1.0; /*no weighting*/
					break;
			case 2:
					wgt = 1.0 / sqrt(dist[i]); /*weighted by inverse distance*/
					break;
			case 3:
					wgt = sqrt(Genergy[j].a); /*weighted by area*/
					break;
			case 4:
					wgt = sqrt(Genergy[j].a / dist[i]); /*weighted by area & inverse distance*/
					break;
			default:
					printf("\nERROR!! wFlag out of range!");
		}
		
		xb = Genergy[j].x - nx; /*relative x-coord*/
		yb = Genergy[j].y - ny; /*relative y-coord*/
		
		A[i] = wgt;
		A[count + i] = xb * wgt;
		A[(2 * count) + i] = yb * wgt;
		A[(3 * count) + i] = xb * yb * wgt;
		A[(4 * count) + i] = xb * xb * wgt;
		A[(5 * count) + i] = yb * yb * wgt;
		
		B[i] = Genergy[j].u * wgt;
	}
	
	/*printf("\nindexing complete");
	/*solve least squares problem using MA44*/
	leastsq(count, 6, 0, A, B, Coeff, 0);
	
	/*Finally evaluate strain energy at the node using the co-efficients*/
	ftemp = (Coeff[0]); /*only require first as relative co-ords were used  ftemp = (Coeff[0]/count);/(double)count */
	
	free(A);
	free(B);
	free(Coeff);
	free(pts);
	free(dist);
	
	return(ftemp);
}

double LstrainV24(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag)
{
    int i,j; /*Incrementor*/
    double r2 = rad * rad; /*radius squared*/
    
    double ftemp,x2,y2,xb,yb; /*variables for squared co-ordinates*/
    double a2;	/*variable for area squared*/
    double wgt;	/*varible for weight function*/
    int *pts;
    pts = malloc(gpoints * sizeof(int)); /*array to store gpoint numbers that lie within radius*/
    double *dist;
    dist = malloc(gpoints * sizeof(double)); /*array to store distance of points from point of interest*/
    
    double min = 0.0;
    int count = 0;	/*initialize count of number of points used to zero*/
    
    for(i=0;i<gpoints;i++)
    {
        ftemp = nx - Genergy[i].x;
        x2 = ftemp * ftemp;	/*squared x-distance from current gauss point to node of interest*/
        ftemp = ny - Genergy[i].y;
        y2 = ftemp * ftemp;	/*squared y-distance from current gauss point to node of interest*/
        
        ftemp = x2 + y2; /*squared distance between gauss point and node*/
        
        if((ftemp < r2) && (r2 >= min))/*if distance less than the radius squared, then add data to arrays*/
        {
            dist[count] = sqrt(ftemp);
            pts[count++] = i;
        }
    }
    
    if(count < 8)  //if not enough points within radius then point must be on an island!
    {
        free(pts);
        free(dist);
        printf("\nNode at [%fl, %fl] considered part of an island in LstrainV24", nx,ny);
        return(0.0); //hence return zero strain energy
    }
    
    
    pts = realloc(pts,count * sizeof(int)); /*resize pts array to number of points inside radius*/
    dist = realloc(dist,count * sizeof(double)); /*resize pts array to number of points inside radius*/
    
    double *A; /*array to store gpoint values*/
    A = malloc( 6 * count * sizeof(double));
    double *B;	/*array to store rhs of equation*/
    B = malloc(count * sizeof(double));
    double *Coeff;
    Coeff = calloc(6,sizeof(double));
    
    for(i=0;i<count;i++)
    {
        j = pts[i];
        
        /*calculate the weight function*/
        switch (wFlag)
        {
            case 1:
                wgt = 1.0; /*no weighting*/
                break;
            case 2:
                wgt = 1.0 / sqrt(dist[i]); /*weighted by inverse distance*/
                break;
            case 3:
                wgt = sqrt(Genergy[j].a); /*weighted by area*/
                break;
            case 4:
                wgt = sqrt(Genergy[j].a / dist[i]); /*weighted by area & inverse distance*/
                break;
            default:
                printf("\nERROR!! wFlag out of range!");
        }
        
        xb = Genergy[j].x - nx; /*relative x-coord*/
        yb = Genergy[j].y - ny; /*relative y-coord*/
        
        A[i] = wgt;
        A[count + i] = xb * wgt;
        A[(2 * count) + i] = yb * wgt;
        A[(3 * count) + i] = xb * yb * wgt;
        A[(4 * count) + i] = xb * xb * wgt;
        A[(5 * count) + i] = yb * yb * wgt;
        
        B[i] = Genergy[j].u * wgt;
    }
    
    /*printf("\nindexing complete");
     /*solve least squares problem using MA44*/
    leastsq(count, 6, 0, A, B, Coeff, 0);
    
    /*Finally evaluate strain energy at the node using the co-efficients*/
    ftemp = (Coeff[0]); /*only require first as relative co-ords were used  ftemp = (Coeff[0]/count);/(double)count */
    
    free(A);
    free(B);
    free(Coeff);
    free(pts);
    free(dist);
    
    return(ftemp);
}
/*Function to calculate Isoparametric Gauss point strain energy values*/	
void GaIsoStrain(int *tnodes, Coord *tcords, double *rhs1, double *rhs2, double alpha, double e11, double v12, double g33, int Gcount, Gstn *Genergy,
					int num, double *Eenergy, double Earea, double wgt)
{
	double j11,j12,j21,j22,Jdet;	/*Variables for Jacobian matrix entries & Determinant*/
	double *Bs;			/*Array to store Strain-Displacement matrix entries*/
	Bs = malloc(8 * sizeof(double));
	double ex1,ey1,txy1,ex2,ey2,txy2,Stn;	/*Varaibles for strains at a gauss point*/
	double gax,gay;			/*Gauss point values*/
	int i,k,temp;			/*Incrementors etc*/
	double ftemp,ftemp2;		/*double temp varaibles*/
	Coord *tdisp1;			/*Temp array to read in primary displacement values*/
	tdisp1 = malloc(4 * sizeof(Coord));
	Coord *tdisp2;			/*Temp array to read in adjoint displacement values*/
	tdisp2 = malloc(4 * sizeof(Coord));
	
	/*first work out global gauss point co-ordinates and store assosiated area*/
	for(i=0;i<4;i++)






	{
		ex1 = 0.0;
		ey1 = 0.0; /*here these are used to sum the global gauss point co-ordinates*/
		
		/*Gauss point Coords in Isoparametric domain*/
		gax = Po[i].x * ga;
		gay = Po[i].y * ga;
		
		for(k=0;k<4;k++)
		{
			/*here txy1 is used to evaluate the shape function at a gauss point*/
			txy1 = (1.0 + (Po[k].x * gax));
			txy1 *= (1.0 + (Po[k].y * gay));
			txy1 *= 0.25;
			
			/*update global gauss point co-ordintae sums*/
			ex1 += txy1 * tcords[k].x;
			ey1 += txy1 * tcords[k].y;
		}
		
		/*store gauss point co-ordinates and associated area ratio*/
		Genergy[Gcount + i].x = ex1;
		Genergy[Gcount + i].y = ey1;
		Genergy[Gcount + i].a = alpha;
	}
	
	double xc =  0.0;
	double yc = 0.0;	/*Element center coords*/
	
	for(k=0;k<4;k++)
	{
		/*update global gauss point co-ordintae sums*/
		xc += tcords[k].x;
		yc += tcords[k].y;
	}
	
	xc *= 0.25;
	yc *= 0.25;
	
	/*Now modify coords of nodes relative to element center*/
	for(i=0;i<4;i++)
	{
		tcords[i].x -= xc;
		tcords[i].y -= yc;
		/*printf("\nMod: tcords[%i].x=%fl, y=%fl",i+1,tcords[i].x,tcords[i].y);*/
	}
	
	/*Populate element displacement arrays (primary & adjoint) */
	for(i=0;i<4;i++)
	{
		temp = tnodes[i] * 2;
		tdisp1[i].x = rhs1[temp];			/*primary X displacement of node i*/
		tdisp1[i].y = rhs1[temp + 1];		/*primary Y displacement of node i*/
		/*printf("\ndisp[%i].x = %fl, y = %fl",i+1,tdisp[i].x,tdisp[i].y);*/
		tdisp2[i].x = rhs2[temp];			/*adjoint X displacement of node i*/
		tdisp2[i].y = rhs2[temp + 1];		/*adjoint Y displacement of node i*/
	}
	
	/*Now calculate strain energies for all four gauss points*/
	for(i=0;i<4;i++)	
	{
		/*First initialize some values to zero*/
		j11 = 0.0;
		j12 = 0.0;
		j21 = 0.0;
		j22 = 0.0;
		
		ex1 = 0.0;
		ey1 = 0.0;
		txy1 = 0.0;
		ex2 = 0.0;
		ey2 = 0.0;
		txy2 = 0.0;
		
		/*Gauss point Coords in Isoparametric domain*/
		gax = Po[i].x * ga;
		gay = Po[i].y * ga;
		
		/*Calcualte jacobian values at current gauss point*/
		for(k=0;k<4;k++)
		{
			ftemp = Po[k].x * (1.0 + (Po[k].y * gay));
			j11 += ftemp * tcords[k].x;
			j12 += ftemp * tcords[k].y;
			ftemp = Po[k].y * (1.0 + (Po[k].x * gax));
			j21 += ftemp * tcords[k].x;
			j22 += ftemp * tcords[k].y;
		}
		
		/*Need to divide by 4, to evaluate properly*/
		j11 *= 0.25;
		j12 *= 0.25;
		j21 *= 0.25;
		j22 *= 0.25;
		
		/*printf("\nj11=%fl, j12=%fl, j21=%fl, j22=%fl",j11,j12,j21,j22);*/
		
		Jdet = (j11 * j22) - (j12 * j21);	/*Deteminant of Jacobian matrix*/
		
		/*Evaluate Strain-Displacement matrix entires*/
		for(k=0;k<4;k++)
		{
			ftemp = 0.25 * Po[k].x * (1.0 + (Po[k].y * gay));
			ftemp2 = 0.25 * Po[k].y * (1.0 + (Po[k].x * gax));
			Bs[k] = (j22 * ftemp) - (j12 * ftemp2);
			Bs[4 + k] = (j11 * ftemp2) - (j21 * ftemp);
			
			/*printf("\nBs[0][%i] = %fl, Bs[1] = %fl",k,Bs[0][k],Bs[1][k]);*/
		}
		
		/*calcualte primary and adjoint strain values at current gauss point*/
		for(k=0;k<4;k++)
		{
			/*primary*/
			ex1 += Bs[k] * tdisp1[k].x;
			ey1 += Bs[4 + k] * tdisp1[k].y;
			txy1 += Bs[4 + k] * tdisp1[k].x;
			txy1 += Bs[k] * tdisp1[k].y;
			/*adjoint*/
			ex2 += Bs[k] * tdisp2[k].x;
			ey2 += Bs[4 + k] * tdisp2[k].y;
			txy2 += Bs[4 + k] * tdisp2[k].x;
			txy2 += Bs[k] * tdisp2[k].y;
		}
		
		/*printf("\nex=%fl, ey=%fl, txy=%fl",ex,ey,txy);
		printf("\nJdet=%fl",Jdet);*/
		
		/*Need to divide by Jdet to get proper values*/
		ex1 /= Jdet;
		ey1 /= Jdet;
		txy1 /= Jdet;
		ex2 /= Jdet;
		ey2 /= Jdet;
		txy2 /= Jdet;
		
		/*Calculate Strain Energy at the gauss point and store*/
		Stn = e11 * ((ex1 * ex2) + (ey1 * ey2));
		Stn += v12 * ( (ex1 * ey2) + (ex2 * ey1) );
		Stn += g33 * txy1 * txy2;
		/*Stn *= 0.5;
		/*printf("\nGstrain %i = %fl",i+1,Stn);*/
		Genergy[Gcount + i].u = Stn;
		Eenergy[num] += wgt * Stn * Jdet; /*update element strain energy*/
	}
	
	free(tdisp1);
	free(tdisp2);
	free(Bs);	
}

/*Function to calculate integration point coords for pentagonal strain energy calculation*/
void Pint(Coord Pcrds[6], Coord Icrds[5])
{
	int i;
	double c1,c2;
	
	for(i=1;i<6;i++)
	{
		/*Calcualte x coordinate as center of triangular section*/
		c1 = Pcrds[i].x;
		c2 = (i < 5) ? Pcrds[i+1].x : Pcrds[1].x;
		Icrds[i-1].x = c1 + c2; /*center is average of the 3 coords -NB Center coord = (0,0)*/
		Icrds[i-1].x *= 0.333333333333;
		
		c1 = Pcrds[i].y;
		c2 = (i < 5) ? Pcrds[i+1].y : Pcrds[1].y;
		Icrds[i-1].y = c1 + c2; /*center is average of the 3 coords*/
		Icrds[i-1].y *= 0.333333333333;
	}
}

/*Function to calculate strain energy at pentagonal intergration points*/
void GaPentStrain(Shape *Pshp, int pnum, Coord *Pcrds, int *tnodes, double *rhs1, double *rhs2, double alpha, double e11, double v12, 
				  double g33, int Gcount, Gstn *Genergy, double cx, double cy, int num, double *Eenergy, double wgt)
{
	int temp,i,j;		/*Incrementors etc*/
	Coord *tdisp1;			/*Temp array to read in primary displacement values*/
	tdisp1 = malloc(5 * sizeof(Coord));
	Coord *tdisp2;			/*Temp array to read in adjoint displacement values*/
	tdisp2 = malloc(5 * sizeof(Coord));
	double ex1,ey1,txy1,ex2,ey2,txy2,Stn,ftemp,ftemp2;	/*Varaibles for strains*/
	double atemp;	/*variable for sub-domain area*/
	int p5 = pnum * 5;	/*constant to access Pshp properly*/
	int p6 = pnum * 6;	/*constant to access Pcrds properly*/
	
	/*first work out global gauss point co-ordinates and store assosiated area*/
	Coord *Icrds;	/*array to store intergration point co-ordinates*/
	Icrds = malloc(5 * sizeof(Coord));
	Coord *pnt;
	pnt = &Pcrds[p6];
	Pint(pnt,Icrds); /*calls function that calculates intergration point co-ordinates*/
	double lx = Pcrds[p6].x;
	double ly = Pcrds[p6].y; /*read in local center co-ordiantes*/
	
	/*for(i=0;i<5;i++)
	{
		printf("\nInt Point %i: x=%fl, y=%fl",i+1,Icrds[i].x,Icrds[i].y);
	}*/
	
	/*now store the co-ordinates and associated area*/
	for(i=0;i<5;i++)
	{
		Genergy[Gcount + i].x = cx + lx + Icrds[i].x;
		Genergy[Gcount + i].y = cy + ly + Icrds[i].y;
		Genergy[Gcount + i].a = alpha;
	}
	
	/*Populate element displacement arrays*/
	for(i=0;i<5;i++)
	{
		temp = tnodes[i] * 2;
		tdisp1[i].x = rhs1[temp];			/*X primary displacement of node i*/
		tdisp1[i].y = rhs1[temp + 1];		/*Y primary displacement of node i*/
		/*printf("\ntdisp[%i].x: %fl, .y: %fl",i,tdisp[i].x,tdisp[i].y);*/
		tdisp2[i].x = rhs2[temp];			/*X adjoint displacement of node i*/
		tdisp2[i].y = rhs2[temp + 1];		/*Y adjoint displacement of node i*/
	}
	
	/*for(i=0;i<5;i++)
	{
		printf("\nfor%i: A=%fl, B=%fl, C=%fl, D=%fl",i+1,Pshp[p5+i].a,Pshp[p5+i].b,Pshp[p5+i].c,Pshp[p5+i].d);
	}*/		
	
	for(i=0;i<5;i++) /*For all intergration points*/
	{
		ex1 = 0.0;
		ey1 = 0.0;
		txy1 = 0.0;
		ex2 = 0.0;
		ey2 = 0.0;
		txy2 = 0.0; /*Initialize strains to zero*/
		
		for(j=0;j<5;j++) /*Calculate strains at the point*/
		{
			ftemp = Pshp[p5+j].d * Icrds[i].y;
			ftemp += Pshp[p5+j].b;
			ex1 += ftemp * tdisp1[j].x;
			ex2 += ftemp * tdisp2[j].x;
			
			ftemp = Pshp[p5+j].d * Icrds[i].x;
			ftemp += Pshp[p5+j].c;
			ey1 += ftemp * tdisp1[j].y;
			ey2 += ftemp * tdisp2[j].y;
			
			ftemp = Pshp[p5+j].d * Icrds[i].y;
			ftemp += Pshp[p5+j].b;
			txy1 += ftemp * tdisp1[j].y;
			txy2 += ftemp * tdisp2[j].y;
			
			ftemp = Pshp[p5+j].d * Icrds[i].x;
			ftemp += Pshp[p5+j].c;
			txy1 += ftemp * tdisp1[j].x;
			txy2 += ftemp * tdisp2[j].x;
		}
		
		/*printf("\nex=%fl, ey=%fl, txy=%fl",ex,ey,txy);
		
		/*Calculate Strain Energy at the integration point*/
		Stn = e11 * ((ex1 * ex2) + (ey1 * ey2));
		Stn += v12 * ( (ex1 * ey2) + (ex2 * ey1) );
		Stn += g33 * txy1 * txy2;
		/*Stn *= 0.5;*/
		
		Genergy[Gcount + i].u = Stn; /*store intergration point strain energy value*/
		/*printf("\nStn%i = %fl",i+1,Stn);*/
		j = (i == 4) ? 0 : i+1;
		atemp = 0.5 * fabs((Pcrds[p6 +1 +i].x * Pcrds[p6 +1 +j].y) + (Pcrds[p6 +1 +i].y * Pcrds[p6 +1 +j].x)); /*sub-domain triangular area*/
		Eenergy[num] += wgt * Stn * atemp; /*update element strain energy*/
	}
	
	free(tdisp1);
	free(tdisp2);
	free(Icrds);
}

/*Function to calculate strain energy values for an element at 4 gauss points*/
void GaINstrain(int *tnodes, double *rhs1, double *rhs2, double alpha, double e11, double v12, double g33, double h, int Gcount,
				Gstn *Genergy,int NumNodes, Coord NodeCoord[NumNodes], int num, double *Eenergy, double Earea, double wgt, double hxz, double hz)
{
	/*Since strain energy is a directionless property we will calculate this as a local property using local X axis*/
	int i,j,temp;	/*incrementors etc*/
	double *Edisp1;	/*Primary Element displacement array*/
	Edisp1 = malloc(24 * sizeof(double));
	double *Edisp2;	/*Adjoint Element displacement array*/
	Edisp2 = malloc(24 * sizeof(double));
	double stnX1,stnY1,stnXY1,stnX2,stnY2,stnXY2,Stn,stnXY1a,stnXY1b,stnXY2a,stnXY2b;	/*Variabes for strain tensors*/
	double gax, gay, ftemp1, ftemp2; /*Variables for gauss point values*/
	double h2 = 0.5 * h;
	double h2xz = 0.5 * hxz;
	double theta, omega, gammer,atemp,dispX,dispZ,Mag,Trans;
	
	/*work out global co-ordinates of element center, from co-ords of node 1 and element edge length*/
	double cx = NodeCoord[tnodes[0]].xz + h2xz;
	double cy = NodeCoord[tnodes[0]].y + h2;
	/*printf("\nK51");
	printf("\nK51");*/
	
	/*first work out global gauss point co-ordinates and store assosiated area*/
	for(i=0;i<4;i++)
	{
		gax = ga * (double)Po[i].x * h2xz;
		gay = ga * (double)Po[i].y * h2;
		Genergy[Gcount + i].x = cx + gax;
		Genergy[Gcount + i].y = cy + gay;
		Genergy[Gcount + i].a = alpha;
	}
	atemp = (hz/h);
	theta = atan(atemp);
	/*printf("\n Theta = %f", theta);*/

	/*Populate element displacement arrays*/
	for(i=0;i<4;i++)
	{
		j = 6 * i;
		temp = tnodes[i] * 6;
		dispX = rhs1[temp];				/*caclulate the angel needed to get the local x displacement*/
		dispZ = rhs1[temp+2];
		/*atemp = (dispX == 0) ? 999999999:(dispZ/dispX);
		omega = atan(atemp);
		gammer = (omega - theta);
		atemp = ((dispZ*dispZ)+(dispX*dispX)); /*Magnitude of the x and z displacements.*
		Mag = (atemp == 0) ? 0:sqrt(atemp);
		Trans = cos(gammer);*/
		
		Edisp1[j] = rhs1[temp];		/*Primary local X displacement of node i*/
		Edisp1[j + 1] = rhs1[temp + 1];	/*Primary Y displacement of node i*/
		Edisp1[j + 2] = rhs1[temp + 2]; /*Primary Z displacement of node i*/
		Edisp1[j + 3] = rhs1[temp + 3];		/*Primary rotX displacement of node i*/
		Edisp1[j + 4] = rhs1[temp + 4];	/*Primary rotY displacement of node i*/
		Edisp1[j + 5] = rhs1[temp + 5]; /*Primary rotZ displacement of node i*/



		Edisp2[j] = rhs2[temp];		/*Primary local X displacement of node i*/
		Edisp2[j + 1] = rhs2[temp + 1];	/*Adjoint Y displacement of node i*/
		Edisp2[j + 2] = rhs1[temp + 2]; /*Adjoint Z displacement of node i*/
		Edisp2[j + 3] = rhs1[temp + 3];		/*Adjoint rotX displacement of node i*/
		Edisp2[j + 4] = rhs1[temp + 4];	/*Adjoint rotY displacement of node i*/
		Edisp2[j + 5] = rhs1[temp + 5]; /*Adjoint rotZ displacement of node i*/
		
	}
	
	/*Evaluate strain energy at each gauss point*/
	for(i=0;i<4;i++)
	{
		stnX1 = 0.0;
		stnY1 = 0.0;
		stnXY1 = 0.0;
		stnX2 = 0.0;
		stnY2 = 0.0;
		stnXY2 = 0.0;	/*Initialize strain to zero for each point*/
		
		gax = ga * (double)Po[i].x;
		gay = ga * (double)Po[i].y;	/*Values at gauss points*/
		
		/*Calculate X,Y and shear strain tensors at the gauss point*/
		for(j=0;j<4;j++)
		{
			temp = j * 6; /* X dof */
			/*compute strain displacement matrix entries*/
			ftemp1 = (double)Po[j].x * (1.0 + ((double)Po[j].y * gay));
			ftemp2 = (double)Po[j].y * (1.0 + ((double)Po[j].x * gax));
			
			/*priamry strains*/
			stnX1 += ftemp1 * Edisp1[temp];
			stnY1 += ftemp2 * Edisp1[temp + 1];
			stnXY1 += ftemp1 * Edisp1[temp + 1];
			stnXY1 += ftemp2 * Edisp1[temp];
			/*printf("\n StnXY1 = %f + %f = %f", ftemp1*Edisp1[temp+1], ftemp2*Edisp1[temp], ftemp2*Edisp1[temp]+ftemp1*Edisp1[temp+1]);*/
			/*printf("\t Ftemp = %f", ftemp2);*/
			
			/*adjoint strains*/
			stnX2 += ftemp1 * Edisp2[temp];
			stnY2 += ftemp2 * Edisp2[temp + 1];
			stnXY2 += ftemp1 * Edisp2[temp + 1];
			stnXY2 += ftemp2 * Edisp2[temp];
		}
		
		/*stnXY1 = fabs(stnXY1a)+fabs(stnXY1b);
		stnXY2 = fabs(stnXY2a)+fabs(stnXY2b);*/
		
		/*printf("\nFor Point %i: ex=%f, ey=%f, txy=%f",i+1,stnX1,stnY1,stnXY1);
		printf("\nFor Point %i: ex=%f, ey=%f, txy=%f",i+1,stnX2,stnY2,stnXY2);*/
		
		/*Calculate Strain Energy at the gauss point*/
		Stn = e11 * ((stnX1 * stnX2) + (stnY1 * stnY2));
		Stn += (v12 * ( (stnX1 * stnY2) + (stnX2 * stnY1) ));
		Stn += g33 * stnXY1 * stnXY2;
		Stn /= 4.0 * Earea; /*compliance tye senstivity, not strain energy*/
		Stn *= alpha; /*multiply by area ratio, Arat = 1.0 for IN elements*/
		Genergy[Gcount + i].u = Stn;	/*store strain energy value for later*/
		Eenergy[num] += wgt * Stn * Earea * 0.25; /*update element strain energy*/
	}
	/*printf("\n");*/
	free(Edisp1);
	free(Edisp2);
}


/*Function to calculate strain energy values for an element at 4 gauss points, for minderlin elements*/
void GaINstrainV2(int *tnodes, double *rhs1, double *rhs2, double alpha, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double h, int Gcount, Gstn *Genergy,int NumNodes, Coord NodeCoord[NumNodes], int num, double *Eenergy, double Earea, double wgt, double hxz, double hz, double *theta, double *Trans, int numply, double plyt,int n)
{
	/*Since strain energy is a directionless property we will calculate this as a local property using local X axis*/
	int i,j,temp,k;	/*incrementors etc*/
	double *UX,*UY,*UZ,*RX,*RY,*RZ;
	UX = malloc(4*sizeof(double));
	UY = malloc(4*sizeof(double));
	UZ = malloc(4*sizeof(double));
	RX = malloc(4*sizeof(double));
	RY = malloc(4*sizeof(double));
	RZ = malloc(4*sizeof(double));
	double Stn;
	Stn = 0;
	double gax, gay, ftemp1, ftemp2; /*Variables for gauss point values*/
	double h2 = 0.5 * h;
	double h2xz = 0.5 * hxz;
	double *K, *KE;
	K = malloc(576*sizeof(double));
	KE = malloc(576*sizeof(double));
	
	/*work out global co-ordinates of element center, from co-ords of node 1 and element edge length*/
	double cx = NodeCoord[tnodes[0]].xz + h2xz;
	double cy = NodeCoord[tnodes[0]].y + h2;
	/*printf("\nK51");
	printf("\nK51");*/
	/*if((n==19)||(n==0)){printf("\n N = %i", n);}*/
	/*first work out global gauss point co-ordinates and store assosiated area*/
	for(i=0;i<4;i++)
	{
		gax = ga * (double)Po[i].x * h2xz;
		gay = ga * (double)Po[i].y * h2;
		Genergy[Gcount + i].x = cx + gax;
		Genergy[Gcount + i].y = cy + gay;
		Genergy[Gcount + i].a = alpha;
	}
	/*Populate element displacement arrays*/
	for(i=0;i<4;i++)
	{
		j = 6 * i;
		temp = tnodes[i] * 6;
		
		UX[i] = rhs1[temp];		/*Primary local X displacement of node i*/
		UY[i] = rhs1[temp + 1];	/*Primary Y displacement of node i*/
		UZ[i] = rhs1[temp + 2]; /*Primary Z displacement of node i*/
		RX[i] = rhs1[temp + 3];		/*Primary rotX displacement of node i*/
		RY[i] = rhs1[temp + 4];	/*Primary rotY displacement of node i*/
		RZ[i] = rhs1[temp + 5]; /*Primary rotZ displacement of node i*/
		/*if((n==19)||(n==0)){printf("\n disp (%f, %f, %f, %f, %f, %f)", UX[i], UY[i], UZ[i], RX[i], RY[i], RZ[i]);}*/
	}
	/*Evaluate strain energy at each gauss point*/
	for(i=0;i<4;i++)
	{
		Stn = 0;		
		gax = ga * (double)Po[i].x;
		gay = ga * (double)Po[i].y;	/*Values at gauss points*/
		/*printf("\n\ngp = %i", i);
		printf("\n\nhxz = %f", hxz);*/			
		/*Calculate the element stiffness matrix around this gauss piont*/
		CompositeKEGauss(i,e1,e2,v12,v21,g12,g23,g13, KE, theta, numply, plyt, h, hxz);
		/*Transfrom the Stiffness matrix to global coordintes*/
		Transform(KE, Trans, K, n);


/*for(k=0;k<24;k++)
{
	printf("\n");
	for(j=0;j<24;j++)
	{
		printf("%0.2f  ", K[j+k*24]);

	}
}*/
		/*if((n==19)||(n==0)){printf("\n%i",i);}*/
		/*Strain energy at the gauss point is U'Kguass U and the mean of all the gauss point strain energy is the element strain energy*/
		for(j=0;j<4;j++)
		{
			k = j*6;
		
			Stn += UX[j] * (UX[0]*K[0+k] + UY[0]*K[24+k] + UZ[0]*K[48+k] + RX[0]*K[72+k] + RY[0]*K[96+k] + RZ[0]*K[120+k] + UX[1]*K[144+k] + UY[1]*K[168+k] + UZ[1]*K[192+k] + RX[1]*K[216+k] + RY[1]*K[240+k] + RZ[1]*K[264+k] + UX[2]*K[288+k] + UY[2]*K[312+k] + UZ[2]*K[336+k] + RX[2]*K[360+k] + RY[2]*K[384+k] + RZ[2]*K[408+k] + UX[3]*K[432+k] + UY[3]*K[456+k] + UZ[3]*K[480+k] + RX[3]*K[504+k] + RY[3]*K[528+k] + RZ[3]*K[552+k]);
			/*if((n==19)||(n==0)){printf("\nStnA = %f", Stn);}*/

			Stn += UY[j] * (UX[0]*K[1+k] + UY[0]*K[25+k] + UZ[0]*K[49+k] + RX[0]*K[73+k] + RY[0]*K[97+k] + RZ[0]*K[121+k] + UX[1]*K[145+k] + UY[1]*K[169+k] + UZ[1]*K[193+k] + RX[1]*K[217+k] + RY[1]*K[241+k] + RZ[1]*K[265+k] + UX[2]*K[289+k] + UY[2]*K[313+k] + UZ[2]*K[337+k] + RX[2]*K[361+k] + RY[2]*K[385+k] + RZ[2]*K[409+k] + UX[3]*K[433+k] + UY[3]*K[457+k] + UZ[3]*K[481+k] + RX[3]*K[505+k] + RY[3]*K[529+k] + RZ[3]*K[553+k]);
			/*if((n==19)||(n==0)){printf("\nStnB = %f", Stn);}*/

			Stn += UZ[j] * (UX[0]*K[2+k] + UY[0]*K[26+k] + UZ[0]*K[50+k] + RX[0]*K[74+k] + RY[0]*K[98+k] + RZ[0]*K[122+k] + UX[1]*K[146+k] + UY[1]*K[170+k] + UZ[1]*K[194+k] + RX[1]*K[218+k] + RY[1]*K[242+k] + RZ[1]*K[266+k] + UX[2]*K[290+k] + UY[2]*K[314+k] + UZ[2]*K[338+k] + RX[2]*K[362+k] + RY[2]*K[386+k] + RZ[2]*K[410+k] + UX[3]*K[434+k] + UY[3]*K[458+k] + UZ[3]*K[482+k] + RX[3]*K[506+k] + RY[3]*K[530+k] + RZ[3]*K[554+k]);
			/*if((n==19)||(n==0)){printf("\nStnC = %f", Stn);	}*/

			Stn += RX[j] * (UX[0]*K[3+k] + UY[0]*K[27+k] + UZ[0]*K[51+k] + RX[0]*K[75+k] + RY[0]*K[99+k] + RZ[0]*K[123+k] + UX[1]*K[147+k] + UY[1]*K[171+k] + UZ[1]*K[195+k] + RX[1]*K[219+k] + RY[1]*K[243+k] + RZ[1]*K[267+k] + UX[2]*K[291+k] + UY[2]*K[315+k] + UZ[2]*K[339+k] + RX[2]*K[363+k] + RY[2]*K[387+k] + RZ[2]*K[411+k] + UX[3]*K[435+k] + UY[3]*K[459+k] + UZ[3]*K[483+k] + RX[3]*K[507+k] + RY[3]*K[531+k] + RZ[3]*K[555+k]);
			/*if((n==19)||(n==0)){printf("\nStnD = %f", Stn);}*/
	
			Stn += RY[j] * (UX[0]*K[4+k] + UY[0]*K[28+k] + UZ[0]*K[52+k] + RX[0]*K[76+k] + RY[0]*K[100+k] + RZ[0]*K[124+k] + UX[1]*K[148+k] + UY[1]*K[172+k] + UZ[1]*K[196+k] + RX[1]*K[220+k] + RY[1]*K[244+k] + RZ[1]*K[268+k] + UX[2]*K[292+k] + UY[2]*K[316+k] + UZ[2]*K[340+k] + RX[2]*K[364+k] + RY[2]*K[388+k] + RZ[2]*K[412+k] + UX[3]*K[436+k] + UY[3]*K[460+k] + UZ[3]*K[484+k] + RX[3]*K[508+k] + RY[3]*K[532+k] + RZ[3]*K[556+k]);
			/*if((n==19)||(n==0)){printf("\nStnE = %f", Stn);}*/
	
			Stn += RZ[j] *(UX[0]*K[5+k] + UY[0]*K[29+k] + UZ[0]*K[53+k] + RX[0]*K[77+k] + RY[0]*K[101+k] + RZ[0]*K[125+k] + UX[1]*K[149+k] + UY[1]*K[173+k] + UZ[1]*K[197+k] + RX[1]*K[221+k] + RY[1]*K[245+k] + RZ[1]*K[269+k] + UX[2]*K[293+k] + UY[2]*K[317+k] + UZ[2]*K[341+k] + RX[2]*K[365+k] + RY[2]*K[389+k] + RZ[2]*K[413+k] + UX[3]*K[437+k] + UY[3]*K[461+k] + UZ[3]*K[485+k] + RX[3]*K[509+k] + RY[3]*K[533+k] + RZ[3]*K[557+k]);
			/*if((n==19)||(n==0)){printf("\nStnF = %f", Stn);}*/
		}
		
		Stn *= alpha; /*multiply by area ratio, Arat = 1.0 for IN elements*/
		/*printf("\nGstrain %i = %fl",i+1,Stn);*/
		Genergy[Gcount + i].u = Stn;	/*store strain energy value for later*/
		Eenergy[num] += wgt * Stn * Earea * 0.25; /*update element strain energy*/
	}
	/*if((n==19)||(n==0)){printf("\n\tStn = %f", Stn);}
	if((n==19)||(n==0)){printf("\n\tEenergy = %f", Eenergy[num]);}*/
	/*printf("\n");*/
	free(UX);
	free(UY);
	free(UZ);
	free(RX);
	free(RY);
	free(RZ);
	free(K);
	free(KE);
}

/*Function to calculate Triangular strain energy value (at center)*/
void GaTriStrain(double alpha, double e11, double v12, double g33, int *tnodes, Coord *tcords, double *rhs1,
				 double *rhs2, int Gcount, Gstn *Genergy,int num, double *Eenergy, double Earea, double wgt)
{
	Coord *Ni;			/*Array for CST shape constants*/
	Ni = malloc(3 * sizeof(Coord));
	Coord *tdisp1;			/*Array to store primary element displacements*/
	tdisp1 = malloc(3 * sizeof(Coord));
	Coord *tdisp2;			/*Array to store adjoint element displacements*/
	tdisp2 = malloc(3 * sizeof(Coord));
	double ex1,ey1,txy1,ex2,ey2,txy2,Stn;	/*Strain variables*/
	int temp,i;
	
	double Tarea = alpha * Earea; /*triangular area*/
	
	/*first work out global gauss point co-ordinates (i.e. element center) and store assosiated area*/
	ex1 = 0.0;
	ey1 = 0.0; /*here these are used for calculating the element center co-ordinates*/
	for(i=0;i<3;i++)
	{
		ex1 += tcords[i].x;
		ey1 += tcords[i].y;
	}
	
	/*store center co-ordinates and element area*/
	Genergy[Gcount].x = 0.333333333333 * ex1;
	Genergy[Gcount].y = 0.333333333333 * ey1;
	Genergy[Gcount].a = alpha;

	/*printf("\nTriArea = %fl",area);*/

	/*Calculate CST shape constants*/
	Ni[0].x = tcords[2].x - tcords[1].x;
	Ni[1].x = tcords[0].x - tcords[2].x;
	Ni[2].x = tcords[1].x - tcords[0].x;
	
	Ni[0].y = tcords[1].y - tcords[2].y;
	Ni[1].y = tcords[2].y - tcords[0].y;
	Ni[2].y = tcords[0].y - tcords[1].y;
	
	/*for(i=0;i<3;i++)
	{
		printf("\nNi[%i].x = %fl, .y = %fl",i,Ni[i].x,Ni[i].y);
	}*/
	
	/*Populate element displacement array*/
	for(i=0;i<3;i++)
	{
		temp = tnodes[i] * 2;
		tdisp1[i].x = rhs1[temp];			/*primary X displacement of node i*/
		tdisp1[i].y = rhs1[temp + 1];		/*primary Y displacement of node i*/
		/*printf("\ntdisp[%i].x: %fl, .y: %fl",i,tdisp[i].x,tdisp[i].y);*/
		tdisp2[i].x = rhs2[temp];			/*adjoint X displacement of node i*/
		tdisp2[i].y = rhs2[temp + 1];		/*adjoint Y displacement of node i*/
	}
	
	ex1 = 0.0;
	ey1 = 0.0;
	txy1 = 0.0;
	ex2 = 0.0;
	ey2 = 0.0;
	txy2 = 0.0; /*Initialize strains to zero*/
	
	/*Now calculate element stain values*/
	for(i=0;i<3;i++)
	{
		/*primary*/
		ex1 += Ni[i].y * tdisp1[i].x;
		ey1 += Ni[i].x * tdisp1[i].y;
		txy1 += Ni[i].x * tdisp1[i].x;
		txy1 += Ni[i].y * tdisp1[i].y;
		/*adjoint*/
		ex2 += Ni[i].y * tdisp2[i].x;
		ey2 += Ni[i].x * tdisp2[i].y;
		txy2 += Ni[i].x * tdisp2[i].x;
		txy2 += Ni[i].y * tdisp2[i].y;
	}
	
	/*printf("\nTri: ex=%fl ey=%fl, txy=%fl",ex,ey,txy);*/
	
	/*Set A strain energy field constant to CST strain energy*/
	Stn = e11 * ((ex1 * ex2) + (ey1 * ey2));
	Stn += v12 * ( (ex1 * ey2) + (ex2 * ey1) );
	Stn += g33 * txy1 * txy2;
	Stn /= 4.0 * Tarea * Tarea; /*compliance tye senstivity, not strain energy*/
	
	/*printf("\nTriStrain = %fl",Stn);*/
	
	/*store the strain energy value*/
	Genergy[Gcount].u = Stn;
	Eenergy[num] += wgt * Stn * Tarea; /*update element strain energy*/
	
	free(Ni);
	free(tdisp1);
	free(tdisp2);
}
