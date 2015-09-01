/*
 *  EMatrix.h
 *  
 *  Created by Peter Dunning on 26/02/2007.
 *
 */
 
#include "ls_types.h"

/*The following 2 functions calculate the values in each 2 x 2 block in the Element Stiffness Matrix 
	These functions integrate the matrix exactly. Trust me it works
	--------------------------------------------------------------------------------------------------*/
	
double Aodd(int i,int j,double e,double g);

double Aeven(int i,int j,double e,double g);

/*---------------------Iso Matrix entry calculation funcions--------------*/
double Ievenv1(int i, int j, double a1, double a2, int flag);

/*New function that calculates isoparametric matrices correctly*/
void IsoMatrix(double *Iso, double c, double d, double g, double a1, double a2, double b1, double b2, double h);

/*---------------------------------IN Element Stiffness Matirx------------------------------*/
void KEMatrix(double *KE, double c, double d, double g);

/*------------------Constant Strain Triangle Stiffness martrix-------------------*/
void Triangle(int Enum, short *NodeStat, int Xi, int Yj, int elemX, int elemY, Elem Number[elemX][elemY], double *Tri, 
				double a1, double b1, double e1, double v1, double g1, double h);

/*Function to calculate a pentagonal element matrix using the least squares method*/
void Pentagon(int pnum, Coord *Pcrds, double *Pent, double e1, double v1, double g1, Shape *Pshp);
