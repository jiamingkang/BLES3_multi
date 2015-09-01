/*
 *  Strain.h
 *  
 *  Created by Peter Dunning on 20/05/2008.
 *
 */
 
#include "ls_types.h"
		
/*Function that calculates the strain energy of a node by a least squares (bi-linear order) filter of near-by gauss points*/
double LstrainV2(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag);

/*Function that calculates the strain energy of a node by a least squares (2nd order) filter of near-by gauss points*/
double LstrainV22(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag);

double LstrainV23(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag);

double LstrainV24(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag);
					
/*Function to calculate Isoparametric Gauss point strain energy values*/	
void GaIsoStrain(int *tnodes, Coord *tcords, double *rhs1, double *rhs2, double alpha, double e11, double v12, double g33, int Gcount, Gstn *Genergy,
				 int num, double *Eenergy, double Earea, double wgt);
					
/*Function to calculate integration point coords for pentagonal strain energy calculation*/
void Pint(Coord Pcrds[6], Coord Icrds[5]);

/*Function to calculate strain energy at pentagonal intergration points*/
void GaPentStrain(Shape *Pshp, int pnum, Coord *Pcrds, int *tnodes, double *rhs1, double *rhs2, double alpha, double e11, double v12, 
					double g33, int Gcount, Gstn *Genergy, double cx, double cy, int num, double *Eenergy, double wgt);
					
/*Function to calculate strain energy values for an element at 4 gauss points*/
void GaINstrain(int *tnodes, double *rhs1, double *rhs2, double alpha, double e11, double v12, double g33, double h, int Gcount,
				Gstn *Genergy,int NumNodes, Coord NodeCoord[NumNodes], int num, double *Eenergy, double Earea, double wgt, double hxz, double hz);

void GaINstrainV2(int *tnodes, double *rhs1, double *rhs2, double alpha, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double h, int Gcount, Gstn *Genergy,int NumNodes, Coord NodeCoord[NumNodes], int num, double *Eenergy, double Earea, double wgt, double hxz, double hz, double *theta, double *Trans, int numply, double plyt,int n);
					
/*Function to calculate Triangular strain energy value (at center)*/
void GaTriStrain(double alpha, double e11, double v12, double g33, int *tnodes, Coord *tcords, double *rhs1,
				 double *rhs2, int Gcount, Gstn *Genergy,int num, double *Eenergy, double Earea, double wgt);
