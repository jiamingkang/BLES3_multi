/*
 *  StrainTest.h
 *  
 *  Created by Peter Dunning on 15/01/2010.
 *
 */
 
#include "ls_types.h"

/*Function that calculates status of elements around a central element, 1 = IN or NIO, 0 = out*/
void LocalStatus(int Xi, int Yj, short StatLoc[3][3], int elemX, int elemY, short *ElemStat, Elem Number[elemX][elemY]);

/*Function that calculates area ratio of elements around a central element*/
void LocalStatus2(int Xi, int Yj, double StatLoc[3][3], int elemX, int elemY, double *alpha, Elem Number[elemX][elemY], double area);

/*Function to select method and variables for testing strain energy methods*/
void StrainTest(int elemX, int elemY, Elem Number[elemX][elemY], short *ElemStat, short *NodeStat, double *alpha, double e11, double v12, double g33,
					double h, Aux *auxNode, double *Nenergy, Shape *Pshp, Coord *Pcrds, double *rhs, int NumNodes, Coord NodeCoord[NumNodes], 
						int Acount, int NIOtotal, int Itotal, int tCount, int qCount, int pCount);
						
/*Function to calculate quadrilateral elememt nodal strain energy values - directly*/	
void IsoStrain2(short GridStat[4], int tnodes[4], Coord tcords[4], double *rhs, double e11, double v12, double g33, short StatLoc[3][3], 
				double StatLoc2[3][3], double *Nenergy, int mFlag);

/*Function to calculate nodal strain energy values for a triangular CST element*/
void TriStrain(double e11, double v12, double g33, int tnodes[5], Coord tcords[5], double *rhs, short StatLoc[3][3], 
				double StatLoc2[3][3], short GridStat[4], double *Nenergy, int mFlag);						

/*Function to calculate nodal strain energy values for an IN element*/
void INstrain(int tnodes[4], short StatLoc[3][3], double StatLoc2[3][3], double *rhs, double e11, double v12, double g33, double h, double *Nenergy, int mFlag);

/*Function to calculate nodal strain energy values for a pentagonal element*/
void PentStrain(Shape *Pshp, int pnum, Coord *Pcrds, int tnodes[5], short StatLoc[3][3], double StatLoc2[3][3], double *rhs, 
					double e11, double v12, double g33, double *Nenergy, short GridStat[4], int mFlag);

/*----------------------------------Functions for intergration point strain energy values----------------------------------*/

/*Function to calculate strain energy values for an IN element at 4 gauss points*/
void GaINstrainTest(int tnodes[4], double *rhs, double area, double e11, double v12, double g33, double h, int Gcount, Gstn *Genergy,
					int NumNodes, Coord NodeCoord[NumNodes]);
					
/*Function to calculate Isoparametric Gauss point strain energy values*/	
void GaIsoStrainTest(int tnodes[4], Coord tcords[4], double *rhs, double area, double e11, double v12, double g33, int Gcount, Gstn *Genergy);

/*Function to calculate Triangular strain energy value (at center)*/
void GaTriStrainTest(double area, double e11, double v12, double g33, int tnodes[5], Coord tcords[5], double *rhs, int Gcount, Gstn *Genergy);

/*Function to calculate strain energy at pentagonal intergration points (center of triangular sub-domains)*/
void GaPentStrainTest(Shape *Pshp, int pnum, Coord *Pcrds, int tnodes[5], double *rhs, double area, 
					double e11, double v12, double g33, int Gcount, Gstn *Genergy, double cx, double cy);
					
/*Function to calculate integration point coords for pentagonal matrix
	(center of triangular sub-domains)*/
void Pint2(int p6, Coord *Pcrds, Coord *Icrds);

/*----------------------------------Functions for center point strain energy values----------------------------------*/
/*NB: intergration point for a triangle is at the center*/

/*Function to calculate strain energy value for an IN element at element center*/
void CeINstrainTest(int tnodes[4], double *rhs, double area, double e11, double v12, double g33, double h, int Gcount, Gstn *Genergy,
					int NumNodes, Coord NodeCoord[NumNodes]);
					
/*Function to calculate Isoparametric Quadrilateral center strain energy value*/	
void CeIsoStrain(int tnodes[4], Coord tcords[4], double *rhs, double area, double e11, double v12, double g33, int Gcount, Gstn *Genergy);

/*Function to calculate strain energy at pentagon center*/
void CePentStrain(Shape *Pshp, int pnum, Coord *Pcrds, int tnodes[5], double *rhs, double area, 
					double e11, double v12, double g33, int Gcount, Gstn *Genergy, double cx, double cy);
					
/*----------------------------------Functions for least squares fitting for strain energy values----------------------------------*/

/*Function that calculates the strain energy of a node by a least squares (bi-linear order) filter of near-by gauss points*/
double Lstrain2(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag);

/*Function that calculates the strain energy of a node by a least squares (2nd order) filter of near-by gauss points*/
double Lstrain22(double nx, double ny, double rad, int gpoints, Gstn *Genergy, int wFlag);
