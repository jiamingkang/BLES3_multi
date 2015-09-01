/*
 *  FixedGrid.h
 *  
 *  Created by Peter Dunning on 26/02/2007.
 *
 */

#include "ls_types.h"

/*Function to calcualte lsf around a rectangular hole*/
void RectHole(int NumNodes, Coord NodeCoord[NumNodes], int NumRect, Coord *Rect, double *lsf, int elemX, int elemY, double hxz[elemX], double maxXZ);

/*Function that calculates an intersection point from determinanats*/
double TwoDpoint(double c, double d, double yb, double y1, double x1);

/*Function that caculates the area of any Polygon*/
double PolyArea(int N, Coord *point);
double PolyArea2(double *lsf, int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], short *NodeStat, short *ElemStat, double *alpha0, double *hxz, double h, double *alpha);
				
/*Function to determine if an auxillary node is within a fixed boundary*/
void AuxBound2(int bsn, Aux *Afix, int NumNodes, int Ntot, Aux *auxNode, int NIOtotal, int *TotFix2, int *FixDofs2);

/*Function to add some auxillary nodes around the boundary*/
void auxAdd(double h, int NIOn, int Xi, int Yj, int elemX, int elemY, short *ElemStat, Elem Number[elemX][elemY], short *NodeStat,
			double *lsf, int *Acount, Aux *auxNode, int NumNodes, Coord NodeCoord[NumNodes], int *NIOnums, double hxz[elemX],double hz[elemX]);

/*Function to count the number of various element types*/			
void elemType(int NumElem, short *ElemStat, int *Itotal, int *Ototal, int *NIOtotal, int *tCount, int *qCount, int *pCount);

/*Function to read IN and auxillary node numbers for an NIO element*/				
void readNodes(int nNod, int NIOnum, int i, int j, int elemX, int elemY, short *NodeStat, Elem Number[elemX][elemY], int *tnodes, Aux *auxNode);
				
/*Function to read IN and auxillary node co-ordinates for an NIO element*/
void readCoord(int nNod, int N2, int i, int j, int elemX, int elemY, short *NodeStat, Elem Number[elemX][elemY], Aux *auxNode, 
				Coord *tcords, int NumNodes, Coord NodeCoord[NumNodes], int *tnodes);

/*Edge calculation function for auxillary node formulation*/
void EdgeCalc3(int Enum, int NIOnum, short *NodeStat, int NumNodes, int i, int j, int elemX, int elemY, Elem Number[elemX][elemY], 
				double h, Coord NodeCoord[NumNodes], Aux *auxNode, double *Edges);
				
/*Function to calculate geometric and center co-ordinates for a pentagonal element*/
void PentCoord(int Enum, double h, double a1, double a2, double b1, double b2, short *NodeStat, int i, int j, 
				int elemX, int elemY, Elem Number[elemX][elemY], int pnum, Coord *Pcrds);
				
/*Function for working out center co-ordinates relative to IN element center*/
double Cent(double a, double b, double h, double A, short p, short q);

/*Function to find the nearest grid node number to a set of co-ordinates*/
int closeNode(double h, double xp, double yp, int NumNodes, Coord NodeCoord[NumNodes]);
