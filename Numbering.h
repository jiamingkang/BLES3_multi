/*
 *  Numbering.h
 *  
 *  Created by Peter Dunning on 26/02/2007.
 *
 */

#include "ls_types.h"

/*Function that numbers all elements and nodes in the FG domain*/
void Numbering(int elemX,int elemY,Elem Number[elemX][elemY]);

/*Node Co-ordinate calculation function*/
void Coordinates(int elemX, int elemY, int NumNodes, double hx, double hy, Elem Number[elemX][elemY], Coord NodeCoord[NumNodes]);

/*Function that orders node numbers into a 2D based on their relative positions*/
void NodeNums2(double h, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int NumNodes, Coord NodeCoord[NumNodes]);
