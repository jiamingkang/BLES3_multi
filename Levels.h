/*
 *  Levels.h
 *  
 *
 *  Created by Peter Dunning on 17/02/2009.
 *
 */
 
#include "ls_types.h"

/*Function to calculate active and mine nodes for the narrow band method*/
void NarBand(int NumNodes, double h, double lBand, int *NumAct, int *NumMine, double *lsf, int *active, int *mine, int NumFix, int *fixed, int *activeHoles, int *NumActHoles);
				
/*Function to calculates boundary length and strain enegy associated with boundary nodes*/
void LamCalc3a(double h, int elemX, int elemY, Elem Number[elemX][elemY], short *ElemStat, short *NodeStat, 
				Aux *auxNode, double *Nenergy, int NumNodes, Coord NodeCoord[NumNodes], int Ntot, Aux *Lbound, double hxz[elemX]);
					
void ArrayCopy(int length, double *original, double *copy);

double LamCalc4(double lambda, double dt, int NumNodes, double *Nenergy, double *lsf, double *lsfHole, int NumFix, int *fixed, 
					int NumAct, int *active, double h, double hbar, double HoleLim, double Voltot, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], short *NodeStat, int elemX, double hxz[elemX]);

/*function to calcualte volume reduction due to hole insertion*/
double LamCalc4a(double lambda, double dt, int NumNodes, double *Nenergy, double *lsf, double *lsfHole, int NumFix, int *fixed, 
					int NumAct, int *active, double h, double hbar, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], short *NodeStat, int elemX, double hxz[elemX]);

/*Function to calcualte volume change for a given lambda value*/
double LamCalc5a(double lambda, double *Vnorm, Aux *Lbound, Aux *auxNode, double h, double alt, double *Nenergy,
					int NIOtotal, int NumNodes, Coord NodeCoord[NumNodes], double *VmaxA, int itt, double maxX, double maxY,
						double Blen_in, double Uave_in, int Ntot, int NumFix, int *fixed, Astn *EmaxA, int Hflag, int elemX, double hxz[elemX]);

/*Function to find lambda value by iteration*/
double LamCalc5b(double delVol, double *Vnorm, Aux *Lbound, int Ntot, Aux *auxNode, double h, Astn *EmaxA, double alt, int NumFix, int *fixed,
					double *Nenergy, int NIOtotal, int NumNodes, Coord NodeCoord[NumNodes], double *VmaxA, int itt, double maxX, double maxY, int Hflag, int elemX, double hxz[elemX]);

/*Non-linear mapping funciton for velcoity to improve convergence*/
void NlMap(int Ntot, double *Vnorm, double Vmax);

/*Function to calculate extension velocities*/
void Vext(int NodeX, int NodeY, int Nodes2[NodeX][NodeY], short *NodeStat, double *Vnorm, double *lsf, int *active, int NumAct, int NumNodes, int Ntot, 
			Coord NodeCoord[NumNodes], double h, Aux *AuxNodes, int sign, double dt, int elemX, double hxz[elemX]);				

/*function that works out Velocity for nodes close to the boundary*/					
void LocalVext2(int NodeX, int NodeY, int Nodes[NodeX][NodeY], short *NodeStat, double *Vnorm, double *Vtemp, double *lsf, double *lsf_temp, short *known, short *trial,
					int NumNodes, Coord NodeCoord[NumNodes], int stop, double h, Aux *AuxNodes, int sign, int elemX, double hxz[elemX]);

/*Function to calcualte gradient using WENO scheme*/
double GradWENO(int Xi, int Yj, int num, double *lsf, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int sign, 
				double Vn, double dt, double h, int elemX, double hxz[elemX]);

/*sub-function for GradWENO*/
double GWsub(double v1,double v2,double v3,double v4,double v5);

/*Function to compute avergae and max difference between two level set functions*/
void PhiComp(int NumNodes, double *lsf1, double *lsf2, int itt, Coord *Delta);

/*Function to smooth level set function using simple linear approach*/
void SmoothPhi(int NumNodes, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], double *lsf, int *fixed, int NumFix);
				
/*Function to re-initalise the lsf as a signed distance function - similar to Vext above*/
void ReInt(int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int NumNodes, double *lsf, double h, 
				int sign, double *lsf_temp, int *fixed, int NumFix, int elemX, double hxy[elemX]);

/*Function to re-initalise the lsf as a signed distance function - similar to Vext above*/
void ReIntlsf2(int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int NumNodes, double *lsf, double h, 
				int sign, double *lsf_temp, int *fixed, int NumFix, int elemX, double hxy[elemX]);
			
/*Function to reinitalise lsf for inital set of trial nodes - similar to LocalVext above*/
void LocalInt(int NodeX, int NodeY, int Nodes[NodeX][NodeY], double *lsf, double *lsf_temp, short *known, short *trial, double h, int sign, int elemX, double hxz[elemX]);

/*function to determine if a number is in a set*/
int InSet(int num, int *array, int length);
