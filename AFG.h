/*
 *  AFG.h
 *  
 *
 *  Created by Peter Dunning on 05/05/2009.
 *
 */
 
#include "ls_types.h"

void AFG_Matrix(double e1, double e2, double v12, double v21, double g12, double g23, double g13,double *KE, double h, short *ElemStat, short *NodeStat, Aux *auxNode, int *irn, int *jcn, double *A, int elemX, int elemY, Elem Number[elemX][elemY], int NumElem, int NumNodes, Coord NodeCoord[NumNodes], double *alpha, short *redun, double hxz[elemX], double *theta, int numply, double plyt, double *Trans);

void AFG_Strain_LS(int elemX, int elemY, Elem Number[elemX][elemY], short *ElemStat, short *NodeStat, double *alpha, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double h, Aux *auxNode, double *Nenergy, double *rhs1, double *rhs2, int NumNodes, Coord NodeCoord[NumNodes], int gpoints, int Acount, int NIOtotal, double *Eenergy, short Hflag, double wgt, double hxz[elemX], double hz[elemX], double *theta, double *Trans, int numply, double plyt);

void CompositeKE(int num, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double *KE, double *theta,int numply, double plyt, int h, double hxz);

void CompositeKEGauss(int GP, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double *KE, double *theta, int numply, double plyt, double h, double hxz);

void ThetaUpdate(int elemX, int elemY,int NumNodes, Elem Number[elemX][elemY],Coord NodeCoord[NumNodes], double *theta, double *lsf, int tow, double *hxz, double h, double maxXZ);

void ThetaUpdate2(int *Lnodes, Coord *LCoord, int num, double *lsf, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], double h, int elemX, double *hxz, double *theta);

void ThetaUpdate3(double *Theta, int elemX, int elemY, Elem Number[elemX][elemY], double *lsf, double h, double *hxz);

void ThetaUpdate4(double *Theta, int elemX, int elemY, Elem Number[elemX][elemY], double *lsf, double *lsfb, int Numlsf, double h, double *hxz);
void ThetaUpdate5(double *Theta, int elemX, int elemY, Elem Number[elemX][elemY], double *lsf, double *lsfb, int Numlsf, double h, double *hxz);

void Transform(double *KE,double *T, double *K, int n);
