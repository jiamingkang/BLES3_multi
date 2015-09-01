/*
 *  FibreSensitivity.h
 *  
 *
 *  Created by Chris Brampton on 06/02/2013.
 *
 */
 
#include "ls_types.h"

/*void DthetaDlsfCalcOLDVERSION(double *DiffThetaLsf, int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double *lsf, double h, double *hxz, int NodeX, int NodeY, int Nodes2[NodeX][NodeY]);*/

void DthetaDlsfCalc(NElemDiff *DTL, int elemX, int elemY, Elem Number[elemX][elemY], double *lsf, double h, double *hxz, int NumNodes, Coord *NodeCoord, int Ntot_new, Aux *AuxNodes, double *NodeDTL, Gstn *NodeDTLg);

void SensitivityCalculation( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], double *rhs1, double *T, int numply, double plyt, double *Theta, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double *ElemThetaStn, double *ElemThetaStnB, double *EF, int itt,int Ntot_new, Aux *AuxNodes, double *NodeThetaStn, Gstn *Nodeseng);

void Nodesensitivity( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes, double *NodeThetaStn, double *NodesenTemp, Gstn *NodeDTLg);

void EnergyFactor(int elemX, int elemY, Elem Number[elemX][elemY], int NodeX, int NodeY, int Nodes2[NodeX][NodeY], double h, double *hxz, double *rhs1, double *rhs0, double *Theta1, double *Theta0, double e1, double e2, double v12, double v21, double g12, double g23, double g13,  double *EF, double *Trans, int numply, double plyt, int itt);

void ScaleDeltaTheta(int elemX, int elemY, double *ElemThetaStn, double *DeltaTheta, double DeltaThetaMax, int itt, double CompA[itt],  Elem Number[elemX][elemY], short *ElemStat, short *ElemStatb, int Numlsf, double *NodeThetaStn, double *NodeDTL, double *DeltaT, int Ntot_new, double *Vnorm, double *SensMaxOld );

void ScaleVnorm(int elemX, int elemY, double *ElemThetaStn , double *DeltaTheta, double DeltaThetaMax, int itt, double CompA[itt],  Elem Number[elemX][elemY], short *ElemStat, short *ElemStatb, int Numlsf, double *NodeThetaStn, double *NodeDTL, double *DeltaT, int Ntot_new , double *Vnorm,double *SensMaxOld);


void DeltaLsfCalc(int NumNodes, double *Deltalsf2, double *DeltaTheta, NElemDiff *DTL, short *ElemStat, double *Vnorm);

void NodeSen(int NumNodes,  NElemDiff *DTL, short *ElemStat, double *NodesenTemp, double *ElemThetaStn);

void CompositeKEGaussDIFF(int GP, double Q11, double Q22, double Q12, double Q16, double Q26, double Q66, double Q44, double Q55, double Q45, double *KE, int numply, double plyt, double h, double hxz, int num);

void CompositeKEDIFF(double Q11, double Q22, double Q12, double Q16, double Q26, double Q66, double Q44, double Q55, double Q45, double *KE, int numply, double plyt, double h, double hxz, int num);

double FibreContScore(int NumNodes, Coord NodeCoord[NumNodes],int elemX, int elemY, Elem Number[elemX][elemY], double h, double *hxz, double Range, double *theta, int itt);

void NodesenP( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes, double *NodeDTL, double *NodesenTemp, Gstn *NodeDTLg, Gstn *Nodeseng, double *Nodesenpp);

void CompositeBmatrix(int GP,  double *Bmatrix, int numply, double plyt, double h, double hxz, int num, int NumNodes, Coord NodeCoord[NumNodes], int *Lnodes, double *detJ);

void MatrixMultiple(double *a, double *b,  double *res, int acol, int arow, int bcol, int brow);

void TransM (double *a, double *b, int mcol, int mrow);

void GPsensitivity( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes, double *NodeDTL, double *NodesenTemp, Gstn *GPsens);
void NodeAuxsens( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes,  Gstn *GPsens, double *Nodesenpp);


void NodeAuxsens2( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes,  double *Nodesenpp, double *NodesenTemp);


void NodeAuxsens3( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h,  int Ntot_new, Aux *AuxNodes,   double *Nodesenpp,double *NodesenTemp);
void NodeAuxsens4( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h,  int Ntot_new, Aux *AuxNodes,   double *Nodesenpp,double *NodesenTemp);