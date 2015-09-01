/*
 *  IFG.h
 *  
 *
 *  Created by Peter Dunning on 13/02/2009.
 *
 */
 // indicator arry for ne_conn(compressed storage)
 #include "ls_types.h"

void IFG_Status(double h,int elemX,int elemY,Elem Numbering[elemX][elemY],int *Acount,int NumElem,int NodeX, int NodeY, int Nodes2[NodeX][NodeY], 
				int NumNodes,Coord NodeCoord[NumNodes],double *lsf,short *NodeStat,short *ElemStat,Aux *auxNode,double ftv_lim, int itt,double hxz[elemX],double hz[elemX],Aux *AuxNodes, Bseg *bptr, int *na_conn_ind, int *na_conn, int *NumAux, int *NumBound);
				
void IFG_Matrix(double e11,double v12,double g33,double *KE,double h,short *ElemStat,short *NodeStat,Aux *auxNode,int *irn, int *jcn, double *A,
					int elemX, int elemY,Elem Number[elemX][elemY],int NumElem,int NumNodes,Coord NodeCoord[NumNodes],Shape *Pshp,Coord *Pcrds,double *alpha);
					
void IFG_Solve(int *irn, int *jcn, double *A, int TotFix, int *FixDofs, double *rhs_in, int NumEntries, int MatrixOrder, short *redun, int NumRhs, int itt);

/*calculate strain energy using least squares of integration points*/						
void IFG_StrainLS(int elemX, int elemY, Elem Number[elemX][elemY], short *ElemStat, short *NodeStat, double *alpha, double e11, double v12, double g33,
					double h, Aux *auxNode, double *Nenergy, Shape *Pshp, Coord *Pcrds, double *rhs1, double *rhs2, int NumNodes, Coord NodeCoord[NumNodes], 
						int gpoints, int Acount, int NIOtotal, double *Eenergy, short Hflag, double wgt);
