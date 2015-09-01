/*
 * FibreSensitivity.c
 *  
 *  Created by Christopher Brampton on 03/02/2012.
 *	Calculates compliance sensitvtiy to fibre angle (defined by the implict function).
 */

#include "AFG.h"
#include "KMatrix.h"
#include "FixedGrid.h"
#include "Levels.h"
#include "ls_types.h"
#include "Strain.h"
#include "FibreSensitivity.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Array storing data for location of the 4 nodes in the square element*/
static Pos PoA[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};
const double gaB = 0.57735026919;	/*Constant for 2-point gauss rule*/

/*Differentiate theta wrt to the level set function*/
void DthetaDlsfCalc(NElemDiff *DTL, int elemX, int elemY, Elem Number[elemX][elemY], double *lsf, double h, double *hxz, int NumNodes, Coord *NodeCoord, int Ntot_new, Aux *AuxNodes, double *NodeDTL, Gstn *NodeDTLg)
{
	int i,j,k,num;							/*interger varribles*/
	int *Lnodes;
	double *Tlsf;					/*Local node data*/
	double *Nx, *Ny;				/*Differnetial of shape functions*/
	double s, n, a, b, X, Y, D, C;	/*local shape funtion co-ordinate values*/
	double dx, dy;					/*Differnetial of lsf wrt to x and y*/
	double mu, dmu;				/*Differnetial of mu wrt to lsf[k] where mu = ((dlsf/dy)/(dlsf/dx))*/
	double datan;					/*Differnetial of theta wrt mu from theta = 90 + arctan(mu)*/
	double v, dv, u, du;		/*Differential varibles used for the quiotent rule application in the claculation of dmu*/
	double dTheta;					/*The Results we are looking for*/
	double ftemp, dtemp;			/*tempoary varibles*/
	int o1, o2, o3;

	b = h;

	Lnodes = malloc(4*sizeof(int));
	Tlsf = malloc(4*sizeof(double));
	Nx = malloc(4*sizeof(double));
	Ny = malloc(4*sizeof(double));
    
    Gstn *GaussValue;
    double gax, gay;
    double h2 = 0.5 * h;
    double h2xz;
    int Gcount=0;
    Gstn *Gtemp;
    GaussValue = malloc(4*elemX*elemY*sizeof(Gstn));
    Gtemp = malloc(4*elemX*elemY*sizeof(Gstn));
	/*Intialise the DTL array entirly to zero*/
	for(i=0;i<NumNodes;i++)
	{
		DTL[i].ai = 0;
		DTL[i].bi = 0;
		DTL[i].ci = 0;
		DTL[i].di = 0;

		DTL[i].an = 0;
		DTL[i].bn = 0;
		DTL[i].cn = 0;
		DTL[i].dn = 0;

		DTL[i].ad = 0;
		DTL[i].bd = 0;
		DTL[i].cd = 0;
		DTL[i].dd = 0;
	}

	/*First calcualte the values dlsf/dy and dlsf/dx for all the elements. These are the basis of our equation theta = 90 + arctan( (dlsf/dy)/(dlsf/dx) )*/
	for(i=0;i<elemX;i++)
	{
		a = hxz[i];
		for(j=0;j<elemY;j++)
		{

			/*Get the local Node data*/
			num = Number[i][j].n - 1;
			Lnodes[0] = Number[i][j].a - 1;
			Lnodes[1] = Number[i][j].b - 1;
			Lnodes[2] = Number[i][j].c - 1;
			Lnodes[3] = Number[i][j].d - 1;

			Tlsf[0] = lsf[Lnodes[0]];
			Tlsf[1] = lsf[Lnodes[1]];
			Tlsf[2] = lsf[Lnodes[2]];
			Tlsf[3] = lsf[Lnodes[3]];
            
            
            
            h2xz = 0.5 * hxz[i];
            /*work out global co-ordinates of element center, from co-ords of node 1 and element edge length*/
            double cx = NodeCoord[Lnodes[0]].xz + h2xz;
            double cy = NodeCoord[Lnodes[0]].y + h2;

			/*if(num==500)
			{
				printf("\n Element = %i", num);
				printf("\n DLT[%i].a  = (%i, %i, %f)", Lnodes[0], DTL[Lnodes[0]].ai, DTL[Lnodes[0]].an, DTL[Lnodes[0]].ad);
				printf("\n DLT[%i].b  = (%i, %i, %f)", Lnodes[1], DTL[Lnodes[1]].bi, DTL[Lnodes[1]].bn, DTL[Lnodes[1]].bd);
				printf("\n DLT[%i].c  = (%i, %i, %f)", Lnodes[2], DTL[Lnodes[2]].ci, DTL[Lnodes[2]].cn, DTL[Lnodes[2]].cd);
				printf("\n DLT[%i].d  = (%i, %i, %f)", Lnodes[3], DTL[Lnodes[3]].di, DTL[Lnodes[3]].dn, DTL[Lnodes[3]].dd);
			}*/


			/*Shape function local cordinates, for this problem a single gauss point at the centre of the element will be adequate*/
			s = 0;
			n = 0;

			dx = 0;
			dy = 0;

			X = ((1)/(4*a*a*b));
			Y = ((1)/(4*a*b*b));
			/*Shape function derivates*/
			Nx[0] = X*(-1+n);
			Nx[1] = X*(1-n);
			Nx[2] = X*(1+n);
			Nx[3] = X*(-1-n);

			Ny[0] = Y*(s-1);
			Ny[1] = Y*(-s-1);
			Ny[2] = Y*(s+1);
			Ny[3] = Y*(-s+1);
            
            
            int x1,x2,x3,x4,y1,y2,y3,y4;
            x1=NodeCoord[Lnodes[0]].xz;
            x2=NodeCoord[Lnodes[1]].xz;
            x3=NodeCoord[Lnodes[2]].xz;
            x4=NodeCoord[Lnodes[3]].xz;
            
            y1=NodeCoord[Lnodes[0]].y;
            y2=NodeCoord[Lnodes[1]].y;
            y3=NodeCoord[Lnodes[2]].y;
            y4=NodeCoord[Lnodes[3]].y;
            
            // printf("\n GP = %i, x = %d,%d,%d,%d, y = %d,%d,%d,%d", GP, x1,x2,x3,x4, y1,y2,y3,y4);
            
            
            /*printf("\n GP = %i, nx = %f, ny = %f", GP, nx, ny);*/
            
            double *Jmatrix, *JmatrixInv, *dndtemp, *res;
            
            Jmatrix=malloc(4*sizeof(double));
            JmatrixInv=malloc(4*sizeof(double));
            dndtemp=malloc(8*sizeof(double));
            res=malloc(8*sizeof(double));
            
            dndtemp[0]=Nx[0];
            dndtemp[1]=Nx[1];
            dndtemp[2]=Nx[2];
            dndtemp[3]=Nx[3];
            dndtemp[4]=Ny[0];
            dndtemp[5]=Ny[1];
            dndtemp[6]=Ny[2];
            dndtemp[7]=Ny[3];
            
            Jmatrix[0]= Nx[0]*x1+Nx[1]*x2+Nx[2]*x3+Nx[3]*x4;
            Jmatrix[1]= Nx[0]*y1+Nx[1]*y2+Nx[2]*y3+Nx[3]*y4;
            Jmatrix[2]= Ny[0]*x1+Ny[1]*x2+Ny[2]*x3+Ny[3]*x4;
            Jmatrix[3]= Ny[0]*y1+Ny[1]*y2+Ny[2]*y3+Ny[3]*y4;
            // invrse J= 1/det(J)      *  [Jmatrix change 0 to 3, and add a '-' in 1 and 2]
            double detJ;
            detJ= Jmatrix[0]*Jmatrix[3]-Jmatrix[1]*Jmatrix[2];
            detJ= fabs(detJ);
            
            JmatrixInv[0]=Jmatrix[3]/detJ;
            JmatrixInv[1]=-1.0*Jmatrix[1]/detJ;
            JmatrixInv[2]=-1.0*Jmatrix[2]/detJ;
            JmatrixInv[3]=Jmatrix[0]/detJ;
            
            
            
            MatrixMultiple(JmatrixInv, dndtemp, res, 2, 2, 4, 2);
            
            Nx[0]=res[0];
            Nx[1]=res[1];
            Nx[2]=res[2];
            Nx[3]=res[3];
            Ny[0]=res[4];
            Ny[1]=res[5];
            Ny[2]=res[6];
            Ny[3]=res[7];
            
            free(res);
            free(Jmatrix);
            free(JmatrixInv);
            free(dndtemp);
            
            
            
            
            
            
			/*if((num==799))
			{
				printf("\nNx[0] = %f", Nx[0]);
				printf("\nNx[1] = %f", Nx[1]);
				printf("\nNx[2] = %f", Nx[2]);
				printf("\nNx[3] = %f", Nx[3]);
			}*/
			
			/*Now use the shape function derivaties to calcualte x and y lsf gradients over the element*/
			for(k=0;k<4;k++)
			{
				dx += Nx[k]*Tlsf[k];
				dy += Ny[k]*Tlsf[k];
				/*if((num==799)){printf("\nNx[%i] * Tlsf[%i] = %f x %f = %f \ndx = %f", k, k, Nx[k], Tlsf[k], Nx[k]*Tlsf[k], dx);}*/
			}
			/*Now we the vaule of each of the lsf functions multiplied by there neighbouring shape function*/
			/*From this we can differentiate theta wrt to the lsf values. This will be do at each node so that the value will be
			the sensitivity of each node lsf value to the element value of theta. This is what we will combine with the sesnitvtiy
			of the element theta value to the golbal compliance to get us the sensitvity of each node lsf value to compliance*/

			/*Perform the cacluation for each node, theta differntiated wrt to each node lsf value*/
			for(k=0;k<4;k++)
			{
				dTheta = 0.0;
				if(k==0){o1=1; o2=2; o3=3;}
				else if(k==1){o1=0; o2=2; o3=3;}
				else if(k==2){o1=0; o2=1; o3=3;}
				else if(k==3){o1=0; o2=1; o3=2;}

				C = Nx[o1]*Tlsf[o1] + Nx[o2]*Tlsf[o2] + Nx[o3]*Tlsf[o3];	/*get the "constant" for the differential term*/
				D = Ny[o1]*Tlsf[o1] + Ny[o2]*Tlsf[o2] + Ny[o3]*Tlsf[o3];	/*which is the other node lsf differentials wrt x and y*/

				u = Ny[k]*Tlsf[k] + D;		v = Nx[k]*Tlsf[k] + C;
				du = Ny[k];				dv = Nx[k];

				/*if((num==799)&&(k==0))
				{
					printf("\n Element = %i", num);
					printf("\ndy: %f	   dx: %f", dy, dx);
					printf("\nlsf: %f, %f, %f, %f",Tlsf[0], Tlsf[1], Tlsf[2], Tlsf[3]);
					printf("\nu: %f, v: %f, du: %f, dv: %f", u, v, du, dv);
				}*/

				if((fabs(dx)>0.000001)&&(fabs(v)>0.000001)) 	/*So long as dlsf/dx is not close to zero we can use the normal formulation*/
				{
					mu = dy/dx;							/*mu = ((dlsf/dx)/(dlsf/dy))*/
					datan = 1.0/(mu*mu + 1.0);				/*Diffential wrt mu of 90 + arctan(mu) = 1/(mu^2+1)*/

					/*Now we can actually do the sum*/
					dtemp = ((v*du) - (u*dv));
					ftemp = (v*v);
					dTheta = (dtemp/ftemp)*datan;		/*Get value of dTheta*/
					/*if((num==500)&&(k==0))
					{
						printf("\nmu: %f datan: %f", mu, datan);
						printf("\ndtemp = %f, ftemp = %f", dtemp, ftemp);
						printf("\ndTheta = %f", dTheta);
					}*/

				}
				else if((fabs(dy)>0.00000000001)&&(fabs(u)>0.000001))  	/*If dlsf/dx is close to zero we can avoid a zero graident by taking advantage of the fact that*/
				{								/* theta = 90 + arctan( (dlsf/dy)/(dlsf/dx) ) = 180 - arctan( (dlsf/dx)/(dlsf/dy) )*/
					mu = dx/dy;							/*mu = ((dlsf/dx)/(dlsf/dy))*/
					datan = (-1.0)/(mu*mu + 1.0);			/*Diffential wrt mu of 90 + arctan(mu) = 1/(mu^2+1)*/
					
					/*Now we can actually do the sum (see how v and u are the other way round this time*/
					dtemp = ((u*dv) - (v*du));
					ftemp = (u*u);
					dTheta = (dtemp/ftemp)*datan;
				}
				else	/*If we get here then both dy and dx are close to zero and the implict function is breaking down so we have a stability problem*/
				{
					printf("\n\n ----------------!!!!!ERRROR Dlsf/Dy and Dlsf/Dx are both close to zero!!!!!!!!!!!!---------------------\n");
					printf(" ----------------!!!!!!!!!!!!!Implicit function is breaking down!!!!!!!!!!!!!!!!!!!---------------------\n");
					printf(" ----------------!!!!!!!!!!!!!!!!!!Model is Losing Stabiltiy!!!!!!!!!!!!!!!!!!!!!!!---------------------\n");
					printf(" ----------------!!!!!!!!!!!!!!!!!!!!!!!!! dTheta[%i] = 0 !!!!!!!!!!!!!!!!!!!!!!!!!---------------------\n", Lnodes[k]);
					printf("\nelem: %i	 dy: %f	   dx: %f",num, dy, dx);
					printf("\nlsf: %f, %f, %f, %f",Tlsf[0], Tlsf[1], Tlsf[2], Tlsf[3]);
					dTheta = 0;
				}
				/*if(num==500){printf("\nLnodes[k] = %i", Lnodes[k]);}
				DTL[Lnodes[k]].ai = 1;	DTL[Lnodes[k]].an = num;  DTL[Lnodes[k]].ad = dTheta;*/

				/*Now Store the DTL array enteries*/
				if(k==0){DTL[Lnodes[k]].ai = 1;	DTL[Lnodes[k]].an = num;  DTL[Lnodes[k]].ad = dTheta;}
				if(k==1){DTL[Lnodes[k]].bi = 1;	DTL[Lnodes[k]].bn = num;  DTL[Lnodes[k]].bd = dTheta; }
				if(k==2){DTL[Lnodes[k]].ci = 1;	DTL[Lnodes[k]].cn = num;  DTL[Lnodes[k]].cd = dTheta; }
				if(k==3){DTL[Lnodes[k]].di = 1;	DTL[Lnodes[k]].dn = num;  DTL[Lnodes[k]].dd = dTheta; }
			
                gax = gaB * (double)PoA[k].x * h2xz;
                gay = gaB * (double)PoA[k].y * h2;
                GaussValue[Gcount].x = cx + gax;
                GaussValue[Gcount].y = cy + gay;
                GaussValue[Gcount].u = dTheta;
                //Gtemp[Gcount].u = dTheta;			/*store gauss point sensitvity for the nodal sensitivty value*/
                
                //printf("dtehat: %d, %.16lf", Gcount, dTheta);
                
                GaussValue[Gcount].a = 0.25*h*hxz[i];	/*Store the gauss point element area contribution*/
                
                Gcount++;
            
            }
            double *Nn;
            Nn = malloc(4*sizeof(double));
            int Gcount2;
            Gcount2=Gcount-4;
			for(k=0;k<4;k++)
            {
                Nn[0]=0.25*(1-gaB*(double)PoA[k].x)*(1-gaB* (double)PoA[k].y);
                Nn[1]=0.25*(1+gaB*(double)PoA[k].x)*(1-gaB* (double)PoA[k].y);
                Nn[2]=0.25*(1+gaB*(double)PoA[k].x)*(1+gaB* (double)PoA[k].y);
                Nn[3]=0.25*(1-gaB*(double)PoA[k].x)*(1+gaB* (double)PoA[k].y);
                
                //PoA[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};
                
                //GaussValue[Gcount2].u=Gtemp[Gcount-4].u*Nn[0]+Gtemp[Gcount-3].u*Nn[1]+Gtemp[Gcount-2].u*Nn[2]+Gtemp[Gcount-1].u*Nn[3];
            
            Gcount2++;
            
            }
            free(Nn);

			/*if(num==500)
			{
				printf("\n Element = %i", num);
				printf("\n DLT[%i].a  = (%i, %i, %f)", Lnodes[0], DTL[Lnodes[0]].ai, DTL[Lnodes[0]].an, DTL[Lnodes[0]].ad);
				printf("\n DLT[%i].b  = (%i, %i, %f)", Lnodes[1], DTL[Lnodes[1]].bi, DTL[Lnodes[1]].bn, DTL[Lnodes[1]].bd);
				printf("\n DLT[%i].c  = (%i, %i, %f)", Lnodes[2], DTL[Lnodes[2]].ci, DTL[Lnodes[2]].cn, DTL[Lnodes[2]].cd);
				printf("\n DLT[%i].d  = (%i, %i, %f)", Lnodes[3], DTL[Lnodes[3]].di, DTL[Lnodes[3]].dn, DTL[Lnodes[3]].dd);
			}*/
		}
	}
    
    
    int kkk;
    
    for(kkk=0;kkk<(4*elemY*elemX); kkk++)
    {
        NodeDTLg[kkk].u=GaussValue[kkk].u;
        NodeDTLg[kkk].x=GaussValue[kkk].x;
        NodeDTLg[kkk].y=GaussValue[kkk].y;
        NodeDTLg[kkk].a=GaussValue[kkk].a;
        
    }

    
    double rad = 2*h;
    double temp2, ftemp2;
    int kk;
    double cx, cy;
    for(kk=0;kk<(Ntot_new);kk++)
    {
        NodeDTL[kk] = 0.0;
    }
    
    /*Now we need to get the nodal values of sensitvity from the elemental values by sending each node through the least squared routine
     Simpler then the standard 2D or 3D code as we have no AUX nodes to concern ourselves with*/
    for(kk=0;kk<Ntot_new;kk++)
    {
        if (kk<NumNodes) {
            
            cx = NodeCoord[kk].xz;
            cy = NodeCoord[kk].y; /*read in co-ordiantes directly*/
            ftemp2 = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            //NodeDTL[kk]= ftemp2; /* multiply smothed strain energy by weight */
        }
        else{
            cx = AuxNodes[kk-NumNodes].x;
            cy = AuxNodes[kk-NumNodes].y; /*read in cveo-ordiantes directly*/
            
            // if (AuxNodes[n-NumNodes].n+1>NumNodes)    // if number is larger than NumNodes, it haven't been calculate
            
            ftemp2 = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            //NodeDTL[kk]= ftemp2; /* multiply smothed strain energy by weight */
            
        }
        /*use least squares to calculate nodal strain energy*/
    }

    
    free(Gtemp);
    
    free(GaussValue);
	
}
/*Calcualtes the Snesitvity of each elements fibre angle to fibre angle */
void SensitivityCalculation( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], double *rhs1, double *Trans, int numply, double plyt, double *Theta, double e1, double e2, double v12, double v21, double g12, double g23, double g13, double *ElemThetaStn, double *ElemThetaStnB, double *EF, int itt, int Ntot_new, Aux *AuxNodes, double *NodeThetaStn, Gstn *Nodeseng)
{
	/*First get the local displacements on every node*/
	FILE *outfile;	/*File varible for output files*/
	char plotname[40];	/*variable to change names of plotting output files*/
	int i, j, k, n, m, p, q, f, num,x,y;	/*interger varibles for the loops*/
	double *UX, *UY, *UZ, *RX, *RY, *RZ;			/*Local displacement values*/
	double *OX, *OY, *OZ, *ORX, *ORY, *ORZ;			/*Local strain factors values*/
	int *Lnodes;									/*local nodes*/
	double temp;
	double C, S, C4, S4, C2S2, CS, C3S, CS3, C2, S2, C3, S3, C2S, CS2; /*Cosine and Sine varibles*/
	double E11, E22, E12, E66, E44, E55;					/*Material Property varibles for strain*/
	double F11, F22, F12, F66, F44, F55;					/*Material Property varibles for stress*/
	double Q11, Q22, Q12, Q16, Q26, Q66, Q44, Q45, Q55;		/*Qmatrix varibles*/
	double P11, P22, P12, P16, P26, P66, P44, P45, P55;		/*Stress matrix material varibles*/
	double Q11d, Q22d, Q12d, Q16d, Q26d, Q66d, Q44d, Q55d, Q45d;	/*Qmatrix Zero degree varibles*/
	double P11d, P22d, P12d, P16d, P26d, P66d, P44d, P55d, P45d;	/*Stress matrix Zero degree material varibles*/
	double Stn, Stn1, Stn2, Stn3, StnB, StnC, StnD;					/*Senstivtiy calculation values*/
	double StnX, StnY, StnZ, StnRx, StnRy;
	double *KE, *SE, *KEb, *SEb;					/*Guass Matircies*/
	double cx, cy;
	double h2 = 0.5 * h;
	double h2xz;
	int Gcount = 0;
	int gpoints = 4*elemX*elemY;
	double gax, gay;
	Gstn *GaussValue;	

	Lnodes = malloc(4*sizeof(int));

	UX = malloc(4*sizeof(double));
	UY = malloc(4*sizeof(double));
	UZ = malloc(4*sizeof(double));
	RX = malloc(4*sizeof(double));
	RY = malloc(4*sizeof(double));
	RZ = malloc(4*sizeof(double));
	
	OX = malloc(4*sizeof(double));
	OY = malloc(4*sizeof(double));
	OZ = malloc(4*sizeof(double));
	ORX = malloc(4*sizeof(double));
	ORY = malloc(4*sizeof(double));
	ORZ = malloc(4*sizeof(double));

	KE = malloc(576*sizeof(double));
	KEb = malloc(576*sizeof(double));
	SE = malloc(576*sizeof(double));
	SEb = malloc(576*sizeof(double));

	GaussValue = malloc(4*elemX*elemY*sizeof(Gstn));

	for(n=0;n<(elemX*elemY);n++)
	{
		ElemThetaStn[n] = 0.0;
	}

	/*In this loop we perform the entire snesitvity anaylisis for every single element, this involves several stages 
	to get all the element properties we need before we actually get start calculating the senstvity*/
	for(n=0;n<elemX;n++)
	{
		for(m=0;m<elemY;m++)
		{

			num = Number[n][m].n-1;		/*Get element number*/
			Lnodes[0] = Number[n][m].a-1;
			Lnodes[1] = Number[n][m].b-1;
			Lnodes[2] = Number[n][m].c-1;
			Lnodes[3] = Number[n][m].d-1;
	
			h2xz = 0.5 * hxz[n];
			/*work out global co-ordinates of element center, from co-ords of node 1 and element edge length*/
			double cx = NodeCoord[Lnodes[0]].xz + h2xz;
			double cy = NodeCoord[Lnodes[0]].y + h2;
			/*---------First we need to get the displacements local to each element from the global displacements----------*/
			p = 9*n;				/* Start position in strorage of Trans matrix for this element*/

			for(i=0;i<4;i++)
			{
				j = 6 * Lnodes[i];		/*Dofs of this node*/
				/*Transform the global displacements (rhs) to the local displacement for each node in the element*/
				UX[i] = rhs1[j]*Trans[0+p] + rhs1[j+1]*Trans[1+p] + rhs1[j+2]*Trans[2+p];		/*local local X displacement of node a*/
				UY[i] = rhs1[j]*Trans[3+p] + rhs1[j+1]*Trans[4+p] + rhs1[j+2]*Trans[5+p];		/*local Y displacement of node a*/
				UZ[i] = rhs1[j]*Trans[6+p] + rhs1[j+1]*Trans[7+p] + rhs1[j+2]*Trans[8+p];		/*local Z displacement of node a*/
				RX[i] = rhs1[j+3]*Trans[0+p] + rhs1[j+4]*Trans[1+p] + rhs1[j+5]*Trans[2+p];		/*local rotX displacement of node a*/
				RY[i] = rhs1[j+3]*Trans[3+p] + rhs1[j+4]*Trans[4+p] + rhs1[j+5]*Trans[5+p];		/*local rotY displacement of node a*/
				RZ[i] = 0.0;	/*rhs1[j+3]*Trans[6+p] + rhs1[j+4]*Trans[7+p] + rhs1[j+5]*Trans[8+p];*/		/*local rotZ displacement of node a*/
                
               // printf("\n disp at element %d is %f", num, UY[i]);
			}
            double *disp;
            disp=malloc(6*4*sizeof(double));
            for(i=0;i<4;i++)
            {
                disp[6*i+0]=UX[i];
                disp[6*i+1]=UY[i];
                disp[6*i+2]=UZ[i];
                disp[6*i+3]=RX[i];
                disp[6*i+4]=RY[i];
                disp[6*i+5]=RZ[i];
            
            }
			/*------------Now we need to change these in to the local stress from the strains-------------------------*/
			/*Note that these are not the the actual stresses, they will become the actual stress when multiplied by the shape
			functions and differnetial matrices included in the Compliance matrix formation (SE). So just like the strains, we are inputing the
			stress equiverlent of the displacment into the sensitvity equation. The sensitvity equations themselves turn these values 
			into the stress and strains in the element. This does work trust me when you look at this again in a few months or years.
			The explination is on the back cover of log book volume 15*/

			/*Get the C material Matrix properties for this element's value of theta*/
			/*Set up rotation matrix*/
           
            if (num==80||num==81||num==40||num==41||num==90||num==89||num==130||num==129||num==399||num==439) {
                printf("\n theta value of element %d: %.16lf", num, Theta[num] );
            }
            C = cos(Theta[num]);
            //printf("C= %f", C);
            S = sin(Theta[num]);
            
            
			//C = cos(3.1415926585/6.0);
            //printf("C= %f", C);
			//S = sin(3.1415926585/6.0);
            
            
            
            //this part is to verify if the sensitivity is correct
            double dT;
            dT=-0.01;
            if (num==80)
            {
            //C = cos(Theta[num]+dT);
            //printf("C= %f", C);
            //S = sin(Theta[num]+dT);
            }
            
            
			C2 = C*C;
			S2 = S*S;
			C3 = C*C*C;
			S3 = S*S*S;
			C2S = C*C*S;
			CS2 = C*S*S;
			C2S2 = C*C*S*S;
			CS = C*S;
			C3S = C*C*C*S;
			CS3 = C*S*S*S;
			C4 = C*C*C*C;
			S4 = S*S*S*S;

			/*Caclcualte material property values*/
			E11 = e1/(1.0-(v12*v21));
			E22 = e2/(1.0-(v12*v21));
			E12 = v21*E11;				/*Could also be E22*v12 as they are the same*/
			E66 = g12;
			E44 = g13;
			E55 = g23;
            
            // Find out the inverse matric of the above one ±±±±!!!!!!!!!!!!!!!!!!
            
            //printf("\n%f,%f,%f,%f,%f,%f", E11,E22,E12,E66,E44,E55);
            
            
            
            /*First we want to get the stress material Matrix*/
             F11 = 1.0/e1;
             F22 = 1.0/e2;
             F12 = -v21/e2;   // Chris made a mistake in here , we use v21 to replace v12, checked with matlab
             F66 = 1.0/g12;
             F44 = 1.0/g13;
             F55 = 1.0/g23;
           // printf("\n%f,%f,%f,%f,%f,%f", F11,F22,F12,F66,F44,F55);
            
           /*
            F11 = 1.0/E11;
            F22 = 1.0/E22;
            F12 = 1.0/E12;
            F66 = 1.0/E66;
            F44 = 1.0/E44;
            F55 = 1.0/E55;
            */

			/*Create C material Matrix*/
			
            
            Q11 = C4*E11 + 2*C2*S2*E12 + S4*E22 + 4*C2*S2*E66;
			Q12 = C2*S2*E11 + S4*E12 + C4*E12 + C2*S2*E22 - 4*C2*S2*E66;
			Q16 = C3*S*E11 + C*S3*E12 - C3*S*E12 - C*S3*E22 - 2*C*S*(C2-S2)*E66;
			Q22 = S4*E11 + 2*C2*S2*E12 + C4*E22 + 4*C2*S2*E66;
			Q26 = C*S3*E11 + C3*S*E12 - C*S3*E12 - C3*S*E22 + 2*C*S*(C2-S2)*E66;
			Q66 = C2S2*E11 - 2*C2S2*E12 + C2S2*E22 + (C4-2*C2S2+S4)*E66;
			Q44 = C2*E44 + S2*E55;
			Q45 = CS*E44 - CS*E55;
			Q55 = S2*E44 + C2*E55;
            
            /* Chris's result, Q12,22,26,45 is not correct, not sure now
            Q11 = C4*E11 + 2*C2*S2*E12 + S4*E22 + 4*C2*S2*E66;
            Q12 = C2*S2*E11 + S4*E12 + C4*E12 + C2*S2*E22 + 2*C2*S2*E66;
            Q16 = C3*S*E11 + C*S3*E12 - C3*S*E12 - C*S3*E22 - 2*C*S*(C2-S2)*E66;
            Q22 = S4*E11 + 2*C2*S2*E12 + C4*E22 + C2*S2*E66;
            Q26 = C*S3*E11 + C3*S*E12 - C*S3*E12 - C3*S*E22 - C*S*(C2-S2)*E66;
            Q66 = C2S2*E11 - 2*C2S2*E12 + C2S2*E22 + (C4-2*C2S2+S4)*E66;
            Q44 = C2*E44 + S2*E55;
            Q45 = CS*E44 + CS*E55;
            Q55 = S2*E44 + C2*E55;
             */
            
            //differential
            /*
            Q11d = -4*C3*S*E11 + 4*(C3*S-C*S3)*E12 + 4*S3*C*E22 + 8*(C3*S-C*S3)*E66;
            
            Q16d = (C4-3*C2*S2)*E11 + (6*C2*S2-S4-C4)*E12 - (3*S2*C2-S4)*E22 + 2*(6*C2*S2-C4-S4)*E66;
            
            Q66d = 2*(C3*S-S3*C)*E11 - 4*(C3*S-S3*C)*E12 + 2*(C3*S-S3*C)*E22 + 8*(S3*C-C3*S)*E66;
            Q44d = 2*S*C*(E55-E44);
            Q55d = 2*S*C*(E44-E55);
            
            
            
            Q12d = 2*(C3*S-C*S3)*E11 + 4*(C*S3 - C3*S)*E12 + 2*(C3*S-C*S3)*E22 + 4*(C3*S-S3*C)*E66;
            
            Q22d = 4*S3*C*E11 + 4*(C3*S-C*S3)*E12 - 4*C3*S*E22 + 2*(C3*S-C*S3)*E66;
            
            Q26d = (3*S2*C2-S4)*E11 + (C4+S4-6*C2*S2)*E12 - (C4-3*S2*C2)*E22 - (C4+S4-6*C2*S2)*E66;
            Q45d = (C2-S2)*(E44+E55);
            */
            //printf("\nQ11 :%f,Q12:%f,Q16:%f,Q22:%f,Q26:%f,Q66:%f,Q44:%f,Q45:%f,Q55:%f",Q11,Q12,Q16,Q22,Q26,Q66,Q44,Q45,Q55);
			
			/*Get the Stress terms by multiplying the stiffness matrix by the displacement*/
			CompositeKEDIFF(Q11, Q22, Q12, Q16, Q26, Q66, Q44, Q55, Q45, KEb, numply, plyt, h, hxz[n],num);
		
			/*Strain * Cmatrix = Stress*/
			for(i=0;i<4;i++)	 /*Calculate it for each node in turn*/
			{
				k= 6*i;		/*k and i are used for indexing purposes*/

				/*UX Term*/
				OX[i] = (UX[0]*KEb[0+k] + UY[0]*KEb[24+k] + UZ[0]*KEb[48+k] + RX[0]*KEb[72+k] + RY[0]*KEb[96+k] + RZ[0]*KEb[120+k] + UX[1]*KEb[144+k] + UY[1]*KEb[168+k] + UZ[1]*KEb[192+k] + RX[1]*KEb[216+k] + RY[1]*KEb[240+k] + RZ[1]*KEb[264+k] + UX[2]*KEb[288+k] + UY[2]*KEb[312+k] + UZ[2]*KEb[336+k] + RX[2]*KEb[360+k] + RY[2]*KEb[384+k] + RZ[2]*KEb[408+k] + UX[3]*KEb[432+k] + UY[3]*KEb[456+k] + UZ[3]*KEb[480+k] + RX[3]*KEb[504+k] + RY[3]*KEb[528+k] + RZ[3]*KEb[552+k]);
				/*UY Term*/																																																					
				OY[i] = (UX[0]*KEb[1+k] + UY[0]*KEb[25+k] + UZ[0]*KEb[49+k] + RX[0]*KEb[73+k] + RY[0]*KEb[97+k] + RZ[0]*KEb[121+k] + UX[1]*KEb[145+k] + UY[1]*KEb[169+k] + UZ[1]*KEb[193+k] + RX[1]*KEb[217+k] + RY[1]*KEb[241+k] + RZ[1]*KEb[265+k] + UX[2]*KEb[289+k] + UY[2]*KEb[313+k] + UZ[2]*KEb[337+k] + RX[2]*KEb[361+k] + RY[2]*KEb[385+k] + RZ[2]*KEb[409+k] + UX[3]*KEb[433+k] + UY[3]*KEb[457+k] + UZ[3]*KEb[481+k] + RX[3]*KEb[505+k] + RY[3]*KEb[529+k] + RZ[3]*KEb[553+k]);
				/*UZ Term*/																																																					
				OZ[i] = (UX[0]*KEb[2+k] + UY[0]*KEb[26+k] + UZ[0]*KEb[50+k] + RX[0]*KEb[74+k] + RY[0]*KEb[98+k] + RZ[0]*KEb[122+k] + UX[1]*KEb[146+k] + UY[1]*KEb[170+k] + UZ[1]*KEb[194+k] + RX[1]*KEb[218+k] + RY[1]*KEb[242+k] + RZ[1]*KEb[266+k] + UX[2]*KEb[290+k] + UY[2]*KEb[314+k] + UZ[2]*KEb[338+k] + RX[2]*KEb[362+k] + RY[2]*KEb[386+k] + RZ[2]*KEb[410+k] + UX[3]*KEb[434+k] + UY[3]*KEb[458+k] + UZ[3]*KEb[482+k] + RX[3]*KEb[506+k] + RY[3]*KEb[530+k] + RZ[3]*KEb[554+k]);
				/*RX Term*/																																																					
				ORX[i] = (UX[0]*KEb[3+k] + UY[0]*KEb[27+k] + UZ[0]*KEb[51+k] + RX[0]*KEb[75+k] + RY[0]*KEb[99+k] + RZ[0]*KEb[123+k] + UX[1]*KEb[147+k] + UY[1]*KEb[171+k] + UZ[1]*KEb[195+k] + RX[1]*KEb[219+k] + RY[1]*KEb[243+k] + RZ[1]*KEb[267+k] + UX[2]*KEb[291+k] + UY[2]*KEb[315+k] + UZ[2]*KEb[339+k] + RX[2]*KEb[363+k] + RY[2]*KEb[387+k] + RZ[2]*KEb[411+k] + UX[3]*KEb[435+k] + UY[3]*KEb[459+k] + UZ[3]*KEb[483+k] + RX[3]*KEb[507+k] + RY[3]*KEb[531+k] + RZ[3]*KEb[555+k]);
				/*RY Term*/																																																					
				ORY[i] = (UX[0]*KEb[4+k] + UY[0]*KEb[28+k] + UZ[0]*KEb[52+k] + RX[0]*KEb[76+k] + RY[0]*KEb[100+k] + RZ[0]*KEb[124+k] + UX[1]*KEb[148+k] + UY[1]*KEb[172+k] + UZ[1]*KEb[196+k] + RX[1]*KEb[220+k] + RY[1]*KEb[244+k] + RZ[1]*KEb[268+k] + UX[2]*KEb[292+k] + UY[2]*KEb[316+k] + UZ[2]*KEb[340+k] + RX[2]*KEb[364+k] + RY[2]*KEb[388+k] + RZ[2]*KEb[412+k] + UX[3]*KEb[436+k] + UY[3]*KEb[460+k] + UZ[3]*KEb[484+k] + RX[3]*KEb[508+k] + RY[3]*KEb[532+k] + RZ[3]*KEb[556+k]);
				/*ORZ Term*/
				ORZ[i] = 0;

				/*OX[i] = UX[i];
				OY[i] = UY[i];
				OZ[i] = UZ[i];
				ORX[i] = RX[i];
				ORY[i] = RY[i];*/
			}

			/*Now we have the stress factors we can set up the KE and SE guasse matrices for the sensitvity calculation*/
			/*if((num==83)||(num==723)/*||(num==1)||(num==761)*)
			{
				printf("\n\n Element = %i (%i,%i)", num,n,m);
				printf("\n Lnodes: (%i, %i, %i, %i)", Lnodes[0], Lnodes[1], Lnodes[2], Lnodes[3]);
				printf("\nNode 0: UX: %f, UY: %f, UZ: %f RX: %f RY: %f RZ: %f", UX[0], UY[0], UZ[0], RX[0], RY[0], RZ[0]);
				printf("\nNode 1: UX: %f, UY: %f, UZ: %f RX: %f RY: %f RZ: %f", UX[1], UY[1], UZ[1], RX[1], RY[1], RZ[1]);
				printf("\nNode 2: UX: %f, UY: %f, UZ: %f RX: %f RY: %f RZ: %f", UX[2], UY[2], UZ[2], RX[2], RY[2], RZ[2]);
				printf("\nNode 3: UX: %f, UY: %f, UZ: %f RX: %f RY: %f RZ: %f", UX[3], UY[3], UZ[3], RX[3], RY[3], RZ[3]);

				printf("\n\n Q11 = %f, Q12 = %f, Q16 = %f", Q11, Q12, Q16);
				printf("\n Q22 = %f, Q26 = %f, Q66 = %f", Q22, Q26, Q66);
				printf("\n Q44 = %f, Q55 = %f, Q45 = %f", Q44, Q55, Q45);

				printf("\n\n Theta[%i] = %f", num, Theta[num]);
				printf("\n C = %f, S = %f", C, S);
				printf("\n C2 = %f, S2 = %f", C2, S2);
				printf("\n C3 = %f, S3 = %f", C3, S3);
				printf("\n C4 = %f, S4 = %f", C4, S4);
				printf("\n CS = %f, C2S2 = %f", CS, C2S2);
				printf("\n C2S = %f, CS2 = %f", C2S, CS2);
				printf("\n C3S = %f, CS3 = %f", C3S, CS3);

				/*printf("\n\n T:(%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)", T[0],T[1],T[2],T[3],T[4],T[5],T[6],T[7],T[8],T[9],T[10],T[11],T[12]);*/
				
				/*printf("\nNode 0: OX: %f, OY: %f, OZ: %f ORX: %f ORY: %f ORZ: %f", OX[0], OY[0], OZ[0], ORX[0], ORY[0], ORZ[0]);
				printf("\nNode 1: OX: %f, OY: %f, OZ: %f ORX: %f ORY: %f ORZ: %f", OX[1], OY[1], OZ[1], ORX[1], ORY[1], ORZ[1]);
				printf("\nNode 2: OX: %f, OY: %f, OZ: %f ORX: %f ORY: %f ORZ: %f", OX[2], OY[2], OZ[2], ORX[2], ORY[2], ORZ[2]);
				printf("\nNode 3: OX: %f, OY: %f, OZ: %f ORX: %f ORY: %f ORZ: %f", OX[3], OY[3], OZ[3], ORX[3], ORY[3], ORZ[3]);
			}*/
       
            
            
			P11 = C4*F11 + 2*C2*S2*F12 + S4*F22 + C2*S2*F66;
			P12 = C2*S2*F11 + S4*F12 + C4*F12 + C2*S2*F22 - C2*S2*F66;
			P16 = 2*C3*S*F11 + 2*C*S3*F12 - 2*C3*S*F12 - 2*C*S3*F22 - C*S*(C2-S2)*F66;
			P22 = S4*F11 + 2*C2*S2*F12 + C4*F22 + C2*S2*F66;
			P26 = 2*C*S3*F11 + 2*C3*S*F12 - 2*C*S3*F12 - 2*C3*S*F22 + CS*(C2-S2)*F66;
			P66 = 4*C2*S2*F11 - 8*C2*S2*F12 + 4*C2*S2*F22 + (C4-2*C2*S2+S4)*F66;
			P44 = C2*F44 + S2*F55;
			P45 = CS*F44 - CS*F55;
			P55 = S2*F44 + C2*F55;
            //printf("\nQ11 :%f,Q12:%f,Q16:%f,Q22:%f,Q26:%f,Q66:%f,Q44:%f,Q45:%f,Q55:%f",P11,P12,P16,P22,P26,P66,P44,P45,P55);
			/*Intialize the Guass matrices*/

			Stn = 0;	/*Intiailise the Stn value*/

			/*Set up differentiation of Q matrix WRT theta in order to get KEzero, repeated for SEzero*/
			
            
            Q11d = -4*C3*S*E11 + 4*(C3*S-C*S3)*E12 + 4*S3*C*E22 + 8*(C3*S-C*S3)*E66;
			Q12d = 2*(C3*S-C*S3)*E11 + 4*(C*S3 - C3*S)*E12 + 2*(C3*S-C*S3)*E22 - 8*(C3*S-S3*C)*E66;
			Q16d = (C4-3*C2*S2)*E11 + (6*C2*S2-S4-C4)*E12 - (3*S2*C2-S4)*E22 + 2*(6*C2*S2-C4-S4)*E66;
			Q22d = 4*S3*C*E11 + 4*(C3*S-C*S3)*E12 - 4*C3*S*E22 + 8*(C3*S-C*S3)*E66;
			Q26d = (3*S2*C2-S4)*E11 + (C4+S4-6*C2*S2)*E12 - (C4-3*S2*C2)*E22 + 2*(C4+S4-6*C2*S2)*E66;
			Q66d = 2*(C3*S-S3*C)*E11 - 4*(C3*S-S3*C)*E12 + 2*(C3*S-S3*C)*E22 + 8*(S3*C-C3*S)*E66;
			Q44d = 2*S*C*(E55-E44);
			Q55d = 2*S*C*(E44-E55);
			Q45d = (C2-S2)*(E44-E55);
            
            // form the dC/dtheta matrix
            double *dQ, *Q;
            int www;
            dQ=malloc(6*6*sizeof(double));
            Q=malloc(6*6*sizeof(double));
            for (www=0; www<36; www++) {
                dQ[www]=0.00000000000000;
                Q[www]=0.00000000000000;
            }
            // this is the sixth row equal to zero
            Q[0]=Q11;
            Q[1]=Q12;
            Q[2]=Q16;
            Q[6]=Q12;
            Q[7]=Q22;
            Q[8]=Q26;
            Q[12]=Q16;
            Q[13]=Q26;
            Q[14]=Q66;
            Q[21]=Q44;
            Q[22]=Q45;
            Q[27]=Q45;
            Q[28]=Q55;
            
            
            dQ[0]=Q11d;
            dQ[1]=Q12d;
            dQ[2]=Q16d;
            dQ[6]=Q12d;
            dQ[7]=Q22d;
            dQ[8]=Q26d;
            dQ[12]=Q16d;
            dQ[13]=Q26d;
            dQ[14]=Q66d;
            dQ[21]=Q44d;
            dQ[22]=Q45d;
            dQ[27]=Q45d;
            dQ[28]=Q55d;
             
             
            ///*  this is the third row equal to 0
            /*
            dQ[0]=Q11d;
            dQ[1]=Q12d;
            dQ[5]=Q16d;
            dQ[6]=Q12d;
            dQ[7]=Q22d;
            dQ[11]=Q26d;
            dQ[30]=Q16d;
            dQ[31]=Q26d;
            dQ[35]=Q66d;
            dQ[21]=Q44d;
            dQ[22]=Q45d;
            dQ[27]=Q45d;
            dQ[28]=Q55d;
            
            */

			P11d = -4*C3*S*F11 + 4*(C3*S-C*S3)*F12 + 4*S3*C*F22 + 2*(C3*S-C*S3)*F66;
			P12d = 2*(C3*S-C*S3)*F11 + 4*(C*S3 - C3*S)*F12 + 2*(C3*S-C*S3)*F22 - 2*(C3*S-S3*C)*F66;
			P16d = 2*(C4-3*C2*S2)*F11 + 2*(6*C2*S2-S4-C4)*F12 - 2*(3*S2*C2-S4)*F22 + (6*C2*S2-C4-S4)*F66;
			P22d = 4*S3*C*F11 + 4*(C3*S-C*S3)*F12 - 4*C3*S*F22 + 2*(C3*S-C*S3)*F66;
			P26d = 2*(3*S2*C2-S4)*F11 + 2*(C4+S4-6*C2*S2)*F12 - 2*(C4-3*S2*C2)*F22 + (C4+S4-6*C2*S2)*F66;
			P66d = 8*(C3*S-S3*C)*F11 - 16*(C3*S-S3*C)*F12 + 8*(C3*S-S3*C)*F22 + 8*(S3*C-C3*S)*F66;
			P44d = 2*S*C*(F55-F44);
			P55d = 2*S*C*(F44-F55);
			P45d = (C2-S2)*(F44-F55);
            
            
            double *dP;
            
            dP=malloc(6*6*sizeof(double));
            for (www=0; www<36; www++) {
                dP[www]=0.00000000000000;
            }
            // this is the sixth row equal to zero
            
            dP[0]=P11d;
            dP[1]=P12d;
            dP[2]=P16d;
            dP[6]=P12d;
            dP[7]=P22d;
            dP[8]=P26d;
            dP[12]=P16d;
            dP[13]=P26d;
            dP[14]=P66d;
            dP[21]=P44d;
            dP[22]=P45d;
            dP[27]=P45d;
            dP[28]=P55d;
            
			/*Now go through each node calculating the gauss matrices and from it the sensitvity contribution of each node*/
			for(i=0;i<4;i++)
			{
				Stn = 0;
				Stn1 = 0;	/*Only need these for each node so reset them. Varible Stn is the one carried over all four nodes*/
				Stn2 = 0;
				Stn3 = 0;
				StnB = 0;
				StnC = 0;
				StnD = 0;
				StnX = 0;
				StnY = 0;
				StnZ = 0;
				StnRx = 0;
				StnRy = 0;

				/*Get Gauss point cordinates*/
				gax = gaB * (double)PoA[i].x * h2xz;
				gay = gaB * (double)PoA[i].y * h2;
				GaussValue[Gcount].x = cx + gax;
				GaussValue[Gcount].y = cy + gay;

				/*Get the non-differential guass matrix fot the strain terms*/				
				CompositeKEGaussDIFF(i, Q11, Q22, Q12, Q16, Q26, Q66, Q44, Q55, Q45, KEb, numply, plyt, h, hxz[n],num);
				/*Get the differential guass matrix for the strain terms*/
				CompositeKEGaussDIFF(i, Q11d, Q22d, Q12d, Q16d, Q26d, Q66d, Q44d, Q55d, Q45d, KE, numply, plyt, h, hxz[n],num);
                //printf("KE: %f", KE[0]);
				/*if((num==0))
				{
					printf("\n");	
					for(x=0;x<24;x++)
					{
						printf("\n");
						for(y=0;y<24;y++)
						{
							printf("%0.2f  ", KE[x+y*24]);

						}
					}	
					printf("\n");			
				}*/
				/*Get the non-differential guass matrix fot the stress terms*/				
				CompositeKEGaussDIFF(i, P11, P22, P12, P16, P26, P66, P44, P55, P45, SEb, numply, plyt, h, hxz[n],num);
				/*Get the guass matrix for stress terms*/
				CompositeKEGaussDIFF(i, P11d, P22d, P12d, P16d, P26d, P66d, P44d, P55d, P45d, SE, numply, plyt, h, hxz[n],5);
                // The term  num in here is useless , no meanning at all.
         
				/*Now calculate the sensitvtiy contribution of the current node*/
				for(j=0;j<4;j++)
				{
					/*if((num==0)){printf("\n(%i,%i)", i, j);}*/

					k= 6*j;		/*k,j and i are used for indexing purposes*/
				/*-----------E'CE-----------*/
					/*UX Term*/
					StnB += UX[j] * (UX[0]*KEb[0+k] + UY[0]*KEb[24+k] + UZ[0]*KEb[48+k] + RX[0]*KEb[72+k] + RY[0]*KEb[96+k] + RZ[0]*KEb[120+k] + UX[1]*KEb[144+k] + UY[1]*KEb[168+k] + UZ[1]*KEb[192+k] + RX[1]*KEb[216+k] + RY[1]*KEb[240+k] + RZ[1]*KEb[264+k] + UX[2]*KEb[288+k] + UY[2]*KEb[312+k] + UZ[2]*KEb[336+k] + RX[2]*KEb[360+k] + RY[2]*KEb[384+k] + RZ[2]*KEb[408+k] + UX[3]*KEb[432+k] + UY[3]*KEb[456+k] + UZ[3]*KEb[480+k] + RX[3]*KEb[504+k] + RY[3]*KEb[528+k] + RZ[3]*KEb[552+k]);
					/*UY Term*/																																																					
					StnB += UY[j] * (UX[0]*KEb[1+k] + UY[0]*KEb[25+k] + UZ[0]*KEb[49+k] + RX[0]*KEb[73+k] + RY[0]*KEb[97+k] + RZ[0]*KEb[121+k] + UX[1]*KEb[145+k] + UY[1]*KEb[169+k] + UZ[1]*KEb[193+k] + RX[1]*KEb[217+k] + RY[1]*KEb[241+k] + RZ[1]*KEb[265+k] + UX[2]*KEb[289+k] + UY[2]*KEb[313+k] + UZ[2]*KEb[337+k] + RX[2]*KEb[361+k] + RY[2]*KEb[385+k] + RZ[2]*KEb[409+k] + UX[3]*KEb[433+k] + UY[3]*KEb[457+k] + UZ[3]*KEb[481+k] + RX[3]*KEb[505+k] + RY[3]*KEb[529+k] + RZ[3]*KEb[553+k]);
					/*UZ Term*/																																																					
					StnB += UZ[j] * (UX[0]*KEb[2+k] + UY[0]*KEb[26+k] + UZ[0]*KEb[50+k] + RX[0]*KEb[74+k] + RY[0]*KEb[98+k] + RZ[0]*KEb[122+k] + UX[1]*KEb[146+k] + UY[1]*KEb[170+k] + UZ[1]*KEb[194+k] + RX[1]*KEb[218+k] + RY[1]*KEb[242+k] + RZ[1]*KEb[266+k] + UX[2]*KEb[290+k] + UY[2]*KEb[314+k] + UZ[2]*KEb[338+k] + RX[2]*KEb[362+k] + RY[2]*KEb[386+k] + RZ[2]*KEb[410+k] + UX[3]*KEb[434+k] + UY[3]*KEb[458+k] + UZ[3]*KEb[482+k] + RX[3]*KEb[506+k] + RY[3]*KEb[530+k] + RZ[3]*KEb[554+k]);
					/*RX Term*/																																																					
					StnB += RX[j] * (UX[0]*KEb[3+k] + UY[0]*KEb[27+k] + UZ[0]*KEb[51+k] + RX[0]*KEb[75+k] + RY[0]*KEb[99+k] + RZ[0]*KEb[123+k] + UX[1]*KEb[147+k] + UY[1]*KEb[171+k] + UZ[1]*KEb[195+k] + RX[1]*KEb[219+k] + RY[1]*KEb[243+k] + RZ[1]*KEb[267+k] + UX[2]*KEb[291+k] + UY[2]*KEb[315+k] + UZ[2]*KEb[339+k] + RX[2]*KEb[363+k] + RY[2]*KEb[387+k] + RZ[2]*KEb[411+k] + UX[3]*KEb[435+k] + UY[3]*KEb[459+k] + UZ[3]*KEb[483+k] + RX[3]*KEb[507+k] + RY[3]*KEb[531+k] + RZ[3]*KEb[555+k]);
					/*RY Term*/																																																					
					StnB += RY[j] * (UX[0]*KEb[4+k] + UY[0]*KEb[28+k] + UZ[0]*KEb[52+k] + RX[0]*KEb[76+k] + RY[0]*KEb[100+k] + RZ[0]*KEb[124+k] + UX[1]*KEb[148+k] + UY[1]*KEb[172+k] + UZ[1]*KEb[196+k] + RX[1]*KEb[220+k] + RY[1]*KEb[244+k] + RZ[1]*KEb[268+k] + UX[2]*KEb[292+k] + UY[2]*KEb[316+k] + UZ[2]*KEb[340+k] + RX[2]*KEb[364+k] + RY[2]*KEb[388+k] + RZ[2]*KEb[412+k] + UX[3]*KEb[436+k] + UY[3]*KEb[460+k] + UZ[3]*KEb[484+k] + RX[3]*KEb[508+k] + RY[3]*KEb[532+k] + RZ[3]*KEb[556+k]);
	
				/*-----------O'SO-----------*/
					/*OX Term*/
					StnC += OX[j] * (OX[0]*SEb[0+k] + OY[0]*SEb[24+k] + OZ[0]*SEb[48+k] + ORX[0]*SEb[72+k] + ORY[0]*SEb[96+k] + RZ[0]*SEb[120+k] + OX[1]*SEb[144+k] + OY[1]*SEb[168+k] + OZ[1]*SEb[192+k] + ORX[1]*SEb[216+k] + ORY[1]*SEb[240+k] + RZ[1]*SEb[264+k] + OX[2]*SEb[288+k] + OY[2]*SEb[312+k] + OZ[2]*SEb[336+k] + ORX[2]*SEb[360+k] + ORY[2]*SEb[384+k] + RZ[2]*SEb[408+k] + OX[3]*SEb[432+k] + OY[3]*SEb[456+k] + OZ[3]*SEb[480+k] + ORX[3]*SEb[504+k] + ORY[3]*SEb[528+k] + RZ[3]*SEb[552+k]);
					/*OY Term*/																																																					
					StnC += OY[j] * (OX[0]*SEb[1+k] + OY[0]*SEb[25+k] + OZ[0]*SEb[49+k] + ORX[0]*SEb[73+k] + ORY[0]*SEb[97+k] + RZ[0]*SEb[121+k] + OX[1]*SEb[145+k] + OY[1]*SEb[169+k] + OZ[1]*SEb[193+k] + ORX[1]*SEb[217+k] + ORY[1]*SEb[241+k] + RZ[1]*SEb[265+k] + OX[2]*SEb[289+k] + OY[2]*SEb[313+k] + OZ[2]*SEb[337+k] + ORX[2]*SEb[361+k] + ORY[2]*SEb[385+k] + RZ[2]*SEb[409+k] + OX[3]*SEb[433+k] + OY[3]*SEb[457+k] + OZ[3]*SEb[481+k] + ORX[3]*SEb[505+k] + ORY[3]*SEb[529+k] + RZ[3]*SEb[553+k]);
					/*OZ Term*/																																																					
					StnC += OZ[j] * (OX[0]*SEb[2+k] + OY[0]*SEb[26+k] + OZ[0]*SEb[50+k] + ORX[0]*SEb[74+k] + ORY[0]*SEb[98+k] + RZ[0]*SEb[122+k] + OX[1]*SEb[146+k] + OY[1]*SEb[170+k] + OZ[1]*SEb[194+k] + ORX[1]*SEb[218+k] + ORY[1]*SEb[242+k] + RZ[1]*SEb[266+k] + OX[2]*SEb[290+k] + OY[2]*SEb[314+k] + OZ[2]*SEb[338+k] + ORX[2]*SEb[362+k] + ORY[2]*SEb[386+k] + RZ[2]*SEb[410+k] + OX[3]*SEb[434+k] + OY[3]*SEb[458+k] + OZ[3]*SEb[482+k] + ORX[3]*SEb[506+k] + ORY[3]*SEb[530+k] + RZ[3]*SEb[554+k]);
					/*ORX Term*/																																																					
					StnC += ORX[j] * (OX[0]*SEb[3+k] + OY[0]*SEb[27+k] + OZ[0]*SEb[51+k] + ORX[0]*SEb[75+k] + ORY[0]*SEb[99+k] + RZ[0]*SEb[123+k] + OX[1]*SEb[147+k] + OY[1]*SEb[171+k] + OZ[1]*SEb[195+k] + ORX[1]*SEb[219+k] + ORY[1]*SEb[243+k] + RZ[1]*SEb[267+k] + OX[2]*SEb[291+k] + OY[2]*SEb[315+k] + OZ[2]*SEb[339+k] + ORX[2]*SEb[363+k] + ORY[2]*SEb[387+k] + RZ[2]*SEb[411+k] + OX[3]*SEb[435+k] + OY[3]*SEb[459+k] + OZ[3]*SEb[483+k] + ORX[3]*SEb[507+k] + ORY[3]*SEb[531+k] + RZ[3]*SEb[555+k]);
					/*ORY Term*/																																																					
					StnC += ORY[j] * (OX[0]*SEb[4+k] + OY[0]*SEb[28+k] + OZ[0]*SEb[52+k] + ORX[0]*SEb[76+k] + ORY[0]*SEb[100+k] + RZ[0]*SEb[124+k] + OX[1]*SEb[148+k] + OY[1]*SEb[172+k] + OZ[1]*SEb[196+k] + ORX[1]*SEb[220+k] + ORY[1]*SEb[244+k] + RZ[1]*SEb[268+k] + OX[2]*SEb[292+k] + OY[2]*SEb[316+k] + OZ[2]*SEb[340+k] + ORX[2]*SEb[364+k] + ORY[2]*SEb[388+k] + RZ[2]*SEb[412+k] + OX[3]*SEb[436+k] + OY[3]*SEb[460+k] + OZ[3]*SEb[484+k] + ORX[3]*SEb[508+k] + ORY[3]*SEb[532+k] + RZ[3]*SEb[556+k]);

				/*-----------E'CO-----------*/
					/*OX Term*/
					StnD += UX[j] * (OX[0]*KEb[0+k] + OY[0]*KEb[24+k] + OZ[0]*KEb[48+k] + ORX[0]*KEb[72+k] + ORY[0]*KEb[96+k] + RZ[0]*KEb[120+k] + OX[1]*KEb[144+k] + OY[1]*KEb[168+k] + OZ[1]*KEb[192+k] + ORX[1]*KEb[216+k] + ORY[1]*KEb[240+k] + RZ[1]*KEb[264+k] + OX[2]*KEb[288+k] + OY[2]*KEb[312+k] + OZ[2]*KEb[336+k] + ORX[2]*KEb[360+k] + ORY[2]*KEb[384+k] + RZ[2]*KEb[408+k] + OX[3]*KEb[432+k] + OY[3]*KEb[456+k] + OZ[3]*KEb[480+k] + ORX[3]*KEb[504+k] + ORY[3]*KEb[528+k] + RZ[3]*KEb[552+k]);
					/*OY Term*/																																																					
					StnD += UY[j] * (OX[0]*KEb[1+k] + OY[0]*KEb[25+k] + OZ[0]*KEb[49+k] + ORX[0]*KEb[73+k] + ORY[0]*KEb[97+k] + RZ[0]*KEb[121+k] + OX[1]*KEb[145+k] + OY[1]*KEb[169+k] + OZ[1]*KEb[193+k] + ORX[1]*KEb[217+k] + ORY[1]*KEb[241+k] + RZ[1]*KEb[265+k] + OX[2]*KEb[289+k] + OY[2]*KEb[313+k] + OZ[2]*KEb[337+k] + ORX[2]*KEb[361+k] + ORY[2]*KEb[385+k] + RZ[2]*KEb[409+k] + OX[3]*KEb[433+k] + OY[3]*KEb[457+k] + OZ[3]*KEb[481+k] + ORX[3]*KEb[505+k] + ORY[3]*KEb[529+k] + RZ[3]*KEb[553+k]);
					/*OZ Term*/																																																					
					StnD += UZ[j] * (OX[0]*KEb[2+k] + OY[0]*KEb[26+k] + OZ[0]*KEb[50+k] + ORX[0]*KEb[74+k] + ORY[0]*KEb[98+k] + RZ[0]*KEb[122+k] + OX[1]*KEb[146+k] + OY[1]*KEb[170+k] + OZ[1]*KEb[194+k] + ORX[1]*KEb[218+k] + ORY[1]*KEb[242+k] + RZ[1]*KEb[266+k] + OX[2]*KEb[290+k] + OY[2]*KEb[314+k] + OZ[2]*KEb[338+k] + ORX[2]*KEb[362+k] + ORY[2]*KEb[386+k] + RZ[2]*KEb[410+k] + OX[3]*KEb[434+k] + OY[3]*KEb[458+k] + OZ[3]*KEb[482+k] + ORX[3]*KEb[506+k] + ORY[3]*KEb[530+k] + RZ[3]*KEb[554+k]);
					/*ORX Term*/																																																					
					StnD += RX[j] * (OX[0]*KEb[3+k] + OY[0]*KEb[27+k] + OZ[0]*KEb[51+k] + ORX[0]*KEb[75+k] + ORY[0]*KEb[99+k] + RZ[0]*KEb[123+k] + OX[1]*KEb[147+k] + OY[1]*KEb[171+k] + OZ[1]*KEb[195+k] + ORX[1]*KEb[219+k] + ORY[1]*KEb[243+k] + RZ[1]*KEb[267+k] + OX[2]*KEb[291+k] + OY[2]*KEb[315+k] + OZ[2]*KEb[339+k] + ORX[2]*KEb[363+k] + ORY[2]*KEb[387+k] + RZ[2]*KEb[411+k] + OX[3]*KEb[435+k] + OY[3]*KEb[459+k] + OZ[3]*KEb[483+k] + ORX[3]*KEb[507+k] + ORY[3]*KEb[531+k] + RZ[3]*KEb[555+k]);
					/*ORY Term*/																																																					
					StnD += RY[j] * (OX[0]*KEb[4+k] + OY[0]*KEb[28+k] + OZ[0]*KEb[52+k] + ORX[0]*KEb[76+k] + ORY[0]*KEb[100+k] + RZ[0]*KEb[124+k] + OX[1]*KEb[148+k] + OY[1]*KEb[172+k] + OZ[1]*KEb[196+k] + ORX[1]*KEb[220+k] + ORY[1]*KEb[244+k] + RZ[1]*KEb[268+k] + OX[2]*KEb[292+k] + OY[2]*KEb[316+k] + OZ[2]*KEb[340+k] + ORX[2]*KEb[364+k] + ORY[2]*KEb[388+k] + RZ[2]*KEb[412+k] + OX[3]*KEb[436+k] + OY[3]*KEb[460+k] + OZ[3]*KEb[484+k] + ORX[3]*KEb[508+k] + ORY[3]*KEb[532+k] + RZ[3]*KEb[556+k]);


				/*-----------E'(dC/dtheta)E----------*/
					/*UX Term*/
					Stn1 += UX[j] * (UX[0]*KE[0+k] + UY[0]*KE[24+k] + UZ[0]*KE[48+k] + RX[0]*KE[72+k] + RY[0]*KE[96+k] + RZ[0]*KE[120+k] + UX[1]*KE[144+k] + UY[1]*KE[168+k] + UZ[1]*KE[192+k] + RX[1]*KE[216+k] + RY[1]*KE[240+k] + RZ[1]*KE[264+k] + UX[2]*KE[288+k] + UY[2]*KE[312+k] + UZ[2]*KE[336+k] + RX[2]*KE[360+k] + RY[2]*KE[384+k] + RZ[2]*KE[408+k] + UX[3]*KE[432+k] + UY[3]*KE[456+k] + UZ[3]*KE[480+k] + RX[3]*KE[504+k] + RY[3]*KE[528+k] + RZ[3]*KE[552+k]);
					/*UY Term*/																																																					
					Stn1 += UY[j] * (UX[0]*KE[1+k] + UY[0]*KE[25+k] + UZ[0]*KE[49+k] + RX[0]*KE[73+k] + RY[0]*KE[97+k] + RZ[0]*KE[121+k] + UX[1]*KE[145+k] + UY[1]*KE[169+k] + UZ[1]*KE[193+k] + RX[1]*KE[217+k] + RY[1]*KE[241+k] + RZ[1]*KE[265+k] + UX[2]*KE[289+k] + UY[2]*KE[313+k] + UZ[2]*KE[337+k] + RX[2]*KE[361+k] + RY[2]*KE[385+k] + RZ[2]*KE[409+k] + UX[3]*KE[433+k] + UY[3]*KE[457+k] + UZ[3]*KE[481+k] + RX[3]*KE[505+k] + RY[3]*KE[529+k] + RZ[3]*KE[553+k]);
					/*UZ Term*/																																																					
					Stn1 += UZ[j] * (UX[0]*KE[2+k] + UY[0]*KE[26+k] + UZ[0]*KE[50+k] + RX[0]*KE[74+k] + RY[0]*KE[98+k] + RZ[0]*KE[122+k] + UX[1]*KE[146+k] + UY[1]*KE[170+k] + UZ[1]*KE[194+k] + RX[1]*KE[218+k] + RY[1]*KE[242+k] + RZ[1]*KE[266+k] + UX[2]*KE[290+k] + UY[2]*KE[314+k] + UZ[2]*KE[338+k] + RX[2]*KE[362+k] + RY[2]*KE[386+k] + RZ[2]*KE[410+k] + UX[3]*KE[434+k] + UY[3]*KE[458+k] + UZ[3]*KE[482+k] + RX[3]*KE[506+k] + RY[3]*KE[530+k] + RZ[3]*KE[554+k]);
					/*RX Term*/																																																					
					Stn1 += RX[j] * (UX[0]*KE[3+k] + UY[0]*KE[27+k] + UZ[0]*KE[51+k] + RX[0]*KE[75+k] + RY[0]*KE[99+k] + RZ[0]*KE[123+k] + UX[1]*KE[147+k] + UY[1]*KE[171+k] + UZ[1]*KE[195+k] + RX[1]*KE[219+k] + RY[1]*KE[243+k] + RZ[1]*KE[267+k] + UX[2]*KE[291+k] + UY[2]*KE[315+k] + UZ[2]*KE[339+k] + RX[2]*KE[363+k] + RY[2]*KE[387+k] + RZ[2]*KE[411+k] + UX[3]*KE[435+k] + UY[3]*KE[459+k] + UZ[3]*KE[483+k] + RX[3]*KE[507+k] + RY[3]*KE[531+k] + RZ[3]*KE[555+k]);
					/*RY Term*/																																																					
					Stn1 += RY[j] * (UX[0]*KE[4+k] + UY[0]*KE[28+k] + UZ[0]*KE[52+k] + RX[0]*KE[76+k] + RY[0]*KE[100+k] + RZ[0]*KE[124+k] + UX[1]*KE[148+k] + UY[1]*KE[172+k] + UZ[1]*KE[196+k] + RX[1]*KE[220+k] + RY[1]*KE[244+k] + RZ[1]*KE[268+k] + UX[2]*KE[292+k] + UY[2]*KE[316+k] + UZ[2]*KE[340+k] + RX[2]*KE[364+k] + RY[2]*KE[388+k] + RZ[2]*KE[412+k] + UX[3]*KE[436+k] + UY[3]*KE[460+k] + UZ[3]*KE[484+k] + RX[3]*KE[508+k] + RY[3]*KE[532+k] + RZ[3]*KE[556+k]);

				/*-----------O'(dS/dtheta)O-----------*/
					/*OX Term*/
					Stn2 += OX[j] * (OX[0]*SE[0+k] + OY[0]*SE[24+k] + OZ[0]*SE[48+k] + ORX[0]*SE[72+k] + ORY[0]*SE[96+k] + RZ[0]*SE[120+k] + OX[1]*SE[144+k] + OY[1]*SE[168+k] + OZ[1]*SE[192+k] + ORX[1]*SE[216+k] + ORY[1]*SE[240+k] + RZ[1]*SE[264+k] + OX[2]*SE[288+k] + OY[2]*SE[312+k] + OZ[2]*SE[336+k] + ORX[2]*SE[360+k] + ORY[2]*SE[384+k] + RZ[2]*SE[408+k] + OX[3]*SE[432+k] + OY[3]*SE[456+k] + OZ[3]*SE[480+k] + ORX[3]*SE[504+k] + ORY[3]*SE[528+k] + RZ[3]*SE[552+k]);
					/*OY Term*/																																																					
					Stn2 += OY[j] * (OX[0]*SE[1+k] + OY[0]*SE[25+k] + OZ[0]*SE[49+k] + ORX[0]*SE[73+k] + ORY[0]*SE[97+k] + RZ[0]*SE[121+k] + OX[1]*SE[145+k] + OY[1]*SE[169+k] + OZ[1]*SE[193+k] + ORX[1]*SE[217+k] + ORY[1]*SE[241+k] + RZ[1]*SE[265+k] + OX[2]*SE[289+k] + OY[2]*SE[313+k] + OZ[2]*SE[337+k] + ORX[2]*SE[361+k] + ORY[2]*SE[385+k] + RZ[2]*SE[409+k] + OX[3]*SE[433+k] + OY[3]*SE[457+k] + OZ[3]*SE[481+k] + ORX[3]*SE[505+k] + ORY[3]*SE[529+k] + RZ[3]*SE[553+k]);
					/*OZ Term*/																																																					
					Stn2 += OZ[j] * (OX[0]*SE[2+k] + OY[0]*SE[26+k] + OZ[0]*SE[50+k] + ORX[0]*SE[74+k] + ORY[0]*SE[98+k] + RZ[0]*SE[122+k] + OX[1]*SE[146+k] + OY[1]*SE[170+k] + OZ[1]*SE[194+k] + ORX[1]*SE[218+k] + ORY[1]*SE[242+k] + RZ[1]*SE[266+k] + OX[2]*SE[290+k] + OY[2]*SE[314+k] + OZ[2]*SE[338+k] + ORX[2]*SE[362+k] + ORY[2]*SE[386+k] + RZ[2]*SE[410+k] + OX[3]*SE[434+k] + OY[3]*SE[458+k] + OZ[3]*SE[482+k] + ORX[3]*SE[506+k] + ORY[3]*SE[530+k] + RZ[3]*SE[554+k]);
					/*ORX Term*/																																																					
					Stn2 += ORX[j] * (OX[0]*SE[3+k] + OY[0]*SE[27+k] + OZ[0]*SE[51+k] + ORX[0]*SE[75+k] + ORY[0]*SE[99+k] + RZ[0]*SE[123+k] + OX[1]*SE[147+k] + OY[1]*SE[171+k] + OZ[1]*SE[195+k] + ORX[1]*SE[219+k] + ORY[1]*SE[243+k] + RZ[1]*SE[267+k] + OX[2]*SE[291+k] + OY[2]*SE[315+k] + OZ[2]*SE[339+k] + ORX[2]*SE[363+k] + ORY[2]*SE[387+k] + RZ[2]*SE[411+k] + OX[3]*SE[435+k] + OY[3]*SE[459+k] + OZ[3]*SE[483+k] + ORX[3]*SE[507+k] + ORY[3]*SE[531+k] + RZ[3]*SE[555+k]);
					/*ORY Term*/																																																					
					Stn2 += ORY[j] * (OX[0]*SE[4+k] + OY[0]*SE[28+k] + OZ[0]*SE[52+k] + ORX[0]*SE[76+k] + ORY[0]*SE[100+k] + RZ[0]*SE[124+k] + OX[1]*SE[148+k] + OY[1]*SE[172+k] + OZ[1]*SE[196+k] + ORX[1]*SE[220+k] + ORY[1]*SE[244+k] + RZ[1]*SE[268+k] + OX[2]*SE[292+k] + OY[2]*SE[316+k] + OZ[2]*SE[340+k] + ORX[2]*SE[364+k] + ORY[2]*SE[388+k] + RZ[2]*SE[412+k] + OX[3]*SE[436+k] + OY[3]*SE[460+k] + OZ[3]*SE[484+k] + ORX[3]*SE[508+k] + ORY[3]*SE[532+k] + RZ[3]*SE[556+k]);

				/*-----------E'(dC/dtheta)O-----------*/
					/*OX Term*/
					Stn3 += UX[j] * (OX[0]*KE[0+k] + OY[0]*KE[24+k] + OZ[0]*KE[48+k] + ORX[0]*KE[72+k] + ORY[0]*KE[96+k] + ORZ[0]*KE[120+k] + OX[1]*KE[144+k] + OY[1]*KE[168+k] + OZ[1]*KE[192+k] + ORX[1]*KE[216+k] + ORY[1]*KE[240+k] + ORZ[1]*KE[264+k] + OX[2]*KE[288+k] + OY[2]*KE[312+k] + OZ[2]*KE[336+k] + ORX[2]*KE[360+k] + ORY[2]*KE[384+k] + ORZ[2]*KE[408+k] + OX[3]*KE[432+k] + OY[3]*KE[456+k] + OZ[3]*KE[480+k] + ORX[3]*KE[504+k] + ORY[3]*KE[528+k] + ORZ[3]*KE[552+k]);
					/*OY Term*/																																																					
					Stn3 += UY[j] * (OX[0]*KE[1+k] + OY[0]*KE[25+k] + OZ[0]*KE[49+k] + ORX[0]*KE[73+k] + ORY[0]*KE[97+k] + ORZ[0]*KE[121+k] + OX[1]*KE[145+k] + OY[1]*KE[169+k] + OZ[1]*KE[193+k] + ORX[1]*KE[217+k] + ORY[1]*KE[241+k] + ORZ[1]*KE[265+k] + OX[2]*KE[289+k] + OY[2]*KE[313+k] + OZ[2]*KE[337+k] + ORX[2]*KE[361+k] + ORY[2]*KE[385+k] + ORZ[2]*KE[409+k] + OX[3]*KE[433+k] + OY[3]*KE[457+k] + OZ[3]*KE[481+k] + ORX[3]*KE[505+k] + ORY[3]*KE[529+k] + ORZ[3]*KE[553+k]);
					/*OZ Term*/																																																					
					Stn3 += UZ[j] * (OX[0]*KE[2+k] + OY[0]*KE[26+k] + OZ[0]*KE[50+k] + ORX[0]*KE[74+k] + ORY[0]*KE[98+k] + ORZ[0]*KE[122+k] + OX[1]*KE[146+k] + OY[1]*KE[170+k] + OZ[1]*KE[194+k] + ORX[1]*KE[218+k] + ORY[1]*KE[242+k] + ORZ[1]*KE[266+k] + OX[2]*KE[290+k] + OY[2]*KE[314+k] + OZ[2]*KE[338+k] + ORX[2]*KE[362+k] + ORY[2]*KE[386+k] + ORZ[2]*KE[410+k] + OX[3]*KE[434+k] + OY[3]*KE[458+k] + OZ[3]*KE[482+k] + ORX[3]*KE[506+k] + ORY[3]*KE[530+k] + ORZ[3]*KE[554+k]);
					/*ORX Term*/																																																					
					Stn3 += RX[j] * (OX[0]*KE[3+k] + OY[0]*KE[27+k] + OZ[0]*KE[51+k] + ORX[0]*KE[75+k] + ORY[0]*KE[99+k] + ORZ[0]*KE[123+k] + OX[1]*KE[147+k] + OY[1]*KE[171+k] + OZ[1]*KE[195+k] + ORX[1]*KE[219+k] + ORY[1]*KE[243+k] + ORZ[1]*KE[267+k] + OX[2]*KE[291+k] + OY[2]*KE[315+k] + OZ[2]*KE[339+k] + ORX[2]*KE[363+k] + ORY[2]*KE[387+k] + ORZ[2]*KE[411+k] + OX[3]*KE[435+k] + OY[3]*KE[459+k] + OZ[3]*KE[483+k] + ORX[3]*KE[507+k] + ORY[3]*KE[531+k] + ORZ[3]*KE[555+k]);
					/*ORY Term*/																																																					
					Stn3 += RY[j] * (OX[0]*KE[4+k] + OY[0]*KE[28+k] + OZ[0]*KE[52+k] + ORX[0]*KE[76+k] + ORY[0]*KE[100+k] + ORZ[0]*KE[124+k] + OX[1]*KE[148+k] + OY[1]*KE[172+k] + OZ[1]*KE[196+k] + ORX[1]*KE[220+k] + ORY[1]*KE[244+k] + ORZ[1]*KE[268+k] + OX[2]*KE[292+k] + OY[2]*KE[316+k] + OZ[2]*KE[340+k] + ORX[2]*KE[364+k] + ORY[2]*KE[388+k] + ORZ[2]*KE[412+k] + OX[3]*KE[436+k] + OY[3]*KE[460+k] + OZ[3]*KE[484+k] + ORX[3]*KE[508+k] + ORY[3]*KE[532+k] + ORZ[3]*KE[556+k]);

					/*Stn1 += StnX+StnY+StnZ+StnRx+StnRy;*/
					/*if((num==83)||(num==723)){
					printf("\n\n num: %i", num);
					printf("\n(%i,%i)", n, m);
					printf("\nNode 0: UX: %f, UY: %f, UZ: %f RX: %f RY: %f RZ: %f", UX[0], UY[0], UZ[0], RX[0], RY[0], RZ[0]);
					printf("\nNode 1: UX: %f, UY: %f, UZ: %f RX: %f RY: %f RZ: %f", UX[1], UY[1], UZ[1], RX[1], RY[1], RZ[1]);
					printf("\nNode 2: UX: %f, UY: %f, UZ: %f RX: %f RY: %f RZ: %f", UX[2], UY[2], UZ[2], RX[2], RY[2], RZ[2]);
					printf("\nNode 3: UX: %f, UY: %f, UZ: %f RX: %f RY: %f RZ: %f", UX[3], UY[3], UZ[3], RX[3], RY[3], RZ[3]);
					/*printf("\n\n T:(%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)", T[0],T[1],T[2],T[3],T[4],T[5],T[6],T[7],T[8],T[9],T[10],T[11],T[12]);*/
					/*printf("\nStn1: %f", Stn1);*/
					/*printf("\nStnX: %f  StnY: %f,  StnZ: %f  StnRx: %f  StnRy:  %f", StnX, StnY, StnZ, StnRx, StnRy);*/
					/*printf("\nKE[%i] = %f", 0+k, KE[0+k]);}*/
					
				}
                
                double *Bmatrix;
                Bmatrix=malloc(24*6*sizeof(double));
                double *strain1, *transtrain1;
                strain1=malloc(6*sizeof(double));
                transtrain1=malloc(6*sizeof(double));
                double *TempMatrix;
                TempMatrix=malloc(6*sizeof(double));
                //initinlize the value , actually it is not neccessary
                int ttt;
                for (ttt=0;ttt<6; ttt++) {
                    strain1[ttt]=0.0;
                    transtrain1[ttt]=0.0;
                    TempMatrix[ttt]=0.0;
                    
                }
                //get the B matrix
                double detJJ;
                detJJ=0.0;
                CompositeBmatrix(i,Bmatrix, numply,  plyt, h,  hxz[n],  num, NumNodes, NodeCoord, Lnodes, &detJJ);
                
                
                
                //Get the K matrix
                double Kmatrix[576];
                double Btransp[144];
                double BD[144];
                TransM (Bmatrix, Btransp, 24, 6);
                MatrixMultiple(Btransp, dQ, BD, 6, 24, 6 , 6);
                MatrixMultiple(BD, Bmatrix, Kmatrix, 6, 24, 24 , 6);
                
                // This K matrix need to times detJ*thickness, then will be the same as Chris' result
                //printf("Kmateix: %f", Kmatrix[0]);
                
                
                
                
                // get the strain
                MatrixMultiple(Bmatrix, disp, strain1, 24, 6, 1, 24);
                
                // get the transpose of strain
                TransM (strain1, transtrain1, 1, 6);
                
                // strain'*dC/dtheta = dtemp
                MatrixMultiple(transtrain1, dQ, TempMatrix, 6, 1, 6, 6);
                
                double *stn11;
                stn11=malloc(1*sizeof(double));
                
                //dtemp*strain=  stn1
                MatrixMultiple(TempMatrix, strain1, stn11, 6, 1, 1, 6);
                // give the value to Stn1...., starin based method done
                
                Stn1=stn11[0];
                // since , the sensitivity need to times det(J)*thickness
                // plyt is the thickness.
                Stn1=Stn1*(detJJ)*plyt;
                
                
                // stress base method beagin, stress= stiffness*strain -----------------------------
                //--------------------------------------------------------------------------------
                //get the stress= C*strain
                double *stress1, *transtress1, *TempMatrix2, *stn22;
                stress1=malloc(6*sizeof(double));
                transtress1=malloc(6*sizeof(double));
                TempMatrix2=malloc(6*sizeof(double));
                stn22=malloc(1*sizeof(double));
                
                for (ttt=0;ttt<6; ttt++) {
                    stress1[ttt]=0.0;
                    transtress1[ttt]=0.0;
                    TempMatrix2[ttt]=0.0;
                    
                }
                
                MatrixMultiple(Q, strain1, stress1, 6, 6, 1, 6);
                
                //transpose of stress   1 is colum, 6 is row
                
                TransM (stress1, transtress1, 1, 6);
                
                //stress'* dS/Dtheta= Tempstress
                MatrixMultiple(transtress1, dP, TempMatrix2, 6, 1, 6, 6);
                
                //dtemp*strain=  stn1
                MatrixMultiple(TempMatrix2, stress1, stn22, 6, 1, 1, 6);
                
                Stn2=stn22[0];
                
                Stn2=Stn2*(detJJ)*plyt;

                
                //printf("\n STn1, %f", Stn1);
                free(stn11);
                free(TempMatrix);
                free(strain1);
                free(transtrain1);
                
                free(stn22);
                free(TempMatrix2);
                free(stress1);
                free(transtress1);
                
            
                
                
               // for (ttt=0;ttt<6; ttt++) {
               //     printf("res= %f",res[ttt]);
              //  }
                
                free(Bmatrix);
                
                
                
                //Stn = (-1.0 + 2.0*EF[num] - EF[num]*EF[num])*Stn1 + EF[num]*EF[num]*Stn2 + (2.0*EF[num]*EF[num] - 2.0*EF[num])*Stn1;
                Stn = (-1.0 + 2.0*0.75 - 0.75*0.75)*Stn1 + 0.75*0.75*Stn2 + (2.0*0.75*0.75 - 2.0*0.75)*Stn1;
                
				//Stn = (-1 + 2*EF[num] - EF[num]*EF[num])*Stn1 + EF[num]*EF[num]*Stn2 + (2*EF[num] - EF[num]*EF[num])*Stn3;	/*Calculate the nodal sensitvity*/
				/*Stn = (-1 + 2*- 0.75*0.75)*Stn1 + 0.75*0.75*Stn2 + (2*0.75 - 0.75*0.75)*Stn3;*/
                
               // Stn = (-1.0 + 2.0*0.75 - 0.75*0.75)*Stn1 + 0.75*0.75*Stn2 + (2.0*0.75*0.75 - 2.0*0.75)*Stn3;
                
				/*ElemThetaStn[num] += Stn*0.25*h*hxz[n];			/*get the mean sensitvity of each gauss point in the element for the element sensitvity*/
				ElemThetaStn[num] += -Stn1*0.25;
                //if (num==80||num==81||num==40||num==41) {
                    printf("\nElement sensititvty of num: %dis equal to %.16lf",num, ElemThetaStn[num]);
               // }
				ElemThetaStnB[num] += (1/12)*(StnB+StnC+StnD)*h*hxz[n];

				GaussValue[Gcount].u = -Stn1;			/*store gauss point sensitvity for the nodal sensitivty value*/
				GaussValue[Gcount].a = 0.25*h*hxz[n];	/*Store the gauss point element area contribution*/
				Gcount++;
                
                
                
			}
            free(disp);
            free(dQ);
            free(dP);
            
		}
	}

	/*Print out files of interesting elemental values*/

	/*free the memeory we no longer need*/
	free(KE);
	free(KEb);
	free(SE);
	free(SEb);
	free(Lnodes);
	free(UX);
	free(UY);
	free(UZ);
	free(RX);
	free(RY);
	free(RZ);
	free(OX);
	free(OY);
	free(OZ);
	free(ORX);
	free(ORY);
	free(ORZ);
	
	double rad = 2.0*h;
	double ftemp;

	for(n=0;n<(Ntot_new);n++)
	{
		NodeThetaStn[n] = 0.0;
	}

	/*Now we need to get the nodal values of sensitvity from the elemental values by sending each node through the least squared routine
	Simpler then the standard 2D or 3D code as we have no AUX nodes to concern ourselves with*/
	for(n=0;n<Ntot_new;n++)
	{
        if (n<NumNodes) {
            
            cx = NodeCoord[n].xz;
            cy = NodeCoord[n].y; /*read in co-ordiantes directly*/
            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
           NodeThetaStn[n]= ftemp; /* multiply smothed strain energy by weight */
        }
        else{
            cx = AuxNodes[n-NumNodes].x;
            cy = AuxNodes[n-NumNodes].y; /*read in co-ordiantes directly*/
            
           // if (AuxNodes[n-NumNodes].n+1>NumNodes)    // if number is larger than NumNodes, it haven't been calculate

            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
           NodeThetaStn[n]= ftemp; /* multiply smothed strain energy by weight */

        }
		/*use least squares to calculate nodal strain energy*/
	}
    
    for(n=0;n<4*elemX*elemY;n++)
    {
        Nodeseng[n].u=GaussValue[n].u;
        Nodeseng[n].a=GaussValue[n].a;
        Nodeseng[n].x=GaussValue[n].x;
        Nodeseng[n].y=GaussValue[n].y;
        
        
        
    }
    
	free(GaussValue);
	
	/*Now send the nodal sensitvities back to the elemental sensitvities*/

	/*sprintf(plotname,"ElemThetaStn%ia.txt",itt); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Element Strain Energy plotting writefile\n");
		}
	else{
	
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;
				fprintf(outfile,"%0.5f\t",ElemThetaStn[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);*/

	for(i=0;i<elemX;i++)
	{
		for(j=0;j<elemY;j++)
		{
			/*Get the element numbers and the nodes attached to it*/
			num = Number[i][j].n-1;
			n = Number[i][j].a-1;
			m = Number[i][j].b-1;
			p = Number[i][j].c-1;
			q = Number[i][j].d-1;
			/*Element sensitvity is the mean of the nodal senstivties*/
			//ElemThetaStn[num] = (NodeThetaStn[n] + NodeThetaStn[m] + NodeThetaStn[p] + NodeThetaStn[q])*0.25;
            if (num==80) {
                //printf("element 80's sens: %.16lf", ElemThetaStn[num]);
            }
            // apparently, the elment sens calcued in here is not correct, the sens at elemnt 80
            // is 0.0052642063600315, while the ture sens with finite different is 0.08658....
		}
	}
	
	/*print out snesitvity files*
	sprintf(plotname,"NodeStn%i.txt",itt); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		fprintf(outfile,"%0.6f\t",NodeThetaStn[Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fprintf(outfile,"%0.6f\t",NodeThetaStn[Number[i][0].b-1]);
		}
		fprintf(outfile,"\n");

		for(j=0;j<elemY;j++)
		{
			fprintf(outfile,"%0.6f\t",NodeThetaStn[Number[0][j].d-1]);
			for(i=0;i<elemX;i++)
			{
				fprintf(outfile,"%0.6f\t",NodeThetaStn[Number[i][j].c-1]);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);
/*--------------------------------------------------------------------------------------------------------------*/
/*--------------Over the many iterations rounding errors actually signficantly shifts the result----------------*/
/*---------Use a Rounding System to cut the error off the end of the number thus resolving our error------------*/
	/*sprintf(plotname,"ElemThetaStn%i.txt",itt); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Element Strain Energy plotting writefile\n");
		}
	else{
	
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;
				fprintf(outfile,"%0.6f\t",ElemThetaStn[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);

	outfile = fopen(plotname, "r");
	if(outfile == NULL){
		printf("\nFailed to open Element Strain Energy plotting writefile\n");
		}
	else{
	
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;
				fscanf(outfile,"%lf\t",&ElemThetaStn[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);
	/*Need to do this twice to get all of them*/
		sprintf(plotname,"ElemThetaStn%i.txt",1/*itt*/); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Element Strain Energy plotting writefile\n");
		}
	else{
	
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;
				fprintf(outfile,"%0.7f\t",ElemThetaStn[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);

	outfile = fopen(plotname, "r");
	if(outfile == NULL){
		printf("\nFailed to open Element Strain Energy plotting writefile\n");
		}
	else{
	
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;
				fscanf(outfile,"%lf\t",&ElemThetaStn[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);

/*--------------------Rounding complete, Numerical errors should be removed now---------------------------------*/
/*--------------------------------------------------------------------------------------------------------------*

	sprintf(plotname,"ElemThetaStn%ib.txt",itt); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Element Strain Energy plotting writefile\n");
		}
	else{
	
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;
				fprintf(outfile,"%0.9f\t",ElemThetaStn[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);
	/**/
	
	//free(NodeThetaStn);
}
void Nodesensitivity( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes, double *NodeDTL, double *NodesenTemp, Gstn *NodeDTLg)
{
   
    double temp;
    double cx, cy;
    double h2 = 0.5 * h;
    double h2xz;
    int Gcount = 0;
    int gpoints = 4*elemX*elemY;
    double gax, gay;
    Gstn *GaussValue;
    int *Lnodes;
    int num;
    
    int n,m,i,j;
    Lnodes = malloc(4*sizeof(int));
    
    GaussValue = malloc(4*elemX*elemY*sizeof(Gstn));

    for(n=0;n<elemX;n++)
    {
        for(m=0;m<elemY;m++)
        {
            
            num = Number[n][m].n-1;		/*Get element number*/
            Lnodes[0] = Number[n][m].a-1;
            Lnodes[1] = Number[n][m].b-1;
            Lnodes[2] = Number[n][m].c-1;
            Lnodes[3] = Number[n][m].d-1;
            
            h2xz = 0.5 * hxz[n];
            /*work out global co-ordinates of element center, from co-ords of node 1 and element edge length*/
            double cx = NodeCoord[Lnodes[0]].xz + h2xz;
            double cy = NodeCoord[Lnodes[0]].y + h2;
            double *Nn;
            Nn = malloc(4*sizeof(double));
            /*Now go through each node calculating the gauss matrices and from it the sensitvity contribution of each node*/
            for(i=0;i<4;i++)
            {
                
                /*Get Gauss point cordinates*/
                gax = gaB * (double)PoA[i].x * h2xz;
                gay = gaB * (double)PoA[i].y * h2;
                GaussValue[Gcount].x = cx + gax;
                GaussValue[Gcount].y = cy + gay;
                Nn[0]=0.25*(1-gaB*(double)PoA[i].x)*(1-gaB* (double)PoA[i].y);
                Nn[1]=0.25*(1+gaB*(double)PoA[i].x)*(1-gaB* (double)PoA[i].y);
                Nn[2]=0.25*(1+gaB*(double)PoA[i].x)*(1+gaB* (double)PoA[i].y);
                Nn[3]=0.25*(1-gaB*(double)PoA[i].x)*(1+gaB* (double)PoA[i].y);
                
                GaussValue[Gcount].u=NodesenTemp[Lnodes[0]]*Nn[0]+NodesenTemp[Lnodes[1]]*Nn[1]+NodesenTemp[Lnodes[2]]*Nn[2]+NodesenTemp[Lnodes[3]]*Nn[3];
                GaussValue[Gcount].a = 0.25*h*hxz[n];	/*Store the gauss point element area contribution*/
                Gcount++;
            }
        }
    }

    double rad = 2.0*h;
    double ftemp;
    
    for(n=0;n<(Ntot_new);n++)
    {
        NodeDTL[n] = 0.0;
    }
    
    /*Now we need to get the nodal values of sensitvity from the elemental values by sending each node through the least squared routine
     Simpler then the standard 2D or 3D code as we have no AUX nodes to concern ourselves with*/
    for(n=0;n<Ntot_new;n++)
    {
        if (n<NumNodes) {
            
            cx = NodeCoord[n].xz;
            cy = NodeCoord[n].y; /*read in co-ordiantes directly*/
            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            NodeDTL[n]= ftemp; /* multiply smothed strain energy by weight */
        }
        else{
            cx = AuxNodes[n-NumNodes].x;
            cy = AuxNodes[n-NumNodes].y; /*read in co-ordiantes directly*/
            
            // if (AuxNodes[n-NumNodes].n+1>NumNodes)    // if number is larger than NumNodes, it haven't been calculate
            
            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            NodeDTL[n]= ftemp; /* multiply smothed strain energy by weight */
            
        }
        
    }
    for(n=0;n<4*elemX*elemY;n++)
    {
        NodeDTLg[n].u=GaussValue[n].u;
         NodeDTLg[n].a=GaussValue[n].a;
         NodeDTLg[n].x=GaussValue[n].x;
         NodeDTLg[n].y=GaussValue[n].y;
        
        
    
    }
    
    free(GaussValue);
    free(Lnodes);

}
void NodesenP( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes, double *NodeDTL, double *NodesenTemp, Gstn *NodeDTLg, Gstn *Nodeseng, double *Nodesenpp)
{
    
    double temp;
    double cx, cy;
    double h2 = 0.5 * h;
    double h2xz;
    int Gcount = 0;
    int gpoints = 4*elemX*elemY;
    double gax, gay;
    Gstn *GaussValue;
    int *Lnodes;
    int num;
    
    int n,m,i,j;
    Lnodes = malloc(4*sizeof(int));
    
    GaussValue = malloc(4*elemX*elemY*sizeof(Gstn));
    
    for(n=0;n<(4*elemX*elemY);n++)
    {
        GaussValue[n].u = NodeDTLg[n].u*Nodeseng[n].u;
        if(n==8||n==9||n==10||n==11||n==88||n==89||n==90||n==91||n==168||n==169||n==170||n==171)
        {
        printf("\n Gauss%d: %f, %f, %f",n,Nodeseng[n].u,NodeDTLg[n].u,GaussValue[n].u);
        }
        GaussValue[n].a= NodeDTLg[n].a;
        GaussValue[n].x= NodeDTLg[n].x;
        GaussValue[n].y= NodeDTLg[n].y;
        Gcount++;
    }
    
    double rad = 2.0*h;
    double ftemp;
    
    for(n=0;n<(Ntot_new);n++)
    {
        Nodesenpp[n] = 0.0;
    }
    
    /*Now we need to get the nodal values of sensitvity from the elemental values by sending each node through the least squared routine
     Simpler then the standard 2D or 3D code as we have no AUX nodes to concern ourselves with*/
    for(n=0;n<Ntot_new;n++)
    {
        if (n<NumNodes) {
            
            cx = NodeCoord[n].xz;
            cy = NodeCoord[n].y; /*read in co-ordiantes directly*/
            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            Nodesenpp[n]= ftemp; /* multiply smothed strain energy by weight */
        }
        else{
            cx = AuxNodes[n-NumNodes].x;
            cy = AuxNodes[n-NumNodes].y; /*read in co-ordiantes directly*/
            
            // if (AuxNodes[n-NumNodes].n+1>NumNodes)    // if number is larger than NumNodes, it haven't been calculate
            
            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            Nodesenpp[n]= ftemp; /* multiply smothed strain energy by weight */
            
        }
        
    }
    
    free(GaussValue);
    free(Lnodes);
    
}

/*Get the energy factor*/
void EnergyFactor(int elemX, int elemY, Elem Number[elemX][elemY], int NodeX, int NodeY, int Nodes2[NodeX][NodeY], double h, double *hxz, double *rhs1, double *rhs0, double *Theta1, double *Theta0, double e1, double e2, double v12, double v21, double g12, double g23, double g13,  double *EF, double *Trans, int numply, double plyt, int itt)
{
	FILE *outfile;	/*File varible for output files*/
	char plotname[40];	/*variable to change names of plotting output files*/
	int i, j, k, n, m, p, q, f, num,x,y;	/*interger varibles for the loops*/
	double *UXE, *UYE, *UZE, *RXE, *RYE, *RZE;			/*Local displacement values this iteration*/
	double *UX0, *UY0, *UZ0, *RX0, *RY0, *RZ0;				/*Local displacement values previous iteration*/
	double *OXE, *OYE, *OZE, *ORXE, *ORYE, *ORZE;			/*Local strain factors values for this iteration*/
	double *OX0, *OY0, *OZ0, *ORX0, *ORY0, *ORZ0;			/*Local strain factors values for previous iteration*/
	double *OX1, *OY1, *OZ1, *ORX1, *ORY1, *ORZ1;			/*Local strain factors values for imbetween iteration*/
	double *CX, *CY, *CZ, *CRX, *CRY, *CRZ;				/*terms for the top of the energy factor equation*/
	double *DX, *DY, *DZ, *DRX, *DRY, *DRZ;				/*terms for the bottom of the energy factor equation*/
	double *KE1, *KE0;
	double C, S, C4, S4, C2S2, CS, C3S, CS3, C2, S2, C3, S3, C2S, CS2; /*Cosine and Sine varibles*/
	double E11, E22, E12, E66, E44, E55;					/*Material Property varibles for strain*/
	double Q11, Q22, Q12, Q16, Q26, Q66, Q44, Q45, Q55;		/*Qmatrix varibles*/
	double Ctemp, Dtemp, Ctemp2, Dtemp2;
	int *Lnodes;

	Lnodes = malloc(4*sizeof(int));

	UXE = malloc(4*sizeof(double));
	UYE = malloc(4*sizeof(double));
	UZE = malloc(4*sizeof(double));
	RXE = malloc(4*sizeof(double));
	RYE = malloc(4*sizeof(double));
	RZE = malloc(4*sizeof(double));
	
	OXE = malloc(4*sizeof(double));
	OYE = malloc(4*sizeof(double));
	OZE = malloc(4*sizeof(double));
	ORXE = malloc(4*sizeof(double));
	ORYE = malloc(4*sizeof(double));
	ORZE = malloc(4*sizeof(double));

	OX1 = malloc(4*sizeof(double));
	OY1 = malloc(4*sizeof(double));
	OZ1 = malloc(4*sizeof(double));
	ORX1 = malloc(4*sizeof(double));
	ORY1 = malloc(4*sizeof(double));
	ORZ1 = malloc(4*sizeof(double));

	UX0 = malloc(4*sizeof(double));
	UY0 = malloc(4*sizeof(double));
	UZ0 = malloc(4*sizeof(double));
	RX0 = malloc(4*sizeof(double));
	RY0 = malloc(4*sizeof(double));
	RZ0 = malloc(4*sizeof(double));
	
	OX0 = malloc(4*sizeof(double));
	OY0 = malloc(4*sizeof(double));
	OZ0 = malloc(4*sizeof(double));
	ORX0 = malloc(4*sizeof(double));
	ORY0 = malloc(4*sizeof(double));
	ORZ0 = malloc(4*sizeof(double));

	CX = malloc(4*sizeof(double));
	CY = malloc(4*sizeof(double));
	CZ = malloc(4*sizeof(double));
	CRX = malloc(4*sizeof(double));
	CRY = malloc(4*sizeof(double));
	CRZ = malloc(4*sizeof(double));

	DX = malloc(4*sizeof(double));
	DY = malloc(4*sizeof(double));
	DZ = malloc(4*sizeof(double));
	DRX = malloc(4*sizeof(double));
	DRY = malloc(4*sizeof(double));
	DRZ = malloc(4*sizeof(double));

	KE0 = malloc(576*sizeof(double));
	KE1 = malloc(576*sizeof(double));
	
	if(itt==0)
	{
		for(i=0;i<(elemX*elemY);i++)
		{
			EF[i] = 0.75;
		}
	}
	else
	{
		/*In this loop we perform the entire snesitvity anaylisis for every single element, this involves several stages 
		to get all the element properties we need before we actually get start calculating the senstvity*/
		for(n=0;n<elemX;n++)
		{
			for(m=0;m<elemY;m++)
			{

				num = Number[n][m].n-1;		/*Get element number*/
				Lnodes[0] = Number[n][m].a-1;
				Lnodes[1] = Number[n][m].b-1;
				Lnodes[2] = Number[n][m].c-1;
				Lnodes[3] = Number[n][m].d-1;

				/*---------First we need to get the displacements local to each element from the global displacements----------*/
				p = 9*n;				/* Start position in strorage of Trans matrix for this element*/

				for(i=0;i<4;i++)	/*Get current iteration's local displacements*/
				{
					j = 6 * Lnodes[i];		/*Dofs of this node*/
					/*Transform the global displacements (rhs) to the local displacement for each node in the element*/
					UXE[i] = rhs1[j]*Trans[0+p] + rhs1[j+1]*Trans[1+p] + rhs1[j+2]*Trans[2+p];		/*local local X displacement of node a*/
					UYE[i] = rhs1[j]*Trans[3+p] + rhs1[j+1]*Trans[4+p] + rhs1[j+2]*Trans[5+p];		/*local Y displacement of node a*/
					UZE[i] = rhs1[j]*Trans[6+p] + rhs1[j+1]*Trans[7+p] + rhs1[j+2]*Trans[8+p];		/*local Z displacement of node a*/
					RXE[i] = rhs1[j+3]*Trans[0+p] + rhs1[j+4]*Trans[1+p] + rhs1[j+5]*Trans[2+p];		/*local rotX displacement of node a*/
					RYE[i] = rhs1[j+3]*Trans[3+p] + rhs1[j+4]*Trans[4+p] + rhs1[j+5]*Trans[5+p];		/*local rotY displacement of node a*/
					RZE[i] = 0;	/*rhs1[j+3]*Trans[6+p] + rhs1[j+4]*Trans[7+p] + rhs1[j+5]*Trans[8+p];*/		/*local rotZ displacement of node a*/
				}
				
				for(i=0;i<4;i++)	/*Get last iteration's local displacements*/
				{
					j = 6 * Lnodes[i];		/*Dofs of this node*/
					/*Transform the global displacements (rhs) to the local displacement for each node in the element*/
					UX0[i] = rhs0[j]*Trans[0+p] + rhs0[j+1]*Trans[1+p] + rhs0[j+2]*Trans[2+p];		/*local local X displacement of node a*/
					UY0[i] = rhs0[j]*Trans[3+p] + rhs0[j+1]*Trans[4+p] + rhs0[j+2]*Trans[5+p];		/*local Y displacement of node a*/
					UZ0[i] = rhs0[j]*Trans[6+p] + rhs0[j+1]*Trans[7+p] + rhs0[j+2]*Trans[8+p];		/*local Z displacement of node a*/
					RX0[i] = rhs0[j+3]*Trans[0+p] + rhs0[j+4]*Trans[1+p] + rhs0[j+5]*Trans[2+p];		/*local rotX displacement of node a*/
					RY0[i] = rhs0[j+3]*Trans[3+p] + rhs0[j+4]*Trans[4+p] + rhs0[j+5]*Trans[5+p];		/*local rotY displacement of node a*/
					RZ0[i] = 0;	/*rhs1[j+3]*Trans[6+p] + rhs1[j+4]*Trans[7+p] + rhs1[j+5]*Trans[8+p];*/		/*local rotZ displacement of node a*/
				}

				/*------------Now we need to change these in to the local stress from the strains-------------------------*/

				/*Get the C material Matrix properties for this element's current value of theta*/
				C = cos(Theta1[num]);
				S = sin(Theta1[num]);
				C2 = C*C;
				S2 = S*S;
				C3 = C*C*C;
				S3 = S*S*S;
				C2S = C*C*S;
				CS2 = C*S*S;
				C2S2 = C*C*S*S;
				CS = C*S;
				C3S = C*C*C*S;
				CS3 = C*S*S*S;
				C4 = C*C*C*C;
				S4 = S*S*S*S;

				/*Caclcualte material property values*/
				E11 = e1/(1-(v12*v21));
				E22 = e2/(1-(v12*v21));
				E12 = v21*E11;				/*Could also be E22*v12 as they are the same*/
				E66 = g12;
				E44 = g13;
				E55 = g23;

				/*Create C material Matrix*/
				Q11 = C4*E11 + 2*C2*S2*E12 + S4*E22 + 4*C2*S2*E66;
				Q12 = C2*S2*E11 + S4*E12 + C4*E12 + C2*S2*E22 - 4*C2*S2*E66;
				Q16 = C3*S*E11 + C*S3*E12 - C3*S*E12 - C*S3*E22 - 2*C*S*(C2-S2)*E66;
				Q22 = S4*E11 + 2*C2*S2*E12 + C4*E22 + 4*C2*S2*E66;
				Q26 = C*S3*E11 + C3*S*E12 - C*S3*E12 - C3*S*E22 + 2*C*S*(C2-S2)*E66;
				Q66 = C2S2*E11 - 2*C2S2*E12 + C2S2*E22 + (C4-2*C2S2+S4)*E66;
				Q44 = C2*E44 + S2*E55;
				Q45 = CS*E44 - CS*E55;
				Q55 = S2*E44 + C2*E55;

				/*Get current iteration stiffness Matrix*/
				CompositeKEDIFF(Q11, Q22, Q12, Q16, Q26, Q66, Q44, Q55, Q45, KE1, numply, plyt, h, hxz[n],num);

				/*Get current iteration Stresses*/
				/*StrainE * CEmatrix = StressE*/
				for(i=0;i<4;i++)	 /*Calculate it for each node in turn*/
				{
					k= 6*i;		/*k and i are used for indexing purposes*/

					/*UXE Term*/
					OXE[i] = (UXE[0]*KE1[0+k] + UYE[0]*KE1[24+k] + UZE[0]*KE1[48+k] + RXE[0]*KE1[72+k] + RYE[0]*KE1[96+k] + RZE[0]*KE1[120+k] + UXE[1]*KE1[144+k] + UYE[1]*KE1[168+k] + UZE[1]*KE1[192+k] + RXE[1]*KE1[216+k] + RYE[1]*KE1[240+k] + RZE[1]*KE1[264+k] + UXE[2]*KE1[288+k] + UYE[2]*KE1[312+k] + UZE[2]*KE1[336+k] + RXE[2]*KE1[360+k] + RYE[2]*KE1[384+k] + RZE[2]*KE1[408+k] + UXE[3]*KE1[432+k] + UYE[3]*KE1[456+k] + UZE[3]*KE1[480+k] + RXE[3]*KE1[504+k] + RYE[3]*KE1[528+k] + RZE[3]*KE1[552+k]);
					/*UYE Term*/																																																					
					OYE[i] = (UXE[0]*KE1[1+k] + UYE[0]*KE1[25+k] + UZE[0]*KE1[49+k] + RXE[0]*KE1[73+k] + RYE[0]*KE1[97+k] + RZE[0]*KE1[121+k] + UXE[1]*KE1[145+k] + UYE[1]*KE1[169+k] + UZE[1]*KE1[193+k] + RXE[1]*KE1[217+k] + RYE[1]*KE1[241+k] + RZE[1]*KE1[265+k] + UXE[2]*KE1[289+k] + UYE[2]*KE1[313+k] + UZE[2]*KE1[337+k] + RXE[2]*KE1[361+k] + RYE[2]*KE1[385+k] + RZE[2]*KE1[409+k] + UXE[3]*KE1[433+k] + UYE[3]*KE1[457+k] + UZE[3]*KE1[481+k] + RXE[3]*KE1[505+k] + RYE[3]*KE1[529+k] + RZE[3]*KE1[553+k]);
					/*UZE Term*/																																																					
					OZE[i] = (UXE[0]*KE1[2+k] + UYE[0]*KE1[26+k] + UZE[0]*KE1[50+k] + RXE[0]*KE1[74+k] + RYE[0]*KE1[98+k] + RZE[0]*KE1[122+k] + UXE[1]*KE1[146+k] + UYE[1]*KE1[170+k] + UZE[1]*KE1[194+k] + RXE[1]*KE1[218+k] + RYE[1]*KE1[242+k] + RZE[1]*KE1[266+k] + UXE[2]*KE1[290+k] + UYE[2]*KE1[314+k] + UZE[2]*KE1[338+k] + RXE[2]*KE1[362+k] + RYE[2]*KE1[386+k] + RZE[2]*KE1[410+k] + UXE[3]*KE1[434+k] + UYE[3]*KE1[458+k] + UZE[3]*KE1[482+k] + RXE[3]*KE1[506+k] + RYE[3]*KE1[530+k] + RZE[3]*KE1[554+k]);
					/*RXE Term*/																																																					
					ORXE[i] = (UXE[0]*KE1[3+k] + UYE[0]*KE1[27+k] + UZE[0]*KE1[51+k] + RXE[0]*KE1[75+k] + RYE[0]*KE1[99+k] + RZE[0]*KE1[123+k] + UXE[1]*KE1[147+k] + UYE[1]*KE1[171+k] + UZE[1]*KE1[195+k] + RXE[1]*KE1[219+k] + RYE[1]*KE1[243+k] + RZE[1]*KE1[267+k] + UXE[2]*KE1[291+k] + UYE[2]*KE1[315+k] + UZE[2]*KE1[339+k] + RXE[2]*KE1[363+k] + RYE[2]*KE1[387+k] + RZE[2]*KE1[411+k] + UXE[3]*KE1[435+k] + UYE[3]*KE1[459+k] + UZE[3]*KE1[483+k] + RXE[3]*KE1[507+k] + RYE[3]*KE1[531+k] + RZE[3]*KE1[555+k]);
					/*RYE Term*/																																																					
					ORYE[i] = (UXE[0]*KE1[4+k] + UYE[0]*KE1[28+k] + UZE[0]*KE1[52+k] + RXE[0]*KE1[76+k] + RYE[0]*KE1[100+k] + RZE[0]*KE1[124+k] + UXE[1]*KE1[148+k] + UYE[1]*KE1[172+k] + UZE[1]*KE1[196+k] + RXE[1]*KE1[220+k] + RYE[1]*KE1[244+k] + RZE[1]*KE1[268+k] + UXE[2]*KE1[292+k] + UYE[2]*KE1[316+k] + UZE[2]*KE1[340+k] + RXE[2]*KE1[364+k] + RYE[2]*KE1[388+k] + RZE[2]*KE1[412+k] + UXE[3]*KE1[436+k] + UYE[3]*KE1[460+k] + UZE[3]*KE1[484+k] + RXE[3]*KE1[508+k] + RYE[3]*KE1[532+k] + RZE[3]*KE1[556+k]);
					/*ORZE Term*/
					ORZE[i] = 0;
				}

				/*Get imbetween iteration Stresses*/
				/*Strain0 * CEmatrix = Stress1*/
				for(i=0;i<4;i++)	 /*Calculate it for each node in turn*/
				{
					k= 6*i;		/*k and i are used for indexing purposes*/
	
					/*UX0 Term*/
					OX1[i] = (UX0[0]*KE1[0+k] + UY0[0]*KE1[24+k] + UZ0[0]*KE1[48+k] + RX0[0]*KE1[72+k] + RY0[0]*KE1[96+k] + RZ0[0]*KE1[120+k] + UX0[1]*KE1[144+k] + UY0[1]*KE1[168+k] + UZ0[1]*KE1[192+k] + RX0[1]*KE1[216+k] + RY0[1]*KE1[240+k] + RZ0[1]*KE1[264+k] + UX0[2]*KE1[288+k] + UY0[2]*KE1[312+k] + UZ0[2]*KE1[336+k] + RX0[2]*KE1[360+k] + RY0[2]*KE1[384+k] + RZ0[2]*KE1[408+k] + UX0[3]*KE1[432+k] + UY0[3]*KE1[456+k] + UZ0[3]*KE1[480+k] + RX0[3]*KE1[504+k] + RY0[3]*KE1[528+k] + RZ0[3]*KE1[552+k]);
					/*UY0 Term*/																																																					
					OY1[i] = (UX0[0]*KE1[1+k] + UY0[0]*KE1[25+k] + UZ0[0]*KE1[49+k] + RX0[0]*KE1[73+k] + RY0[0]*KE1[97+k] + RZ0[0]*KE1[121+k] + UX0[1]*KE1[145+k] + UY0[1]*KE1[169+k] + UZ0[1]*KE1[193+k] + RX0[1]*KE1[217+k] + RY0[1]*KE1[241+k] + RZ0[1]*KE1[265+k] + UX0[2]*KE1[289+k] + UY0[2]*KE1[313+k] + UZ0[2]*KE1[337+k] + RX0[2]*KE1[361+k] + RY0[2]*KE1[385+k] + RZ0[2]*KE1[409+k] + UX0[3]*KE1[433+k] + UY0[3]*KE1[457+k] + UZ0[3]*KE1[481+k] + RX0[3]*KE1[505+k] + RY0[3]*KE1[529+k] + RZ0[3]*KE1[553+k]);
					/*UZ0 Term*/																																																					
					OZ1[i] = (UX0[0]*KE1[2+k] + UY0[0]*KE1[26+k] + UZ0[0]*KE1[50+k] + RX0[0]*KE1[74+k] + RY0[0]*KE1[98+k] + RZ0[0]*KE1[122+k] + UX0[1]*KE1[146+k] + UY0[1]*KE1[170+k] + UZ0[1]*KE1[194+k] + RX0[1]*KE1[218+k] + RY0[1]*KE1[242+k] + RZ0[1]*KE1[266+k] + UX0[2]*KE1[290+k] + UY0[2]*KE1[314+k] + UZ0[2]*KE1[338+k] + RX0[2]*KE1[362+k] + RY0[2]*KE1[386+k] + RZ0[2]*KE1[410+k] + UX0[3]*KE1[434+k] + UY0[3]*KE1[458+k] + UZ0[3]*KE1[482+k] + RX0[3]*KE1[506+k] + RY0[3]*KE1[530+k] + RZ0[3]*KE1[554+k]);
					/*RX0 Term*/																																																					
					ORX1[i] = (UX0[0]*KE1[3+k] + UY0[0]*KE1[27+k] + UZ0[0]*KE1[51+k] + RX0[0]*KE1[75+k] + RY0[0]*KE1[99+k] + RZ0[0]*KE1[123+k] + UX0[1]*KE1[147+k] + UY0[1]*KE1[171+k] + UZ0[1]*KE1[195+k] + RX0[1]*KE1[219+k] + RY0[1]*KE1[243+k] + RZ0[1]*KE1[267+k] + UX0[2]*KE1[291+k] + UY0[2]*KE1[315+k] + UZ0[2]*KE1[339+k] + RX0[2]*KE1[363+k] + RY0[2]*KE1[387+k] + RZ0[2]*KE1[411+k] + UX0[3]*KE1[435+k] + UY0[3]*KE1[459+k] + UZ0[3]*KE1[483+k] + RX0[3]*KE1[507+k] + RY0[3]*KE1[531+k] + RZ0[3]*KE1[555+k]);
					/*RY0 Term*/																																																					
					ORY1[i] = (UX0[0]*KE1[4+k] + UY0[0]*KE1[28+k] + UZ0[0]*KE1[52+k] + RX0[0]*KE1[76+k] + RY0[0]*KE1[100+k] + RZ0[0]*KE1[124+k] + UX0[1]*KE1[148+k] + UY0[1]*KE1[172+k] + UZ0[1]*KE1[196+k] + RX0[1]*KE1[220+k] + RY0[1]*KE1[244+k] + RZ0[1]*KE1[268+k] + UX0[2]*KE1[292+k] + UY0[2]*KE1[316+k] + UZ0[2]*KE1[340+k] + RX0[2]*KE1[364+k] + RY0[2]*KE1[388+k] + RZ0[2]*KE1[412+k] + UX0[3]*KE1[436+k] + UY0[3]*KE1[460+k] + UZ0[3]*KE1[484+k] + RX0[3]*KE1[508+k] + RY0[3]*KE1[532+k] + RZ0[3]*KE1[556+k]);
					/*ORZ0 Term*/
					ORZ1[i] = 0;
				}

				/*Get the C material Matrix properties for this element's previous itteration value of theta*/
				C = cos(Theta0[num]);
				S = sin(Theta0[num]);
				C2 = C*C;
				S2 = S*S;
				C3 = C*C*C;
				S3 = S*S*S;
				C2S = C*C*S;
				CS2 = C*S*S;
				C2S2 = C*C*S*S;
				CS = C*S;
				C3S = C*C*C*S;
				CS3 = C*S*S*S;
				C4 = C*C*C*C;
				S4 = S*S*S*S;

				/*Create C material Matrix*/
				Q11 = C4*E11 + 2*C2*S2*E12 + S4*E22 + 4*C2*S2*E66;
				Q12 = C2*S2*E11 + S4*E12 + C4*E12 + C2*S2*E22 - 4*C2*S2*E66;
				Q16 = C3*S*E11 + C*S3*E12 - C3*S*E12 - C*S3*E22 - 2*C*S*(C2-S2)*E66;
				Q22 = S4*E11 + 2*C2*S2*E12 + C4*E22 + 4*C2*S2*E66;
				Q26 = C*S3*E11 + C3*S*E12 - C*S3*E12 - C3*S*E22 + 2*C*S*(C2-S2)*E66;
				Q66 = C2S2*E11 - 2*C2S2*E12 + C2S2*E22 + (C4-2*C2S2+S4)*E66;
				Q44 = C2*E44 + S2*E55;
				Q45 = CS*E44 - CS*E55;
				Q55 = S2*E44 + C2*E55;
			
				/*Get current iteration stiffness Matrix*/
				CompositeKEDIFF(Q11, Q22, Q12, Q16, Q26, Q66, Q44, Q55, Q45, KE0, numply, plyt, h, hxz[n],num);

				/*Get Previous iteration Stresses*/
				/*Strain0 * C0matrix = Stress0*/
				for(i=0;i<4;i++)	 /*Calculate it for each node in turn*/
				{
					k= 6*i;		/*k and i are used for indexing purposes*/

					/*UX0 Term*/
					OX0[i] = (UX0[0]*KE0[0+k] + UY0[0]*KE0[24+k] + UZ0[0]*KE0[48+k] + RX0[0]*KE0[72+k] + RY0[0]*KE0[96+k] + RZ0[0]*KE0[120+k] + UX0[1]*KE0[144+k] + UY0[1]*KE0[168+k] + UZ0[1]*KE0[192+k] + RX0[1]*KE0[216+k] + RY0[1]*KE0[240+k] + RZ0[1]*KE0[264+k] + UX0[2]*KE0[288+k] + UY0[2]*KE0[312+k] + UZ0[2]*KE0[336+k] + RX0[2]*KE0[360+k] + RY0[2]*KE0[384+k] + RZ0[2]*KE0[408+k] + UX0[3]*KE0[432+k] + UY0[3]*KE0[456+k] + UZ0[3]*KE0[480+k] + RX0[3]*KE0[504+k] + RY0[3]*KE0[528+k] + RZ0[3]*KE0[552+k]);
					/*UY0 Term*/																																																					
					OY0[i] = (UX0[0]*KE0[1+k] + UY0[0]*KE0[25+k] + UZ0[0]*KE0[49+k] + RX0[0]*KE0[73+k] + RY0[0]*KE0[97+k] + RZ0[0]*KE0[121+k] + UX0[1]*KE0[145+k] + UY0[1]*KE0[169+k] + UZ0[1]*KE0[193+k] + RX0[1]*KE0[217+k] + RY0[1]*KE0[241+k] + RZ0[1]*KE0[265+k] + UX0[2]*KE0[289+k] + UY0[2]*KE0[313+k] + UZ0[2]*KE0[337+k] + RX0[2]*KE0[361+k] + RY0[2]*KE0[385+k] + RZ0[2]*KE0[409+k] + UX0[3]*KE0[433+k] + UY0[3]*KE0[457+k] + UZ0[3]*KE0[481+k] + RX0[3]*KE0[505+k] + RY0[3]*KE0[529+k] + RZ0[3]*KE0[553+k]);
					/*UZ0 Term*/																																																					
					OZ0[i] = (UX0[0]*KE0[2+k] + UY0[0]*KE0[26+k] + UZ0[0]*KE0[50+k] + RX0[0]*KE0[74+k] + RY0[0]*KE0[98+k] + RZ0[0]*KE0[122+k] + UX0[1]*KE0[146+k] + UY0[1]*KE0[170+k] + UZ0[1]*KE0[194+k] + RX0[1]*KE0[218+k] + RY0[1]*KE0[242+k] + RZ0[1]*KE0[266+k] + UX0[2]*KE0[290+k] + UY0[2]*KE0[314+k] + UZ0[2]*KE0[338+k] + RX0[2]*KE0[362+k] + RY0[2]*KE0[386+k] + RZ0[2]*KE0[410+k] + UX0[3]*KE0[434+k] + UY0[3]*KE0[458+k] + UZ0[3]*KE0[482+k] + RX0[3]*KE0[506+k] + RY0[3]*KE0[530+k] + RZ0[3]*KE0[554+k]);
					/*RX0 Term*/																																																					
					ORX0[i] = (UX0[0]*KE0[3+k] + UY0[0]*KE0[27+k] + UZ0[0]*KE0[51+k] + RX0[0]*KE0[75+k] + RY0[0]*KE0[99+k] + RZ0[0]*KE0[123+k] + UX0[1]*KE0[147+k] + UY0[1]*KE0[171+k] + UZ0[1]*KE0[195+k] + RX0[1]*KE0[219+k] + RY0[1]*KE0[243+k] + RZ0[1]*KE0[267+k] + UX0[2]*KE0[291+k] + UY0[2]*KE0[315+k] + UZ0[2]*KE0[339+k] + RX0[2]*KE0[363+k] + RY0[2]*KE0[387+k] + RZ0[2]*KE0[411+k] + UX0[3]*KE0[435+k] + UY0[3]*KE0[459+k] + UZ0[3]*KE0[483+k] + RX0[3]*KE0[507+k] + RY0[3]*KE0[531+k] + RZ0[3]*KE0[555+k]);
					/*RY0 Term*/																																																					
					ORY0[i] = (UX0[0]*KE0[4+k] + UY0[0]*KE0[28+k] + UZ0[0]*KE0[52+k] + RX0[0]*KE0[76+k] + RY0[0]*KE0[100+k] + RZ0[0]*KE0[124+k] + UX0[1]*KE0[148+k] + UY0[1]*KE0[172+k] + UZ0[1]*KE0[196+k] + RX0[1]*KE0[220+k] + RY0[1]*KE0[244+k] + RZ0[1]*KE0[268+k] + UX0[2]*KE0[292+k] + UY0[2]*KE0[316+k] + UZ0[2]*KE0[340+k] + RX0[2]*KE0[364+k] + RY0[2]*KE0[388+k] + RZ0[2]*KE0[412+k] + UX0[3]*KE0[436+k] + UY0[3]*KE0[460+k] + UZ0[3]*KE0[484+k] + RX0[3]*KE0[508+k] + RY0[3]*KE0[532+k] + RZ0[3]*KE0[556+k]);
					/*ORZ0 Term*/
					ORZ0[i] = 0;

				}

				/*Calculate the difference between the stresses (Oe-O1) and (O0-O1)*/
				for(i=0;i<4;i++)
				{
					CX[i] = OXE[i] - OX1[i];
					CY[i] = OYE[i] - OY1[i];
					CZ[i] = OZE[i] - OZ1[i];
					CRX[i] = ORXE[i] - ORX1[i];
					CRY[i] = ORYE[i] - ORY1[i];
					CRZ[i] = ORZE[i] - ORZ1[i];
	
					DX[i] = OX0[i] - OX1[i];
					DY[i] = OY0[i] - OY1[i];
					DZ[i] = OZ0[i] - OZ1[i];
					DRX[i] = ORX0[i] - ORX1[i];
					DRY[i] = ORY0[i] - ORY1[i];
					DRZ[i] = ORZ0[i] - ORZ1[i];
				}

				/*Calculate the normal moduli of the vector differences stress differences |Oe-O1| and |O0-O1|*/
				Dtemp = 0;
				Ctemp = 0;
				for(i=0;i<4;i++)
				{
					Ctemp += (CX[i]*CX[i]) + (CY[i]*CY[i]) + (CZ[i]*CZ[i]) + (CRX[i]*CRX[i]) + (CRY[i]*CRY[i]) + (CRZ[i]*CRZ[i]);
					Dtemp += (DX[i]*DX[i]) + (DY[i]*DY[i]) + (DZ[i]*DZ[i]) + (DRX[i]*DRX[i]) + (DRY[i]*DRY[i]) + (DRZ[i]*DRZ[i]);
				}
				Ctemp2 = sqrt(Ctemp);
				Dtemp2 = sqrt(Dtemp);
			
				EF[num] = Ctemp2/Dtemp2;

				if(EF[num]>0.9999999999999){ EF[num] = 0.9999;}
				else if (EF[num]<0.000000001){ EF[num] = 0;}
			}
		}
	}

	free(Lnodes);
	free(UXE);
	free(UYE);
	free(UZE);
	free(RXE);
	free(RYE);
	free(RZE);
	free(OXE);
	free(OYE);
	free(OZE);
	free(ORXE);
	free(ORYE);
	free(ORZE);
	free(OX1);
	free(OY1);
	free(OZ1);
	free(ORX1);
	free(ORY1);
	free(ORZ1);
	free(UX0);
	free(UY0);
	free(UZ0);
	free(RX0);
	free(RY0);
	free(RZ0);
	free(OX0);
	free(OY0);
	free(OZ0);
	free(ORX0);
	free(ORY0);
	free(ORZ0);
	free(CX);
	free(CY);
	free(CZ);
	free(CRX);
	free(CRY);
	free(CRZ);
	free(DX);
	free(DY);
	free(DZ);
	free(DRX);
	free(DRY);
	free(DRZ);
	free(KE0);
	free(KE1);

}
/*we didn't use this part........jiaming*/
void GlobalLineSens(int elemX, int elemY, Elem Number[elemX][elemY], short *ElemStat, double *ElemThetaStn, double *ElemThetaStnC, int itt)
{

	FILE *outfile;	/*File varible for output files*/
	char plotname[50];	/*variable to change names of plotting output files*/
	int i, j, k, num, num2, countFS, countCS;
	int count = 0;
	int Ncut, Nfull;
	int NumElem = elemX*elemY;
	SET *CutSet;
	SET *FullSet;
	double MinD, dtemp, etemp, xtemp, ytemp;
	FullSet = malloc(NumElem*sizeof(SET));
	CutSet = malloc(NumElem*sizeof(SET));
	double X, Y, X2, Y2;

	k = (elemY<elemX) ? elemY:elemX;
	int  *Neighbours;
	double dist;
	Neighbours = malloc(k*sizeof(int));

	for(j=0;j<k;j++)
	{
		Neighbours[j] = -1;
	}

	countFS = 0;
	countCS = 0;
	/*First we need to establish which elements are cut and which aren't and put them in a set*/
	for(i=0;i<elemX;i++)
	{
		for(j=0;j<elemY;j++)
		{
			num = Number[i][j].n-1;
			ElemThetaStnC[num] = ElemStat[num] * ElemThetaStn[num]; /*Whilst we are here intailise ElemThetaStnC to ElemThetaStn value if uncut, and 0  if cut*/
			/*printf("\nETS[%i] = %f", num, ElemThetaStn[num]);*/
			if(ElemStat[num]==0)
			{
				FullSet[countFS].r = num;
				FullSet[countFS].x = i;
				FullSet[countFS].y = j;
				countFS++;
			}
			else if(ElemStat[num]==1)
			{
				CutSet[countCS].r = num;
				CutSet[countCS].x = i;
				CutSet[countCS].y = j;
				countCS++;
			}
		}
	}
	/*printf("\ncountCS = %i", countCS);
	printf("\ncountFS = %i", countFS);*/

	Nfull = countFS;
	Ncut = countCS;

	/*sprintf(plotname,"ElemThetaStn%iC.txt",itt); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Element Strain Energy plotting writefile\n");
		}
	else{
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;
				fprintf(outfile,"%f\t",ElemThetaStnC[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);*/

	count = 1.0;

	for(i=0;i<Nfull;i++)	/*For each of the uncut elements*/
	{
		MinD = 100000000;	/*Intialise MinD to a large number*/
		count = 1.0;

		X = FullSet[i].x;		/*Get the key values for this uncut element*/
		Y = FullSet[i].y;
		num = FullSet[i].r;

		for(j=0;j<Ncut;j++)	/*Look at each of the cut elements and see which one is largest*/
		{
			X2 = CutSet[j].x;		/*Get the key values for this cut element*/
			Y2 = CutSet[j].y;
			num2 = CutSet[j].r;

			xtemp = X2-X;							/*Calculate the distance between the two elements*/
			ytemp = Y2-Y;
			etemp = xtemp*xtemp + ytemp*ytemp;
			dtemp = sqrt(etemp);

			if(dtemp<MinD)					/*If we have found a new minimum then get rid of all old neighbours and set up a new set*/
			{
				for(k=0;k<count;k++)
				{
					Neighbours[k] = -1;
				}
				Neighbours[0] = num2;
				count = 1.0;
				MinD = dtemp;
				dist = MinD;
			}
			else if(dtemp==MinD)		/*If this element is equal distance away as the current neighbour then add it to the set*/
			{
				Neighbours[count] = num2;
				count++;
			}
		/*printf("\ncount = %f", count);
		printf("\nMinD = %f\tdtemp = %f", MinD, dtemp);*/

		}
		/*printf("\nnum = %i", num);
		printf("\ncount = %f", count);*/
		/*Now at the sensitviy of the uncut element to the sensitvtiy of it's nearest cut Neighbour(s). May apply a weighting system here later.*/
		for(k=0;k<count;k++)
		{
			dtemp = count;
			/*printf("\n1/dtemp = 1/%f = %f", dtemp, 1/dtemp);*/
			num2 = Neighbours[k];
			etemp = ElemThetaStn[num];
			/*printf("\netemp = ETS[%i] = %f", num, ElemThetaStn[num]);*/
			ElemThetaStnC[num2] += etemp/(dtemp*dist*dist);		/*Divide contribution of uncut element's to the cut element's sesntivity by the number of elements it contributing to*/
			/*printf("\nETSC[%i] = %f / %f * %f = %f = %f", num2, etemp, dtemp, dist,etemp/dtemp ,ElemThetaStnC[num2]);/**/
		}
	}
	/*printf("\ncount = %i", count);*/


	/*We should have the new Global "ElemThetaStnC" sensitvities so print out and the free memory values*
	
	/**/
	free(FullSet);
	free(CutSet);
	free(Neighbours);

}
/*Calculate the change in theta from the sensitvities*/
void ScaleDeltaTheta(int elemX, int elemY, double *ElemThetaStn , double *DeltaTheta, double DeltaThetaMax, int itt, double CompA[itt],  Elem Number[elemX][elemY], short *ElemStat, short *ElemStatb, int Numlsf, double *NodeThetaStn, double *NodeDTL, double *DeltaT, int Ntot_new , double *Vnorm,double *SensMaxOld)
{
	FILE *outfile;	/*File varible for output files*/
	char plotname[40];	/*variable to change names of plotting output files*/
	int i,j,num,p,k;
	double Ctemp, Ctemp2, MaxC, MinC, SensMax, SensTemp2, SensTemp, SensSum, SensMaxTemp,SensMaxTemp2;

	//SensMaxTemp = *SensMaxOld; /*Set SensMax to 0*/
    SensMax=0;
	SensSum = 0; /*Set SensMax to 0*/
	/*First find the maximum sensitvity values*/
	for(i=0; i<(elemX+1)*(elemY+1);i++)
	{
		for(p=0;p<Numlsf;p++)
		{
			SensMax = (fabs(Vnorm[i])>fabs(SensMax)) ? fabs(Vnorm[i]):fabs(SensMax);
			/*k = i + (elemX*elemY)*p; 
			SensMax = (fabs(ElemThetaStn[i]*ElemStat[k])>fabs(SensMax)) ? fabs(ElemThetaStn[i]):fabs(SensMax);*/
		}
	}
    for(i=0; i<(elemX*elemY);i++)
    {
        for(p=0;p<Numlsf;p++)
        {
           // SensMax = (fabs(ElemThetaStn[i]*ElemStatb[i])>fabs(SensMax)) ? fabs(ElemThetaStn[i]):fabs(SensMax);
        }
    }
    
	/*SensMax = SensSum/(elemX*elemY);*/
    //*SensMaxOld=SensMax;
    //SensMaxTemp2=SensMax;
    //if (SensMax-SensMaxTemp<0.0000000001)
    //{
   //     SensMax=SensMaxTemp;
   // }
    // set the sensmax is equal to the sensmax in his own itteration instead the old one which is =0.
    

    
    //SensMax = (fabs(SensMaxTemp)>fabs(SensMax)) ? fabs(SensMaxTemp):fabs(SensMax);
	printf("\n SensMax = %f \n", SensMax);
    
	/*Now if more that 10 itterations have past find the minmum and Maximum compliance from the last 10 itterations*/
	if(itt<10)
	{
		Ctemp = 1.0;
	}
	else
	{
		MaxC = 0;				/*Set to small value*/
		MinC = 100000000000;	/*Set to very large value*/
		/*Get max and min Compliance*/
		for(i=0; i<(10);i++)
		{
			j = itt-i;
			MaxC = (CompA[j]>MaxC) ? CompA[j]:MaxC;
			MinC = (CompA[j]<MinC) ? CompA[j]:MinC;
		}
		/*printf("\n MaxC = %f", MaxC);
		printf("\n MinC = %f", MinC);*/
		
		Ctemp2 = (MaxC - MinC)/MaxC;	/*Now calculate Ctemp from these values*/
		Ctemp = sqrt(sqrt(sqrt(sqrt(Ctemp2))));			/*Sqaure route relaxes the positive feed back of scaling this off compliance*/
		Ctemp = 1.0;
	}

	/*Now scale delta theta for each element to the compliance and maximum sensitvtiy*/

    //DeltaThetaMax=0.1;
    double Maxthetachange;
    //Maxthetachange=3.0*(3.14159265/180.0);
    
    //*DeltaT = Ctemp*Maxthetachange/fabs(SensMax);
    //*DeltaT = Ctemp*DeltaThetaMax/(1+fabs(SensMax));
    *DeltaT = Ctemp*DeltaThetaMax/(fabs(SensMax));		/*Get desired change in theta.*/
    
    
    
		/*printf("\tDeltaTheta = %f", DeltaTheta[i]);*/


	/*Print out the result*
	sprintf(plotname,"DeltaTheta%i.txt",itt);
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Element Strain Energy plotting writefile\n");
	}
	else{

		for(j=elemY-1;j>=0;j--)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;	
				fprintf(outfile,"%f\t",DeltaTheta[num]); 
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);

	/*For ease of calculation and viewing that was calculated in degrees, but computer porgrams always use radiums so transform it*/

}

/*Use the nodal sensitvity of lsf value to fibre angle in each element and the element sensitvity to get the change in lsf at each node*/
void DeltaLsfCalc(int NumNodes, double *Deltalsf2, double *DeltaTheta, NElemDiff *DTL, short *ElemStat, double *Vnorm)
{
	int i, count,num_kk, index;
	double dtemp;
    int kk;
	for(i=0;i<NumNodes;i++)
	{
		count = 0;
		dtemp = 0;
        kk=0;
		index = DTL[i].ai;			/*Is their an element next to the node in this position*/
		if(index == 1)
		{
			num_kk = DTL[i].an;
					/*Update total lsf change*/
            kk+=ElemStat[num_kk];
            //Deltalsf2[i] = Vnorm[i]*ElemStat[num_kk];
            
            count += 1;
		}

		index = DTL[i].bi;			/*Is their an element next to the node in this position*/
		if(index == 1)
		{
			num_kk = DTL[i].bn;
			/*Update total lsf change*/
			count += 1;
            kk+=ElemStat[num_kk];
            
		}

		index = DTL[i].ci;			/*Is their an element next to the node in this position*/
		if(index == 1)
		{
			num_kk = DTL[i].cn;
			
			count += 1;
            kk+=ElemStat[num_kk];
		}

		index = DTL[i].di;			/*Is their an element next to the node in this position*/
		if(index == 1)
		{
			num_kk = DTL[i].dn;
			
			count += 1;
            kk+=ElemStat[num_kk];
		}
        /*
            if (Vnorm[i]>=0.00000001)
            {
                Deltalsf2[i] = log(Vnorm[i]+1.0);
            }
            else if (Vnorm[i]<=-0.00000001) {
                Vnorm[i]=-Vnorm[i]+1.0;
                Deltalsf2[i] = -log(Vnorm[i]);
            }
            else {
            
                Deltalsf2[i] = Vnorm[i];

            }*/
        //if (kk!=0) {
             Deltalsf2[i] = Vnorm[i];
       // }
       // else{
        //    Deltalsf2[i] = 0.0;
      //  }
        /*/count*//*Change in lsf is the mean of all the elment lsf values*/
        
        if(kk!=0)
        {
            if (i==83||i==82||i==84||i==85) {
                 //printf("\nVorm: %d, %.16lf", i, Vnorm[i]);
            }
       
        }
	}
}

void NodeSen(int NumNodes,  NElemDiff *DTL, short *ElemStat,  double *NodesenTemp, double *ElemThetaStn)
{
    int i, count,num, index;
    double dtemp;

    for(i=0;i<NumNodes;i++)
    {
        count = 0;
        dtemp = 0;
        index = DTL[i].ai;			/*Is their an element next to the node in this position*/
        if(index == 1)
        {
            num = DTL[i].an;
            //dtemp += DTL[i].ad*ElemThetaStn[num]*ElemStat[num];		/*Update total lsf change*/
            dtemp += DTL[i].ad*ElemThetaStn[num];
            printf("\nDTL at node %d is : %.16lf, %d", i, DTL[i].ad, num);
            count += 1;
        }
        
        index = DTL[i].bi;			/*Is their an element next to the node in this position*/
        if(index == 1)
        {
            num = DTL[i].bn;
            dtemp += DTL[i].bd*ElemThetaStn[num];		/*Update total lsf change*/
            count += 1;
            printf("\nDTL at node %d is : %.16lf, %d", i, DTL[i].bd, num);
        }
        
        index = DTL[i].ci;			/*Is their an element next to the node in this position*/
        if(index == 1)
        {
            num = DTL[i].cn;
            dtemp += DTL[i].cd*ElemThetaStn[num];		/*Update total lsf change*/
            count += 1;
            printf("\nDTL at node %d is : %.16lf, %d", i, DTL[i].cd, num);
        }
        
        index = DTL[i].di;			/*Is their an element next to the node in this position*/
        if(index == 1)
        {
            num = DTL[i].dn;
            dtemp += DTL[i].dd*ElemThetaStn[num] ;		/*Update total lsf change*/
            count += 1;
            printf("\nDTL at node %d is : %.16lf, %d", i, DTL[i].dd, num);
        }
        
        NodesenTemp[i] = dtemp;
        if (i==82||i==80||i==133||i==224||i==362||i==400) {
            printf("\n NodesenTemp %d: %.16lf", i, NodesenTemp[i]);
        }
        
    }
    
    for (i=0; i<NumNodes; i++) {
        printf("\nnodesens: %d, %f", i, NodesenTemp[i]);
    }
}

/*Gets the guass stiffness matirx for sensitvity evaluation*/
void CompositeKEGaussDIFF(int GP, double Q11, double Q22, double Q12, double Q16, double Q26, double Q66, double Q44, double Q55, double Q45, double *KE, int numply, double plyt, double h, double hxz, int num)
{
/*NOTE For the sake of the dgemm functions set up here is column strong (arrays go down thecolumns first, then along the rows). The fina out put however will be the other way around so that we can use it in assmeble*/
int n,m,i,j,p,q;
double dn1dx,dn1dy,dn2dx,dn2dy,dn3dx,dn3dy,dn4dx,dn4dy,nx,ny,n1j,n2j,n3j,n4j;
/*Intialize all the matrices we'll be using*/
double *A, *D, *B, *Cmatrix, *F, *Qmatrix;
/*Material property matrices*/
A = malloc(9*sizeof(double));
D = malloc(9*sizeof(double));
B = malloc(9*sizeof(double));
F = malloc(4*sizeof(double));
Cmatrix = malloc(4*sizeof(double));
Qmatrix = malloc(9*sizeof(double));
/*Shape function Matrices*/
double *BB,*BM,*BS;
BM = malloc(24*sizeof(double));
BB = malloc(36*sizeof(double));
BS = malloc(24*sizeof(double));

/*StiffnessMatrixes*/
double *Km, *Kp;
Km = malloc(64*sizeof(double));
Kp = malloc(144*sizeof(double));

	for(n=0;n<3;n++)
	{
		for(m=0;m<3;m++)
		{
			A[n+3*m] = 0;
			D[n+3*m] = 0;
			B[n+3*m] = 0;
		}
	}
	for(p=0;p<8;p++)
	{
		for(q=0;q<8;q++)
		{
			Km[p+8*q] = 0;
		}
	}
	for(p=0;p<12;p++)
	{
		for(q=0;q<12;q++)
		{
			Kp[p+12*q] = 0;
		}
	}				
		
	for(q=0;q<numply;q++)
	{

		Qmatrix[0] = Q11;
		Qmatrix[1] = Q12;
		Qmatrix[2] = Q16;
		Qmatrix[3] = Q12;
		Qmatrix[4] = Q22;
		Qmatrix[5] = Q26;
		Qmatrix[6] = Q16;
		Qmatrix[7] = Q26;
		Qmatrix[8] = Q66;

		Cmatrix[0] = Q44;
		Cmatrix[1] = Q45;
		Cmatrix[2] = Q45;
		Cmatrix[3] = Q55;

		for(n=0;n<3;n++)
		{
			for(m=0;m<3;m++)
			{
				A[n+3*m] += plyt*Qmatrix[n+3*m];
				B[n+3*m] += -1*0.5*plyt*plyt*Qmatrix[n+3*m];
				D[n+3*m] += (plyt*plyt*plyt*Qmatrix[n+3*m])/12.0;
				if((n<2)&&(m<2)){F[n+2*m] = (5.0/6.0)*plyt*Cmatrix[n+2*m];}
			}
		}				
	}
	/*if((num == 0))

	{
	printf("\nCmatrix\n");
	printf("%f\t %f\n", Cmatrix[0], Cmatrix[2]);
	printf("%f\t %f\n", Cmatrix[1], Cmatrix[3]);


	printf("\nQmatrix\n");
	printf("%f\t %f\t %f\n", Qmatrix[0], Qmatrix[3], Qmatrix[6]);
	printf("%f\t %f\t %f\n", Qmatrix[1], Qmatrix[4], Qmatrix[7]);
	printf("%f\t %f\t %f\n", Qmatrix[2], Qmatrix[5], Qmatrix[8]);


	/*printf("\nAmatrix\n");
	printf("%f\t %f\t %f\n", A[0], A[3], A[6]);
	printf("%f\t %f\t %f\n", A[1], A[4], A[7]);
	printf("%f\t %f\t %f\n", A[2], A[5], A[8]);

	printf("\nBmatrix\n");
	printf("%f\t %f\t %f\n", B[0], B[3], B[6]);
	printf("%f\t %f\t %f\n", B[1], B[4], B[7]);
	printf("%f\t %f\t %f\n", B[2], B[5], B[8]);


	printf("\nDmatrix\n");
	printf("%f\t %f\t %f\n", D[0], D[3], D[6]);
	printf("%f\t %f\t %f\n", D[1], D[4], D[7]);
	printf("%f\t %f\t %f\n", D[2], D[5], D[8]);

	printf("\nFmatrix\n");
	printf("%f\t %f\n", F[0], F[2]);
	printf("%f\t %f\n", F[1], F[3]);
	}*/

	free(Cmatrix);
	free(Qmatrix);

        if((GP==0)||(GP==3)){nx = -0.5773502692;}
		else{nx = 0.5773502692;}

		if((GP==0)||(GP==1)){ny = -0.5773502692;}
		else{ny = 0.5773502692;}

		/*printf("\n GP = %i, nx = %f, ny = %f", GP, nx, ny);*/


		n1j = 0.125*(1-nx)*(1-ny);
		n2j = 0.125*(1+nx)*(1-ny);
		n3j = 0.125*(1+nx)*(1+ny);
		n4j = 0.125*(1-nx)*(1+ny);

		dn1dx = -0.25*(1-ny)/hxz;
		dn2dx = 0.25*(1-ny)/hxz;
		dn3dx = 0.25*(1+ny)/hxz;
		dn4dx = -0.25*(1+ny)/hxz;

		dn1dy = -0.25*(1-nx)/h;
		dn2dy = -0.25*(1+nx)/h;
		dn3dy = 0.25*(1+nx)/h;
		dn4dy = 0.25*(1-nx)/h;

		/*printf("\n dn1dx = %f", dn1dx);
		printf("\n dn1dx = %f", dn2dx);
		printf("\n dn1dx = %f", dn3dx);
		printf("\n dn1dx = %f", dn4dx);

		printf("\n dn1dx = %f", dn1dy);
		printf("\n dn1dx = %f", dn2dy);
		printf("\n dn1dx = %f", dn3dy);

		printf("\n dn1dx = %f", dn4dy);


		printf("\n n1j = %f", n1j);
		printf("\n n2j = %f", n2j);
		printf("\n n3j = %f", n3j);

		printf("\n n4j = %f\n", n4j);*/

		BM[0] = dn1dx;
		BM[3] = 0;
		BM[6] = dn2dx;
		BM[9] = 0;
		BM[12] = dn3dx;
		BM[15] = 0;
		BM[18] = dn4dx;
		BM[21] = 0;
		
		BM[1] = 0;
		BM[4] = dn1dy;
		BM[7] = 0;
		BM[10] = dn2dy;
		BM[13] = 0;
		BM[16] = dn3dy;
		BM[19] = 0;
		BM[22] = dn4dy;

		BM[2] = dn1dy;
		BM[5] = dn1dx;
		BM[8] = dn2dy;
		BM[11] = dn2dx;
		BM[14] = dn3dy;
		BM[17] = dn3dx;
		BM[20] = dn4dy;
		BM[23] = dn4dx;

		BB[0] = 0;
		BB[3] = 0;
		BB[6] = -dn1dx;
		BB[9] = 0;
		BB[12] = 0;
		BB[15] = -dn2dx;
		BB[18] = 0;
		BB[21] = 0;
		BB[24] = -dn3dx;
		BB[27] = 0;
		BB[30] = 0;
		BB[33] = -dn4dx;

		BB[1] = 0;
		BB[4] = dn1dy;
		BB[7] = 0;
		BB[10] = 0;
		BB[13] = dn2dy;
		BB[16] = 0;
		BB[19] = 0;
		BB[22] = dn3dy;
		BB[25] = 0;
		BB[28] = 0;
		BB[31] = dn4dy;
		BB[34] = 0;

		BB[2] = 0;
		BB[5] = dn1dx;
		BB[8] = -dn1dy;
		BB[11] = 0;
		BB[14] = dn2dx;
		BB[17] = -dn2dy;
		BB[20] = 0;
		BB[23] = dn3dx;
		BB[26] = -dn3dy;
		BB[29] = 0;
		BB[32] = dn4dx;
		BB[35] = -dn4dy;

		BS[0] = dn1dx;
		BS[2] = 0;
		BS[4] = n1j;
		BS[6] = dn2dx;
		BS[8] = 0;
		BS[10] = n2j;
		BS[12] = dn3dx;
		BS[14] = 0;
		BS[16] = n3j;
		BS[18] = dn4dx;
		BS[20] = 0;
		BS[22] = n4j;

		BS[1] = dn1dy;
		BS[3] = -n1j;
		BS[5] = 0;
		BS[7] = dn2dy;
		BS[9] = -n2j;
		BS[11] = 0;
		BS[13] = dn3dy;
		BS[15] = -n3j;
		BS[17] = 0;
		BS[19] = dn4dy;
		BS[21] = -n4j;
		BS[23] = 0;

/*if((num == 361)||(num==378))

{
printf("BMshape\n");

printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[0], BM[3], BM[6], BM[9], BM[12], BM[15], BM[18], BM[21]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[1], BM[4], BM[7], BM[10], BM[13], BM[16], BM[19], BM[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[2], BM[5], BM[8], BM[11], BM[14], BM[17], BM[20], BM[23]);
}

/*printf("BBshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[0], BB[3], BB[6], BB[9], BB[12], BB[15], BB[18], BB[21], BB[24], BB[27], BB[30], BB[33]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[1], BB[4], BB[7], BB[10], BB[13], BB[16], BB[19], BB[22], BB[25], BB[28], BB[31], BB[34]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[2], BB[5], BB[8], BB[11], BB[14], BB[17], BB[20], BB[23], BB[26], BB[29], BB[32], BB[35]);

printf("BSshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[0], BS[2], BS[4], BS[6], BS[8], BS[10], BS[12], BS[14], BS[16], BS[18], BS[20], BS[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[1], BS[3], BS[5], BS[7], BS[9], BS[11], BS[13], BS[15], BS[17], BS[19], BS[21], BS[23]);*/



		/*Multiply the shape functions and the material matrices together using dgemm*/
/*printf("\nA");

printf("\nA");*/

/*Membrain terms*/		
Km[0] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[0] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[1] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[2]);
Km[1] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[0] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[1] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[2]);
Km[2] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[0] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[1] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[2]);
Km[3] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[0] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[1] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[2]);
Km[4] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[0] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[1] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[2]);
Km[5] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[0] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[1] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[2]);
Km[6] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[0] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[1] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[2]);
Km[7] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[0] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[1] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[2]);

Km[8] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[3] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[4] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[5]);
Km[9] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[3] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[4] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[5]);
Km[10] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[3] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[4] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[5]);
Km[11] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[3] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[4] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[5]);
Km[12] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[3] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[4] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[5]);
Km[13] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[3] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[4] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[5]);
Km[14] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[3] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[4] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[5]);
Km[15] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[3] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[4] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[5]);

Km[16] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[6] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[7] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[8]);
Km[17] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[6] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[7] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[8]);
Km[18] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[6] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[7] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[8]);
Km[19] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[6] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[7] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[8]);
Km[20] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[6] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[7] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[8]);
Km[21] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[6] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[7] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[8]);
Km[22] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[6] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[7] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[8]);
Km[23] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[6] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[7] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[8]);

Km[24] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[9] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[10] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[11]);
Km[25] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[9] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[10] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[11]);
Km[26] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[9] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[10] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[11]);
Km[27] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[9] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[10] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[11]);
Km[28] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[9] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[10] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[11]);
Km[29] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[9] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[10] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[11]);
Km[30] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[9] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[10] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[11]);
Km[31] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[9] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[10] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[11]);
		
Km[32] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[12] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[13] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[14]);
Km[33] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[12] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[13] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[14]);
Km[34] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[12] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[13] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[14]);
Km[35] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[12] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[13] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[14]);
Km[36] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[12] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[13] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[14]);
Km[37] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[12] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[13] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[14]);
Km[38] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[12] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[13] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[14]);
Km[39] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[12] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[13] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[14]);

Km[40] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[15] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[16] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[17]);
Km[41] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[15] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[16] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[17]);
Km[42] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[15] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[16] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[17]);
Km[43] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[15] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[16] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[17]);
Km[44] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[15] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[16] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[17]);
Km[45] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[15] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[16] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[17]);
Km[46] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[15] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[16] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[17]);
Km[47] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[15] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[16] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[17]);

Km[48] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[18] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[19] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[20]);
Km[49] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[18] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[19] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[20]);
Km[50] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[18] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[19] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[20]);
Km[51] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[18] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[19] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[20]);
Km[52] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[18] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[19] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[20]);
Km[53] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[18] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[19] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[20]);
Km[54] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[18] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[19] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[20]);
Km[55] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[18] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[19] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[20]);

Km[56] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[21] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[22] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[23]);
Km[57] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[21] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[22] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[23]);
Km[58] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[21] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[22] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[23]);
Km[59] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[21] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[22] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[23]);
Km[60] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[21] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[22] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[23]);
Km[61] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[21] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[22] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[23]);
Km[62] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[21] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[22] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[23]);
Km[63] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[21] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[22] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[23]);


/*Bending terms*/
Kp[0] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[0] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[1] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[2]);
Kp[1] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[0] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[1] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[2]);
Kp[2] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[0] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[1] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[2]);
Kp[3] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[0] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[1] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[2]);
Kp[4] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[0] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[1] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[2]);
Kp[5] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[0] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[1] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[2]);
Kp[6] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[0] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[1] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[2]);
Kp[7] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[0] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[1] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[2]);
Kp[8] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[0] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[1] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[2]);
Kp[9] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[0] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[1] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[2]);
Kp[10] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[0] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[1] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[2]);
Kp[11] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[0] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[1] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[2]);

Kp[12] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[3] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[4] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[5]);
Kp[13] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[3] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[4] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[5]);
Kp[14] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[3] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[4] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[5]);
Kp[15] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[3] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[4] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[5]);
Kp[16] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[3] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[4] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[5]);
Kp[17] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[3] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[4] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[5]);
Kp[18] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[3] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[4] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[5]);
Kp[19] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[3] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[4] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[5]);
Kp[20] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[3] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[4] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[5]);
Kp[21] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[3] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[4] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[5]);
Kp[22] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[3] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[4] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[5]);
Kp[23] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[3] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[4] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[5]);

Kp[24] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[6] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[7] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[8]);
Kp[25] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[6] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[7] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[8]);
Kp[26] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[6] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[7] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[8]);
Kp[27] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[6] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[7] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[8]);
Kp[28] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[6] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[7] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[8]);
Kp[29] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[6] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[7] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[8]);
Kp[30] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[6] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[7] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[8]);
Kp[31] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[6] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[7] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[8]);
Kp[32] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[6] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[7] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[8]);
Kp[33] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[6] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[7] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[8]);
Kp[34] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[6] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[7] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[8]);
Kp[35] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[6] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[7] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[8]);

Kp[36] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[9] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[10] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[11]);
Kp[37] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[9] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[10] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[11]);
Kp[38] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[9] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[10] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[11]);
Kp[39] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[9] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[10] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[11]);
Kp[40] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[9] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[10] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[11]);
Kp[41] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[9] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[10] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[11]);
Kp[42] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[9] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[10] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[11]);
Kp[43] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[9] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[10] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[11]);
Kp[44] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[9] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[10] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[11]);
Kp[45] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[9] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[10] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[11]);
Kp[46] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[9] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[10] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[11]);
Kp[47] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[9] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[10] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[11]);

Kp[48] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[12] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[13] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[14]);
Kp[49] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[12] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[13] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[14]);
Kp[50] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[12] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[13] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[14]);
Kp[51] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[12] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[13] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[14]);
Kp[52] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[12] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[13] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[14]);
Kp[53] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[12] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[13] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[14]);
Kp[54] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[12] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[13] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[14]);
Kp[55] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[12] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[13] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[14]);
Kp[56] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[12] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[13] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[14]);
Kp[57] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[12] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[13] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[14]);
Kp[58] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[12] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[13] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[14]);
Kp[59] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[12] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[13] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[14]);

Kp[60] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[15] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[16] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[17]);
Kp[61] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[15] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[16] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[17]);
Kp[62] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[15] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[16] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[17]);
Kp[63] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[15] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[16] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[17]);
Kp[64] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[15] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[16] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[17]);
Kp[65] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[15] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[16] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[17]);
Kp[66] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[15] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[16] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[17]);
Kp[67] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[15] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[16] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[17]);
Kp[68] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[15] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[16] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[17]);
Kp[69] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[15] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[16] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[17]);
Kp[70] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[15] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[16] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[17]);
Kp[71] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[15] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[16] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[17]);

Kp[72] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[18] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[19] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[20]);
Kp[73] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[18] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[19] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[20]);
Kp[74] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[18] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[19] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[20]);
Kp[75] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[18] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[19] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[20]);
Kp[76] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[18] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[19] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[20]);
Kp[77] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[18] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[19] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[20]);
Kp[78] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[18] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[19] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[20]);
Kp[79] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[18] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[19] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[20]);
Kp[80] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[18] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[19] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[20]);
Kp[81] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[18] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[19] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[20]);
Kp[82] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[18] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[19] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[20]);
Kp[83] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[18] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[19] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[20]);

Kp[84] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[21] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[22] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[23]);
Kp[85] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[21] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[22] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[23]);
Kp[86] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[21] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[22] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[23]);
Kp[87] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[21] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[22] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[23]);
Kp[88] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[21] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[22] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[23]);
Kp[89] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[21] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[22] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[23]);
Kp[90] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[21] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[22] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[23]);
Kp[91] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[21] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[22] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[23]);
Kp[92] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[21] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[22] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[23]);
Kp[93] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[21] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[22] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[23]);
Kp[94] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[21] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[22] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[23]);
Kp[95] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[21] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[22] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[23]);

Kp[96] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[24] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[25] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[26]);
Kp[97] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[24] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[25] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[26]);
Kp[98] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[24] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[25] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[26]);
Kp[99] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[24] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[25] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[26]);
Kp[100] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[24] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[25] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[26]);
Kp[101] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[24] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[25] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[26]);
Kp[102] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[24] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[25] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[26]);
Kp[103] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[24] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[25] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[26]);
Kp[104] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[24] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[25] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[26]);
Kp[105] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[24] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[25] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[26]);
Kp[106] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[24] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[25] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[26]);
Kp[107] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[24] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[25] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[26]);

Kp[108] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[27] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[28] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[29]);
Kp[109] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[27] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[28] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[29]);
Kp[110] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[27] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[28] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[29]);
Kp[111] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[27] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[28] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[29]);
Kp[112] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[27] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[28] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[29]);
Kp[113] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[27] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[28] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[29]);
Kp[114] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[27] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[28] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[29]);
Kp[115] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[27] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[28] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[29]);
Kp[116] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[27] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[28] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[29]);
Kp[117] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[27] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[28] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[29]);
Kp[118] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[27] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[28] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[29]);
Kp[119] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[27] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[28] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[29]);

Kp[120] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[30] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[31] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[32]);
Kp[121] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[30] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[31] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[32]);
Kp[122] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[30] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[31] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[32]);
Kp[123] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[30] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[31] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[32]);
Kp[124] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[30] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[31] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[32]);
Kp[125] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[30] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[31] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[32]);
Kp[126] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[30] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[31] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[32]);
Kp[127] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[30] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[31] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[32]);
Kp[128] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[30] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[31] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[32]);
Kp[129] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[30] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[31] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[32]);
Kp[130] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[30] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[31] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[32]);
Kp[131] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[30] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[31] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[32]);

Kp[132] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[33] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[34] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[35]);
Kp[133] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[33] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[34] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[35]);
Kp[134] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[33] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[34] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[35]);
Kp[135] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[33] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[34] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[35]);
Kp[136] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[33] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[34] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[35]);
Kp[137] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[33] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[34] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[35]);
Kp[138] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[33] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[34] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[35]);
Kp[139] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[33] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[34] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[35]);
Kp[140] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[33] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[34] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[35]);
Kp[141] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[33] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[34] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[35]);
Kp[142] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[33] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[34] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[35]);
Kp[143] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[33] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[34] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[35]);

    
    
free(BB);
free(BM);
free(A);
free(D);
free(B);

/*Shear terms go through single guass intergration only*/

    
    nx = 0;
ny = 0;

n1j = 0.25*(1-nx)*(1-ny);
n2j = 0.25*(1+nx)*(1-ny);
n3j = 0.25*(1+nx)*(1+ny);
n4j = 0.25*(1-nx)*(1+ny);

dn1dx = -0.5*(1-ny)/hxz;
dn2dx = 0.5*(1-ny)/hxz;
dn3dx = 0.5*(1+ny)/hxz;
dn4dx = -0.5*(1+ny)/hxz;

dn1dy = -0.5*(1-nx)/h;
dn2dy = -0.5*(1+nx)/h;
dn3dy = 0.5*(1+nx)/h;
dn4dy = 0.5*(1-nx)/h;


BS[0] = dn1dx;
BS[2] = 0;
BS[4] = n1j;
BS[6] = dn2dx;
BS[8] = 0;
BS[10] = n2j;
BS[12] = dn3dx;
BS[14] = 0;
BS[16] = n3j;
BS[18] = dn4dx;
BS[20] = 0;
BS[22] = n4j;

BS[1] = dn1dy;
BS[3] = -n1j;
BS[5] = 0;
BS[7] = dn2dy;
BS[9] = -n2j;
BS[11] = 0;
BS[13] = dn3dy;
BS[15] = -n3j;
BS[17] = 0;
BS[19] = dn4dy;
BS[21] = -n4j;
BS[23] = 0;

/*printf("BSshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[0], BS[2], BS[4], BS[6], BS[8], BS[10], BS[12], BS[14], BS[16], BS[18], BS[20], BS[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[1], BS[3], BS[5], BS[7], BS[9], BS[11], BS[13], BS[15], BS[17], BS[19], BS[21], BS[23]);*/

/*Shear terms*/
    
    
Kp[0] += ((BS[0]*F[0]+BS[1]*F[1])*BS[0] + (BS[0]*F[2]+BS[1]*F[3])*BS[1]);
Kp[1] += ((BS[2]*F[0]+BS[3]*F[1])*BS[0] + (BS[2]*F[2]+BS[3]*F[3])*BS[1]);
Kp[2] += ((BS[4]*F[0]+BS[5]*F[1])*BS[0] + (BS[4]*F[2]+BS[5]*F[3])*BS[1]);
Kp[3] += ((BS[6]*F[0]+BS[7]*F[1])*BS[0] + (BS[6]*F[2]+BS[7]*F[3])*BS[1]);
Kp[4] += ((BS[8]*F[0]+BS[9]*F[1])*BS[0] + (BS[8]*F[2]+BS[9]*F[3])*BS[1]);
Kp[5] += ((BS[10]*F[0]+BS[11]*F[1])*BS[0] + (BS[10]*F[2]+BS[11]*F[3])*BS[1]);
Kp[6] += ((BS[12]*F[0]+BS[13]*F[1])*BS[0] + (BS[12]*F[2]+BS[13]*F[3])*BS[1]);
Kp[7] += ((BS[14]*F[0]+BS[15]*F[1])*BS[0] + (BS[14]*F[2]+BS[15]*F[3])*BS[1]);
Kp[8] += ((BS[16]*F[0]+BS[17]*F[1])*BS[0] + (BS[16]*F[2]+BS[17]*F[3])*BS[1]);
Kp[9] += ((BS[18]*F[0]+BS[19]*F[1])*BS[0] + (BS[18]*F[2]+BS[19]*F[3])*BS[1]);
Kp[10] += ((BS[20]*F[0]+BS[21]*F[1])*BS[0] + (BS[20]*F[2]+BS[21]*F[3])*BS[1]);
Kp[11] += ((BS[22]*F[0]+BS[23]*F[1])*BS[0] + (BS[22]*F[2]+BS[23]*F[3])*BS[1]);

Kp[12] += ((BS[0]*F[0]+BS[1]*F[1])*BS[2] + (BS[0]*F[2]+BS[1]*F[3])*BS[3]);
Kp[13] += ((BS[2]*F[0]+BS[3]*F[1])*BS[2] + (BS[2]*F[2]+BS[3]*F[3])*BS[3]);
Kp[14] += ((BS[4]*F[0]+BS[5]*F[1])*BS[2] + (BS[4]*F[2]+BS[5]*F[3])*BS[3]);
Kp[15] += ((BS[6]*F[0]+BS[7]*F[1])*BS[2] + (BS[6]*F[2]+BS[7]*F[3])*BS[3]);
Kp[16] += ((BS[8]*F[0]+BS[9]*F[1])*BS[2] + (BS[8]*F[2]+BS[9]*F[3])*BS[3]);
Kp[17] += ((BS[10]*F[0]+BS[11]*F[1])*BS[2] + (BS[10]*F[2]+BS[11]*F[3])*BS[3]);
Kp[18] += ((BS[12]*F[0]+BS[13]*F[1])*BS[2] + (BS[12]*F[2]+BS[13]*F[3])*BS[3]);
Kp[19] += ((BS[14]*F[0]+BS[15]*F[1])*BS[2] + (BS[14]*F[2]+BS[15]*F[3])*BS[3]);
Kp[20] += ((BS[16]*F[0]+BS[17]*F[1])*BS[2] + (BS[16]*F[2]+BS[17]*F[3])*BS[3]);
Kp[21] += ((BS[18]*F[0]+BS[19]*F[1])*BS[2] + (BS[18]*F[2]+BS[19]*F[3])*BS[3]);
Kp[22] += ((BS[20]*F[0]+BS[21]*F[1])*BS[2] + (BS[20]*F[2]+BS[21]*F[3])*BS[3]);
Kp[23] += ((BS[22]*F[0]+BS[23]*F[1])*BS[2] + (BS[22]*F[2]+BS[23]*F[3])*BS[3]);

Kp[24] += ((BS[0]*F[0]+BS[1]*F[1])*BS[4] + (BS[0]*F[2]+BS[1]*F[3])*BS[5]);
Kp[25] += ((BS[2]*F[0]+BS[3]*F[1])*BS[4] + (BS[2]*F[2]+BS[3]*F[3])*BS[5]);
Kp[26] += ((BS[4]*F[0]+BS[5]*F[1])*BS[4] + (BS[4]*F[2]+BS[5]*F[3])*BS[5]);
Kp[27] += ((BS[6]*F[0]+BS[7]*F[1])*BS[4] + (BS[6]*F[2]+BS[7]*F[3])*BS[5]);
Kp[28] += ((BS[8]*F[0]+BS[9]*F[1])*BS[4] + (BS[8]*F[2]+BS[9]*F[3])*BS[5]);
Kp[29] += ((BS[10]*F[0]+BS[11]*F[1])*BS[4] + (BS[10]*F[2]+BS[11]*F[3])*BS[5]);
Kp[30] += ((BS[12]*F[0]+BS[13]*F[1])*BS[4] + (BS[12]*F[2]+BS[13]*F[3])*BS[5]);
Kp[31] += ((BS[14]*F[0]+BS[15]*F[1])*BS[4] + (BS[14]*F[2]+BS[15]*F[3])*BS[5]);

Kp[32] += ((BS[16]*F[0]+BS[17]*F[1])*BS[4] + (BS[16]*F[2]+BS[17]*F[3])*BS[5]);
Kp[33] += ((BS[18]*F[0]+BS[19]*F[1])*BS[4] + (BS[18]*F[2]+BS[19]*F[3])*BS[5]);
Kp[34] += ((BS[20]*F[0]+BS[21]*F[1])*BS[4] + (BS[20]*F[2]+BS[21]*F[3])*BS[5]);
Kp[35] += ((BS[22]*F[0]+BS[23]*F[1])*BS[4] + (BS[22]*F[2]+BS[23]*F[3])*BS[5]);

Kp[36] += ((BS[0]*F[0]+BS[1]*F[1])*BS[6] + (BS[0]*F[2]+BS[1]*F[3])*BS[7]);
Kp[37] += ((BS[2]*F[0]+BS[3]*F[1])*BS[6] + (BS[2]*F[2]+BS[3]*F[3])*BS[7]);
Kp[38] += ((BS[4]*F[0]+BS[5]*F[1])*BS[6] + (BS[4]*F[2]+BS[5]*F[3])*BS[7]);
Kp[39] += ((BS[6]*F[0]+BS[7]*F[1])*BS[6] + (BS[6]*F[2]+BS[7]*F[3])*BS[7]);
Kp[40] += ((BS[8]*F[0]+BS[9]*F[1])*BS[6] + (BS[8]*F[2]+BS[9]*F[3])*BS[7]);
Kp[41] += ((BS[10]*F[0]+BS[11]*F[1])*BS[6] + (BS[10]*F[2]+BS[11]*F[3])*BS[7]);
Kp[42] += ((BS[12]*F[0]+BS[13]*F[1])*BS[6] + (BS[12]*F[2]+BS[13]*F[3])*BS[7]);
Kp[43] += ((BS[14]*F[0]+BS[15]*F[1])*BS[6] + (BS[14]*F[2]+BS[15]*F[3])*BS[7]);
Kp[44] += ((BS[16]*F[0]+BS[17]*F[1])*BS[6] + (BS[16]*F[2]+BS[17]*F[3])*BS[7]);
Kp[45] += ((BS[18]*F[0]+BS[19]*F[1])*BS[6] + (BS[18]*F[2]+BS[19]*F[3])*BS[7]);
Kp[46] += ((BS[20]*F[0]+BS[21]*F[1])*BS[6] + (BS[20]*F[2]+BS[21]*F[3])*BS[7]);
Kp[47] += ((BS[22]*F[0]+BS[23]*F[1])*BS[6] + (BS[22]*F[2]+BS[23]*F[3])*BS[7]);

Kp[48] += ((BS[0]*F[0]+BS[1]*F[1])*BS[8] + (BS[0]*F[2]+BS[1]*F[3])*BS[9]);
Kp[49] += ((BS[2]*F[0]+BS[3]*F[1])*BS[8] + (BS[2]*F[2]+BS[3]*F[3])*BS[9]);
Kp[50] += ((BS[4]*F[0]+BS[5]*F[1])*BS[8] + (BS[4]*F[2]+BS[5]*F[3])*BS[9]);
Kp[51] += ((BS[6]*F[0]+BS[7]*F[1])*BS[8] + (BS[6]*F[2]+BS[7]*F[3])*BS[9]);
Kp[52] += ((BS[8]*F[0]+BS[9]*F[1])*BS[8] + (BS[8]*F[2]+BS[9]*F[3])*BS[9]);
Kp[53] += ((BS[10]*F[0]+BS[11]*F[1])*BS[8] + (BS[10]*F[2]+BS[11]*F[3])*BS[9]);
Kp[54] += ((BS[12]*F[0]+BS[13]*F[1])*BS[8] + (BS[12]*F[2]+BS[13]*F[3])*BS[9]);
Kp[55] += ((BS[14]*F[0]+BS[15]*F[1])*BS[8] + (BS[14]*F[2]+BS[15]*F[3])*BS[9]);
Kp[56] += ((BS[16]*F[0]+BS[17]*F[1])*BS[8] + (BS[16]*F[2]+BS[17]*F[3])*BS[9]);
Kp[57] += ((BS[18]*F[0]+BS[19]*F[1])*BS[8] + (BS[18]*F[2]+BS[19]*F[3])*BS[9]);
Kp[58] += ((BS[20]*F[0]+BS[21]*F[1])*BS[8] + (BS[20]*F[2]+BS[21]*F[3])*BS[9]);
Kp[59] += ((BS[22]*F[0]+BS[23]*F[1])*BS[8] + (BS[22]*F[2]+BS[23]*F[3])*BS[9]);

Kp[60] += ((BS[0]*F[0]+BS[1]*F[1])*BS[10] + (BS[0]*F[2]+BS[1]*F[3])*BS[11]);
Kp[61] += ((BS[2]*F[0]+BS[3]*F[1])*BS[10] + (BS[2]*F[2]+BS[3]*F[3])*BS[11]);
Kp[62] += ((BS[4]*F[0]+BS[5]*F[1])*BS[10] + (BS[4]*F[2]+BS[5]*F[3])*BS[11]);
Kp[63] += ((BS[6]*F[0]+BS[7]*F[1])*BS[10] + (BS[6]*F[2]+BS[7]*F[3])*BS[11]);
Kp[64] += ((BS[8]*F[0]+BS[9]*F[1])*BS[10] + (BS[8]*F[2]+BS[9]*F[3])*BS[11]);
Kp[65] += ((BS[10]*F[0]+BS[11]*F[1])*BS[10] + (BS[10]*F[2]+BS[11]*F[3])*BS[11]);
Kp[66] += ((BS[12]*F[0]+BS[13]*F[1])*BS[10] + (BS[12]*F[2]+BS[13]*F[3])*BS[11]);
Kp[67] += ((BS[14]*F[0]+BS[15]*F[1])*BS[10] + (BS[14]*F[2]+BS[15]*F[3])*BS[11]);
Kp[68] += ((BS[16]*F[0]+BS[17]*F[1])*BS[10] + (BS[16]*F[2]+BS[17]*F[3])*BS[11]);
Kp[69] += ((BS[18]*F[0]+BS[19]*F[1])*BS[10] + (BS[18]*F[2]+BS[19]*F[3])*BS[11]);
Kp[70] += ((BS[20]*F[0]+BS[21]*F[1])*BS[10] + (BS[20]*F[2]+BS[21]*F[3])*BS[11]);

Kp[71] += ((BS[22]*F[0]+BS[23]*F[1])*BS[10] + (BS[22]*F[2]+BS[23]*F[3])*BS[11]);

Kp[72] += ((BS[0]*F[0]+BS[1]*F[1])*BS[12] + (BS[0]*F[2]+BS[1]*F[3])*BS[13]);
Kp[73] += ((BS[2]*F[0]+BS[3]*F[1])*BS[12] + (BS[2]*F[2]+BS[3]*F[3])*BS[13]);
Kp[74] += ((BS[4]*F[0]+BS[5]*F[1])*BS[12] + (BS[4]*F[2]+BS[5]*F[3])*BS[13]);
Kp[75] += ((BS[6]*F[0]+BS[7]*F[1])*BS[12] + (BS[6]*F[2]+BS[7]*F[3])*BS[13]);
Kp[76] += ((BS[8]*F[0]+BS[9]*F[1])*BS[12] + (BS[8]*F[2]+BS[9]*F[3])*BS[13]);
Kp[77] += ((BS[10]*F[0]+BS[11]*F[1])*BS[12] + (BS[10]*F[2]+BS[11]*F[3])*BS[13]);
Kp[78] += ((BS[12]*F[0]+BS[13]*F[1])*BS[12] + (BS[12]*F[2]+BS[13]*F[3])*BS[13]);
Kp[79] += ((BS[14]*F[0]+BS[15]*F[1])*BS[12] + (BS[14]*F[2]+BS[15]*F[3])*BS[13]);
Kp[80] += ((BS[16]*F[0]+BS[17]*F[1])*BS[12] + (BS[16]*F[2]+BS[17]*F[3])*BS[13]);
Kp[81] += ((BS[18]*F[0]+BS[19]*F[1])*BS[12] + (BS[18]*F[2]+BS[19]*F[3])*BS[13]);
Kp[82] += ((BS[20]*F[0]+BS[21]*F[1])*BS[12] + (BS[20]*F[2]+BS[21]*F[3])*BS[13]);
Kp[83] += ((BS[22]*F[0]+BS[23]*F[1])*BS[12] + (BS[22]*F[2]+BS[23]*F[3])*BS[13]);

Kp[84] += ((BS[0]*F[0]+BS[1]*F[1])*BS[14] + (BS[0]*F[2]+BS[1]*F[3])*BS[15]);
Kp[85] += ((BS[2]*F[0]+BS[3]*F[1])*BS[14] + (BS[2]*F[2]+BS[3]*F[3])*BS[15]);
Kp[86] += ((BS[4]*F[0]+BS[5]*F[1])*BS[14] + (BS[4]*F[2]+BS[5]*F[3])*BS[15]);
Kp[87] += ((BS[6]*F[0]+BS[7]*F[1])*BS[14] + (BS[6]*F[2]+BS[7]*F[3])*BS[15]);
Kp[88] += ((BS[8]*F[0]+BS[9]*F[1])*BS[14] + (BS[8]*F[2]+BS[9]*F[3])*BS[15]);
Kp[89] += ((BS[10]*F[0]+BS[11]*F[1])*BS[14] + (BS[10]*F[2]+BS[11]*F[3])*BS[15]);
Kp[90] += ((BS[12]*F[0]+BS[13]*F[1])*BS[14] + (BS[12]*F[2]+BS[13]*F[3])*BS[15]);
Kp[91] += ((BS[14]*F[0]+BS[15]*F[1])*BS[14] + (BS[14]*F[2]+BS[15]*F[3])*BS[15]);
Kp[92] += ((BS[16]*F[0]+BS[17]*F[1])*BS[14] + (BS[16]*F[2]+BS[17]*F[3])*BS[15]);
Kp[93] += ((BS[18]*F[0]+BS[19]*F[1])*BS[14] + (BS[18]*F[2]+BS[19]*F[3])*BS[15]);
Kp[94] += ((BS[20]*F[0]+BS[21]*F[1])*BS[14] + (BS[20]*F[2]+BS[21]*F[3])*BS[15]);
Kp[95] += ((BS[22]*F[0]+BS[23]*F[1])*BS[14] + (BS[22]*F[2]+BS[23]*F[3])*BS[15]);

Kp[96] += ((BS[0]*F[0]+BS[1]*F[1])*BS[16] + (BS[0]*F[2]+BS[1]*F[3])*BS[17]);
Kp[97] += ((BS[2]*F[0]+BS[3]*F[1])*BS[16] + (BS[2]*F[2]+BS[3]*F[3])*BS[17]);
Kp[98] += ((BS[4]*F[0]+BS[5]*F[1])*BS[16] + (BS[4]*F[2]+BS[5]*F[3])*BS[17]);
Kp[99] += ((BS[6]*F[0]+BS[7]*F[1])*BS[16] + (BS[6]*F[2]+BS[7]*F[3])*BS[17]);
Kp[100] += ((BS[8]*F[0]+BS[9]*F[1])*BS[16] + (BS[8]*F[2]+BS[9]*F[3])*BS[17]);
Kp[101] += ((BS[10]*F[0]+BS[11]*F[1])*BS[16] + (BS[10]*F[2]+BS[11]*F[3])*BS[17]);
Kp[102] += ((BS[12]*F[0]+BS[13]*F[1])*BS[16] + (BS[12]*F[2]+BS[13]*F[3])*BS[17]);
Kp[103] += ((BS[14]*F[0]+BS[15]*F[1])*BS[16] + (BS[14]*F[2]+BS[15]*F[3])*BS[17]);
Kp[104] += ((BS[16]*F[0]+BS[17]*F[1])*BS[16] + (BS[16]*F[2]+BS[17]*F[3])*BS[17]);
Kp[105] += ((BS[18]*F[0]+BS[19]*F[1])*BS[16] + (BS[18]*F[2]+BS[19]*F[3])*BS[17]);
Kp[106] += ((BS[20]*F[0]+BS[21]*F[1])*BS[16] + (BS[20]*F[2]+BS[21]*F[3])*BS[17]);
Kp[107] += ((BS[22]*F[0]+BS[23]*F[1])*BS[16] + (BS[22]*F[2]+BS[23]*F[3])*BS[17]);

Kp[108] += ((BS[0]*F[0]+BS[1]*F[1])*BS[18] + (BS[0]*F[2]+BS[1]*F[3])*BS[19]);
Kp[109] += ((BS[2]*F[0]+BS[3]*F[1])*BS[18] + (BS[2]*F[2]+BS[3]*F[3])*BS[19]);
Kp[110] += ((BS[4]*F[0]+BS[5]*F[1])*BS[18] + (BS[4]*F[2]+BS[5]*F[3])*BS[19]);
Kp[111] += ((BS[6]*F[0]+BS[7]*F[1])*BS[18] + (BS[6]*F[2]+BS[7]*F[3])*BS[19]);
Kp[112] += ((BS[8]*F[0]+BS[9]*F[1])*BS[18] + (BS[8]*F[2]+BS[9]*F[3])*BS[19]);
Kp[113] += ((BS[10]*F[0]+BS[11]*F[1])*BS[18] + (BS[10]*F[2]+BS[11]*F[3])*BS[19]);
Kp[114] += ((BS[12]*F[0]+BS[13]*F[1])*BS[18] + (BS[12]*F[2]+BS[13]*F[3])*BS[19]);
Kp[115] += ((BS[14]*F[0]+BS[15]*F[1])*BS[18] + (BS[14]*F[2]+BS[15]*F[3])*BS[19]);
Kp[116] += ((BS[16]*F[0]+BS[17]*F[1])*BS[18] + (BS[16]*F[2]+BS[17]*F[3])*BS[19]);
Kp[117] += ((BS[18]*F[0]+BS[19]*F[1])*BS[18] + (BS[18]*F[2]+BS[19]*F[3])*BS[19]);
Kp[118] += ((BS[20]*F[0]+BS[21]*F[1])*BS[18] + (BS[20]*F[2]+BS[21]*F[3])*BS[19]);
Kp[119] += ((BS[22]*F[0]+BS[23]*F[1])*BS[18] + (BS[22]*F[2]+BS[23]*F[3])*BS[19]);

Kp[120] += ((BS[0]*F[0]+BS[1]*F[1])*BS[20] + (BS[0]*F[2]+BS[1]*F[3])*BS[21]);
Kp[121] += ((BS[2]*F[0]+BS[3]*F[1])*BS[20] + (BS[2]*F[2]+BS[3]*F[3])*BS[21]);
Kp[122] += ((BS[4]*F[0]+BS[5]*F[1])*BS[20] + (BS[4]*F[2]+BS[5]*F[3])*BS[21]);
Kp[123] += ((BS[6]*F[0]+BS[7]*F[1])*BS[20] + (BS[6]*F[2]+BS[7]*F[3])*BS[21]);
Kp[124] += ((BS[8]*F[0]+BS[9]*F[1])*BS[20] + (BS[8]*F[2]+BS[9]*F[3])*BS[21]);
Kp[125] += ((BS[10]*F[0]+BS[11]*F[1])*BS[20] + (BS[10]*F[2]+BS[11]*F[3])*BS[21]);
Kp[126] += ((BS[12]*F[0]+BS[13]*F[1])*BS[20] + (BS[12]*F[2]+BS[13]*F[3])*BS[21]);
Kp[127] += ((BS[14]*F[0]+BS[15]*F[1])*BS[20] + (BS[14]*F[2]+BS[15]*F[3])*BS[21]);
Kp[128] += ((BS[16]*F[0]+BS[17]*F[1])*BS[20] + (BS[16]*F[2]+BS[17]*F[3])*BS[21]);
Kp[129] += ((BS[18]*F[0]+BS[19]*F[1])*BS[20] + (BS[18]*F[2]+BS[19]*F[3])*BS[21]);
Kp[130] += ((BS[20]*F[0]+BS[21]*F[1])*BS[20] + (BS[20]*F[2]+BS[21]*F[3])*BS[21]);
Kp[131] += ((BS[22]*F[0]+BS[23]*F[1])*BS[20] + (BS[22]*F[2]+BS[23]*F[3])*BS[21]);

Kp[132] += ((BS[0]*F[0]+BS[1]*F[1])*BS[22] + (BS[0]*F[2]+BS[1]*F[3])*BS[23]);
Kp[133] += ((BS[2]*F[0]+BS[3]*F[1])*BS[22] + (BS[2]*F[2]+BS[3]*F[3])*BS[23]);
Kp[134] += ((BS[4]*F[0]+BS[5]*F[1])*BS[22] + (BS[4]*F[2]+BS[5]*F[3])*BS[23]);
Kp[135] += ((BS[6]*F[0]+BS[7]*F[1])*BS[22] + (BS[6]*F[2]+BS[7]*F[3])*BS[23]);
Kp[136] += ((BS[8]*F[0]+BS[9]*F[1])*BS[22] + (BS[8]*F[2]+BS[9]*F[3])*BS[23]);
Kp[137] += ((BS[10]*F[0]+BS[11]*F[1])*BS[22] + (BS[10]*F[2]+BS[11]*F[3])*BS[23]);
Kp[138] += ((BS[12]*F[0]+BS[13]*F[1])*BS[22] + (BS[12]*F[2]+BS[13]*F[3])*BS[23]);
Kp[139] += ((BS[14]*F[0]+BS[15]*F[1])*BS[22] + (BS[14]*F[2]+BS[15]*F[3])*BS[23]);
Kp[140] += ((BS[16]*F[0]+BS[17]*F[1])*BS[22] + (BS[16]*F[2]+BS[17]*F[3])*BS[23]);
Kp[141] += ((BS[18]*F[0]+BS[19]*F[1])*BS[22] + (BS[18]*F[2]+BS[19]*F[3])*BS[23]);
Kp[142] += ((BS[20]*F[0]+BS[21]*F[1])*BS[22] + (BS[20]*F[2]+BS[21]*F[3])*BS[23]);
Kp[143] += ((BS[22]*F[0]+BS[23]*F[1])*BS[22] + (BS[22]*F[2]+BS[23]*F[3])*BS[23]);

free(F);
free(BS);

/*if((num == 361)||(num==378))
{
printf("\n Membrain Matrix");
for(p=0;p<8;p++)
{
	printf("\n");
	for(j=0;j<8;j++)
	{
		printf("%0.3f\t", Km[j*8+p]);
	}
}
printf("\n");
}*/

/*Check the reuslts so far*/
/*printf("\n Membrain Matrix");
for(p=0;p<8;p++)
{
	printf("\n");
	for(j=0;j<8;j++)
	{
		printf("%0.3f\t", Km[j*8+p]);
	}
}
printf("\n");
printf("\n Plate Matrix");
for(i=0;i<12;i++)
{
	printf("\n");
	for(j=0;j<12;j++)
	{
		printf("%0.3f\t", Kp[j*12+i]);
	}
}*/

/*Assemble Full Matirx (Try to transpose it so that it is correct in C matrix form by not hugely important since it should be symetrical*/
for(i=0;i<576;i++)
{
	KE[i] = 0.0000000000000000;	/*Intialize every value to zero*/
}

/*Fill in the values*/
for(j=0;j<4;j++)
{
	p =24*6*j;
	KE[p] = Km[0+2*j];
	KE[1+p] = Km[8+2*j];
	KE[6+p] = Km[16+2*j];
	KE[7+p] = Km[24+2*j];
	KE[12+p] = Km[32+2*j];
	KE[13+p] = Km[40+2*j];
	KE[18+p] = Km[48+2*j];
	KE[19+p] = Km[56+2*j];
	KE[24+p] = Km[1+2*j];

	KE[25+p] = Km[9+2*j];
	KE[30+p] = Km[17+2*j];
	KE[31+p] = Km[25+2*j];
	KE[36+p] = Km[33+2*j];
	KE[37+p] = Km[41+2*j];
	KE[42+p] = Km[49+2*j];
	KE[43+p] = Km[57+2*j];

	q =24*6*j+24*2;
	KE[2+q] = Kp[0+3*j];
	KE[3+q] = Kp[12+3*j];
	KE[4+q] = Kp[24+3*j];
	KE[8+q] = Kp[36+3*j];
	KE[9+q] = Kp[48+3*j];
	KE[10+q] = Kp[60+3*j];
	KE[14+q] = Kp[72+3*j];
	KE[15+q] = Kp[84+3*j];
	KE[16+q] = Kp[96+3*j];
	KE[20+q] = Kp[108+3*j];

	KE[21+q] = Kp[120+3*j];
	KE[22+q] = Kp[132+3*j];
	q =24*6*j+24*3;
	KE[2+q] = Kp[1+3*j];
	KE[3+q] = Kp[13+3*j];
	KE[4+q] = Kp[25+3*j];
	KE[8+q] = Kp[37+3*j];
	KE[9+q] = Kp[49+3*j];
	KE[10+q] = Kp[61+3*j];
	KE[14+q] = Kp[73+3*j];
	KE[15+q] = Kp[85+3*j];
	KE[16+q] = Kp[97+3*j];
	KE[20+q] = Kp[109+3*j];
	KE[21+q] = Kp[121+3*j];
	KE[22+q] = Kp[133+3*j];
	q =24*6*j+24*4;
	KE[2+q] = Kp[2+3*j];
	KE[3+q] = Kp[14+3*j];
	KE[4+q] = Kp[26+3*j];
	KE[8+q] = Kp[38+3*j];
	KE[9+q] = Kp[50+3*j];
	KE[10+q] = Kp[62+3*j];
	KE[14+q] = Kp[74+3*j];
	KE[15+q] = Kp[86+3*j];
	KE[16+q] = Kp[98+3*j];
	KE[20+q] = Kp[110+3*j];
	KE[21+q] = Kp[122+3*j];
	KE[22+q] = Kp[134+3*j];

}

/*printf("\n Element Matrix");
for(i=0;i<24;i++)
{
	printf("\n");
	for(j=0;j<24;j++)
	{
		printf("%0.2f  ", KE[j+i*24]);

	}
}*/

free(Km);
free(Kp);
/*printf("\n Temp Membrain Matrix");*/
/*for(p=0;p<12;p++)
{
	printf("\n");
	for(j=0;j<3;j++)
	{
		printf("%0.3f\t", tempKp[j*12+p]);
	}
}*/
/*printf("\nB");
printf("\nB");	*/	

}

/*Gets the stiffness matirx for Stress calculation*/ 
void CompositeKEDIFF(double Q11, double Q22, double Q12, double Q16, double Q26, double Q66, double Q44, double Q55, double Q45, double *KE, int numply, double plyt, double h, double hxz, int num)
{
/*NOTE For the sake of the dgemm functions set up here is column strong (arrays go down thecolumns first, then along the rows). The fina out put however will be the other way around so that we can use it in assmeble*/
int n,m,i,j,p,q,GP;
double dn1dx,dn1dy,dn2dx,dn2dy,dn3dx,dn3dy,dn4dx,dn4dy,nx,ny,n1j,n2j,n3j,n4j;
/*Intialize all the matrices we'll be using*/
double *A, *D, *B, *Cmatrix, *F, *Qmatrix;
/*Material property matrices*/
A = malloc(9*sizeof(double));
D = malloc(9*sizeof(double));
B = malloc(9*sizeof(double));
F = malloc(4*sizeof(double));
Cmatrix = malloc(4*sizeof(double));
Qmatrix = malloc(9*sizeof(double));
/*Shape function Matrices*/
double *BB,*BM,*BS;
BM = malloc(24*sizeof(double));
BB = malloc(36*sizeof(double));
BS = malloc(24*sizeof(double));

/*StiffnessMatrixes*/
double *Km, *Kp;
Km = malloc(64*sizeof(double));
Kp = malloc(144*sizeof(double));

	for(n=0;n<3;n++)
	{
		for(m=0;m<3;m++)
		{
			A[n+3*m] = 0;
			D[n+3*m] = 0;
			B[n+3*m] = 0;
		}
	}
	for(p=0;p<8;p++)
	{
		for(q=0;q<8;q++)
		{
			Km[p+8*q] = 0;
		}
	}
	for(p=0;p<12;p++)
	{
		for(q=0;q<12;q++)
		{
			Kp[p+12*q] = 0;
		}
	}				
		
	for(q=0;q<numply;q++)
	{

		Qmatrix[0] = Q11;
		Qmatrix[1] = Q12;
		Qmatrix[2] = Q16;
		Qmatrix[3] = Q12;
		Qmatrix[4] = Q22;
		Qmatrix[5] = Q26;
		Qmatrix[6] = Q16;
		Qmatrix[7] = Q26;
		Qmatrix[8] = Q66;

		Cmatrix[0] = Q44;
		Cmatrix[1] = Q45;
		Cmatrix[2] = Q45;
		Cmatrix[3] = Q55;

		for(n=0;n<3;n++)
		{
			for(m=0;m<3;m++)
			{
				A[n+3*m] += plyt*Qmatrix[n+3*m];
				B[n+3*m] += -1*0.5*plyt*plyt*Qmatrix[n+3*m];
				D[n+3*m] += (plyt*plyt*plyt*Qmatrix[n+3*m])/12.0;
				if((n<2)&&(m<2)){F[n+2*m] = (5.0/6.0)*plyt*Cmatrix[n+2*m];}
			}
		}				
	}
	/*if((num == 0))

	{
	printf("\nCmatrix\n");
	printf("%f\t %f\n", Cmatrix[0], Cmatrix[2]);
	printf("%f\t %f\n", Cmatrix[1], Cmatrix[3]);


	printf("\nQmatrix\n");
	printf("%f\t %f\t %f\n", Qmatrix[0], Qmatrix[3], Qmatrix[6]);
	printf("%f\t %f\t %f\n", Qmatrix[1], Qmatrix[4], Qmatrix[7]);
	printf("%f\t %f\t %f\n", Qmatrix[2], Qmatrix[5], Qmatrix[8]);


	/*printf("\nAmatrix\n");
	printf("%f\t %f\t %f\n", A[0], A[3], A[6]);
	printf("%f\t %f\t %f\n", A[1], A[4], A[7]);
	printf("%f\t %f\t %f\n", A[2], A[5], A[8]);

	printf("\nBmatrix\n");
	printf("%f\t %f\t %f\n", B[0], B[3], B[6]);
	printf("%f\t %f\t %f\n", B[1], B[4], B[7]);
	printf("%f\t %f\t %f\n", B[2], B[5], B[8]);


	printf("\nDmatrix\n");
	printf("%f\t %f\t %f\n", D[0], D[3], D[6]);
	printf("%f\t %f\t %f\n", D[1], D[4], D[7]);
	printf("%f\t %f\t %f\n", D[2], D[5], D[8]);

	printf("\nFmatrix\n");
	printf("%f\t %f\n", F[0], F[2]);
	printf("%f\t %f\n", F[1], F[3]);
	}*/

	free(Cmatrix);
	free(Qmatrix);

for(GP=0;GP<4;GP++)
{
		if((GP==0)||(GP==3)){nx = -0.5773502692;}
		else{nx = 0.5773502692;}

		if((GP==0)||(GP==1)){ny = -0.5773502692;}
		else{ny = 0.5773502692;}

		/*printf("\n GP = %i, nx = %f, ny = %f", GP, nx, ny);*/


		n1j = 0.125*(1-nx)*(1-ny);
		n2j = 0.125*(1+nx)*(1-ny);
		n3j = 0.125*(1+nx)*(1+ny);
		n4j = 0.125*(1-nx)*(1+ny);

		dn1dx = -0.25*(1-ny)/hxz;
		dn2dx = 0.25*(1-ny)/hxz;
		dn3dx = 0.25*(1+ny)/hxz;
		dn4dx = -0.25*(1+ny)/hxz;

		dn1dy = -0.25*(1-nx)/h;
		dn2dy = -0.25*(1+nx)/h;
		dn3dy = 0.25*(1+nx)/h;
		dn4dy = 0.25*(1-nx)/h;

		/*printf("\n dn1dx = %f", dn1dx);
		printf("\n dn1dx = %f", dn2dx);
		printf("\n dn1dx = %f", dn3dx);
		printf("\n dn1dx = %f", dn4dx);

		printf("\n dn1dx = %f", dn1dy);
		printf("\n dn1dx = %f", dn2dy);
		printf("\n dn1dx = %f", dn3dy);

		printf("\n dn1dx = %f", dn4dy);


		printf("\n n1j = %f", n1j);
		printf("\n n2j = %f", n2j);
		printf("\n n3j = %f", n3j);

		printf("\n n4j = %f\n", n4j);*/

		BM[0] = dn1dx;
		BM[3] = 0;
		BM[6] = dn2dx;
		BM[9] = 0;
		BM[12] = dn3dx;
		BM[15] = 0;
		BM[18] = dn4dx;
		BM[21] = 0;
		
		BM[1] = 0;
		BM[4] = dn1dy;
		BM[7] = 0;
		BM[10] = dn2dy;
		BM[13] = 0;
		BM[16] = dn3dy;
		BM[19] = 0;
		BM[22] = dn4dy;

		BM[2] = dn1dy;
		BM[5] = dn1dx;
		BM[8] = dn2dy;
		BM[11] = dn2dx;
		BM[14] = dn3dy;
		BM[17] = dn3dx;
		BM[20] = dn4dy;
		BM[23] = dn4dx;

		BB[0] = 0;
		BB[3] = 0;
		BB[6] = -dn1dx;
		BB[9] = 0;
		BB[12] = 0;
		BB[15] = -dn2dx;
		BB[18] = 0;
		BB[21] = 0;
		BB[24] = -dn3dx;
		BB[27] = 0;
		BB[30] = 0;
		BB[33] = -dn4dx;

		BB[1] = 0;
		BB[4] = dn1dy;
		BB[7] = 0;
		BB[10] = 0;
		BB[13] = dn2dy;
		BB[16] = 0;
		BB[19] = 0;
		BB[22] = dn3dy;
		BB[25] = 0;
		BB[28] = 0;
		BB[31] = dn4dy;
		BB[34] = 0;

		BB[2] = 0;
		BB[5] = dn1dx;
		BB[8] = -dn1dy;
		BB[11] = 0;
		BB[14] = dn2dx;
		BB[17] = -dn2dy;
		BB[20] = 0;
		BB[23] = dn3dx;
		BB[26] = -dn3dy;
		BB[29] = 0;
		BB[32] = dn4dx;
		BB[35] = -dn4dy;

		BS[0] = dn1dx;
		BS[2] = 0;
		BS[4] = n1j;
		BS[6] = dn2dx;
		BS[8] = 0;
		BS[10] = n2j;
		BS[12] = dn3dx;
		BS[14] = 0;
		BS[16] = n3j;
		BS[18] = dn4dx;
		BS[20] = 0;
		BS[22] = n4j;

		BS[1] = dn1dy;
		BS[3] = -n1j;
		BS[5] = 0;
		BS[7] = dn2dy;
		BS[9] = -n2j;
		BS[11] = 0;
		BS[13] = dn3dy;
		BS[15] = -n3j;
		BS[17] = 0;
		BS[19] = dn4dy;
		BS[21] = -n4j;
		BS[23] = 0;

/*if((num == 361)||(num==378))

{
printf("BMshape\n");

printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[0], BM[3], BM[6], BM[9], BM[12], BM[15], BM[18], BM[21]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[1], BM[4], BM[7], BM[10], BM[13], BM[16], BM[19], BM[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BM[2], BM[5], BM[8], BM[11], BM[14], BM[17], BM[20], BM[23]);
}

/*printf("BBshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[0], BB[3], BB[6], BB[9], BB[12], BB[15], BB[18], BB[21], BB[24], BB[27], BB[30], BB[33]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[1], BB[4], BB[7], BB[10], BB[13], BB[16], BB[19], BB[22], BB[25], BB[28], BB[31], BB[34]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BB[2], BB[5], BB[8], BB[11], BB[14], BB[17], BB[20], BB[23], BB[26], BB[29], BB[32], BB[35]);

printf("BSshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[0], BS[2], BS[4], BS[6], BS[8], BS[10], BS[12], BS[14], BS[16], BS[18], BS[20], BS[22]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[1], BS[3], BS[5], BS[7], BS[9], BS[11], BS[13], BS[15], BS[17], BS[19], BS[21], BS[23]);*/



		/*Multiply the shape functions and the material matrices together using dgemm*/
/*printf("\nA");

printf("\nA");*/

/*Membrain terms*/		
Km[0] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[0] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[1] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[2]);
Km[1] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[0] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[1] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[2]);
Km[2] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[0] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[1] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[2]);
Km[3] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[0] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[1] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[2]);
Km[4] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[0] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[1] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[2]);
Km[5] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[0] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[1] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[2]);
Km[6] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[0] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[1] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[2]);
Km[7] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[0] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[1] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[2]);

Km[8] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[3] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[4] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[5]);
Km[9] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[3] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[4] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[5]);
Km[10] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[3] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[4] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[5]);
Km[11] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[3] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[4] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[5]);
Km[12] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[3] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[4] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[5]);
Km[13] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[3] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[4] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[5]);
Km[14] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[3] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[4] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[5]);
Km[15] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[3] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[4] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[5]);

Km[16] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[6] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[7] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[8]);
Km[17] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[6] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[7] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[8]);
Km[18] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[6] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[7] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[8]);
Km[19] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[6] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[7] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[8]);
Km[20] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[6] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[7] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[8]);
Km[21] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[6] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[7] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[8]);
Km[22] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[6] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[7] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[8]);
Km[23] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[6] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[7] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[8]);

Km[24] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[9] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[10] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[11]);
Km[25] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[9] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[10] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[11]);
Km[26] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[9] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[10] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[11]);
Km[27] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[9] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[10] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[11]);
Km[28] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[9] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[10] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[11]);
Km[29] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[9] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[10] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[11]);
Km[30] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[9] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[10] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[11]);
Km[31] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[9] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[10] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[11]);
		
Km[32] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[12] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[13] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[14]);
Km[33] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[12] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[13] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[14]);
Km[34] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[12] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[13] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[14]);
Km[35] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[12] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[13] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[14]);
Km[36] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[12] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[13] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[14]);
Km[37] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[12] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[13] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[14]);
Km[38] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[12] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[13] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[14]);
Km[39] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[12] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[13] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[14]);

Km[40] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[15] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[16] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[17]);
Km[41] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[15] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[16] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[17]);
Km[42] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[15] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[16] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[17]);
Km[43] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[15] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[16] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[17]);
Km[44] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[15] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[16] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[17]);
Km[45] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[15] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[16] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[17]);
Km[46] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[15] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[16] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[17]);
Km[47] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[15] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[16] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[17]);

Km[48] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[18] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[19] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[20]);
Km[49] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[18] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[19] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[20]);
Km[50] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[18] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[19] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[20]);
Km[51] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[18] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[19] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[20]);
Km[52] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[18] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[19] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[20]);
Km[53] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[18] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[19] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[20]);
Km[54] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[18] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[19] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[20]);
Km[55] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[18] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[19] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[20]);

Km[56] += ((BM[0]*A[0]+BM[1]*A[1]+BM[2]*A[2])*BM[21] + (BM[0]*A[3]+BM[1]*A[4]+BM[2]*A[5])*BM[22] + (BM[0]*A[6]+BM[1]*A[7]+BM[2]*A[8])*BM[23]);
Km[57] += ((BM[3]*A[0]+BM[4]*A[1]+BM[5]*A[2])*BM[21] + (BM[3]*A[3]+BM[4]*A[4]+BM[5]*A[5])*BM[22] + (BM[3]*A[6]+BM[4]*A[7]+BM[5]*A[8])*BM[23]);
Km[58] += ((BM[6]*A[0]+BM[7]*A[1]+BM[8]*A[2])*BM[21] + (BM[6]*A[3]+BM[7]*A[4]+BM[8]*A[5])*BM[22] + (BM[6]*A[6]+BM[7]*A[7]+BM[8]*A[8])*BM[23]);
Km[59] += ((BM[9]*A[0]+BM[10]*A[1]+BM[11]*A[2])*BM[21] + (BM[9]*A[3]+BM[10]*A[4]+BM[11]*A[5])*BM[22] + (BM[9]*A[6]+BM[10]*A[7]+BM[11]*A[8])*BM[23]);
Km[60] += ((BM[12]*A[0]+BM[13]*A[1]+BM[14]*A[2])*BM[21] + (BM[12]*A[3]+BM[13]*A[4]+BM[14]*A[5])*BM[22] + (BM[12]*A[6]+BM[13]*A[7]+BM[14]*A[8])*BM[23]);
Km[61] += ((BM[15]*A[0]+BM[16]*A[1]+BM[17]*A[2])*BM[21] + (BM[15]*A[3]+BM[16]*A[4]+BM[17]*A[5])*BM[22] + (BM[15]*A[6]+BM[16]*A[7]+BM[17]*A[8])*BM[23]);
Km[62] += ((BM[18]*A[0]+BM[19]*A[1]+BM[20]*A[2])*BM[21] + (BM[18]*A[3]+BM[19]*A[4]+BM[20]*A[5])*BM[22] + (BM[18]*A[6]+BM[19]*A[7]+BM[20]*A[8])*BM[23]);
Km[63] += ((BM[21]*A[0]+BM[22]*A[1]+BM[23]*A[2])*BM[21] + (BM[21]*A[3]+BM[22]*A[4]+BM[23]*A[5])*BM[22] + (BM[21]*A[6]+BM[22]*A[7]+BM[23]*A[8])*BM[23]);


/*Bending terms*/
Kp[0] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[0] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[1] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[2]);
Kp[1] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[0] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[1] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[2]);
Kp[2] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[0] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[1] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[2]);
Kp[3] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[0] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[1] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[2]);
Kp[4] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[0] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[1] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[2]);
Kp[5] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[0] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[1] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[2]);
Kp[6] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[0] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[1] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[2]);
Kp[7] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[0] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[1] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[2]);
Kp[8] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[0] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[1] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[2]);
Kp[9] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[0] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[1] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[2]);
Kp[10] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[0] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[1] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[2]);
Kp[11] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[0] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[1] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[2]);

Kp[12] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[3] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[4] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[5]);
Kp[13] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[3] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[4] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[5]);
Kp[14] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[3] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[4] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[5]);
Kp[15] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[3] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[4] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[5]);
Kp[16] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[3] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[4] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[5]);
Kp[17] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[3] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[4] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[5]);
Kp[18] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[3] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[4] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[5]);
Kp[19] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[3] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[4] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[5]);
Kp[20] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[3] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[4] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[5]);
Kp[21] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[3] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[4] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[5]);
Kp[22] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[3] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[4] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[5]);
Kp[23] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[3] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[4] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[5]);

Kp[24] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[6] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[7] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[8]);
Kp[25] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[6] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[7] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[8]);
Kp[26] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[6] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[7] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[8]);
Kp[27] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[6] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[7] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[8]);
Kp[28] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[6] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[7] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[8]);
Kp[29] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[6] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[7] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[8]);
Kp[30] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[6] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[7] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[8]);
Kp[31] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[6] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[7] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[8]);
Kp[32] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[6] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[7] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[8]);
Kp[33] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[6] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[7] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[8]);
Kp[34] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[6] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[7] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[8]);
Kp[35] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[6] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[7] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[8]);

Kp[36] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[9] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[10] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[11]);
Kp[37] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[9] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[10] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[11]);
Kp[38] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[9] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[10] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[11]);
Kp[39] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[9] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[10] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[11]);
Kp[40] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[9] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[10] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[11]);
Kp[41] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[9] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[10] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[11]);
Kp[42] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[9] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[10] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[11]);
Kp[43] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[9] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[10] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[11]);
Kp[44] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[9] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[10] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[11]);
Kp[45] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[9] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[10] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[11]);
Kp[46] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[9] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[10] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[11]);
Kp[47] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[9] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[10] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[11]);

Kp[48] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[12] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[13] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[14]);
Kp[49] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[12] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[13] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[14]);
Kp[50] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[12] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[13] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[14]);
Kp[51] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[12] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[13] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[14]);
Kp[52] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[12] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[13] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[14]);
Kp[53] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[12] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[13] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[14]);
Kp[54] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[12] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[13] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[14]);
Kp[55] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[12] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[13] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[14]);
Kp[56] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[12] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[13] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[14]);
Kp[57] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[12] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[13] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[14]);
Kp[58] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[12] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[13] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[14]);
Kp[59] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[12] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[13] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[14]);

Kp[60] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[15] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[16] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[17]);
Kp[61] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[15] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[16] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[17]);
Kp[62] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[15] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[16] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[17]);
Kp[63] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[15] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[16] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[17]);
Kp[64] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[15] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[16] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[17]);
Kp[65] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[15] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[16] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[17]);
Kp[66] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[15] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[16] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[17]);
Kp[67] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[15] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[16] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[17]);
Kp[68] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[15] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[16] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[17]);
Kp[69] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[15] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[16] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[17]);
Kp[70] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[15] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[16] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[17]);
Kp[71] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[15] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[16] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[17]);

Kp[72] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[18] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[19] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[20]);
Kp[73] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[18] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[19] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[20]);
Kp[74] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[18] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[19] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[20]);
Kp[75] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[18] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[19] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[20]);
Kp[76] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[18] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[19] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[20]);
Kp[77] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[18] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[19] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[20]);
Kp[78] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[18] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[19] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[20]);
Kp[79] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[18] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[19] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[20]);
Kp[80] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[18] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[19] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[20]);
Kp[81] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[18] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[19] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[20]);
Kp[82] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[18] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[19] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[20]);
Kp[83] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[18] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[19] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[20]);

Kp[84] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[21] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[22] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[23]);
Kp[85] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[21] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[22] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[23]);
Kp[86] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[21] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[22] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[23]);
Kp[87] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[21] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[22] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[23]);
Kp[88] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[21] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[22] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[23]);
Kp[89] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[21] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[22] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[23]);
Kp[90] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[21] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[22] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[23]);
Kp[91] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[21] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[22] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[23]);
Kp[92] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[21] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[22] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[23]);
Kp[93] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[21] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[22] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[23]);
Kp[94] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[21] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[22] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[23]);
Kp[95] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[21] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[22] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[23]);

Kp[96] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[24] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[25] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[26]);
Kp[97] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[24] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[25] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[26]);
Kp[98] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[24] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[25] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[26]);
Kp[99] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[24] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[25] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[26]);
Kp[100] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[24] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[25] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[26]);
Kp[101] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[24] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[25] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[26]);
Kp[102] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[24] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[25] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[26]);
Kp[103] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[24] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[25] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[26]);
Kp[104] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[24] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[25] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[26]);
Kp[105] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[24] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[25] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[26]);
Kp[106] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[24] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[25] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[26]);
Kp[107] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[24] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[25] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[26]);

Kp[108] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[27] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[28] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[29]);
Kp[109] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[27] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[28] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[29]);
Kp[110] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[27] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[28] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[29]);
Kp[111] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[27] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[28] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[29]);
Kp[112] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[27] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[28] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[29]);
Kp[113] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[27] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[28] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[29]);
Kp[114] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[27] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[28] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[29]);
Kp[115] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[27] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[28] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[29]);
Kp[116] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[27] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[28] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[29]);
Kp[117] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[27] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[28] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[29]);
Kp[118] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[27] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[28] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[29]);
Kp[119] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[27] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[28] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[29]);

Kp[120] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[30] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[31] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[32]);
Kp[121] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[30] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[31] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[32]);
Kp[122] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[30] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[31] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[32]);
Kp[123] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[30] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[31] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[32]);
Kp[124] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[30] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[31] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[32]);
Kp[125] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[30] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[31] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[32]);
Kp[126] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[30] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[31] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[32]);
Kp[127] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[30] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[31] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[32]);
Kp[128] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[30] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[31] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[32]);
Kp[129] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[30] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[31] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[32]);
Kp[130] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[30] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[31] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[32]);
Kp[131] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[30] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[31] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[32]);

Kp[132] += ((BB[0]*D[0]+BB[1]*D[1]+BB[2]*D[2])*BB[33] + (BB[0]*D[3]+BB[1]*D[4]+BB[2]*D[5])*BB[34] + (BB[0]*D[6]+BB[1]*D[7]+BB[2]*D[8])*BB[35]);
Kp[133] += ((BB[3]*D[0]+BB[4]*D[1]+BB[5]*D[2])*BB[33] + (BB[3]*D[3]+BB[4]*D[4]+BB[5]*D[5])*BB[34] + (BB[3]*D[6]+BB[4]*D[7]+BB[5]*D[8])*BB[35]);
Kp[134] += ((BB[6]*D[0]+BB[7]*D[1]+BB[8]*D[2])*BB[33] + (BB[6]*D[3]+BB[7]*D[4]+BB[8]*D[5])*BB[34] + (BB[6]*D[6]+BB[7]*D[7]+BB[8]*D[8])*BB[35]);
Kp[135] += ((BB[9]*D[0]+BB[10]*D[1]+BB[11]*D[2])*BB[33] + (BB[9]*D[3]+BB[10]*D[4]+BB[11]*D[5])*BB[34] + (BB[9]*D[6]+BB[10]*D[7]+BB[11]*D[8])*BB[35]);
Kp[136] += ((BB[12]*D[0]+BB[13]*D[1]+BB[14]*D[2])*BB[33] + (BB[12]*D[3]+BB[13]*D[4]+BB[14]*D[5])*BB[34] + (BB[12]*D[6]+BB[13]*D[7]+BB[14]*D[8])*BB[35]);
Kp[137] += ((BB[15]*D[0]+BB[16]*D[1]+BB[17]*D[2])*BB[33] + (BB[15]*D[3]+BB[16]*D[4]+BB[17]*D[5])*BB[34] + (BB[15]*D[6]+BB[16]*D[7]+BB[17]*D[8])*BB[35]);
Kp[138] += ((BB[18]*D[0]+BB[19]*D[1]+BB[20]*D[2])*BB[33] + (BB[18]*D[3]+BB[19]*D[4]+BB[20]*D[5])*BB[34] + (BB[18]*D[6]+BB[19]*D[7]+BB[20]*D[8])*BB[35]);
Kp[139] += ((BB[21]*D[0]+BB[22]*D[1]+BB[23]*D[2])*BB[33] + (BB[21]*D[3]+BB[22]*D[4]+BB[23]*D[5])*BB[34] + (BB[21]*D[6]+BB[22]*D[7]+BB[23]*D[8])*BB[35]);
Kp[140] += ((BB[24]*D[0]+BB[25]*D[1]+BB[26]*D[2])*BB[33] + (BB[24]*D[3]+BB[25]*D[4]+BB[26]*D[5])*BB[34] + (BB[24]*D[6]+BB[25]*D[7]+BB[26]*D[8])*BB[35]);
Kp[141] += ((BB[27]*D[0]+BB[28]*D[1]+BB[29]*D[2])*BB[33] + (BB[27]*D[3]+BB[28]*D[4]+BB[29]*D[5])*BB[34] + (BB[27]*D[6]+BB[28]*D[7]+BB[29]*D[8])*BB[35]);
Kp[142] += ((BB[30]*D[0]+BB[31]*D[1]+BB[32]*D[2])*BB[33] + (BB[30]*D[3]+BB[31]*D[4]+BB[32]*D[5])*BB[34] + (BB[30]*D[6]+BB[31]*D[7]+BB[32]*D[8])*BB[35]);
Kp[143] += ((BB[33]*D[0]+BB[34]*D[1]+BB[35]*D[2])*BB[33] + (BB[33]*D[3]+BB[34]*D[4]+BB[35]*D[5])*BB[34] + (BB[33]*D[6]+BB[34]*D[7]+BB[35]*D[8])*BB[35]);

//////shear term
/*
    Kp[0] += ((BS[0]*F[0]+BS[1]*F[1])*BS[0] + (BS[0]*F[2]+BS[1]*F[3])*BS[1]);
    Kp[1] += ((BS[2]*F[0]+BS[3]*F[1])*BS[0] + (BS[2]*F[2]+BS[3]*F[3])*BS[1]);
    Kp[2] += ((BS[4]*F[0]+BS[5]*F[1])*BS[0] + (BS[4]*F[2]+BS[5]*F[3])*BS[1]);
    Kp[3] += ((BS[6]*F[0]+BS[7]*F[1])*BS[0] + (BS[6]*F[2]+BS[7]*F[3])*BS[1]);
    Kp[4] += ((BS[8]*F[0]+BS[9]*F[1])*BS[0] + (BS[8]*F[2]+BS[9]*F[3])*BS[1]);
    Kp[5] += ((BS[10]*F[0]+BS[11]*F[1])*BS[0] + (BS[10]*F[2]+BS[11]*F[3])*BS[1]);
    Kp[6] += ((BS[12]*F[0]+BS[13]*F[1])*BS[0] + (BS[12]*F[2]+BS[13]*F[3])*BS[1]);
    Kp[7] += ((BS[14]*F[0]+BS[15]*F[1])*BS[0] + (BS[14]*F[2]+BS[15]*F[3])*BS[1]);
    Kp[8] += ((BS[16]*F[0]+BS[17]*F[1])*BS[0] + (BS[16]*F[2]+BS[17]*F[3])*BS[1]);
    Kp[9] += ((BS[18]*F[0]+BS[19]*F[1])*BS[0] + (BS[18]*F[2]+BS[19]*F[3])*BS[1]);
    Kp[10] += ((BS[20]*F[0]+BS[21]*F[1])*BS[0] + (BS[20]*F[2]+BS[21]*F[3])*BS[1]);
    Kp[11] += ((BS[22]*F[0]+BS[23]*F[1])*BS[0] + (BS[22]*F[2]+BS[23]*F[3])*BS[1]);
    
    Kp[12] += ((BS[0]*F[0]+BS[1]*F[1])*BS[2] + (BS[0]*F[2]+BS[1]*F[3])*BS[3]);
    Kp[13] += ((BS[2]*F[0]+BS[3]*F[1])*BS[2] + (BS[2]*F[2]+BS[3]*F[3])*BS[3]);
    Kp[14] += ((BS[4]*F[0]+BS[5]*F[1])*BS[2] + (BS[4]*F[2]+BS[5]*F[3])*BS[3]);
    Kp[15] += ((BS[6]*F[0]+BS[7]*F[1])*BS[2] + (BS[6]*F[2]+BS[7]*F[3])*BS[3]);
    Kp[16] += ((BS[8]*F[0]+BS[9]*F[1])*BS[2] + (BS[8]*F[2]+BS[9]*F[3])*BS[3]);
    Kp[17] += ((BS[10]*F[0]+BS[11]*F[1])*BS[2] + (BS[10]*F[2]+BS[11]*F[3])*BS[3]);
    Kp[18] += ((BS[12]*F[0]+BS[13]*F[1])*BS[2] + (BS[12]*F[2]+BS[13]*F[3])*BS[3]);
    Kp[19] += ((BS[14]*F[0]+BS[15]*F[1])*BS[2] + (BS[14]*F[2]+BS[15]*F[3])*BS[3]);
    Kp[20] += ((BS[16]*F[0]+BS[17]*F[1])*BS[2] + (BS[16]*F[2]+BS[17]*F[3])*BS[3]);
    Kp[21] += ((BS[18]*F[0]+BS[19]*F[1])*BS[2] + (BS[18]*F[2]+BS[19]*F[3])*BS[3]);
    Kp[22] += ((BS[20]*F[0]+BS[21]*F[1])*BS[2] + (BS[20]*F[2]+BS[21]*F[3])*BS[3]);
    Kp[23] += ((BS[22]*F[0]+BS[23]*F[1])*BS[2] + (BS[22]*F[2]+BS[23]*F[3])*BS[3]);
    
    Kp[24] += ((BS[0]*F[0]+BS[1]*F[1])*BS[4] + (BS[0]*F[2]+BS[1]*F[3])*BS[5]);
    Kp[25] += ((BS[2]*F[0]+BS[3]*F[1])*BS[4] + (BS[2]*F[2]+BS[3]*F[3])*BS[5]);
    Kp[26] += ((BS[4]*F[0]+BS[5]*F[1])*BS[4] + (BS[4]*F[2]+BS[5]*F[3])*BS[5]);
    Kp[27] += ((BS[6]*F[0]+BS[7]*F[1])*BS[4] + (BS[6]*F[2]+BS[7]*F[3])*BS[5]);
    Kp[28] += ((BS[8]*F[0]+BS[9]*F[1])*BS[4] + (BS[8]*F[2]+BS[9]*F[3])*BS[5]);
    Kp[29] += ((BS[10]*F[0]+BS[11]*F[1])*BS[4] + (BS[10]*F[2]+BS[11]*F[3])*BS[5]);
    Kp[30] += ((BS[12]*F[0]+BS[13]*F[1])*BS[4] + (BS[12]*F[2]+BS[13]*F[3])*BS[5]);
    Kp[31] += ((BS[14]*F[0]+BS[15]*F[1])*BS[4] + (BS[14]*F[2]+BS[15]*F[3])*BS[5]);
    
    Kp[32] += ((BS[16]*F[0]+BS[17]*F[1])*BS[4] + (BS[16]*F[2]+BS[17]*F[3])*BS[5]);
    Kp[33] += ((BS[18]*F[0]+BS[19]*F[1])*BS[4] + (BS[18]*F[2]+BS[19]*F[3])*BS[5]);
    Kp[34] += ((BS[20]*F[0]+BS[21]*F[1])*BS[4] + (BS[20]*F[2]+BS[21]*F[3])*BS[5]);
    Kp[35] += ((BS[22]*F[0]+BS[23]*F[1])*BS[4] + (BS[22]*F[2]+BS[23]*F[3])*BS[5]);
    
    Kp[36] += ((BS[0]*F[0]+BS[1]*F[1])*BS[6] + (BS[0]*F[2]+BS[1]*F[3])*BS[7]);
    Kp[37] += ((BS[2]*F[0]+BS[3]*F[1])*BS[6] + (BS[2]*F[2]+BS[3]*F[3])*BS[7]);
    Kp[38] += ((BS[4]*F[0]+BS[5]*F[1])*BS[6] + (BS[4]*F[2]+BS[5]*F[3])*BS[7]);
    Kp[39] += ((BS[6]*F[0]+BS[7]*F[1])*BS[6] + (BS[6]*F[2]+BS[7]*F[3])*BS[7]);
    Kp[40] += ((BS[8]*F[0]+BS[9]*F[1])*BS[6] + (BS[8]*F[2]+BS[9]*F[3])*BS[7]);
    Kp[41] += ((BS[10]*F[0]+BS[11]*F[1])*BS[6] + (BS[10]*F[2]+BS[11]*F[3])*BS[7]);
    Kp[42] += ((BS[12]*F[0]+BS[13]*F[1])*BS[6] + (BS[12]*F[2]+BS[13]*F[3])*BS[7]);
    Kp[43] += ((BS[14]*F[0]+BS[15]*F[1])*BS[6] + (BS[14]*F[2]+BS[15]*F[3])*BS[7]);
    Kp[44] += ((BS[16]*F[0]+BS[17]*F[1])*BS[6] + (BS[16]*F[2]+BS[17]*F[3])*BS[7]);
    Kp[45] += ((BS[18]*F[0]+BS[19]*F[1])*BS[6] + (BS[18]*F[2]+BS[19]*F[3])*BS[7]);
    Kp[46] += ((BS[20]*F[0]+BS[21]*F[1])*BS[6] + (BS[20]*F[2]+BS[21]*F[3])*BS[7]);
    Kp[47] += ((BS[22]*F[0]+BS[23]*F[1])*BS[6] + (BS[22]*F[2]+BS[23]*F[3])*BS[7]);
    
    Kp[48] += ((BS[0]*F[0]+BS[1]*F[1])*BS[8] + (BS[0]*F[2]+BS[1]*F[3])*BS[9]);
    Kp[49] += ((BS[2]*F[0]+BS[3]*F[1])*BS[8] + (BS[2]*F[2]+BS[3]*F[3])*BS[9]);
    Kp[50] += ((BS[4]*F[0]+BS[5]*F[1])*BS[8] + (BS[4]*F[2]+BS[5]*F[3])*BS[9]);
    Kp[51] += ((BS[6]*F[0]+BS[7]*F[1])*BS[8] + (BS[6]*F[2]+BS[7]*F[3])*BS[9]);
    Kp[52] += ((BS[8]*F[0]+BS[9]*F[1])*BS[8] + (BS[8]*F[2]+BS[9]*F[3])*BS[9]);
    Kp[53] += ((BS[10]*F[0]+BS[11]*F[1])*BS[8] + (BS[10]*F[2]+BS[11]*F[3])*BS[9]);
    Kp[54] += ((BS[12]*F[0]+BS[13]*F[1])*BS[8] + (BS[12]*F[2]+BS[13]*F[3])*BS[9]);
    Kp[55] += ((BS[14]*F[0]+BS[15]*F[1])*BS[8] + (BS[14]*F[2]+BS[15]*F[3])*BS[9]);
    Kp[56] += ((BS[16]*F[0]+BS[17]*F[1])*BS[8] + (BS[16]*F[2]+BS[17]*F[3])*BS[9]);
    Kp[57] += ((BS[18]*F[0]+BS[19]*F[1])*BS[8] + (BS[18]*F[2]+BS[19]*F[3])*BS[9]);
    Kp[58] += ((BS[20]*F[0]+BS[21]*F[1])*BS[8] + (BS[20]*F[2]+BS[21]*F[3])*BS[9]);
    Kp[59] += ((BS[22]*F[0]+BS[23]*F[1])*BS[8] + (BS[22]*F[2]+BS[23]*F[3])*BS[9]);
    
    Kp[60] += ((BS[0]*F[0]+BS[1]*F[1])*BS[10] + (BS[0]*F[2]+BS[1]*F[3])*BS[11]);
    Kp[61] += ((BS[2]*F[0]+BS[3]*F[1])*BS[10] + (BS[2]*F[2]+BS[3]*F[3])*BS[11]);
    Kp[62] += ((BS[4]*F[0]+BS[5]*F[1])*BS[10] + (BS[4]*F[2]+BS[5]*F[3])*BS[11]);
    Kp[63] += ((BS[6]*F[0]+BS[7]*F[1])*BS[10] + (BS[6]*F[2]+BS[7]*F[3])*BS[11]);
    Kp[64] += ((BS[8]*F[0]+BS[9]*F[1])*BS[10] + (BS[8]*F[2]+BS[9]*F[3])*BS[11]);
    Kp[65] += ((BS[10]*F[0]+BS[11]*F[1])*BS[10] + (BS[10]*F[2]+BS[11]*F[3])*BS[11]);
    Kp[66] += ((BS[12]*F[0]+BS[13]*F[1])*BS[10] + (BS[12]*F[2]+BS[13]*F[3])*BS[11]);
    Kp[67] += ((BS[14]*F[0]+BS[15]*F[1])*BS[10] + (BS[14]*F[2]+BS[15]*F[3])*BS[11]);
    Kp[68] += ((BS[16]*F[0]+BS[17]*F[1])*BS[10] + (BS[16]*F[2]+BS[17]*F[3])*BS[11]);
    Kp[69] += ((BS[18]*F[0]+BS[19]*F[1])*BS[10] + (BS[18]*F[2]+BS[19]*F[3])*BS[11]);
    Kp[70] += ((BS[20]*F[0]+BS[21]*F[1])*BS[10] + (BS[20]*F[2]+BS[21]*F[3])*BS[11]);
    
    Kp[71] += ((BS[22]*F[0]+BS[23]*F[1])*BS[10] + (BS[22]*F[2]+BS[23]*F[3])*BS[11]);
    
    Kp[72] += ((BS[0]*F[0]+BS[1]*F[1])*BS[12] + (BS[0]*F[2]+BS[1]*F[3])*BS[13]);
    Kp[73] += ((BS[2]*F[0]+BS[3]*F[1])*BS[12] + (BS[2]*F[2]+BS[3]*F[3])*BS[13]);
    Kp[74] += ((BS[4]*F[0]+BS[5]*F[1])*BS[12] + (BS[4]*F[2]+BS[5]*F[3])*BS[13]);
    Kp[75] += ((BS[6]*F[0]+BS[7]*F[1])*BS[12] + (BS[6]*F[2]+BS[7]*F[3])*BS[13]);
    Kp[76] += ((BS[8]*F[0]+BS[9]*F[1])*BS[12] + (BS[8]*F[2]+BS[9]*F[3])*BS[13]);
    Kp[77] += ((BS[10]*F[0]+BS[11]*F[1])*BS[12] + (BS[10]*F[2]+BS[11]*F[3])*BS[13]);
    Kp[78] += ((BS[12]*F[0]+BS[13]*F[1])*BS[12] + (BS[12]*F[2]+BS[13]*F[3])*BS[13]);
    Kp[79] += ((BS[14]*F[0]+BS[15]*F[1])*BS[12] + (BS[14]*F[2]+BS[15]*F[3])*BS[13]);
    Kp[80] += ((BS[16]*F[0]+BS[17]*F[1])*BS[12] + (BS[16]*F[2]+BS[17]*F[3])*BS[13]);
    Kp[81] += ((BS[18]*F[0]+BS[19]*F[1])*BS[12] + (BS[18]*F[2]+BS[19]*F[3])*BS[13]);
    Kp[82] += ((BS[20]*F[0]+BS[21]*F[1])*BS[12] + (BS[20]*F[2]+BS[21]*F[3])*BS[13]);
    Kp[83] += ((BS[22]*F[0]+BS[23]*F[1])*BS[12] + (BS[22]*F[2]+BS[23]*F[3])*BS[13]);
    
    Kp[84] += ((BS[0]*F[0]+BS[1]*F[1])*BS[14] + (BS[0]*F[2]+BS[1]*F[3])*BS[15]);
    Kp[85] += ((BS[2]*F[0]+BS[3]*F[1])*BS[14] + (BS[2]*F[2]+BS[3]*F[3])*BS[15]);
    Kp[86] += ((BS[4]*F[0]+BS[5]*F[1])*BS[14] + (BS[4]*F[2]+BS[5]*F[3])*BS[15]);
    Kp[87] += ((BS[6]*F[0]+BS[7]*F[1])*BS[14] + (BS[6]*F[2]+BS[7]*F[3])*BS[15]);
    Kp[88] += ((BS[8]*F[0]+BS[9]*F[1])*BS[14] + (BS[8]*F[2]+BS[9]*F[3])*BS[15]);
    Kp[89] += ((BS[10]*F[0]+BS[11]*F[1])*BS[14] + (BS[10]*F[2]+BS[11]*F[3])*BS[15]);
    Kp[90] += ((BS[12]*F[0]+BS[13]*F[1])*BS[14] + (BS[12]*F[2]+BS[13]*F[3])*BS[15]);
    Kp[91] += ((BS[14]*F[0]+BS[15]*F[1])*BS[14] + (BS[14]*F[2]+BS[15]*F[3])*BS[15]);
    Kp[92] += ((BS[16]*F[0]+BS[17]*F[1])*BS[14] + (BS[16]*F[2]+BS[17]*F[3])*BS[15]);
    Kp[93] += ((BS[18]*F[0]+BS[19]*F[1])*BS[14] + (BS[18]*F[2]+BS[19]*F[3])*BS[15]);
    Kp[94] += ((BS[20]*F[0]+BS[21]*F[1])*BS[14] + (BS[20]*F[2]+BS[21]*F[3])*BS[15]);
    Kp[95] += ((BS[22]*F[0]+BS[23]*F[1])*BS[14] + (BS[22]*F[2]+BS[23]*F[3])*BS[15]);
    
    Kp[96] += ((BS[0]*F[0]+BS[1]*F[1])*BS[16] + (BS[0]*F[2]+BS[1]*F[3])*BS[17]);
    Kp[97] += ((BS[2]*F[0]+BS[3]*F[1])*BS[16] + (BS[2]*F[2]+BS[3]*F[3])*BS[17]);
    Kp[98] += ((BS[4]*F[0]+BS[5]*F[1])*BS[16] + (BS[4]*F[2]+BS[5]*F[3])*BS[17]);
    Kp[99] += ((BS[6]*F[0]+BS[7]*F[1])*BS[16] + (BS[6]*F[2]+BS[7]*F[3])*BS[17]);
    Kp[100] += ((BS[8]*F[0]+BS[9]*F[1])*BS[16] + (BS[8]*F[2]+BS[9]*F[3])*BS[17]);
    Kp[101] += ((BS[10]*F[0]+BS[11]*F[1])*BS[16] + (BS[10]*F[2]+BS[11]*F[3])*BS[17]);
    Kp[102] += ((BS[12]*F[0]+BS[13]*F[1])*BS[16] + (BS[12]*F[2]+BS[13]*F[3])*BS[17]);
    Kp[103] += ((BS[14]*F[0]+BS[15]*F[1])*BS[16] + (BS[14]*F[2]+BS[15]*F[3])*BS[17]);
    Kp[104] += ((BS[16]*F[0]+BS[17]*F[1])*BS[16] + (BS[16]*F[2]+BS[17]*F[3])*BS[17]);
    Kp[105] += ((BS[18]*F[0]+BS[19]*F[1])*BS[16] + (BS[18]*F[2]+BS[19]*F[3])*BS[17]);
    Kp[106] += ((BS[20]*F[0]+BS[21]*F[1])*BS[16] + (BS[20]*F[2]+BS[21]*F[3])*BS[17]);
    Kp[107] += ((BS[22]*F[0]+BS[23]*F[1])*BS[16] + (BS[22]*F[2]+BS[23]*F[3])*BS[17]);
    
    Kp[108] += ((BS[0]*F[0]+BS[1]*F[1])*BS[18] + (BS[0]*F[2]+BS[1]*F[3])*BS[19]);
    Kp[109] += ((BS[2]*F[0]+BS[3]*F[1])*BS[18] + (BS[2]*F[2]+BS[3]*F[3])*BS[19]);
    Kp[110] += ((BS[4]*F[0]+BS[5]*F[1])*BS[18] + (BS[4]*F[2]+BS[5]*F[3])*BS[19]);
    Kp[111] += ((BS[6]*F[0]+BS[7]*F[1])*BS[18] + (BS[6]*F[2]+BS[7]*F[3])*BS[19]);
    Kp[112] += ((BS[8]*F[0]+BS[9]*F[1])*BS[18] + (BS[8]*F[2]+BS[9]*F[3])*BS[19]);
    Kp[113] += ((BS[10]*F[0]+BS[11]*F[1])*BS[18] + (BS[10]*F[2]+BS[11]*F[3])*BS[19]);
    Kp[114] += ((BS[12]*F[0]+BS[13]*F[1])*BS[18] + (BS[12]*F[2]+BS[13]*F[3])*BS[19]);
    Kp[115] += ((BS[14]*F[0]+BS[15]*F[1])*BS[18] + (BS[14]*F[2]+BS[15]*F[3])*BS[19]);
    Kp[116] += ((BS[16]*F[0]+BS[17]*F[1])*BS[18] + (BS[16]*F[2]+BS[17]*F[3])*BS[19]);
    Kp[117] += ((BS[18]*F[0]+BS[19]*F[1])*BS[18] + (BS[18]*F[2]+BS[19]*F[3])*BS[19]);
    Kp[118] += ((BS[20]*F[0]+BS[21]*F[1])*BS[18] + (BS[20]*F[2]+BS[21]*F[3])*BS[19]);
    Kp[119] += ((BS[22]*F[0]+BS[23]*F[1])*BS[18] + (BS[22]*F[2]+BS[23]*F[3])*BS[19]);
    
    Kp[120] += ((BS[0]*F[0]+BS[1]*F[1])*BS[20] + (BS[0]*F[2]+BS[1]*F[3])*BS[21]);
    Kp[121] += ((BS[2]*F[0]+BS[3]*F[1])*BS[20] + (BS[2]*F[2]+BS[3]*F[3])*BS[21]);
    Kp[122] += ((BS[4]*F[0]+BS[5]*F[1])*BS[20] + (BS[4]*F[2]+BS[5]*F[3])*BS[21]);
    Kp[123] += ((BS[6]*F[0]+BS[7]*F[1])*BS[20] + (BS[6]*F[2]+BS[7]*F[3])*BS[21]);
    Kp[124] += ((BS[8]*F[0]+BS[9]*F[1])*BS[20] + (BS[8]*F[2]+BS[9]*F[3])*BS[21]);
    Kp[125] += ((BS[10]*F[0]+BS[11]*F[1])*BS[20] + (BS[10]*F[2]+BS[11]*F[3])*BS[21]);
    Kp[126] += ((BS[12]*F[0]+BS[13]*F[1])*BS[20] + (BS[12]*F[2]+BS[13]*F[3])*BS[21]);
    Kp[127] += ((BS[14]*F[0]+BS[15]*F[1])*BS[20] + (BS[14]*F[2]+BS[15]*F[3])*BS[21]);
    Kp[128] += ((BS[16]*F[0]+BS[17]*F[1])*BS[20] + (BS[16]*F[2]+BS[17]*F[3])*BS[21]);
    Kp[129] += ((BS[18]*F[0]+BS[19]*F[1])*BS[20] + (BS[18]*F[2]+BS[19]*F[3])*BS[21]);
    Kp[130] += ((BS[20]*F[0]+BS[21]*F[1])*BS[20] + (BS[20]*F[2]+BS[21]*F[3])*BS[21]);
    Kp[131] += ((BS[22]*F[0]+BS[23]*F[1])*BS[20] + (BS[22]*F[2]+BS[23]*F[3])*BS[21]);
    
    Kp[132] += ((BS[0]*F[0]+BS[1]*F[1])*BS[22] + (BS[0]*F[2]+BS[1]*F[3])*BS[23]);
    Kp[133] += ((BS[2]*F[0]+BS[3]*F[1])*BS[22] + (BS[2]*F[2]+BS[3]*F[3])*BS[23]);
    Kp[134] += ((BS[4]*F[0]+BS[5]*F[1])*BS[22] + (BS[4]*F[2]+BS[5]*F[3])*BS[23]);
    Kp[135] += ((BS[6]*F[0]+BS[7]*F[1])*BS[22] + (BS[6]*F[2]+BS[7]*F[3])*BS[23]);
    Kp[136] += ((BS[8]*F[0]+BS[9]*F[1])*BS[22] + (BS[8]*F[2]+BS[9]*F[3])*BS[23]);
    Kp[137] += ((BS[10]*F[0]+BS[11]*F[1])*BS[22] + (BS[10]*F[2]+BS[11]*F[3])*BS[23]);
    Kp[138] += ((BS[12]*F[0]+BS[13]*F[1])*BS[22] + (BS[12]*F[2]+BS[13]*F[3])*BS[23]);
    Kp[139] += ((BS[14]*F[0]+BS[15]*F[1])*BS[22] + (BS[14]*F[2]+BS[15]*F[3])*BS[23]);
    Kp[140] += ((BS[16]*F[0]+BS[17]*F[1])*BS[22] + (BS[16]*F[2]+BS[17]*F[3])*BS[23]);
    Kp[141] += ((BS[18]*F[0]+BS[19]*F[1])*BS[22] + (BS[18]*F[2]+BS[19]*F[3])*BS[23]);
    Kp[142] += ((BS[20]*F[0]+BS[21]*F[1])*BS[22] + (BS[20]*F[2]+BS[21]*F[3])*BS[23]);
    Kp[143] += ((BS[22]*F[0]+BS[23]*F[1])*BS[22] + (BS[22]*F[2]+BS[23]*F[3])*BS[23]);
 */

}
free(BB);
free(BM);
free(A);
free(D);
free(B);
    
   // free(F);
    //free(BS);
/*Shear terms go through single guass intergration only*/

    nx = 0;
ny = 0;

n1j = 0.25*(1-nx)*(1-ny);
n2j = 0.25*(1+nx)*(1-ny);
n3j = 0.25*(1+nx)*(1+ny);
n4j = 0.25*(1-nx)*(1+ny);

dn1dx = -0.5*(1-ny)/hxz;
dn2dx = 0.5*(1-ny)/hxz;
dn3dx = 0.5*(1+ny)/hxz;
dn4dx = -0.5*(1+ny)/hxz;

dn1dy = -0.5*(1-nx)/h;
dn2dy = -0.5*(1+nx)/h;
dn3dy = 0.5*(1+nx)/h;
dn4dy = 0.5*(1-nx)/h;


BS[0] = dn1dx;
BS[2] = 0;
BS[4] = n1j;
BS[6] = dn2dx;
BS[8] = 0;
BS[10] = n2j;
BS[12] = dn3dx;
BS[14] = 0;
BS[16] = n3j;
BS[18] = dn4dx;
BS[20] = 0;
BS[22] = n4j;

BS[1] = dn1dy;
BS[3] = -n1j;
BS[5] = 0;
BS[7] = dn2dy;
BS[9] = -n2j;
BS[11] = 0;
BS[13] = dn3dy;
BS[15] = -n3j;
BS[17] = 0;
BS[19] = dn4dy;
BS[21] = -n4j;
BS[23] = 0;
 
 
 Kp[0] += ((BS[0]*F[0]+BS[1]*F[1])*BS[0] + (BS[0]*F[2]+BS[1]*F[3])*BS[1]);
 Kp[1] += ((BS[2]*F[0]+BS[3]*F[1])*BS[0] + (BS[2]*F[2]+BS[3]*F[3])*BS[1]);
 Kp[2] += ((BS[4]*F[0]+BS[5]*F[1])*BS[0] + (BS[4]*F[2]+BS[5]*F[3])*BS[1]);
 Kp[3] += ((BS[6]*F[0]+BS[7]*F[1])*BS[0] + (BS[6]*F[2]+BS[7]*F[3])*BS[1]);
 Kp[4] += ((BS[8]*F[0]+BS[9]*F[1])*BS[0] + (BS[8]*F[2]+BS[9]*F[3])*BS[1]);
 Kp[5] += ((BS[10]*F[0]+BS[11]*F[1])*BS[0] + (BS[10]*F[2]+BS[11]*F[3])*BS[1]);
 Kp[6] += ((BS[12]*F[0]+BS[13]*F[1])*BS[0] + (BS[12]*F[2]+BS[13]*F[3])*BS[1]);
 Kp[7] += ((BS[14]*F[0]+BS[15]*F[1])*BS[0] + (BS[14]*F[2]+BS[15]*F[3])*BS[1]);
 Kp[8] += ((BS[16]*F[0]+BS[17]*F[1])*BS[0] + (BS[16]*F[2]+BS[17]*F[3])*BS[1]);
 Kp[9] += ((BS[18]*F[0]+BS[19]*F[1])*BS[0] + (BS[18]*F[2]+BS[19]*F[3])*BS[1]);
 Kp[10] += ((BS[20]*F[0]+BS[21]*F[1])*BS[0] + (BS[20]*F[2]+BS[21]*F[3])*BS[1]);
 Kp[11] += ((BS[22]*F[0]+BS[23]*F[1])*BS[0] + (BS[22]*F[2]+BS[23]*F[3])*BS[1]);
 
 Kp[12] += ((BS[0]*F[0]+BS[1]*F[1])*BS[2] + (BS[0]*F[2]+BS[1]*F[3])*BS[3]);
 Kp[13] += ((BS[2]*F[0]+BS[3]*F[1])*BS[2] + (BS[2]*F[2]+BS[3]*F[3])*BS[3]);
 Kp[14] += ((BS[4]*F[0]+BS[5]*F[1])*BS[2] + (BS[4]*F[2]+BS[5]*F[3])*BS[3]);
 Kp[15] += ((BS[6]*F[0]+BS[7]*F[1])*BS[2] + (BS[6]*F[2]+BS[7]*F[3])*BS[3]);
 Kp[16] += ((BS[8]*F[0]+BS[9]*F[1])*BS[2] + (BS[8]*F[2]+BS[9]*F[3])*BS[3]);
 Kp[17] += ((BS[10]*F[0]+BS[11]*F[1])*BS[2] + (BS[10]*F[2]+BS[11]*F[3])*BS[3]);
 Kp[18] += ((BS[12]*F[0]+BS[13]*F[1])*BS[2] + (BS[12]*F[2]+BS[13]*F[3])*BS[3]);
 Kp[19] += ((BS[14]*F[0]+BS[15]*F[1])*BS[2] + (BS[14]*F[2]+BS[15]*F[3])*BS[3]);
 Kp[20] += ((BS[16]*F[0]+BS[17]*F[1])*BS[2] + (BS[16]*F[2]+BS[17]*F[3])*BS[3]);
 Kp[21] += ((BS[18]*F[0]+BS[19]*F[1])*BS[2] + (BS[18]*F[2]+BS[19]*F[3])*BS[3]);
 Kp[22] += ((BS[20]*F[0]+BS[21]*F[1])*BS[2] + (BS[20]*F[2]+BS[21]*F[3])*BS[3]);
 Kp[23] += ((BS[22]*F[0]+BS[23]*F[1])*BS[2] + (BS[22]*F[2]+BS[23]*F[3])*BS[3]);
 
 Kp[24] += ((BS[0]*F[0]+BS[1]*F[1])*BS[4] + (BS[0]*F[2]+BS[1]*F[3])*BS[5]);
 Kp[25] += ((BS[2]*F[0]+BS[3]*F[1])*BS[4] + (BS[2]*F[2]+BS[3]*F[3])*BS[5]);
 Kp[26] += ((BS[4]*F[0]+BS[5]*F[1])*BS[4] + (BS[4]*F[2]+BS[5]*F[3])*BS[5]);
 Kp[27] += ((BS[6]*F[0]+BS[7]*F[1])*BS[4] + (BS[6]*F[2]+BS[7]*F[3])*BS[5]);
 Kp[28] += ((BS[8]*F[0]+BS[9]*F[1])*BS[4] + (BS[8]*F[2]+BS[9]*F[3])*BS[5]);
 Kp[29] += ((BS[10]*F[0]+BS[11]*F[1])*BS[4] + (BS[10]*F[2]+BS[11]*F[3])*BS[5]);
 Kp[30] += ((BS[12]*F[0]+BS[13]*F[1])*BS[4] + (BS[12]*F[2]+BS[13]*F[3])*BS[5]);
 Kp[31] += ((BS[14]*F[0]+BS[15]*F[1])*BS[4] + (BS[14]*F[2]+BS[15]*F[3])*BS[5]);
 
 Kp[32] += ((BS[16]*F[0]+BS[17]*F[1])*BS[4] + (BS[16]*F[2]+BS[17]*F[3])*BS[5]);
 Kp[33] += ((BS[18]*F[0]+BS[19]*F[1])*BS[4] + (BS[18]*F[2]+BS[19]*F[3])*BS[5]);
 Kp[34] += ((BS[20]*F[0]+BS[21]*F[1])*BS[4] + (BS[20]*F[2]+BS[21]*F[3])*BS[5]);
 Kp[35] += ((BS[22]*F[0]+BS[23]*F[1])*BS[4] + (BS[22]*F[2]+BS[23]*F[3])*BS[5]);
 
 Kp[36] += ((BS[0]*F[0]+BS[1]*F[1])*BS[6] + (BS[0]*F[2]+BS[1]*F[3])*BS[7]);
 Kp[37] += ((BS[2]*F[0]+BS[3]*F[1])*BS[6] + (BS[2]*F[2]+BS[3]*F[3])*BS[7]);
 Kp[38] += ((BS[4]*F[0]+BS[5]*F[1])*BS[6] + (BS[4]*F[2]+BS[5]*F[3])*BS[7]);
 Kp[39] += ((BS[6]*F[0]+BS[7]*F[1])*BS[6] + (BS[6]*F[2]+BS[7]*F[3])*BS[7]);
 Kp[40] += ((BS[8]*F[0]+BS[9]*F[1])*BS[6] + (BS[8]*F[2]+BS[9]*F[3])*BS[7]);
 Kp[41] += ((BS[10]*F[0]+BS[11]*F[1])*BS[6] + (BS[10]*F[2]+BS[11]*F[3])*BS[7]);
 Kp[42] += ((BS[12]*F[0]+BS[13]*F[1])*BS[6] + (BS[12]*F[2]+BS[13]*F[3])*BS[7]);
 Kp[43] += ((BS[14]*F[0]+BS[15]*F[1])*BS[6] + (BS[14]*F[2]+BS[15]*F[3])*BS[7]);
 Kp[44] += ((BS[16]*F[0]+BS[17]*F[1])*BS[6] + (BS[16]*F[2]+BS[17]*F[3])*BS[7]);
 Kp[45] += ((BS[18]*F[0]+BS[19]*F[1])*BS[6] + (BS[18]*F[2]+BS[19]*F[3])*BS[7]);
 Kp[46] += ((BS[20]*F[0]+BS[21]*F[1])*BS[6] + (BS[20]*F[2]+BS[21]*F[3])*BS[7]);
 Kp[47] += ((BS[22]*F[0]+BS[23]*F[1])*BS[6] + (BS[22]*F[2]+BS[23]*F[3])*BS[7]);
 
 Kp[48] += ((BS[0]*F[0]+BS[1]*F[1])*BS[8] + (BS[0]*F[2]+BS[1]*F[3])*BS[9]);
 Kp[49] += ((BS[2]*F[0]+BS[3]*F[1])*BS[8] + (BS[2]*F[2]+BS[3]*F[3])*BS[9]);
 Kp[50] += ((BS[4]*F[0]+BS[5]*F[1])*BS[8] + (BS[4]*F[2]+BS[5]*F[3])*BS[9]);
 Kp[51] += ((BS[6]*F[0]+BS[7]*F[1])*BS[8] + (BS[6]*F[2]+BS[7]*F[3])*BS[9]);
 Kp[52] += ((BS[8]*F[0]+BS[9]*F[1])*BS[8] + (BS[8]*F[2]+BS[9]*F[3])*BS[9]);
 Kp[53] += ((BS[10]*F[0]+BS[11]*F[1])*BS[8] + (BS[10]*F[2]+BS[11]*F[3])*BS[9]);
 Kp[54] += ((BS[12]*F[0]+BS[13]*F[1])*BS[8] + (BS[12]*F[2]+BS[13]*F[3])*BS[9]);
 Kp[55] += ((BS[14]*F[0]+BS[15]*F[1])*BS[8] + (BS[14]*F[2]+BS[15]*F[3])*BS[9]);
 Kp[56] += ((BS[16]*F[0]+BS[17]*F[1])*BS[8] + (BS[16]*F[2]+BS[17]*F[3])*BS[9]);
 Kp[57] += ((BS[18]*F[0]+BS[19]*F[1])*BS[8] + (BS[18]*F[2]+BS[19]*F[3])*BS[9]);
 Kp[58] += ((BS[20]*F[0]+BS[21]*F[1])*BS[8] + (BS[20]*F[2]+BS[21]*F[3])*BS[9]);
 Kp[59] += ((BS[22]*F[0]+BS[23]*F[1])*BS[8] + (BS[22]*F[2]+BS[23]*F[3])*BS[9]);
 
 Kp[60] += ((BS[0]*F[0]+BS[1]*F[1])*BS[10] + (BS[0]*F[2]+BS[1]*F[3])*BS[11]);
 Kp[61] += ((BS[2]*F[0]+BS[3]*F[1])*BS[10] + (BS[2]*F[2]+BS[3]*F[3])*BS[11]);
 Kp[62] += ((BS[4]*F[0]+BS[5]*F[1])*BS[10] + (BS[4]*F[2]+BS[5]*F[3])*BS[11]);
 Kp[63] += ((BS[6]*F[0]+BS[7]*F[1])*BS[10] + (BS[6]*F[2]+BS[7]*F[3])*BS[11]);
 Kp[64] += ((BS[8]*F[0]+BS[9]*F[1])*BS[10] + (BS[8]*F[2]+BS[9]*F[3])*BS[11]);
 Kp[65] += ((BS[10]*F[0]+BS[11]*F[1])*BS[10] + (BS[10]*F[2]+BS[11]*F[3])*BS[11]);
 Kp[66] += ((BS[12]*F[0]+BS[13]*F[1])*BS[10] + (BS[12]*F[2]+BS[13]*F[3])*BS[11]);
 Kp[67] += ((BS[14]*F[0]+BS[15]*F[1])*BS[10] + (BS[14]*F[2]+BS[15]*F[3])*BS[11]);
 Kp[68] += ((BS[16]*F[0]+BS[17]*F[1])*BS[10] + (BS[16]*F[2]+BS[17]*F[3])*BS[11]);
 Kp[69] += ((BS[18]*F[0]+BS[19]*F[1])*BS[10] + (BS[18]*F[2]+BS[19]*F[3])*BS[11]);
 Kp[70] += ((BS[20]*F[0]+BS[21]*F[1])*BS[10] + (BS[20]*F[2]+BS[21]*F[3])*BS[11]);
 
 Kp[71] += ((BS[22]*F[0]+BS[23]*F[1])*BS[10] + (BS[22]*F[2]+BS[23]*F[3])*BS[11]);
 
 Kp[72] += ((BS[0]*F[0]+BS[1]*F[1])*BS[12] + (BS[0]*F[2]+BS[1]*F[3])*BS[13]);
 Kp[73] += ((BS[2]*F[0]+BS[3]*F[1])*BS[12] + (BS[2]*F[2]+BS[3]*F[3])*BS[13]);
 Kp[74] += ((BS[4]*F[0]+BS[5]*F[1])*BS[12] + (BS[4]*F[2]+BS[5]*F[3])*BS[13]);
 Kp[75] += ((BS[6]*F[0]+BS[7]*F[1])*BS[12] + (BS[6]*F[2]+BS[7]*F[3])*BS[13]);
 Kp[76] += ((BS[8]*F[0]+BS[9]*F[1])*BS[12] + (BS[8]*F[2]+BS[9]*F[3])*BS[13]);
 Kp[77] += ((BS[10]*F[0]+BS[11]*F[1])*BS[12] + (BS[10]*F[2]+BS[11]*F[3])*BS[13]);
 Kp[78] += ((BS[12]*F[0]+BS[13]*F[1])*BS[12] + (BS[12]*F[2]+BS[13]*F[3])*BS[13]);
 Kp[79] += ((BS[14]*F[0]+BS[15]*F[1])*BS[12] + (BS[14]*F[2]+BS[15]*F[3])*BS[13]);
 Kp[80] += ((BS[16]*F[0]+BS[17]*F[1])*BS[12] + (BS[16]*F[2]+BS[17]*F[3])*BS[13]);
 Kp[81] += ((BS[18]*F[0]+BS[19]*F[1])*BS[12] + (BS[18]*F[2]+BS[19]*F[3])*BS[13]);
 Kp[82] += ((BS[20]*F[0]+BS[21]*F[1])*BS[12] + (BS[20]*F[2]+BS[21]*F[3])*BS[13]);
 Kp[83] += ((BS[22]*F[0]+BS[23]*F[1])*BS[12] + (BS[22]*F[2]+BS[23]*F[3])*BS[13]);
 
 Kp[84] += ((BS[0]*F[0]+BS[1]*F[1])*BS[14] + (BS[0]*F[2]+BS[1]*F[3])*BS[15]);
 Kp[85] += ((BS[2]*F[0]+BS[3]*F[1])*BS[14] + (BS[2]*F[2]+BS[3]*F[3])*BS[15]);
 Kp[86] += ((BS[4]*F[0]+BS[5]*F[1])*BS[14] + (BS[4]*F[2]+BS[5]*F[3])*BS[15]);
 Kp[87] += ((BS[6]*F[0]+BS[7]*F[1])*BS[14] + (BS[6]*F[2]+BS[7]*F[3])*BS[15]);
 Kp[88] += ((BS[8]*F[0]+BS[9]*F[1])*BS[14] + (BS[8]*F[2]+BS[9]*F[3])*BS[15]);
 Kp[89] += ((BS[10]*F[0]+BS[11]*F[1])*BS[14] + (BS[10]*F[2]+BS[11]*F[3])*BS[15]);
 Kp[90] += ((BS[12]*F[0]+BS[13]*F[1])*BS[14] + (BS[12]*F[2]+BS[13]*F[3])*BS[15]);
 Kp[91] += ((BS[14]*F[0]+BS[15]*F[1])*BS[14] + (BS[14]*F[2]+BS[15]*F[3])*BS[15]);
 Kp[92] += ((BS[16]*F[0]+BS[17]*F[1])*BS[14] + (BS[16]*F[2]+BS[17]*F[3])*BS[15]);
 Kp[93] += ((BS[18]*F[0]+BS[19]*F[1])*BS[14] + (BS[18]*F[2]+BS[19]*F[3])*BS[15]);
 Kp[94] += ((BS[20]*F[0]+BS[21]*F[1])*BS[14] + (BS[20]*F[2]+BS[21]*F[3])*BS[15]);
 Kp[95] += ((BS[22]*F[0]+BS[23]*F[1])*BS[14] + (BS[22]*F[2]+BS[23]*F[3])*BS[15]);
 
 Kp[96] += ((BS[0]*F[0]+BS[1]*F[1])*BS[16] + (BS[0]*F[2]+BS[1]*F[3])*BS[17]);
 Kp[97] += ((BS[2]*F[0]+BS[3]*F[1])*BS[16] + (BS[2]*F[2]+BS[3]*F[3])*BS[17]);
 Kp[98] += ((BS[4]*F[0]+BS[5]*F[1])*BS[16] + (BS[4]*F[2]+BS[5]*F[3])*BS[17]);
 Kp[99] += ((BS[6]*F[0]+BS[7]*F[1])*BS[16] + (BS[6]*F[2]+BS[7]*F[3])*BS[17]);
 Kp[100] += ((BS[8]*F[0]+BS[9]*F[1])*BS[16] + (BS[8]*F[2]+BS[9]*F[3])*BS[17]);
 Kp[101] += ((BS[10]*F[0]+BS[11]*F[1])*BS[16] + (BS[10]*F[2]+BS[11]*F[3])*BS[17]);
 Kp[102] += ((BS[12]*F[0]+BS[13]*F[1])*BS[16] + (BS[12]*F[2]+BS[13]*F[3])*BS[17]);
 Kp[103] += ((BS[14]*F[0]+BS[15]*F[1])*BS[16] + (BS[14]*F[2]+BS[15]*F[3])*BS[17]);
 Kp[104] += ((BS[16]*F[0]+BS[17]*F[1])*BS[16] + (BS[16]*F[2]+BS[17]*F[3])*BS[17]);
 Kp[105] += ((BS[18]*F[0]+BS[19]*F[1])*BS[16] + (BS[18]*F[2]+BS[19]*F[3])*BS[17]);
 Kp[106] += ((BS[20]*F[0]+BS[21]*F[1])*BS[16] + (BS[20]*F[2]+BS[21]*F[3])*BS[17]);
 Kp[107] += ((BS[22]*F[0]+BS[23]*F[1])*BS[16] + (BS[22]*F[2]+BS[23]*F[3])*BS[17]);
 
 Kp[108] += ((BS[0]*F[0]+BS[1]*F[1])*BS[18] + (BS[0]*F[2]+BS[1]*F[3])*BS[19]);
 Kp[109] += ((BS[2]*F[0]+BS[3]*F[1])*BS[18] + (BS[2]*F[2]+BS[3]*F[3])*BS[19]);
 Kp[110] += ((BS[4]*F[0]+BS[5]*F[1])*BS[18] + (BS[4]*F[2]+BS[5]*F[3])*BS[19]);
 Kp[111] += ((BS[6]*F[0]+BS[7]*F[1])*BS[18] + (BS[6]*F[2]+BS[7]*F[3])*BS[19]);
 Kp[112] += ((BS[8]*F[0]+BS[9]*F[1])*BS[18] + (BS[8]*F[2]+BS[9]*F[3])*BS[19]);
 Kp[113] += ((BS[10]*F[0]+BS[11]*F[1])*BS[18] + (BS[10]*F[2]+BS[11]*F[3])*BS[19]);
 Kp[114] += ((BS[12]*F[0]+BS[13]*F[1])*BS[18] + (BS[12]*F[2]+BS[13]*F[3])*BS[19]);
 Kp[115] += ((BS[14]*F[0]+BS[15]*F[1])*BS[18] + (BS[14]*F[2]+BS[15]*F[3])*BS[19]);
 Kp[116] += ((BS[16]*F[0]+BS[17]*F[1])*BS[18] + (BS[16]*F[2]+BS[17]*F[3])*BS[19]);
 Kp[117] += ((BS[18]*F[0]+BS[19]*F[1])*BS[18] + (BS[18]*F[2]+BS[19]*F[3])*BS[19]);
 Kp[118] += ((BS[20]*F[0]+BS[21]*F[1])*BS[18] + (BS[20]*F[2]+BS[21]*F[3])*BS[19]);
 Kp[119] += ((BS[22]*F[0]+BS[23]*F[1])*BS[18] + (BS[22]*F[2]+BS[23]*F[3])*BS[19]);
 
 Kp[120] += ((BS[0]*F[0]+BS[1]*F[1])*BS[20] + (BS[0]*F[2]+BS[1]*F[3])*BS[21]);
 Kp[121] += ((BS[2]*F[0]+BS[3]*F[1])*BS[20] + (BS[2]*F[2]+BS[3]*F[3])*BS[21]);
 Kp[122] += ((BS[4]*F[0]+BS[5]*F[1])*BS[20] + (BS[4]*F[2]+BS[5]*F[3])*BS[21]);
 Kp[123] += ((BS[6]*F[0]+BS[7]*F[1])*BS[20] + (BS[6]*F[2]+BS[7]*F[3])*BS[21]);
 Kp[124] += ((BS[8]*F[0]+BS[9]*F[1])*BS[20] + (BS[8]*F[2]+BS[9]*F[3])*BS[21]);
 Kp[125] += ((BS[10]*F[0]+BS[11]*F[1])*BS[20] + (BS[10]*F[2]+BS[11]*F[3])*BS[21]);
 Kp[126] += ((BS[12]*F[0]+BS[13]*F[1])*BS[20] + (BS[12]*F[2]+BS[13]*F[3])*BS[21]);
 Kp[127] += ((BS[14]*F[0]+BS[15]*F[1])*BS[20] + (BS[14]*F[2]+BS[15]*F[3])*BS[21]);
 Kp[128] += ((BS[16]*F[0]+BS[17]*F[1])*BS[20] + (BS[16]*F[2]+BS[17]*F[3])*BS[21]);
 Kp[129] += ((BS[18]*F[0]+BS[19]*F[1])*BS[20] + (BS[18]*F[2]+BS[19]*F[3])*BS[21]);
 Kp[130] += ((BS[20]*F[0]+BS[21]*F[1])*BS[20] + (BS[20]*F[2]+BS[21]*F[3])*BS[21]);
 Kp[131] += ((BS[22]*F[0]+BS[23]*F[1])*BS[20] + (BS[22]*F[2]+BS[23]*F[3])*BS[21]);
 
 Kp[132] += ((BS[0]*F[0]+BS[1]*F[1])*BS[22] + (BS[0]*F[2]+BS[1]*F[3])*BS[23]);
 Kp[133] += ((BS[2]*F[0]+BS[3]*F[1])*BS[22] + (BS[2]*F[2]+BS[3]*F[3])*BS[23]);
 Kp[134] += ((BS[4]*F[0]+BS[5]*F[1])*BS[22] + (BS[4]*F[2]+BS[5]*F[3])*BS[23]);
 Kp[135] += ((BS[6]*F[0]+BS[7]*F[1])*BS[22] + (BS[6]*F[2]+BS[7]*F[3])*BS[23]);
 Kp[136] += ((BS[8]*F[0]+BS[9]*F[1])*BS[22] + (BS[8]*F[2]+BS[9]*F[3])*BS[23]);
 Kp[137] += ((BS[10]*F[0]+BS[11]*F[1])*BS[22] + (BS[10]*F[2]+BS[11]*F[3])*BS[23]);
 Kp[138] += ((BS[12]*F[0]+BS[13]*F[1])*BS[22] + (BS[12]*F[2]+BS[13]*F[3])*BS[23]);
 Kp[139] += ((BS[14]*F[0]+BS[15]*F[1])*BS[22] + (BS[14]*F[2]+BS[15]*F[3])*BS[23]);
 Kp[140] += ((BS[16]*F[0]+BS[17]*F[1])*BS[22] + (BS[16]*F[2]+BS[17]*F[3])*BS[23]);
 Kp[141] += ((BS[18]*F[0]+BS[19]*F[1])*BS[22] + (BS[18]*F[2]+BS[19]*F[3])*BS[23]);
 Kp[142] += ((BS[20]*F[0]+BS[21]*F[1])*BS[22] + (BS[20]*F[2]+BS[21]*F[3])*BS[23]);
 Kp[143] += ((BS[22]*F[0]+BS[23]*F[1])*BS[22] + (BS[22]*F[2]+BS[23]*F[3])*BS[23]);
 
 
 free(F);
 free(BS);
 
 

/*printf("BSshape\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[0], BS[2], BS[4], BS[6], BS[8], BS[10], BS[12], BS[14], BS[16], BS[18], BS[20], BS[22]);

printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", BS[1], BS[3], BS[5], BS[7], BS[9], BS[11], BS[13], BS[15], BS[17], BS[19], BS[21], BS[23]);*/

/*Shear terms*/

//
    ///
    ///
    


/*if((num == 361)||(num==378))
{

printf("\n Membrain Matrix");
for(p=0;p<8;p++)
{
	printf("\n");

	for(j=0;j<8;j++)
	{
		printf("%0.3f\t", Km[j*8+p]);
	}

}
printf("\n");
}*/

/*Check the reuslts so far*/
/*printf("\n Membrain Matrix");
for(p=0;p<8;p++)
{
	printf("\n");

	for(j=0;j<8;j++)
	{
		printf("%0.3f\t", Km[j*8+p]);
	}

}
printf("\n");
printf("\n Plate Matrix");
for(i=0;i<12;i++)

{
	printf("\n");
	for(j=0;j<12;j++)
	{

		printf("%0.3f\t", Kp[j*12+i]);
	}
}*/

/*Assemble Full Matirx (Try to transpose it so that it is correct in C matrix form by not hugely important since it should be symetrical*/
for(i=0;i<576;i++)
{
	KE[i] = 0.0000000000000000;	/*Intialize every value to zero*/
}

/*Fill in the values*/
for(j=0;j<4;j++)
{
	p =24*6*j;
	KE[p] = Km[0+2*j];
	KE[1+p] = Km[8+2*j];
	KE[6+p] = Km[16+2*j];
	KE[7+p] = Km[24+2*j];
	KE[12+p] = Km[32+2*j];
	KE[13+p] = Km[40+2*j];
	KE[18+p] = Km[48+2*j];
	KE[19+p] = Km[56+2*j];
	KE[24+p] = Km[1+2*j];

	KE[25+p] = Km[9+2*j];
	KE[30+p] = Km[17+2*j];
	KE[31+p] = Km[25+2*j];
	KE[36+p] = Km[33+2*j];
	KE[37+p] = Km[41+2*j];
	KE[42+p] = Km[49+2*j];
	KE[43+p] = Km[57+2*j];

	q =24*6*j+24*2;
	KE[2+q] = Kp[0+3*j];
	KE[3+q] = Kp[12+3*j];
	KE[4+q] = Kp[24+3*j];
	KE[8+q] = Kp[36+3*j];
	KE[9+q] = Kp[48+3*j];
	KE[10+q] = Kp[60+3*j];
	KE[14+q] = Kp[72+3*j];
	KE[15+q] = Kp[84+3*j];
	KE[16+q] = Kp[96+3*j];
	KE[20+q] = Kp[108+3*j];

	KE[21+q] = Kp[120+3*j];
	KE[22+q] = Kp[132+3*j];
	q =24*6*j+24*3;
	KE[2+q] = Kp[1+3*j];
	KE[3+q] = Kp[13+3*j];
	KE[4+q] = Kp[25+3*j];
	KE[8+q] = Kp[37+3*j];
	KE[9+q] = Kp[49+3*j];
	KE[10+q] = Kp[61+3*j];
	KE[14+q] = Kp[73+3*j];
	KE[15+q] = Kp[85+3*j];
	KE[16+q] = Kp[97+3*j];
	KE[20+q] = Kp[109+3*j];
	KE[21+q] = Kp[121+3*j];
	KE[22+q] = Kp[133+3*j];
	q =24*6*j+24*4;
	KE[2+q] = Kp[2+3*j];
	KE[3+q] = Kp[14+3*j];
	KE[4+q] = Kp[26+3*j];
	KE[8+q] = Kp[38+3*j];
	KE[9+q] = Kp[50+3*j];
	KE[10+q] = Kp[62+3*j];
	KE[14+q] = Kp[74+3*j];
	KE[15+q] = Kp[86+3*j];
	KE[16+q] = Kp[98+3*j];
	KE[20+q] = Kp[110+3*j];
	KE[21+q] = Kp[122+3*j];
	KE[22+q] = Kp[134+3*j];

}

/*printf("\n Element Matrix");
for(i=0;i<24;i++)
{
	printf("\n");

	for(j=0;j<24;j++)
	{
		printf("%0.2f  ", KE[j+i*24]);


	}
}*/

free(Km);
free(Kp);
/*printf("\n Temp Membrain Matrix");*/
/*for(p=0;p<12;p++)
{
	printf("\n");

	for(j=0;j<3;j++)
	{
		printf("%0.3f\t", tempKp[j*12+p]);
	}

}*/
/*printf("\nB");
printf("\nB");	*/	

}

/*Function to calculate fiber continuity score*/
double FibreContScore(int NumNodes, Coord NodeCoord[NumNodes],int elemX, int elemY, Elem Number[elemX][elemY], double h, double *hxz, double Range, double *theta, int itt)
{
	FILE *outfile;	/*File varible for output files*/
	char plotname[40];	/*variable to change names of plotting output files*/
	int i,j,n,m,num,num2;
	double Cx,Cy,Cx2,Cy2,count,Score,DeltaT;
	double Xtemp, Ytemp, dtemp, dist;
	double FCS = 0.0;
	double Tcount = elemX*elemY;
	double *FCSE;
	FCSE = malloc(elemX*elemY*sizeof(double));
	
	Range *= 0.01745329252;  /*Define in Radiums*/

	for(i=0;i<elemX;i++)
	{
		for(j=0;j<elemY;j++)
		{
			count = 0.0;
			Score = 0.0;

			num = Number[i][j].n-1;		/*Get current element number and co-ordinates*/
			n = Number[i][j].a-1;
			Cx = NodeCoord[n].xz + (hxz[i]/2);
			Cy = NodeCoord[n].y + (h/2);
			
			/*Now for each negihbouring element find out if change in theta is less then the range*/

			if(i!=0)
			{
				count++;
				num2 = Number[i-1][j].n-1;		/*Get current neighbour element number and co-ordinates*/
				m = Number[i-1][j].a-1;
				Cx2 = NodeCoord[m].xz + (hxz[i-1]/2);
				Cy2 = NodeCoord[m].y + (h/2);

				Xtemp = Cx-Cx2;				/*Get distance between the elements*/
				Ytemp = Cy-Cy2;	
				dtemp = Xtemp*Xtemp + Ytemp*Ytemp;
				dist = sqrt(dtemp);

				DeltaT = theta[num]-theta[num2];
				/*printf("\n1 DeltaT = %f = %f - %f", DeltaT, theta[num], theta[num2]);*/

				if(fabs(DeltaT)<fabs(Range))
				{
					Score += 1.0/*dist*/;
				}

			}

			if(i!=(elemX-1))
			{
				count++;
				num2 = Number[i+1][j].n-1;		/*Get current neighbour element number and co-ordinates*/
				m = Number[i+1][j].a-1;
				Cx2 = NodeCoord[m].xz + (hxz[i+1]/2);
				Cy2 = NodeCoord[m].y + (h/2);

				Xtemp = Cx-Cx2;				/*Get distance between the elements*/
				Ytemp = Cy-Cy2;	
				dtemp = Xtemp*Xtemp + Ytemp*Ytemp;
				dist = sqrt(dtemp);

				DeltaT = theta[num]-theta[num2];
				/*printf("\n2 DeltaT = %f = %f - %f", DeltaT, theta[num], theta[num2]);*/

				if(fabs(DeltaT)<fabs(Range))
				{
					Score += 1.0/*dist*/;
				}

			}

			if(j!=0)
			{
				count++;
				num2 = Number[i][j-1].n-1;		/*Get current neighbour element number and co-ordinates*/
				m = Number[i][j-1].a-1;
				Cx2 = NodeCoord[m].xz + (hxz[i]/2);
				Cy2 = NodeCoord[m].y + (h/2);

				Xtemp = Cx-Cx2;				/*Get distance between the elements*/
				Ytemp = Cy-Cy2;	
				dtemp = Xtemp*Xtemp + Ytemp*Ytemp;
				dist = sqrt(dtemp);

				DeltaT = theta[num]-theta[num2];
				/*printf("\n3 DeltaT = %f = %f - %f", DeltaT, theta[num], theta[num2]);*/

				if(fabs(DeltaT)<fabs(Range))
				{
					Score += 1.0/*dist*/;
				}

			}

			if(j!=(elemY-1))
			{
				count++;
				num2 = Number[i][j+1].n-1;		/*Get current neighbour element number and co-ordinates*/
				m = Number[i][j+1].a-1;
				Cx2 = NodeCoord[m].xz + (hxz[i]/2);
				Cy2 = NodeCoord[m].y + (h/2);

				Xtemp = Cx-Cx2;				/*Get distance between the elements*/
				Ytemp = Cy-Cy2;	
				dtemp = Xtemp*Xtemp + Ytemp*Ytemp;
				dist = sqrt(dtemp);

				DeltaT = theta[num]-theta[num2];
				/*printf("\n4 DeltaT = %f = %f - %f", DeltaT, theta[num], theta[num2]);*/

				if(fabs(DeltaT)<fabs(Range))
				{
					Score += 1.0/*dist*/;
				}

			}

			if((i!=0)&&(j!=0))
			{
				count++;
				num2 = Number[i-1][j-1].n-1;		/*Get current neighbour element number and co-ordinates*/
				m = Number[i-1][j-1].a-1;
				Cx2 = NodeCoord[m].xz + (hxz[i-1]/2);
				Cy2 = NodeCoord[m].y + (h/2);

				Xtemp = Cx-Cx2;				/*Get distance between the elements*/
				Ytemp = Cy-Cy2;	
				dtemp = Xtemp*Xtemp + Ytemp*Ytemp;
				dist = sqrt(dtemp);

				DeltaT = theta[num]-theta[num2];
				/*printf("\n5 DeltaT = %f = %f - %f", DeltaT, theta[num], theta[num2]);*/

				if(fabs(DeltaT)<fabs(Range))
				{
					Score += 1.0/*dist*/;
				}

			}

			if((i!=elemX-1)&&(j!=0))
			{
				count++;
				num2 = Number[i+1][j-1].n-1;		/*Get current neighbour element number and co-ordinates*/
				m = Number[i+1][j-1].a-1;
				Cx2 = NodeCoord[m].xz + (hxz[i+1]/2);
				Cy2 = NodeCoord[m].y + (h/2);

				Xtemp = Cx-Cx2;				/*Get distance between the elements*/
				Ytemp = Cy-Cy2;	
				dtemp = Xtemp*Xtemp + Ytemp*Ytemp;
				dist = sqrt(dtemp);

				DeltaT = theta[num]-theta[num2];
				/*printf("\n6 DeltaT = %f = %f - %f", DeltaT, theta[num], theta[num2]);*/

				if(fabs(DeltaT)<fabs(Range))
				{
					Score += 1.0/*dist*/;
				}

			}

			if((i!=0)&&(j!=elemY-1))
			{
				count++;
				num2 = Number[i-1][j+1].n-1;		/*Get current neighbour element number and co-ordinates*/
				m = Number[i-1][j+1].a-1;
				Cx2 = NodeCoord[m].xz + (hxz[i-1]/2);
				Cy2 = NodeCoord[m].y + (h/2);

				Xtemp = Cx-Cx2;				/*Get distance between the elements*/
				Ytemp = Cy-Cy2;	
				dtemp = Xtemp*Xtemp + Ytemp*Ytemp;
				dist = sqrt(dtemp);

				DeltaT = theta[num]-theta[num2];
				/*printf("\n7 DeltaT = %f = %f - %f", DeltaT, theta[num], theta[num2]);*/

				if(fabs(DeltaT)<fabs(Range))
				{
					Score += 1.0/*dist*/;
				}

			}

			if((i!=elemX-1)&&(j!=elemY-1))
			{
				count++;
				num2 = Number[i+1][j+1].n-1;		/*Get current neighbour element number and co-ordinates*/
				m = Number[i+1][j+1].a-1;
				Cx2 = NodeCoord[m].xz + (hxz[i+1]/2);
				Cy2 = NodeCoord[m].y + (h/2);

				Xtemp = Cx-Cx2;				/*Get distance between the elements*/
				Ytemp = Cy-Cy2;	
				dtemp = Xtemp*Xtemp + Ytemp*Ytemp;
				dist = sqrt(dtemp);

				DeltaT = theta[num]-theta[num2];
				/*printf("\n8 DeltaT = %f = %f - %f", DeltaT, theta[num], theta[num2]);*/

				if(fabs(DeltaT)<fabs(Range))
				{
					Score += 1.0/*dist*/;
				}

			}


			/*Once all the Elements Have been looked at then we need to calculate the over all FSC*/
			FCS += (Score/count);
			FCSE[num] = (Score);
			/*printf("\nScore = %f  count = %f  FCS = %f",Score, count, FCS);*/
		}
	}

	/*sprintf(plotname,"FCSPlot%i.txt",itt); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				fprintf(outfile,"%lf\t",FCSE[Number[i][j].n-1]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);*/

FCS *= 100;
FCS /= Tcount;
free(FCSE);
return(FCS);
}






//GP indicate the number of gauss point, GP=0 is the first one
void CompositeBmatrix(int GP,  double *Bmatrix, int numply, double plyt, double h, double hxz, int num, int NumNodes, Coord NodeCoord[NumNodes], int *Lnodes, double *detJJ)
{
    /*NOTE For the sake of the dgemm functions set up here is column strong (arrays go down thecolumns first, then along the rows). The fina out put however will be the other way around so that we can use it in assmeble*/
    int i;
    double dn1dx,dn1dy,dn2dx,dn2dy,dn3dx,dn3dy,dn4dx,dn4dy,nx,ny,n1j,n2j,n3j,n4j;
    /*Intialize all the matrices we'll be using*/
    //double  *Bmatrix;  //needs to be define outside
    /*Material property matrices*/
    double detJ;
   // Bmatrix = malloc(144*sizeof(double)); // need to be define outside

    /*Shape function Matrices*/
    double *BB,*BM,*BS;
    BM = malloc(24*sizeof(double));
    BB = malloc(36*sizeof(double));
    BS = malloc(24*sizeof(double));
    
   // double cx = NodeCoord[Lnodes[0]].xz + h2xz;
    //double cy = NodeCoord[Lnodes[0]].y + h2;
    
    
    for(i=0;i<144;i++)
    {
        Bmatrix[i] = 0.0000000000000000;	/*Intialize every value to zero*/
    }
    

        if((GP==0)||(GP==3)){nx = -0.5773502692;}
        else{nx = 0.5773502692;}
        
        if((GP==0)||(GP==1)){ny = -0.5773502692;}
        else{ny = 0.5773502692;}
    
    int x1,x2,x3,x4,y1,y2,y3,y4;
    x1=NodeCoord[Lnodes[0]].xz;
    x2=NodeCoord[Lnodes[1]].xz;
    x3=NodeCoord[Lnodes[2]].xz;
    x4=NodeCoord[Lnodes[3]].xz;
    
    y1=NodeCoord[Lnodes[0]].y;
    y2=NodeCoord[Lnodes[1]].y;
    y3=NodeCoord[Lnodes[2]].y;
    y4=NodeCoord[Lnodes[3]].y;
    
   // printf("\n GP = %i, x = %d,%d,%d,%d, y = %d,%d,%d,%d", GP, x1,x2,x3,x4, y1,y2,y3,y4);
    
 
        /*printf("\n GP = %i, nx = %f, ny = %f", GP, nx, ny);*/
        
        
        n1j = 0.25*(1-nx)*(1-ny);
        n2j = 0.25*(1+nx)*(1-ny);
        n3j = 0.25*(1+nx)*(1+ny);
        n4j = 0.25*(1-nx)*(1+ny);
        
        dn1dx = -0.25*(1-ny)/hxz;
        dn2dx = 0.25*(1-ny)/hxz;
        dn3dx = 0.25*(1+ny)/hxz;
        dn4dx = -0.25*(1+ny)/hxz;
        
        dn1dy = -0.25*(1-nx)/h;
        dn2dy = -0.25*(1+nx)/h;
        dn3dy = 0.25*(1+nx)/h;
        dn4dy = 0.25*(1-nx)/h;
   
    double *Jmatrix, *JmatrixInv, *dndtemp, *res;
    
    Jmatrix=malloc(4*sizeof(double));
    JmatrixInv=malloc(4*sizeof(double));
    dndtemp=malloc(8*sizeof(double));
    res=malloc(8*sizeof(double));
    
    dndtemp[0]=dn1dx;
    dndtemp[1]=dn2dx;
    dndtemp[2]=dn3dx;
    dndtemp[3]=dn4dx;
    dndtemp[4]=dn1dy;
    dndtemp[5]=dn2dy;
    dndtemp[6]=dn3dy;
    dndtemp[7]=dn4dy;
    
    Jmatrix[0]= dn1dx*x1+dn2dx*x2+dn3dx*x3+dn4dx*x4;
    Jmatrix[1]= dn1dx*y1+dn2dx*y2+dn3dx*y3+dn4dx*y4;
    Jmatrix[2]= dn1dy*x1+dn2dy*x2+dn3dy*x3+dn4dy*x4;
    Jmatrix[3]= dn1dy*y1+dn2dy*y2+dn3dy*y3+dn4dy*y4;
    double Jmax1,Jmax2,Jmax3;
    Jmax1=Jmatrix[1];
    Jmax2=Jmatrix[2];
    Jmax3=Jmatrix[3];
    
    // invrse J= 1/det(J)      *  [Jmatrix change 0 to 3, and add a '-' in 1 and 2]
    
    detJ= Jmatrix[0]*Jmatrix[3]-Jmatrix[1]*Jmatrix[2];
    detJ= fabs(detJ);
    // pass the determinate value of Jacobian matrix to the main program(global vairiable)
    *detJJ=detJ;
    
    
    JmatrixInv[0]=Jmatrix[3]/detJ;
    JmatrixInv[1]=-1.0*Jmatrix[1]/detJ;
    JmatrixInv[2]=-1.0*Jmatrix[2]/detJ;
    JmatrixInv[3]=Jmatrix[0]/detJ;
    //printf("Jinv: %f, %f, %f, %f",JmatrixInv[0],JmatrixInv[1],JmatrixInv[2],JmatrixInv[3]);
    

    
    MatrixMultiple(JmatrixInv, dndtemp, res, 2, 2, 4, 2);
    
    dn1dx=res[0];
    dn2dx=res[1];
    dn3dx=res[2];
    dn4dx=res[3];
    dn1dy=res[4];
    dn2dy=res[5];
    dn3dy=res[6];
    dn4dy=res[7];
    
    free(res);
    free(Jmatrix);
    free(JmatrixInv);
    free(dndtemp);
    
    BM[0] = dn1dx;
    BM[3] = 0;
    BM[6] = dn2dx;
    BM[9] = 0;
    BM[12] = dn3dx;
    BM[15] = 0;
    BM[18] = dn4dx;
    BM[21] = 0;
    
    BM[1] = 0;
    BM[4] = dn1dy;
    BM[7] = 0;
    BM[10] = dn2dy;
    BM[13] = 0;
    BM[16] = dn3dy;
    BM[19] = 0;
    BM[22] = dn4dy;
    
    BM[2] = dn1dy;
    BM[5] = dn1dx;
    BM[8] = dn2dy;
    BM[11] = dn2dx;
    BM[14] = dn3dy;
    BM[17] = dn3dx;
    BM[20] = dn4dy;
    BM[23] = dn4dx;
    
    BB[0] = 0;
    BB[3] = 0;
    BB[6] = -dn1dx;
    BB[9] = 0;
    BB[12] = 0;
    BB[15] = -dn2dx;
    BB[18] = 0;
    BB[21] = 0;
    BB[24] = -dn3dx;
    BB[27] = 0;
    BB[30] = 0;
    BB[33] = -dn4dx;
    
    BB[1] = 0;
    BB[4] = dn1dy;
    BB[7] = 0;
    BB[10] = 0;
    BB[13] = dn2dy;
    BB[16] = 0;
    BB[19] = 0;
    BB[22] = dn3dy;
    BB[25] = 0;
    BB[28] = 0;
    BB[31] = dn4dy;
    BB[34] = 0;
    
    BB[2] = 0;
    BB[5] = dn1dx;
    BB[8] = -dn1dy;
    BB[11] = 0;
    BB[14] = dn2dx;
    BB[17] = -dn2dy;
    BB[20] = 0;
    BB[23] = dn3dx;
    BB[26] = -dn3dy;
    BB[29] = 0;
    BB[32] = dn4dx;
    BB[35] = -dn4dy;
    
    BS[0] = dn1dx;
    BS[2] = 0;
    BS[4] = n1j;
    BS[6] = dn2dx;
    BS[8] = 0;
    BS[10] = n2j;
    BS[12] = dn3dx;
    BS[14] = 0;
    BS[16] = n3j;
    BS[18] = dn4dx;
    BS[20] = 0;
    BS[22] = n4j;
    
    BS[1] = dn1dy;
    BS[3] = -n1j;
    BS[5] = 0;
    BS[7] = dn2dy;
    BS[9] = -n2j;
    BS[11] = 0;
    BS[13] = dn3dy;
    BS[15] = -n3j;
    BS[17] = 0;
    BS[19] = dn4dy;
    BS[21] = -n4j;
    BS[23] = 0;
    


        
        int kk;
        for (kk=0; kk<4; kk++) {
        
        
        
        Bmatrix[0+6*kk]=     BM[0+6*kk];
        Bmatrix[24+6*kk]=    BM[1+6*kk];
        Bmatrix[48+6*kk]=    BM[2+6*kk];
        Bmatrix[1+6*kk]=     BM[3+6*kk];
        Bmatrix[25+6*kk]=    BM[4+6*kk];
        Bmatrix[49+6*kk]=    BM[5+6*kk];
        
        Bmatrix[75+6*kk]=    BB[0+9*kk];
        Bmatrix[97+6*kk]=    BB[1+9*kk];
        Bmatrix[121+6*kk]=   BB[2+9*kk];
        Bmatrix[76+6*kk]=    BB[3+9*kk];
        Bmatrix[98+6*kk]=    BB[4+9*kk];
        Bmatrix[122+6*kk]=   BB[5+9*kk];
        Bmatrix[77+6*kk]=    BB[6+9*kk];
        Bmatrix[99+6*kk]=    BB[7+9*kk];
        Bmatrix[123+6*kk]=   BB[8+9*kk];
        
            
            
    //shear term
            Bmatrix[97+6*kk]+=    BS[0+6*kk];
            Bmatrix[121+6*kk]+=   BS[1+6*kk];
            Bmatrix[98+6*kk]+=    BS[2+6*kk];
            Bmatrix[122+6*kk]+=   BS[3+6*kk];
            Bmatrix[99+6*kk]+=    BS[4+6*kk];
            Bmatrix[123+6*kk]+=   BS[5+6*kk];
            
        
        }
    
    free(BB);
    free(BM);
    
    /////// shear term, use single gauss point
/*
    nx = 0;
    ny = 0;
    
    n1j = 0.25*(1-nx)*(1-ny);
    n2j = 0.25*(1+nx)*(1-ny);
    n3j = 0.25*(1+nx)*(1+ny);
    n4j = 0.25*(1-nx)*(1+ny);
    
    dn1dx = -0.25*(1-ny)/hxz;
    dn2dx = 0.25*(1-ny)/hxz;
    dn3dx = 0.25*(1+ny)/hxz;
    dn4dx = -0.5*(1+ny)/hxz;
    
    dn1dy = -0.25*(1-nx)/h;
    dn2dy = -0.25*(1+nx)/h;
    dn3dy = 0.25*(1+nx)/h;
    dn4dy = 0.25*(1-nx)/h;
    
    
    BS[0] = dn1dx;
    BS[2] = 0;
    BS[4] = n1j;
    BS[6] = dn2dx;
    BS[8] = 0;
    BS[10] = n2j;
    BS[12] = dn3dx;
    BS[14] = 0;
    BS[16] = n3j;
    BS[18] = dn4dx;
    BS[20] = 0;
    BS[22] = n4j;
    
    BS[1] = dn1dy;
    BS[3] = -n1j;
    BS[5] = 0;
    BS[7] = dn2dy;
    BS[9] = -n2j;
    BS[11] = 0;
    BS[13] = dn3dy;
    BS[15] = -n3j;
    BS[17] = 0;
    BS[19] = dn4dy;
    BS[21] = -n4j;
    BS[23] = 0;


    for (kk=0; kk<4; kk++) {
        
        i=kk;

        
    }
    */
    
    
    free(BS);

        
}


// c indicate the row number of matrix a, d is for matrix d

void MatrixMultiple(double *a, double *b, double *res, int acol, int arow, int bcol, int brow)

{
    //res= malloc(arow*bcol*sizeof(double));
    //res=malloc(acol*brow);
    double temp;
    int i,j,k;
    for (i=0; i<arow; i++)
    {
        for(j=0; j<bcol; j++)
        {
            
            temp=0.0;
            for(k=0;k<acol;k++)
            {
                temp +=a[i*acol+k]*b[(k)*bcol+j];
            
            }
            
            
            res[i*bcol+j]=temp;
        
        
        }
        
    }
}


void TransM(double *a, double *b, int mcol, int mrow)
{
    double *temp;
    temp=malloc(mcol*mrow*sizeof(double));
    int i,j,k;
    for (i=0; i<mrow; i++) {
        
        for (j=0;j<mcol; j++) {
            
            temp[i*mcol+j]=a[mcol*j+i];
            
        }
        
    }
    for (k=0; k<mcol*mrow; k++) {
        b[k]=temp[k];
    }
    free(temp);
}




void GPsensitivity( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes, double *NodeDTL, double *NodesenTemp, Gstn *GPsens)
{
    
    double temp;
    double cx, cy;
    double h2 = 0.5 * h;
    double h2xz;
    int Gcount = 0;
    int gpoints = 4*elemX*elemY;
    double gax, gay;
    Gstn *GaussValue;
    int *Lnodes;
    int num;
    
    int n,m,i,j;
    Lnodes = malloc(4*sizeof(int));
    
    GaussValue = malloc(4*elemX*elemY*sizeof(Gstn));
    
    for(n=0;n<elemX;n++)
    {
        for(m=0;m<elemY;m++)
        {
            
            num = Number[n][m].n-1;		/*Get element number*/
            Lnodes[0] = Number[n][m].a-1;
            Lnodes[1] = Number[n][m].b-1;
            Lnodes[2] = Number[n][m].c-1;
            Lnodes[3] = Number[n][m].d-1;
            
            h2xz = 0.5 * hxz[n];
            /*work out global co-ordinates of element center, from co-ords of node 1 and element edge length*/
            double cx = NodeCoord[Lnodes[0]].xz + h2xz;
            double cy = NodeCoord[Lnodes[0]].y + h2;
            double *Nn;
            Nn = malloc(4*sizeof(double));
            /*Now go through each node calculating the gauss matrices and from it the sensitvity contribution of each node*/
            for(i=0;i<4;i++)
            {
                
                /*Get Gauss point cordinates*/
                gax = gaB * (double)PoA[i].x * h2xz;
                gay = gaB * (double)PoA[i].y * h2;
                GaussValue[Gcount].x = cx + gax;
                GaussValue[Gcount].y = cy + gay;
                Nn[0]=0.25*(1-gaB*(double)PoA[i].x)*(1-gaB* (double)PoA[i].y);
                Nn[1]=0.25*(1+gaB*(double)PoA[i].x)*(1-gaB* (double)PoA[i].y);
                Nn[2]=0.25*(1+gaB*(double)PoA[i].x)*(1+gaB* (double)PoA[i].y);
                Nn[3]=0.25*(1-gaB*(double)PoA[i].x)*(1+gaB* (double)PoA[i].y);
                
                GaussValue[Gcount].u=NodesenTemp[Lnodes[0]]*Nn[0]+NodesenTemp[Lnodes[1]]*Nn[1]+NodesenTemp[Lnodes[2]]*Nn[2]+NodesenTemp[Lnodes[3]]*Nn[3];
                GaussValue[Gcount].a = 0.25*h*hxz[n];	/*Store the gauss point element area contribution*/
                printf("\nSens at gauss point: %d, %f", Gcount, GaussValue[Gcount].u);
                Gcount++;
            }
        }
    }
    for(n=0;n<4*elemX*elemY;n++)
    {
        GPsens[n].u=GaussValue[n].u;
        GPsens[n].a=GaussValue[n].a;
        GPsens[n].x=GaussValue[n].x;
        GPsens[n].y=GaussValue[n].y;
        
    }
    
    
    free(GaussValue);
    free(Lnodes);
    
}



void NodeAuxsens( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes,  Gstn *GPsens, double *Nodesenpp)
{
    
    double temp;
    double cx, cy;
    double h2 = 0.5 * h;
    double h2xz;
    int Gcount = 0;
    int gpoints = 4*elemX*elemY;
    double gax, gay;
    Gstn *GaussValue;
    int *Lnodes;
    int num;
    
    int n,m,i,j;
    Lnodes = malloc(4*sizeof(int));
    
    GaussValue = malloc(4*elemX*elemY*sizeof(Gstn));
    
    for(n=0;n<(4*elemX*elemY);n++)
    {
        GaussValue[n].u = GPsens[n].u;
        if(n==8||n==9||n==10||n==11||n==88||n==89||n==90||n==91||n==168||n==169||n==170||n==171)
        {
           // printf("\n Gauss%d:  %f",n,GaussValue[n].u);
        }
        GaussValue[n].a= GPsens[n].a;
        GaussValue[n].x= GPsens[n].x;
        GaussValue[n].y= GPsens[n].y;
        Gcount++;
    }
    
    double rad = 2.0*h;
    double ftemp;
    
    for(n=0;n<(Ntot_new);n++)
    {
        Nodesenpp[n] = 0.0;
    }
    
    /*Now we need to get the nodal values of sensitvity from the elemental values by sending each node through the least squared routine
     Simpler then the standard 2D or 3D code as we have no AUX nodes to concern ourselves with*/
    for(n=0;n<Ntot_new;n++)
    {
        if (n<NumNodes) {
            
            cx = NodeCoord[n].xz;
            cy = NodeCoord[n].y; /*read in co-ordiantes directly*/
            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            Nodesenpp[n]= ftemp; /* multiply smothed strain energy by weight */
            printf("\nMesh nodes: no.%d, [%f,%f]", n, cx,cy);
        }
        else{
            cx = AuxNodes[n-NumNodes].x;
            cy = AuxNodes[n-NumNodes].y; /*read in co-ordiantes directly*/
            printf("\nAuxnode no.%d, [%f,%f]", n, cx,cy);
            // if (AuxNodes[n-NumNodes].n+1>NumNodes)    // if number is larger than NumNodes, it haven't been calculate
            
            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            Nodesenpp[n]= ftemp; /* multiply smothed strain energy by weight */
            
        }
        
    }
    
    free(GaussValue);
    free(Lnodes);
    
}

void NodeAuxsens2( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h, double *hxz, int Ntot_new, Aux *AuxNodes,  double *Nodesenpp, double *NodesenTemp)
{
    
    double temp;
    double cx, cy;
    double h2 = 0.5 * h;
    double h2xz;
    int Gcount = 0;
    int gpoints = 4*elemX*elemY;
    double gax, gay;
    Gstn *GaussValue;

    int num;
    
    int n,m,i,j;

    
    GaussValue = malloc((elemX+1)*(elemY+1)*sizeof(Gstn));
    
    for(n=0;n<(elemX+1)*(elemY+1);n++)
    {
            double cx = NodeCoord[n].xz;
            double cy = NodeCoord[n].y;
        
        GaussValue[Gcount].x = cx;
        GaussValue[Gcount].y = cy;

                
        GaussValue[Gcount].u=NodesenTemp[n];
        //GaussValue[Gcount].a = 0.25*h*hxz[n];	/*Store the gauss point element area contribution*/
        GaussValue[Gcount].a = 1.0;
        Gcount++;
    }
    
    double rad = 3.0*h;
    double ftemp;
    
    for(n=0;n<(Ntot_new);n++)
    {
        Nodesenpp[n] = 0.0;
    }
    
    /*Now we need to get the nodal values of sensitvity from the elemental values by sending each node through the least squared routine
     Simpler then the standard 2D or 3D code as we have no AUX nodes to concern ourselves with*/
    for(n=0;n<Ntot_new;n++)
    {
        if (n<NumNodes) {
            
            Nodesenpp[n]= NodesenTemp[n]; /* multiply smothed strain energy by weight */
            
        }
        else{
            cx = AuxNodes[n-NumNodes].x;
            cy = AuxNodes[n-NumNodes].y; /*read in co-ordiantes directly*/
            
            // if (AuxNodes[n-NumNodes].n+1>NumNodes)    // if number is larger than NumNodes, it haven't been calculate
            
            ftemp = LstrainV23(cx, cy, rad, Gcount, GaussValue,4); /*weighted by alpha / distance*/
            Nodesenpp[n]= ftemp; /* multiply smothed strain energy by weight */
            
        }
        
    }
    
    free(GaussValue);


}


void NodeAuxsens3( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h,  int Ntot_new, Aux *AuxNodes,   double *Nodesenpp,double *NodesenTemp)
{
    double nx, ny;

    int *Lnodes;

    double x2,y2;
    int n,i;
    Lnodes = malloc(4*sizeof(int));
    
    double min = 0.0;
    
    
    double rad = 1.0*h;
    
    double r2 = rad * rad;
    
    
    double ftemp;
    
    for(n=0;n<(Ntot_new);n++)
    {
        Nodesenpp[n] = 0.0;
    }
    
    /*Now we need to get the nodal values of sensitvity from the elemental values by sending each node through the least squared routine
     Simpler then the standard 2D or 3D code as we have no AUX nodes to concern ourselves with*/
    for(n=0;n<Ntot_new;n++)
    {
        if (n<NumNodes) {
            
            nx = NodeCoord[n].xz;
            ny = NodeCoord[n].y; /*read in co-ordiantes directly*/
            Nodesenpp[n]= NodesenTemp[n]; /* multiply smothed strain energy by weight */
        }
        else{
            
            nx = AuxNodes[n-NumNodes].x;
            ny = AuxNodes[n-NumNodes].y; /*read in co-ordiantes directly*/
            printf("\nAuxnode no.%d, [%f,%f]", n, nx,ny);
            int count=0;
            int NumofNode[2];
            for(i=0;i<NumNodes;i++)
            {
                double cx = NodeCoord[i].xz;
                double cy = NodeCoord[i].y;
                
                
                ftemp = nx - cx;
                x2 = ftemp * ftemp;	/*squared x-distance from current gauss point to node of interest*/
                ftemp = ny - cy;
                y2 = ftemp * ftemp;	/*squared y-distance from current gauss point to node of interest*/
                
                ftemp = x2 + y2; /*squared distance between gauss point and node*/
                
                if((ftemp < r2) && (r2 >= min))/*if distance less than the radius squared, then add data to arrays*/
                {
                    NumofNode[count]=i;
                    count++;
                }
            }
            
    
            Nodesenpp[n]= NodesenTemp[NumofNode[0]]+NodesenTemp[NumofNode[1]]; /* multiply smothed strain energy by weight */
            
        }
        
    }
    

    free(Lnodes);
    
}



// get sensitivity at auxillary point by using shape function, actually it is linear interploation
// now this function only work with square elemets, which the length and width are both 1
//

void NodeAuxsens4( int elemX, int elemY, Elem Number[elemX][elemY], int NumNodes, Coord NodeCoord[NumNodes], double h,  int Ntot_new, Aux *AuxNodes,   double *Nodesenpp,double *NodesenTemp)
{
    double nx, ny;
    
    int *Lnodes;
    
    double x2,y2;
    int n,i;
    Lnodes = malloc(4*sizeof(int));
    
    double min = 0.0;
    
    
    double rad = 1.0*h;
    
    double r2 = rad * rad;
    
    
    double ftemp;
    
    for(n=0;n<(Ntot_new);n++)
    {
        Nodesenpp[n] = 0.0;
    }
    
    /*Now we need to get the nodal values of sensitvity from the elemental values by sending each node through the least squared routine
     Simpler then the standard 2D or 3D code as we have no AUX nodes to concern ourselves with*/
    for(n=0;n<Ntot_new;n++)
    {
        if (n<NumNodes) {
            
            nx = NodeCoord[n].xz;
            ny = NodeCoord[n].y; /*read in co-ordiantes directly*/
            Nodesenpp[n]= NodesenTemp[n]; /* multiply smothed strain energy by weight */
        }
        else{
            
            nx = AuxNodes[n-NumNodes].x;
            ny = AuxNodes[n-NumNodes].y; /*read in co-ordiantes directly*/
            printf("\nAuxnode no.%d, [%f,%f]", n, nx,ny);
            int count=0;
            int NumofNode[2];
            for(i=0;i<NumNodes;i++)
            {
                double cx = NodeCoord[i].xz;
                double cy = NodeCoord[i].y;
                
                
                ftemp = nx - cx;
                x2 = ftemp * ftemp;	/*squared x-distance from current gauss point to node of interest*/
                ftemp = ny - cy;
                y2 = ftemp * ftemp;	/*squared y-distance from current gauss point to node of interest*/
                
                ftemp = x2 + y2; /*squared distance between gauss point and node*/
                
                if((ftemp < r2) && (r2 >= min))/*if distance less than the radius squared, then add data to arrays*/
                {
                    NumofNode[count]=i;
                    count++;
                }
            }
            
            // node in x or y
            
            double x1,y1,x2,y2;
            x1=NodeCoord[NumofNode[0]].xz;
            x2=NodeCoord[NumofNode[1]].xz;
            y1=NodeCoord[NumofNode[0]].y;
            y2=NodeCoord[NumofNode[1]].y;
            
            
            if (abs(x1-x2)<=0.00000001 ) {
                
                Nodesenpp[n]= NodesenTemp[NumofNode[0]]+ (ny-y1)*(NodesenTemp[NumofNode[1]]-NodesenTemp[NumofNode[0]])/(y2-y1);
            }
            else
            {
                Nodesenpp[n]= NodesenTemp[NumofNode[0]]+ (nx-x1)*(NodesenTemp[NumofNode[1]]-NodesenTemp[NumofNode[0]])/(x2-x1);
            }
           // Nodesenpp[n]= NodesenTemp[NumofNode[0]]+NodesenTemp[NumofNode[1]]; /* multiply smothed strain energy by weight */
            
        }
        
    }
    
    
    free(Lnodes);
    
}


void ScaleVnorm(int elemX, int elemY, double *ElemThetaStn , double *DeltaTheta, double DeltaThetaMax, int itt, double CompA[itt],  Elem Number[elemX][elemY], short *ElemStat, short *ElemStatb, int Numlsf, double *NodeThetaStn, double *NodeDTL, double *DeltaT, int Ntot_new , double *Vnorm,double *SensMaxOld)
{
    FILE *outfile;	/*File varible for output files*/
    char plotname[40];	/*variable to change names of plotting output files*/
    int i,j,num,p,k;
    double Ctemp, Ctemp2, MaxC, MinC, SensMax, SensTemp2, SensTemp, SensSum, SensMaxTemp,SensMaxTemp2;
    
    //SensMaxTemp = *SensMaxOld; /*Set SensMax to 0*/
    SensMax=*SensMaxOld;
    SensSum = 0; /*Set SensMax to 0*/
    /*First find the maximum sensitvity values*/
    for(i=0; i<(elemX+1)*(elemY+1);i++)
    {
    
            SensMax = (fabs(Vnorm[i])>fabs(SensMax)) ? fabs(Vnorm[i]):fabs(SensMax);
            /*k = i + (elemX*elemY)*p;
             SensMax = (fabs(ElemThetaStn[i]*ElemStat[k])>fabs(SensMax)) ? fabs(ElemThetaStn[i]):fabs(SensMax);*/

    }
    
    *SensMaxOld=SensMax;
    
    
    /*SensMax = SensSum/(elemX*elemY);*/
    //*SensMaxOld=SensMax;
    //SensMaxTemp2=SensMax;
    //if (SensMax-SensMaxTemp<0.0000000001)
    //{
    //     SensMax=SensMaxTemp;
    // }
    // set the sensmax is equal to the sensmax in his own itteration instead the old one which is =0.
    
    
    
    //SensMax = (fabs(SensMaxTemp)>fabs(SensMax)) ? fabs(SensMaxTemp):fabs(SensMax);
    printf("\n SensMax = %f \n", SensMax);
    
    /*Now if more that 10 itterations have past find the minmum and Maximum compliance from the last 10 itterations*/
    if(itt<10)
    {
        Ctemp = 1.0;
    }
    else
    {
        MaxC = 0;				/*Set to small value*/
        MinC = 100000000000;	/*Set to very large value*/
        /*Get max and min Compliance*/
        for(i=0; i<(10);i++)
        {
            j = itt-i;
            MaxC = (CompA[j]>MaxC) ? CompA[j]:MaxC;
            MinC = (CompA[j]<MinC) ? CompA[j]:MinC;
        }
        /*printf("\n MaxC = %f", MaxC);
         printf("\n MinC = %f", MinC);*/
        
        Ctemp2 = (MaxC - MinC)/MaxC;	/*Now calculate Ctemp from these values*/
        Ctemp = sqrt(sqrt(sqrt(sqrt(Ctemp2))));			/*Sqaure route relaxes the positive feed back of scaling this off compliance*/
        Ctemp = 1.0;
    }
    
    /*Now scale delta theta for each element to the compliance and maximum sensitvtiy*/
    
    //DeltaThetaMax=0.1;
    double Maxthetachange;
    //Maxthetachange=3.0*(3.14159265/180.0);
    
    //*DeltaT = Ctemp*Maxthetachange/fabs(SensMax);
    //*DeltaT = Ctemp*DeltaThetaMax/(1+fabs(SensMax));
    *DeltaT = Ctemp*DeltaThetaMax/(fabs(SensMax));		/*Get desired change in theta.*/
    
    if (fabs(SensMax)<0.1)
    {
        *DeltaT=1.0;
    
    }
    
    /*printf("\tDeltaTheta = %f", DeltaTheta[i]);*/
    
    
    /*Print out the result*
     sprintf(plotname,"DeltaTheta%i.txt",itt);
     outfile = fopen(plotname, "w");
     if(outfile == NULL){
     printf("\nFailed to open Element Strain Energy plotting writefile\n");
     }
     else{
     
     for(j=elemY-1;j>=0;j--)
     {
     for(i=0;i<elemX;i++)
     {
     num = Number[i][j].n - 1;	
     fprintf(outfile,"%f\t",DeltaTheta[num]); 
     }
     fprintf(outfile,"\n");
     }
     }
     fclose(outfile);
     
     /*For ease of calculation and viewing that was calculated in degrees, but computer porgrams always use radiums so transform it*/
    
}
