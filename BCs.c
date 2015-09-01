/*
 *  BCs.c
 *  
 *  Created by Peter Dunning on 22/06/2010.
 *  functions to apply boundary conditions and construct required force vectors
 *
 */

#include "BCs.h"
#include "FixedGrid.h"
#include "ls_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*function that finds all zero dispacement boundary conditions and assembles into one array */
void FixDisp(double h, int NumNodes, Coord NodeCoord[NumNodes], int *NumFix, int *fixed, int *TotFix, int *FixDofs, 
				int *LAfix, Aux *Afix, FILE *infile)
{
	int i,j,k,dtemp,temp; /*incrementors etc. */	
	double xtemp, ytemp;	/*Temporary variables to determine closest node to set of co-ordinates & direction*/
	int NumFixtemp = 0;	/*temporary variable to count number of nodes with fixed level set values, initially zero*/
	
	/*--------------------------This section applies zero displacement boundary conditions to nodes--------------------------------*/
	int bn; /*Number of Zero displacement Boundary conditions */
	fscanf(infile, "%i", &bn);
	printf("\n\nNumber of Boundary conditions on nodes = %i",bn);
	int *ZeroDisplace; /*Array Containing info for rows and columns to be set to zero*/
	ZeroDisplace = malloc(bn * sizeof(int));
	
	for(i=0;i<bn;i++)
	{
		fscanf(infile, "%lf", &xtemp);
		fscanf(infile, "%lf", &ytemp);
		fscanf(infile, "%i", &dtemp);
		
		temp = closeNode(h,xtemp,ytemp,NumNodes,NodeCoord);	/*find closest node to zero-displacement co-ordinates*/
		/*printf("\nZero Boundary Condition %i applied at Node %i in direction %i",n+1,temp+1,dtemp);*/
		if(*NumFix != 0)
		{
			fixed[NumFixtemp++] = temp; /*fix level set values for this node throughout the optimisation*/
		}
		/*Now fix the specified dof*/
		if(dtemp !=0)
		{
			temp *= 6;
			temp += dtemp;
			ZeroDisplace[i] = temp;
		}
	}

	/*--------------------------This section determines dof's to be fixed relating to a fixed (rectangular) area----------------------------*/
	short *dofFix;
	short *nodeFixtemp;
	dofFix = calloc(6 * NumNodes,sizeof(short)); /*Initialize array to store dof data for fixed areas*/
	/*Coord *Boxs;
	Boxs = malloc(4 * sizeof(double)); /*Array to store min and max co-ords of a rectangular area*/
	double BX0,BX1,BY0,BY1;
	
	int bsn; /*Number of areas where boundary conditions are to be applied*/
	fscanf(infile, "%i", &bsn); /*Read in number of fixed (rectangular) areas*/
	
	printf("\nNumber of Fixed Areas = %i",bsn);
	for(i=0;i<bsn;i++)
	{
		/*Read in minimum (Boxs[0]) and maximum (Boxs[1]) co-ordinates of current fixed area*/
		fscanf(infile, "%lf", &BX0);
		fscanf(infile, "%lf", &BX1);
		fscanf(infile, "%lf", &BY0);
		fscanf(infile, "%lf", &BY1);

		/*printf("\nBoxs[0].x = %lf", BX0);
		printf("\nBoxs[1].x = %lf", BX1);
		printf("\nBoxs[0].y = %lf", BY0);
		printf("\nBoxs[1].y = %lf\n", BY1);*/
		
		nodeFixtemp = malloc(NumNodes*sizeof(short));	/*Tempary array to store node data for current fixed area*/	
		
		for(j=0;j<NumNodes;j++) /*For all nodes - determine if the node lies within the fixed rectangular area*/
		{
			nodeFixtemp[j] = 0;

			if( ((NodeCoord[j].x - BX1) < 0.000001) && ((NodeCoord[j].x - BX0) > -0.000001) ) /*If within x bounds*/
			{
				if( ((NodeCoord[j].y - BY1) < 0.000001) && ((NodeCoord[j].y - BY0) > -0.000001) ) /*If also within y bounds*/
				{
					nodeFixtemp[j] = 1; /*Node is to be fixed*/
					
					if(*NumFix != 0)
					{
						fixed[NumFixtemp++] = j; /*fix level set values for this node throughout the optimisation*/
					}
					/*printf("\nNode %i to be fixed",j+1);*/
				}
			}
			/*printf("\nnodeFixtemp[%i] = %i", j, nodeFixtemp[j]);*/
		}
		
		/*if(*LAfix != 0) /*If using BFG, copy fixed area data*
		{
			temp = 2 * i; /*to store data properly*
			Afix[temp].x = Boxs[0].x;
			Afix[temp].y = Boxs[0].y;
			Afix[temp+1].x = Boxs[1].x;
			Afix[temp+1].y = Boxs[1].y;
		}*/
		
		for(k=0;k<6;k++) /*for each direction (x,y)*/
		{
			fscanf(infile, "%i", &temp); /*read in binary flag to determine if current direction is to fixed*/
			
			if(*LAfix != 0)	/*If using BFG, copy fixed direction flags*/
			{
				Afix[(2 * i) + k].n = temp;
			}
			
			if(temp == 1) /*If current direction is to be fixed*/
			{
				for(j=0;j<NumNodes;j++) /* for all nodes in domain update assossiated dof's to be fixed*/
				{
					if(nodeFixtemp[j] != 0)
					{
						temp = j * 6;
						temp += k;
						dofFix[temp] = 1;
						/*printf("\nZero Boundary Condition applied at Node %i in direction %i",j+1,dtemp);*/
					}
				}
			}
		}
		free(nodeFixtemp);	/*Free up some memory*/
	}

	dtemp = 6 * NumNodes; /*Number of dofs*/
	int LzeroSurface = 0; /*Variable for length of Zero displacement boundary conditions related to areas*/
	for(i=0;i<dtemp;i++) /*Now sum all entries in dofFix array to determine the number of dofs to be fixed*/
	{
		LzeroSurface += dofFix[i];
	}

	printf("\nNumber of Area fixed dofs = %i", LzeroSurface);

	int *ZeroSurface; /*Initialize array that will contain info on the dof numbers that are fixed*/
	ZeroSurface = malloc(LzeroSurface * sizeof(int));

	temp = 0;
	for(i=0;i<dtemp;i++)
	{
		if(dofFix[i] == 1) /*If current dof is to be fixed add it to the ZeroSurface array*/
		{
			ZeroSurface[temp++] = i+1;
		}
	}

	/*----------------------------------Now combine all fixed dofs into one array (in assending sequence)----------------------------------*/
	*TotFix = bn + LzeroSurface + 6; /*Total number of fixed dof, +2 for boundary condition application function later*/
	printf("\nTotal Number of fixed dofs = %i", *TotFix-6);
	int Fmax,Fmin;	/*Varaibles for subsequent sequencing algorithm*/
	FixDofs[0] = 0;
	FixDofs[*TotFix-5] = 100 * dtemp; /*set this to a large number to avoid problems with auxillary nodes*/

	for(i=(*TotFix-6);i>=0;i--)
	{
		Fmax = FixDofs[i+1];	/*Set max to previous entry*/
		Fmin = 0;
		
		/*Two loops to find maximum dof less than the current lowest in FixDofs*/
		for(j=0;j<bn;j++)
		{
			if((ZeroDisplace[j] < Fmax) && (ZeroDisplace[j] > Fmin))
			{
				Fmin = ZeroDisplace[j];
			}
		}
		
		for(j=0;j<LzeroSurface;j++)
		{
			if((ZeroSurface[j] < Fmax) && (ZeroSurface[j] > Fmin))
			{
				Fmin = ZeroSurface[j];
			}
		}
		
		FixDofs[i] = Fmin;
	}
	
	/*--------------------- Cleaning up etc. -------------------*/
	if(*NumFix != 0) /*if lsf values are fixed at bcs*/
	{
		*NumFix = NumFixtemp; 
	}
	
	if(*LAfix != 0)	/*If using IFG, update length of Afix array*/
	{
		*LAfix = 2 * bsn; /*twice number of fixed areas*/
	}
	
	/*Free up memory used in this function*/
	free(dofFix);
	free(ZeroDisplace);
	free(ZeroSurface);
	/*free(Boxs);*/
	/*Print FixDofs to screen
	for(i=1;i<*TotFix;i++)
	{
		printf("\nFixDofs[%i] = %i",i,FixDofs[i]);
	}*/

}

/*function that assembles the load vector for Objective (1), minimisation of deterministic strain energy*/
void LoadVector(double h, int NumNodes, Coord NodeCoord[NumNodes], int *NumForces, int *rhs_ind, double *rhs_val,
					int NumPforce, double *Pforce, int NumAforce, double *Aforce)
{
	int i,j,k,temp,node; /*incrementors etc. */	
	double ftemp, ftemp2, dtemp, etemp, gtemp, htemp, itemp, xtemp, ytemp;	/*Temporary variables to determine closest node to set of co-ordinates & direction*/
	double cosT, sinT; /*variables for sin and cos of loading direction*/

	/*------------------------------------------This section applies single point forces at nodes------------------------------------*/
	/*First define force array, i.e. rhs of equation*/
	double *rhs;
	rhs = calloc(6 * NumNodes,sizeof(double)); /*rhs defined using dynamic memory allocation*/

	for(i=0;i<NumPforce;i++)
	{
		j = i * 8; /*place indicator for Pforce array*/

		node = closeNode(h,Pforce[j],Pforce[j+1],NumNodes,NodeCoord);	/*find closest node to applied force co-ordinates*/
		node *= 6; /*dof in x-direction*/	
		ftemp = Pforce[j+2]; /*magnitude of force in X*/
		dtemp = Pforce[j+3]; /*direction of load theta in Y*/
		etemp = Pforce[j+4]; /*direction of load alpha in Z*/
		gtemp = Pforce[j+5]; /*direction of load alpha in ROTX*/
		htemp = Pforce[j+6]; /*direction of load alpha in ROTY*/
		itemp = Pforce[j+7]; /*direction of load alpha in ROTZ*/
		
		/*cosT = cos(dtemp);
		sinT = sin(dtemp); /*work out sin and cos of loading direction*
		
		if(fabs(cosT) < 0.001) /*load is applied on the vertical (y) direction*
		{
			rhs[node+1] = (sinT < 0.0) ? -ftemp : ftemp; /*apply load (with correct sign)*
		}	
		else if(fabs(sinT) < 0.001) /*load is applied on the horizontal (x) direction*
		{
			rhs[node] = (cosT < 0.0) ? -ftemp : ftemp; /*apply load (with correct sign)*
		}
		else
		{*/
			rhs[node] = ftemp; /*compute & apply x-component*/
			rhs[node+1] = dtemp; /*compute & apply y-component*/
			rhs[node+2] = etemp; /*compute & apply z-component*/
			rhs[node+3] = gtemp; /*compute & apply rotx-component*/
			rhs[node+4] = htemp; /*compute & apply roty-component*/
			rhs[node+5] = itemp; /*compute & apply rotz-component*/

		/*}*/
	}
	
	/*/////////////// THE REST OF BCs.c IS CURRENTLY UNUSED SO WE ARE LEAVING IT OUT OF THE Curved Plate MODDELLING FOR NOW/////////////////////////////*/
	/*-----------------------------------------This section applies Unifrom Forces to an Area-----------------------------------------*/
	/*int Xnum1, Xnum2, Ynum1, Ynum2; /*variables to track end node numbers*
	short *nodeFixtemp;
	
	for(i=0;i<NumAforce;i++)
	{
		k = i * 8; /*place indicator for Aforce array*
		nodeFixtemp = calloc(NumNodes,sizeof(short));	/*Tempary array to store node data for current area*
		
		ftemp = 100000000.0; /*initalize minimum x-cord*
		ftemp2 = 100000000.0; /*initalize minimum y-cord*
		xtemp = -1.0;	 /*initalize maximum x-cord*
		ytemp = -1.0;	 /*initalize maximum y-cord*
		
		/*printf("\nXmin=%lf, Xmax=%lf, Ymin=%lf, Ymax=%lf",Aforce[k],Aforce[k+1],Aforce[k+2],Aforce[k+3]);*
		
		for(j=0;j<NumNodes;j++) /*For all nodes - determine if the node lies within the fixed rectangular area*
		{
			if( ((NodeCoord[j].x - Aforce[k+1]) < 0.0000001) && ((NodeCoord[j].x - Aforce[k]) > -0.0000001) )  /*If within x bounds*
			{
				if( ((NodeCoord[j].y - Aforce[k+3]) < 0.0000001) && ((NodeCoord[j].y - Aforce[k+2]) > -0.000001) )  /*If also within y bounds*
				{
					nodeFixtemp[j] = 1; /*Node has a force applied*
					/*printf("\nForce applied to Node %i",j+1);*
					
					/*section to find max and min node co-ords in area*
					Xnum1 = (NodeCoord[j].x < ftemp) ? j : Xnum1; /*update min x node number if required*
					ftemp = (NodeCoord[j].x < ftemp) ? NodeCoord[j].x : ftemp; /*check minimum x-cord*
					Ynum1 = (NodeCoord[j].y < ftemp2) ? j : Ynum1; /*update min y node number if required*
					ftemp2 = (NodeCoord[j].y < ftemp2) ? NodeCoord[j].y : ftemp2; /*check minimum y-cord*
					Xnum2 = (NodeCoord[j].x > xtemp) ? j : Xnum2; /*update max x node number if required*
					xtemp = (NodeCoord[j].x > xtemp) ? NodeCoord[j].x : xtemp; /*check maximum x-cord*
					Ynum2 = (NodeCoord[j].y > ytemp) ? j : Ynum2; /*update max y node number if required*
					ytemp = (NodeCoord[j].y > ytemp) ? NodeCoord[j].y : ytemp; /*check maximum y-cord*
				}
			}
		}
		
		/*define the end nodes*
		if(fabs(ftemp - xtemp) < h) /*if difference between max and min x-cords is < h, then use y numbers instead*
		{
			/*first check y co-ord difference isn't < h*
			if(fabs(ftemp2 - ytemp) >= h)
			{
				Xnum1 = Ynum1;
				Xnum2 = Ynum2; /*store end node numbers as Xnum's *
			}
			else
			{
				printf("\nERROR! Force area smaller than h x h !!!");
				printf("\nXnum1=%i, Xnum2=%i",Xnum1,Xnum2);
			}
		}

		ftemp = Aforce[k+4]; /*magnitude of force per unit length*
		ftemp *= h; /*mutiply force per unit length, by element edge length*
		dtemp = Aforce[k+6]; /*direction of load*
	
		cosT = cos(dtemp);
		sinT = sin(dtemp); /*work out sin and cos of loading direction*
		
		for(j=0;j<NumNodes;j++) /*for all nodes*
		{
			if(nodeFixtemp[j] != 0) /*If current node is effected by the force add force to rhs array*
			{
				ytemp = ((j == Xnum1) || (j == Xnum2)) ? 0.5 : 1.0; /*factor to scale loads on end nodes*
				temp = 2 * j; /*dof indicator*/
				
				/*----add loads to first column of rhs----*
				ftemp2 = ftemp * ytemp; /*mutiply force by scale factor*
				
				if(fabs(cosT) < 0.001) /*load is applied on the vertical (y) direction*
				{
					rhs[temp+1] = (sinT < 0.0) ? -ftemp2 : ftemp2; /*apply load (with correct sign)*
				}	
				else if(fabs(sinT) < 0.001) /*load is applied on the horizontal (x) direction*
				{
					rhs[temp] = (cosT < 0.0) ? -ftemp2 : ftemp2; /*apply load (with correct sign)*
				}
				else
				{
					rhs[temp] = cosT * ftemp2; /*compute & apply x-component*
					rhs[temp+1] = sinT * ftemp2; /*compute & apply y-component*
				}
			}
		}
		
		free(nodeFixtemp);
	}*/

	/*Print non-zero rhs entries to screen & reduce storage*/
	temp = 0; /*now used to count non-zero rhs entries*/
	dtemp = 6 * NumNodes; /*number of dofs*/
	for(j=0;j<dtemp;j++)
	{
		if(fabs(rhs[j]) > 0.000001)
		{
			rhs_ind[temp] = j;
			rhs_val[temp++] = rhs[j];
			printf("\nForce at dof %i = %lf",j+1,rhs[j]);
		}
	}
		
	*NumForces = temp; /*store total number of entries in rhs array*/
	printf("\n\nForces applied, %i entries in rhs array\n",temp);
	free(rhs);
}

/*function to compute load vectors and weights required to calculate expected strain energy (compliance)*/
void UncLoads(double h, int NumNodes, Coord NodeCoord[NumNodes], int *NumUn, int *unc_ind, double *unc_val, 
				int *NumRhs, double *unc_wgt, int NumPforce, double *Pforce, int NumAforce, double *Aforce)
{
	int i,j,k,temp,node,dtemp,num,NumF,pnt; /*incrementors etc. */	
	double ftemp, ftemp2, xtemp, ytemp;	/*Temporary variables to determine closest node to set of co-ordinates & direction*/
	double cosT, sinT; /*variables for sin and cos of loading direction*/
	double mu, sig, muT, sigT; /*variables for means and standard devs of uncertain magnitude and direction*/
	double ex, ex2; /*varibales for expeonential components of pseudo loads*/
	double wgt1, wgt2; /*variables for load case weights*/
	
	int NumDof = 2 * NumNodes; /*number of dofs*/

	for(i=0;i<NumPforce;i++)
	{
		j = i * 6; /*place indicator for Pforce array*/
		/* read in force co-ordinates means and standard devs*/
		xtemp = Pforce[j];
		ytemp = Pforce[j+1];	/*location*/
		mu = Pforce[j+2];		/*mean magnitude*/
		sig = Pforce[j+3];		/*stand dev, magnitude*/
		muT = Pforce[j+4];		/*mean direction*/
		sigT = Pforce[j+5];		/*stand dev, direction*/
		
		/*---check mean direction to see if it's parallel with the grid---*/
		cosT = cos(muT);
		sinT = sin(muT);
		if(fabs(cosT) < 0.001)
		{
			/*force applied in y-direction*/
			dtemp = 2;
		}
		else if(fabs(sinT) < 0.001)
		{
			/*force applied in x-direction*/
			dtemp = 1;
		}
		else
		{
			/*force has both x & y components*/
			dtemp = 3;
		}
		
		/*calculate exponential values*/
		ftemp = -1.0 * sigT * sigT;
		ex = exp(ftemp);
		ex2 = ex * ex;
		/*printf("\nex=%f, ex2=%f",ex,ex2);*/
		
		/*calculate weights for this load*/
		wgt1 = (mu * mu) * (ex2 - ex);
		wgt1 += sig * sig * ex2;	/*weight for mean loading direction*/
		
		wgt2 = (mu * mu) + (sig * sig);
		wgt2 *= 0.5 * (1.0 - ex2);	/*weight for horizontal and vertical load directions*/
		
		/*printf("\nx=%lf, y=%lf, mu=%lf, sig=%lf, muT=%lf, sigT=%lf",xtemp,ytemp,mu,sig,muT,sigT);*/

		temp = closeNode(h,xtemp,ytemp,NumNodes,NodeCoord);	/*find closest node to applied force co-ordinates*/
		temp *= 2; /* X dof */
		
		/*----Add loads to first column of rhs----*/
		ftemp = mu * sqrt(ex);
		if(dtemp == 1)
		{
			unc_ind[*NumUn] = temp;
			unc_val[*NumUn] = (cosT < 0.0) ? -ftemp : ftemp; /*X component*/
			pnt = *NumUn;
			*NumUn = ++pnt;
		}
		else if(dtemp == 2)
		{
			unc_ind[*NumUn] = temp + 1;
			unc_val[*NumUn] = (sinT < 0.0) ? -ftemp : ftemp; /*Y component*/
			pnt = *NumUn;
			*NumUn = ++pnt;
		}
		else if(dtemp == 3)
		{
			unc_ind[*NumUn] = temp;
			unc_val[*NumUn] = ftemp * cos(muT); /*X component*/
			pnt = *NumUn;
			*NumUn = ++pnt;
			unc_ind[*NumUn] = temp + 1;
			unc_val[*NumUn] = ftemp * sin(muT); /*Y component*/
			pnt = *NumUn;
			*NumUn = ++pnt;
		}
			
		/*----Next load condition at mean direction----*/		
		if(fabs(wgt1) > 0.0001)
		{
			if(dtemp == 1)
			{	
				unc_ind[*NumUn] = temp + (*NumRhs * NumDof);
				unc_val[*NumUn] = 1.0; /*unit load in X direction*/
				pnt = *NumUn;
				*NumUn = ++pnt;
				unc_wgt[*NumRhs] = wgt1 + wgt2; /*combined weight*/
			}
			else if(dtemp == 2)
			{	
				unc_ind[*NumUn] = temp + 1 + (*NumRhs * NumDof);
				printf("\nunc_ind = %i, NumUn = %i",unc_ind[*NumUn],*NumUn);
				unc_val[*NumUn] = 1.0; /*unit load in X direction*/
				pnt = *NumUn;
				*NumUn = ++pnt;
				unc_wgt[*NumRhs] = wgt1 + wgt2; /*combined weight*/
			}
			else if(dtemp == 3)
			{
				
				unc_ind[*NumUn] = temp + (*NumRhs * NumDof);
				unc_val[*NumUn] = cosT; /*X component of a unit load*/
				pnt = *NumUn;
				*NumUn = ++pnt;
				
				unc_ind[*NumUn] = temp + 1 + (*NumRhs * NumDof);
				unc_val[*NumUn] = sinT; /*Y component of a unit load*/
				pnt = *NumUn;
				*NumUn = ++pnt;
				unc_wgt[*NumRhs] = wgt1;		
			}
			/*update number of rhs count*/
			pnt = *NumRhs;
			*NumRhs = ++pnt;
		
			/*additional loads for uncertain loading direction*/
			if(fabs(wgt2) > 0.0001)
			{
				if((dtemp == 1) || (dtemp == 3))
				{
					unc_ind[*NumUn] = temp + 1 + (*NumRhs * NumDof);
					unc_val[*NumUn] = 1.0; /*Y unit load*/
					pnt = *NumUn;
					*NumUn = ++pnt;
					unc_wgt[*NumRhs] = wgt2;
					/*update number of rhs count*/
					pnt = *NumRhs;
					*NumRhs = ++pnt;
				}
				if((dtemp == 2) || (dtemp == 3))
				{
					unc_ind[*NumUn] = temp + (*NumRhs * NumDof);
					unc_val[*NumUn] = 1.0; /*X unit load*/
					pnt = *NumUn;
					*NumUn = ++pnt;
					unc_wgt[*NumRhs] = wgt2;
					/*update number of rhs count*/
					pnt = *NumRhs;
					*NumRhs = ++pnt;
				}
			}
			
		}
	}

	/*-----------------------------------------This section applies Unifrom Forces to an Area-----------------------------------------*/
	int Xnum1, Xnum2, Ynum1, Ynum2; /*variables to track end node numbers*/
	short *nodeFixtemp;
	
	for(i=0;i<NumAforce;i++)
	{
		k = i * 8; /*place indicator for Aforce array*/
		
		mu = Aforce[k+4];	/*mean magnitude - force / unit length*/
		sig = Aforce[k+5];	/*stand dev, magnitude*/
		muT = Aforce[k+6];	/*mean direction*/
		sigT = Aforce[k+7]; /*stand dev, direction*/
		
		mu *= h;
		sig *= h; /*need to multiply magnitude values b edge length to get the total force / element */	
		
		/*---check mean direction to see if it's parallel with the grid---*/
		cosT = cos(muT);
		sinT = sin(muT);
		if(fabs(cosT) < 0.001)
		{
			/*force applied in y-direction*/
			dtemp = 2;
		}
		else if(fabs(sinT) < 0.001)
		{
			/*force applied in x-direction*/
			dtemp = 1;
		}
		else
		{
			/*force has both x & y components*/
			dtemp = 3;
		}
		
		/*work out extra number of required rhs' */
		num = 3; /*initalize to maximum*/
		if(fabs(sigT) < 0.00001)
		{
			num -= 2; /*minus 2 if direction is certain*/
			if(fabs(sig) < 0.00001)
			{
				num--; /*minus 1 if magnitude is certain*/
			}
		}
		else if(dtemp != 3)
		{
			num--; /*else if mean direction orthogonal to grid, then minus 1*/
		}
		
		/*calculate exponential values*/
		ftemp = -1.0 * sigT * sigT;
		ex = exp(ftemp);
		ex2 = ex * ex;
		/*printf("\nex=%f, ex2=%f",ex,ex2);*/
		
		/*calculate weights for this load*/
		wgt1 = (mu * mu) * (ex2 - ex);
		wgt1 += sig * sig * ex2;	/*weight for mean loading direction*/
		
		wgt2 = (mu * mu) + (sig * sig);
		wgt2 *= 0.5 * (1.0 - ex2);	/*weight for horizontal and vertical load directions*/
		
		nodeFixtemp = calloc(NumNodes,sizeof(short));	/*Tempary array to store node data for current area*/	
		
		ftemp = 1000000000.0; /*initalize minimum x-cord*/
		ftemp2 = 1000000000.0; /*initalize minimum y-cord*/
		xtemp = -1.0;	 /*initalize maximum x-cord*/
		ytemp = -1.0;	 /*initalize maximum y-cord*/
		
		for(j=0;j<NumNodes;j++) /*For all nodes - determine if the node lies within the fixed rectangular area*/
		{
			if( ((NodeCoord[j].x - Aforce[k+1]) < 0.000001) && ((NodeCoord[j].x - Aforce[k]) > -0.000001) )  /*If within x bounds*/
			{
				if( ((NodeCoord[j].y - Aforce[k+3]) < 0.000001) && ((NodeCoord[j].y - Aforce[k+2]) > -0.000001) )  /*If also within y bounds*/
				{
					nodeFixtemp[j] = 1; /*Node has a force applied*/
					
					/*section to find max and min node co-ords in area*/
					Xnum1 = (NodeCoord[j].x < ftemp) ? j : Xnum1; /*update min x node number if required*/
					ftemp = (NodeCoord[j].x < ftemp) ? NodeCoord[j].x : ftemp; /*check minimum x-cord*/
					Ynum1 = (NodeCoord[j].y < ftemp2) ? j : Ynum1; /*update min y node number if required*/
					ftemp2 = (NodeCoord[j].y < ftemp2) ? NodeCoord[j].y : ftemp2; /*check minimum y-cord*/
					Xnum2 = (NodeCoord[j].x > xtemp) ? j : Xnum2; /*update max x node number if required*/
					xtemp = (NodeCoord[j].x > xtemp) ? NodeCoord[j].x : xtemp; /*check maximum x-cord*/
					Ynum2 = (NodeCoord[j].y > ytemp) ? j : Ynum2; /*update max y node number if required*/
					ytemp = (NodeCoord[j].y > ytemp) ? NodeCoord[j].y : ytemp; /*check maximum y-cord*/
					
					/*printf("\nForce applied to Node %i",j+1);*/
				}
			}
		}
		
		/*printf("\nftemp=%lf, ftemp2=%lf, xtemp=%lf, ytemp=%lf",ftemp,ftemp2,xtemp,ytemp);
		printf("\nXnum1=%i, Xnum2=%i, Ynum1=%i, Ynum2=%i",Xnum1,Xnum2,Ynum1,Ynum2);*/
		
		/*define the end nodes*/
		if(fabs(ftemp - xtemp) < h) /*if difference between max and min x-cords is < h, then use y numbers instead*/
		{
			/*first check y co-ord difference isn't < h*/
			if(fabs(ftemp2 - ytemp) >= h)
			{
				Xnum1 = Ynum1;
				Xnum2 = Ynum2; /*store end node numbers as Xnum's */
			}
			else
			{
				printf("\nERROR! Force area smaller than h x h !!!");
			}
		}
		/*printf("\nXnum1=%i, Xnum2=%i",Xnum1,Xnum2);*/
		
		for(j=0;j<NumNodes;j++) /*for all nodes*/
		{
			if(nodeFixtemp[j] != 0) /*If current node is effected by the force add force to rhs array*/
			{
				ytemp = ((j == Xnum1) || (j == Xnum2)) ? 0.5 : 1.0; /*factor to scale loads on end nodes*/
				temp = 2 * j; /*dof indicator*/
				
					/*----add loads to first column of rhs----*/
					ftemp = mu * sqrt(ex) * ytemp;
					if(dtemp == 1)
					{
						unc_ind[*NumUn] = temp;
						unc_val[*NumUn] = (cosT < 0.0) ? -ftemp : ftemp; /*X component*/
						pnt = *NumUn;
						*NumUn = ++pnt;
					}
					else if(dtemp == 2)
					{
						unc_ind[*NumUn] = temp + 1;
						unc_val[*NumUn] = (sinT < 0.0) ? -ftemp : ftemp; /*Y component*/
						pnt = *NumUn;
						*NumUn = ++pnt;
					}
					else if(dtemp == 3)
					{
						unc_ind[*NumUn] = temp;
						unc_val[*NumUn] = ftemp * cos(muT); /*X component*/
						pnt = *NumUn;
						*NumUn = ++pnt;
						unc_ind[*NumUn] = temp + 1;
						unc_val[*NumUn] = ftemp * sin(muT); /*Y component*/
						pnt = *NumUn;
						*NumUn = ++pnt;
					}
						
					/*----next load condition at mean direction----*/
					if(fabs(wgt1) > 0.0001)
					{
						if(dtemp == 1)
						{
							unc_ind[*NumUn] = temp + (*NumRhs * NumDof);
							unc_val[*NumUn] = ytemp; /*X unit load*/
							pnt = *NumUn;
							*NumUn = ++pnt;
							unc_wgt[*NumRhs] = wgt1 + wgt2; /*combined weight*/
						}
						else if(dtemp == 2)
						{
							unc_ind[*NumUn] = temp + 1 + (*NumRhs * NumDof);
							unc_val[*NumUn] = ytemp; /*Y unit load*/
							pnt = *NumUn;
							*NumUn = ++pnt;
							unc_wgt[*NumRhs] = wgt1 + wgt2; /*combined weight*/
						}
						else if(dtemp == 3)
						{
							unc_ind[*NumUn] = temp + (*NumRhs * NumDof);
							unc_val[*NumUn] = ytemp * cos(muT); /*X component of unit load*/
							pnt = *NumUn;
							*NumUn = ++pnt;
							unc_ind[*NumUn] = temp + 1 + (*NumRhs * NumDof);
							unc_val[*NumUn] = ytemp * sin(muT); /*Y component of unit laod*/
							pnt = *NumUn;
							*NumUn = ++pnt;
							unc_wgt[*NumRhs] = wgt1;
						}
						
						/*update number of rhs count*/
						pnt = *NumRhs;
						*NumRhs = ++pnt;
					
						/*additional loads for uncertain loading direction*/
						if(fabs(wgt2) > 0.0001)
						{
							if((dtemp == 1) || (dtemp == 3))
							{
								unc_ind[*NumUn] = temp + 1 + (*NumRhs * NumDof);
								unc_val[*NumUn] = ytemp; /*Y unit load*/
								pnt = *NumUn;
								*NumUn = ++pnt;
								unc_wgt[*NumRhs] = wgt2;
								/*update number of rhs count*/
								pnt = *NumRhs;
								*NumRhs = ++pnt;
							}
							if((dtemp == 2) || (dtemp == 3))
							{
								unc_ind[*NumUn] = temp + (*NumRhs * NumDof);
								unc_val[*NumUn] = ytemp; /*X unit load*/
								pnt = *NumUn;
								*NumUn = ++pnt;
								unc_wgt[*NumRhs] = wgt2;
								/*update number of rhs count*/
								pnt = *NumRhs;
								*NumRhs = ++pnt;
							}
						}
						/*need to reset number of rhs' for next node, to ensure dofs for next node are correct*/
						pnt = *NumRhs;
						*NumRhs = pnt-num;
					}
			}
		}
		
		/*finally add back the extra rhs' required for this load*/
		pnt = *NumRhs;
		*NumRhs = pnt+num;
		
		free(nodeFixtemp);	/*Free some memory*/
	}
}

/*function to compute "Co-Variance" matrix for variance minimisation*/
void VarMat(double h, int NumNodes, Coord NodeCoord[NumNodes], int *NumVar, int *var_ind, double *var_val, 
				int NumPforce, double *Pforce, int NumAforce, double *Aforce)
{
	int i,j,k,temp;	/*incrementors etc.*/
	double mu,sig,muT,sigT; /*statistical data*/
	double cosT,sinT,cos2,sin2;	/*tri values*/
	
	printf("\n\nWarning!! Variance min is only approximate for loading direction uncertainty");
	
	/*--------------First compute values for single point laods-------------*/
	for(i=0;i<NumPforce;i++)
	{
		j = i * 6; /*place indicator for Pforce array*/
		/* read in force means and standard devs*/
		mu = Pforce[j+2] * Pforce[j+2];		/*mean magnitude, squared*/
		sig = Pforce[j+3] * Pforce[j+3];	/*stand dev, magnitude, squared*/
		muT = Pforce[j+4];					/*mean direction*/
		sigT = Pforce[j+5] * Pforce[j+5];	/*stand dev, direction, squared*/
		
		/*---compute tri function values and squares---*/
		cosT = cos(muT);
		sinT = sin(muT);
		cos2 = cosT * cosT;
		sin2 = sinT * sinT;
		
		var_ind[*NumVar] = closeNode(h,Pforce[j],Pforce[j+1],NumNodes,NodeCoord);	/*find closest node to applied force co-ordinates*/
		
		/*----Compute "Co-Variance" matrix enrites for cuurent node----*/
		j = 3 * *NumVar;	/*place indicator for var_val array*/
		var_val[j] = 2.0 * ( (sin2 * mu * sigT) + (cos2 * sig) );
		var_val[j+1] = 2.0 * ( (cos2 * mu * sigT) + (sin2 * sig) );
		var_val[j+2] = 2.0 * ( cosT * sinT * ( sig - (mu * sigT) ) ); /*times 2 for actual sensitivities*/
		temp = *NumVar;
		*NumVar = ++temp;
	}
	
	/*--------------Now compute values for distributed laods-------------*/
	int Xnum1, Xnum2, Ynum1, Ynum2; /*variables to track end node numbers*/
	double ftemp,ftemp2,xtemp,ytemp;
	double fact1,fact2,fact3; /*factors for co-variance matrix*/
	short *nodeFixtemp;
	
	for(i=0;i<NumAforce;i++)
	{
		k = i * 8; /*place indicator for Aforce array*/
		
		mu = Aforce[k+4] * Aforce[k+4];		/*mean magnitude - force / unit length, squared*/
		sig = Aforce[k+5] * Aforce[k+5];	/*stand dev, magnitude, squared*/
		muT = Aforce[k+6];					/*mean direction*/
		sigT = Aforce[k+7] * Aforce[k+7];	/*stand dev, direction, squared*/
		
		mu *= h * h;
		sig *= h * h; /*need to multiply magnitude values by edge length (squared) to get the total force / element */	
		
		/*---compute tri function values and squares---*/
		cosT = cos(muT);
		sinT = sin(muT);
		cos2 = cosT * cosT;
		sin2 = sinT * sinT;
		
		/*compute co-variance entries from statistical data*/
		fact1 = 2.0 * ( (sin2 * mu * sigT) + (cos2 * sig) );
		fact2 = 2.0 * ( (cos2 * mu * sigT) + (sin2 * sig) );
		fact3 = 2.0 * ( cosT * sinT * ( sig - (mu * sigT) ) ); /* times 2 for adjoint forces */
		
		nodeFixtemp = calloc(NumNodes,sizeof(short));	/*Tempary array to store node data for current area*/	
		
		ftemp = 1000000000.0; /*initalize minimum x-cord*/
		ftemp2 = 1000000000.0; /*initalize minimum y-cord*/
		xtemp = -1.0;	 /*initalize maximum x-cord*/
		ytemp = -1.0;	 /*initalize maximum y-cord*/
		
		for(j=0;j<NumNodes;j++) /*For all nodes - determine if the node lies within the fixed rectangular area*/
		{
			if( ((NodeCoord[j].x - Aforce[k+1]) < 0.000001) && ((NodeCoord[j].x - Aforce[k]) > -0.000001) )  /*If within x bounds*/
			{
				if( ((NodeCoord[j].y - Aforce[k+3]) < 0.000001) && ((NodeCoord[j].y - Aforce[k+2]) > -0.000001) )  /*If also within y bounds*/
				{
					nodeFixtemp[j] = 1; /*Node has a force applied*/
					
					/*section to find max and min node co-ords in area*/
					Xnum1 = (NodeCoord[j].x < ftemp) ? j : Xnum1; /*update min x node number if required*/
					ftemp = (NodeCoord[j].x < ftemp) ? NodeCoord[j].x : ftemp; /*check minimum x-cord*/
					Ynum1 = (NodeCoord[j].y < ftemp2) ? j : Ynum1; /*update min y node number if required*/
					ftemp2 = (NodeCoord[j].y < ftemp2) ? NodeCoord[j].y : ftemp2; /*check minimum y-cord*/
					Xnum2 = (NodeCoord[j].x > xtemp) ? j : Xnum2; /*update max x node number if required*/
					xtemp = (NodeCoord[j].x > xtemp) ? NodeCoord[j].x : xtemp; /*check maximum x-cord*/
					Ynum2 = (NodeCoord[j].y > ytemp) ? j : Ynum2; /*update max y node number if required*/
					ytemp = (NodeCoord[j].y > ytemp) ? NodeCoord[j].y : ytemp; /*check maximum y-cord*/
					
					/*printf("\nForce applied to Node %i",j+1);*/
				}
			}
		}
		
		/*define the end nodes*/
		if(fabs(ftemp - xtemp) < h) /*if difference between max and min x-cords is < h, then use y numbers instead*/
		{
			/*first check y co-ord difference isn't < h*/
			if(fabs(ftemp2 - ytemp) >= h)
			{
				Xnum1 = Ynum1;
				Xnum2 = Ynum2; /*store end node numbers as Xnum's */
			}
			else
			{
				printf("\nERROR! Force area smaller than h x h !!!");
			}
		}
		/*printf("\nXnum1=%i, Xnum2=%i",Xnum1,Xnum2);*/
		
		for(j=0;j<NumNodes;j++) /*for all nodes*/
		{
			if(nodeFixtemp[j] != 0) /*If current node is effected by the force add force to rhs array*/
			{
				ytemp = ((j == Xnum1) || (j == Xnum2)) ? 0.5 : 1.0; /*factor to scale loads on end nodes*/				
				var_ind[*NumVar] = j;	/*store node affected by distributed load*/
		
				/*----Compute "Co-Variance" matrix enrites for curent node----*/
				k = 3 * *NumVar;	/*place indicator for var_val array*/
				var_val[k] = ytemp * fact1;
				var_val[k+1] = ytemp * fact2;
				var_val[k+2] = ytemp * fact3;
				/*finally update the number of nodes with uncertain loads*/
				temp = *NumVar;
				*NumVar = ++temp;
			}
		}

		free(nodeFixtemp);	/*Free some memory*/
	}
}
