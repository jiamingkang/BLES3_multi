/*
 *  BCs.h
 *  
 *  Created by Peter Dunning on 22/06/2010.
 *  functions to apply boundary conditions and construct required force vectors
 *
 */
 
#include "ls_types.h"
#include <stdio.h>

/*function that finds all zero dispacement boundary conditions and assembles into one array */
void FixDisp(double h, int NumNodes, Coord NodeCoord[NumNodes], int *NumFix, int *fixed, int *TotFix, int *FixDofs, 
				int *LAfix, Aux *Afix, FILE *infile);
				
/*function that assembles the load vector for Objective (1), minimisation of deterministic strain energy*/
void LoadVector(double h, int NumNodes, Coord NodeCoord[NumNodes], int *NumForces, int *rhs_ind, double *rhs_val,
				int NumPforce, double *Pforce, int NumAforce, double *Aforce);

/*function to compute load vectors and weights required to calculate expected strain energy (compliance)*/
void UncLoads(double h, int NumNodes, Coord NodeCoord[NumNodes], int *NumUn, int *unc_ind, double *unc_val, 
				int *NumRhs, double *unc_wgt, int NumPforce, double *Pforce, int NumAforce, double *Aforce);
				
/*function to compute "Co-Variance" matrix for variance minimisation*/
void VarMat(double h, int NumNodes, Coord NodeCoord[NumNodes], int *NumVar, int *var_ind, double *var_val, 
				int NumPforce, double *Pforce, int NumAforce, double *Aforce);
