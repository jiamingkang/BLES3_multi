/*
 *  IFG.c
 *  
 *	Created by Peter Dunning
 *
 *	Version 1.2	at 19/01/2009
 *	Modified for BLES program
 *  Fully working for displacement calculation
 *	FUlly working for Strain Energy calculation
 *	Based on FG2D code
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "IFG.h"
#include "FixedGrid.h"
#include "EMatrix.h"
#include "KMatrix.h"
#include "ls_types.h"
#include "Solve.h"
#include "Strain.h"

void IFG_Status(double h,int elemX,int elemY,Elem Number[elemX][elemY],int *Acount,int NumElem,int NodeX, int NodeY, int Nodes2[NodeX][NodeY],
				int NumNodes,Coord NodeCoord[NumNodes],double *lsf,short *NodeStat,short *ElemStat,Aux *auxNode,double ftv_lim,int itt,double hxz[elemX],double hz[elemX], Aux *AuxNodes, Bseg *bptr, int *na_conn_ind, int *na_conn, int *NumAux, int *NumBound)
{

	printf("\nIN IFG STATUS\n");
	int n,m,o;	/*Incrementors*/
	int temp,txtemp,num; /*Temperary variables*/
	int count,begin; /*Incrementors*/
	double ftemp,ftemp2,htemp;
	double lsf1, lsf2;
	FILE *outfile;	/*File varible for output files*/
	char plotname[40];	/*variable to change names of plotting output files*/
	
	/*-----Calculate node and elements status, add auxillary nodes-----*/
	/*ftemp = (ftv_lim > 0.0) ? 0.001 : 0.001; /*lower limit for a boundary node (1% or 0.1% element edge length)*/

	/*Determine Node status*/
	for(n=0;n<NumNodes;n++)
	{
		/*If signed distance function is 0.0 (NB: Inside within the fit to vertex limit), then node lies on the boundary*/
		if((lsf[n] > -0.000001) && (lsf[n] < 0.000001))
		{
			NodeStat[n] = 2;		
		}
		/*Otherwise if lsf = -ve node is OUT, +ve then node is IN*/
		else
		{
			NodeStat[n] = (lsf[n] > 0.0) ? 1 : 0;
		}
	}


	/*New bit for Element status determination and local node data*/
	int n4; /*varible for 4 x num*/
	short sumO, sumB, sumI;	/*Variables for summing IN, OUT and Boundary nodes*/
	int Bcount = 0; /*Variable to count number of boundary elements*/
	int Lnodes[4];	/*Array to temp store local node data*/
	
	double *fnodes;
	double lsf_min,lsf_val; /*varibales for min node lsf value*/
	
	/*Need to check that each element is only cut once by the implicit boundary*/
	int Ncount = 0; /*variable to count no. nodes that need changing to status = 2*/
	int *chNodes;
	chNodes = malloc(NumNodes * sizeof(int)); /*array to store node numbers that will be changed to status = 2*/
	fnodes = malloc(4 * sizeof(double)); /*array to store actual distances from boundary using interpolation*/
					
	for(m=0;m<elemY;m++)
	{
		for(n=0;n<elemX;n++)	/*For all elements*/
		{
			num = Number[n][m].n - 1;	/*Element number, -1 for C array storage*/
			n4 = 4 * num;
			/*Read in local node numbers, -1 for C array storage*/
			Lnodes[0] = Number[n][m].a - 1;
			Lnodes[1] = Number[n][m].b - 1;
			Lnodes[2] = Number[n][m].c - 1;
			Lnodes[3] = Number[n][m].d - 1;

			for(o=0;o<4;o++) {
				fnodes[o] = 100.0; /*initialize to a large number*/
				}
			
			count = 0;
			for(o=0;o<4;o++)
			{
				temp = (o == 3) ? 0 : (o + 1); /*find next node round*/
				/*If edge is cut (1 node is IN other is OUT) then update count*/
				if((NodeStat[Lnodes[o]] + NodeStat[Lnodes[temp]]) == 1)
				{
						count++;
						lsf1 = lsf[Lnodes[o]];
						lsf2 = lsf[Lnodes[temp]];
						ftemp = fabs(lsf1 / (lsf1 - lsf2));
						ftemp2 = fabs(lsf2 / (lsf1 - lsf2));
						fnodes[o] = (ftemp < fnodes[o]) ? ftemp : fnodes[o];
						fnodes[temp]  = (ftemp2 < fnodes[temp]) ? ftemp2 : fnodes[temp];
						/*printf("\nElement = %i", num+1);
						printf("\nCount = %i", count);*/
					}
			}
			
			if(count > 2) /*if more than 2 edges are cut then need to move boundary so that element is only cut once*/
			{
				printf("\nElement %i has %i edges cut/n",num+1,count);
				temp = -1;
				lsf_min = 10.0; /*set inital minimum to a large number*/
				for(o=0;o<4;o++)
				{
					/*printf("\tfnodes[%i] = %f",o,fnodes[o]);*/
					if( ((NodeStat[Lnodes[o]] != 2)) && (fnodes[o] < lsf_min) )
					{
						lsf_min = fnodes[o];
						temp = 1;
						printf("\nLnodes[%i] = %i", o, Lnodes[o]);
						printf("\nNodeStat[%i] = %i", Lnodes[o], NodeStat[Lnodes[o]]);
						printf("\nfnodes[%i] = %f", o, fnodes[o]);
						printf("\nlsf_min = %f\n", lsf_min);
					}
				}
				
				if(temp == 1)
				{
					for(o=0;o<4;o++)
					{
						if( ((NodeStat[Lnodes[o]] != 2)) && ((fnodes[o] - lsf_min) < 0.000001) )
						{
							chNodes[Ncount++] = Lnodes[o]; /*add node to array change array*/
							/*printf("\nNode %i, at [%lf, %lf], changed to 2 from %i",Lnodes[o]+1,NodeCoord[Lnodes[o]].x,NodeCoord[Lnodes[o]].y, NodeStat[Lnodes[o]]);*/
						}
					}
				}
			}
		}
	}

	free(fnodes);

		for(n=0;n<Ncount;n++)
	{
		temp = chNodes[n];
		NodeStat[temp] = 2; /*Make nodes with smallest lsf values into boundary nodes*/
		/*printf("\nNodeStat[%i] = %i", temp, NodeStat[temp]);*/
	}
	
	free(chNodes);
					
	for(m=0;m<elemY;m++)
	{
		for(n=0;n<elemX;n++)	/*For all elements*/
		{
			num = Number[n][m].n - 1;	/*Element number, -1 for C array storage*/
			/*Read in local node numbers, -1 for C array storage*/
			Lnodes[0] = Number[n][m].a - 1;
			Lnodes[1] = Number[n][m].b - 1;
			Lnodes[2] = Number[n][m].c - 1;
			Lnodes[3] = Number[n][m].d - 1;
			
            //jhbhjbjhbhj
			/*Set node sums back to zero*/
			sumO = 0;
			sumI = 0;
			sumB = 0;
			
			/*Find number of OUT and Boundary nodes for the element*/
			for(o=0;o<4;o++)
			{
				sumO += (NodeStat[Lnodes[o]] == 0) ? 1 : 0;
				sumI += (NodeStat[Lnodes[o]] == 1) ? 1 : 0;
				sumB += (NodeStat[Lnodes[o]] == 2) ? 1 : 0;
			}
			
			if(sumO == 0)	/*If No nodes are OUT then Element is IN*/
			{
				ElemStat[num] = 0;/*4;*/
			}
				
			else if (sumI == 0)	/*If no nodes are IN then Element is OUT*/
			{
				ElemStat[num] = 0;
			}
			
			else	/*element is cut*/
			{
				ElemStat[num] = 1;
			}
		}
	}

	/*printf("\nElement Status Determined");	
	
	printf("\nBcount = %i",Bcount);*/
	
    // Discretizr the boundary///////////////////////////////////////////////////////////////////////////
    double tol=1.0e-6;
    int *na_count = calloc(NumNodes, sizeof(int)); // count number of aux nodes connected to an element

    int ncnt;
    int count2;
    int Acount_b;
    double Ax, Ay;
    int temp2;
    count = 0;	// initialize count of auxiallry nodes
    count2 = 0; // initialize count of boundary segments
    temp=0;
    int An[4]; // Array to store intersection node numbers

    
    
    // For all elements
    for(m=0;m<elemY;m++)
    {
        for(n=0;n<elemX;n++)
        {
            num = Number[n][m].n-1; //Element number
            if(ElemStat[num] != 0)
            {
                Acount_b = 0;
                // Read in local node numbers
                Lnodes[0] = Number[n][m].a-1;
                Lnodes[1] = Number[n][m].b-1;
                Lnodes[2] = Number[n][m].c-1;
                Lnodes[3] = Number[n][m].d-1;
                
                // Look at each edge in turn to determine if cut, or part of the boundary
                for(o=0;o<4;o++)
                {
                    temp = (o == 3) ? 0 : (o + 1); // find next node round
                    sumI = NodeStat[Lnodes[o]] + NodeStat[Lnodes[temp]]; // sum node status
                    
                    // If edge is cut (1 node is IN other is OUT) then update count
                    if(sumI == 1)
                    {
                        // Compute co-ordintes of interection point (interpolation)
                        lsf1 = lsf[Lnodes[o]];
                        lsf2 = lsf[Lnodes[temp]];
                        ftemp = (lsf1 + lsf2) / (lsf1 - lsf2);
                        ftemp += 1.0;
                        ftemp *= 0.5 * h;
                        
                        if( (ftemp > h) || (ftemp < 0.0) ) {
                            printf("ERROR! Boundary distance = %lf, near node %i",ftemp,Lnodes[o]);
                        }
                        
                        // co-ordinates of the aux node, depends on the cut edge
                        switch(o) {
                            case 0:
                                Ax = NodeCoord[Lnodes[o]].x + ftemp;
                                Ay = NodeCoord[Lnodes[o]].y;
                                break;
                            case 1:
                                Ax = NodeCoord[Lnodes[o]].x;
                                Ay = NodeCoord[Lnodes[o]].y + ftemp;
                                break;
                            case 2:
                                Ax = NodeCoord[Lnodes[o]].x - ftemp;
                                Ay = NodeCoord[Lnodes[o]].y;
                                break;
                            case 3:
                                Ax = NodeCoord[Lnodes[o]].x;
                                Ay = NodeCoord[Lnodes[o]].y - ftemp;
                                break;
                        }
                        
                        // Add the auxillary boundary node
                        // need to check that the auxillary node doesn't already exist
                        An[Acount_b] = -1;
                        int end = na_count[Lnodes[o]];
                        if(end > 0)
                        {
                            ncnt = 4*Lnodes[o];
                            while(end > 0)
                            {
                                if( (fabs(Ax - AuxNodes[na_conn[ncnt]].x) < tol)
                                   && (fabs(Ay - AuxNodes[na_conn[ncnt]].y) < tol) )
                                {
                                    An[Acount_b] = na_conn[ncnt];
                                    end = 0;
                                }
                                ncnt++;
                                end--;
                            }
                        }
                        
                        // If auxillary node does not exist, then create it
                        if(An[Acount_b] == -1)
                        {
                            // add grid node -> aux node connectivity (for both nodes)
                            ncnt = na_count[ Lnodes[o] ]; // number connected aux nodes so far
                            ncnt += Lnodes[o]*4; // point to correct location
                            na_conn[ncnt] = count; // add connectivity
                            na_count[ Lnodes[o] ]++; // increase count
                            
                            ncnt = na_count[ Lnodes[temp] ]; // number connected aux nodes so far
                            ncnt += Lnodes[temp]*4; // point to correct location
                            na_conn[ncnt] = count; // add connectivity
                            na_count[ Lnodes[temp] ]++;  // increase count
                            
                            AuxNodes[count].x = Ax;
                            AuxNodes[count].y = Ay;
                            AuxNodes[count].n = count+NumNodes+1;
                            An[Acount_b] = count;
                            count++; // increase Aux node count
                        }
                        
                        Acount_b++; //update count of boundary intersection points for this element
                    }
                    
                    // If both nodes are on the boundary
                    else if(sumI == 4)
                    {
                        // Add the elemet edge to boundary segment data
                        bptr[count2].n1 = Lnodes[o];    // node number of point 1
                        bptr[count2].n2= Lnodes[temp]; // node number of point 2
                        bptr[count2++].e = num;			// associated elememt number
                    }
                }
                
                // For cut elements determine boundary segment(s)
                
                // if there are two edges cut, then a boundary segment must cross both
                if(Acount_b == 2)
                {
                    bptr[count2].n1 = NumNodes + An[0];	// node number of point 1
                    bptr[count2].n2 = NumNodes + An[1];	// node number of point 2
                    bptr[count2++].e = num;				// associated elememt number
                }
                
                // if there is only one cut edge, then the boundary must also cross an element node
                else if(Acount_b == 1)
                {
                    // find a node that is on the boundary & has an OUT neighbour
                    for(o=0;o<4;o++)
                    {
                        // If node is on boudnary, check its neighbours
                        if(NodeStat[Lnodes[o]] == 2)
                        {
                            temp = (o == 3) ? 0 : (o + 1);   // find next node round
                            temp2 = (o == 0) ? 3 : (o - 1);  // find previous node
                            
                            // If a neighbour is an OUT node, then add boundary segment
                            if( (NodeStat[Lnodes[temp]] == 0) || (NodeStat[Lnodes[temp2]] == 0) )
                            {
                                bptr[count2].n1 = NumNodes + An[0];	// node number of point 1
                                bptr[count2].n2 = Lnodes[o];		// node number of point 2
                                bptr[count2++].e = num;				// associated elememt number
                            }
                        }
                    }
                }
                
                // if there are four cut edges, then determine which Aux node pairs form the boundary
                else if(Acount_b == 4)
                {
                    // first determine lsf value at element centre - sign is only of interest
                    ftemp = 0.0;
                    for(o=0;o<4;o++) {
                        ftemp += lsf[Lnodes[o]];
                    }
                    
                    // Now look at status of node 1
                    temp = NodeStat[Lnodes[0]];
                    
                    // Node pairs that correspond to boundary segments can easily be determined
                    if( ( (temp == 1) && (ftemp > 0.0) ) || ( (temp == 0) && (ftemp < 0.0) ) )
                    {
                        bptr[count2].n1 = NumNodes + An[0];		// node number of point 1
                        bptr[count2].n2 = NumNodes + An[1];		// node number of point 2
                        bptr[count2++].e = num;					// associated elememt number
                        
                        bptr[count2].n1 = NumNodes + An[2];		// node number of point 1
                        bptr[count2].n2 = NumNodes + An[3];		// node number of point 2
                        bptr[count2++].e = num;					// associated elememt number
                    }
                    
                    else
                    {
                        bptr[count2].n1 =  NumNodes + An[0];    // node number of point 1
                        bptr[count2].n2 =  NumNodes + An[3];	// node number of point 2
                        bptr[count2++].e = num;					// associated elememt number
                        
                        bptr[count2].n1 =  NumNodes + An[1];    // node number of point 1
                        bptr[count2].n2 =  NumNodes + An[2];	// node number of point 2
                        bptr[count2++].e = num;					// associated elememt number
                    }
                    
                    // update element status to indicate if centre of element is IN or OUT
                    //  4 = centre IN, 5 = centre OUT
                    ElemStat[num] = (ftemp > 0.0) ? 4:5;
                    
                }
                
                // if no edges are cut and element is not IN
                //  then boundary segment must cross diagonal
                else if( (Acount_b == 0) && (ElemStat[num] != 10) )
                {
                    // find the two boudnary nodes
                    for(o=0;o<4;o++)
                    {
                        if( NodeStat[Lnodes[o]] == 2)
                        {
                            An[Acount_b] = Lnodes[o];
                            Acount_b++;
                        }
                    }
                    bptr[count2].n1 = An[0];    // node number of point 1
                    bptr[count2].n2 = An[1];	// node number of point 2
                    bptr[count2++].e = num;		// associated elememt number
                }
            }
        }
    }
    
    
    	// read back the totals for NumAux & NumBound

    *NumAux=count;
    *NumBound=count2;
   // AuxNodes = realloc(AuxNodes, (count * sizeof(Aux)));
    bptr = realloc(bptr, (count2 * sizeof(Bseg)));
    // re-size auxillary node array and boundary segment data array to min size
    
    na_conn = realloc(na_conn, (2*count*sizeof(int)));
    
    
    count = 0;
    for(n=0;n<NumNodes;n++)
    {
        temp2 = na_count[n];
        na_conn_ind[n]=count;
        if(temp2 > 0)
        {
            temp = 4*n; // point ot start of connected aux nodes
            while(temp2 > 0)
            {
                na_conn[count++] = na_conn[temp++];
                temp2--;
            }
        }
    }
    na_conn_ind[NumNodes] = count; // end point
    free(na_count);
    free(na_conn);

    
    
	int *NIOnums;
	NIOnums = calloc(Bcount,sizeof(int)); /*Array to store element numbers of NIO elements*/
	count = 0;
	free(NIOnums);
	//free(Lnodes);
	
}


void IFG_Solve(int *irn, int *jcn, double *A, int TotFix, int *FixDofs, double *rhs_in, int NumEntries, int MatrixOrder, short *redun, int NumRhs, int itt)
{

    
    
    
	FILE *outfile;	/*File varible for output files*/
	char plotname[40];
	int n,m,o,txtemp,count,temp,temp2;
	/*define some new index arrays by dynamic memory allocation*/
	int *irn_in, *jcn_in;
	double *A_in;
	irn_in = calloc(NumEntries,sizeof(int));
	jcn_in = calloc(NumEntries,sizeof(int));
	A_in = calloc(NumEntries,sizeof(double));
	
	short include; /*variable to determine if an entry is to be included*/
	
	int *FixDofs2;
	FixDofs2 = malloc(MatrixOrder * sizeof(int)); /*initally set to maximum possible*/
	FixDofs2[0] = 0; /*set inital value*/
	
	/*need to remove redundant dofs*/
	count = 1; /*variable to count along FixDofs*/
	txtemp = 1; /*variable to count along FixDofs2*/
	
	for(n=0;n<MatrixOrder;n++)
	{
		if(FixDofs[count] == n+1) /*if dof already fixed then simply add to new array*/
		{
			FixDofs2[txtemp++] = n+1;
			/*printf("\ndof %i fixed by IFG_Solve",n);
			printf("\t FixDof[%i] = %i", count, FixDofs[count]);*/
			count++; /*look for next fixed dof*/
		}
		
		else if(redun[n] == 0)
		{
			FixDofs2[txtemp++] = n+1;
			/*printf("\ndof %i redundant in IFG_Solve",n);*/
		}
	}
	
	printf("\nFixed or removed dofs = %i",txtemp-1);
	FixDofs2[txtemp] = 10 * NumRhs * MatrixOrder; /*last entry in array to above maximum*/
	int TotFix2 = txtemp+1; /*set length of new array*/
	/*printf("\ncheck = %i, %i... %i, %i",FixDofs2[0],FixDofs2[1],FixDofs2[TotFix2-2],FixDofs2[TotFix2-1]);*/
	txtemp=0;

	for(n=0;n<NumEntries;n++)
	{
		include = 1;
		
		/*Find the number of fixed dofs upto the row number, of fix if necessary*/
		for(m=0;m<(TotFix2-1);m++)
		{	
			/*printf("\nirn=%i, up=%i, dn=%i",irn[n],FixDofs[m+1],FixDofs[m]);*/
			
			if(irn[n] == FixDofs2[m])
			{
				include = 0;
				break;
			}
			
			if((irn[n] < FixDofs2[m+1]) && (irn[n] > FixDofs2[m]))
			{
				break;
			}
			
		}
		
		/*Find the number of fixed dofs upto the column number, of fix if necessary*/
		if(include == 1)
		{
			for(o=0;o<(TotFix2-1);o++)
			{	
				if(jcn[n] == FixDofs2[o])
				{
					include = 0;
					break;
				}
				
				else if((jcn[n] < FixDofs2[o+1]) && (jcn[n] > FixDofs2[o]))
				{
					break;
				}
				
			}
		}
		
		if(include == 1)
		{
			A_in[txtemp] = A[n];
			irn_in[txtemp] = irn[n] - m;
			jcn_in[txtemp] = jcn[n] - o;
			txtemp++;
		}
	}

	/*Update the number of index entries*/
	NumEntries = txtemp;
	/*printf("\nNumEntries Reduced to %i",NumEntries);*/

	/*Reallocate memory for index arrays*/
	irn_in = realloc(irn_in,NumEntries * sizeof(int));
	jcn_in = realloc(jcn_in,NumEntries * sizeof(int));
	A_in = realloc(A_in,NumEntries * sizeof(double));

	/*Now resize the Force (rhs) array*/
	double *rhs;
	int mtemp = MatrixOrder - TotFix2 + 2; /*column length of the rhs matix sent to MA57 solver*/
	rhs = calloc((mtemp * NumRhs),sizeof(double)); /*rhs defined using dynamic memory allocation*/

	count = 1;
	m=0;
	for(n=0;n<MatrixOrder;n++)
	{
		txtemp = FixDofs2[count]-1;
		if(n != txtemp)		/*If dof not fixed copy accross the value*/
		{
			for(o=0;o<NumRhs;o++)
			{
				temp = o * mtemp; /*to get to correct place in rhs matrix*/
				temp2 = o * MatrixOrder; /*to get to correct place in rhs_in matrix*/
				rhs[m+temp] = rhs_in[n+temp2];
				/*if(rhs[m+temp] != 0){printf("\nrhs[%i] = %f\n", m+temp, rhs[m+temp]);}*/
			}
			m++;
		}
		
		else {				/*Else start looking at next fixed dof*/
			count++;
		}
	}
	

	printf("\nNumber of degrees of freedom to be solved = %i",MatrixOrder);

	/* Write index information file*/ 
	outfile = fopen("Index.txt", "w");
	if(outfile == NULL){
		printf("\nFailed to open Index writefile\n");
		}
	else{
	fprintf(outfile,"entry\tirn\tjrn\tK\t\n"); /*column headings*/
	/*fprintf(outfile,"%i %i",MatrixOrder, NumEntries);*/
		for(n=0;n<NumEntries;n++)
		{
			fprintf(outfile,"%i",(n+1));
			fprintf(outfile,"\t");
			fprintf(outfile," %i\t",irn_in[n]);
			fprintf(outfile,"%i\t",jcn_in[n]);
			fprintf(outfile,"%lf\t",A_in[n]);

			fprintf(outfile,"\n");
		}
		for(n=0;n<mtemp;n++)
		{
			fprintf(outfile,"%lf\n ",rhs[n]);
		}
	}

	fclose(outfile);

	printf("\nIndex Information file Written\n\n");
	/*double X,Y;
	X = 7/0;*/
	/*------------------------------------------------------------------------------------------------------------
	/
	/		Section 8: Solve the Finite Element Equation for Nodal Displacements
	/
	/			Inputs: * Global Stiffness Matrix Order and Force Vector (rhs) array from Section 7
	/					* Number of Index Entries and Index Array pointers from Section 7
	/			
	/			Outputs: * Nodal Displacement Array and data file
	/
	/-------------------------------------------------------------------------------------------------------------*/

	/*Get User to input amount of info output to be displayed on screen by the solver*/
	int pinfo = 1;	/*variable for solution info printing 0 -> 4*/

	/*Solve the Equation*/
	solve(mtemp, NumEntries, irn_in, jcn_in, A_in, NumRhs, rhs, pinfo);
	/*solve2(mtemp, NumEntries, irn_in, jcn_in, A_in, rhs);*/

	/*Free Index array memeory*/
	free(A_in);
	free(irn_in);
	free(jcn_in);

	/*Re-insert zero displacements into the displacement array*/
	count = 1;
	m=0;
	txtemp = FixDofs2[count]-1;
	for(n=0;n<MatrixOrder;n++)
	{
		if(n != txtemp)		/*If dof not fixed copy accross the value*/
		{
			for(o=0;o<NumRhs;o++)
			{
				temp = o * MatrixOrder; /*to get to correct place in rhs_in matrix*/
				temp2 = o * mtemp; /*to get to correct place in rhs matrix*/
				rhs_in[n+temp] = rhs[m+temp2];
			}
			m++;
		}
		
		else {				/*Else dof is fixed! i.e zero displacement -> then look for next fixed dof*/
			for(o=0;o<NumRhs;o++)
			{
				temp = o * MatrixOrder;
				rhs_in[n+temp] = 0.0;
			}
			txtemp = FixDofs2[++count]-1;
		}
	}

	/*Free old rhs array*/
	free(rhs);
	free(FixDofs2);

	/* Write Displacement information file*/ 
	int nnods = div(MatrixOrder,2).quot;
	/*printf("\nnnods = %i",nnods);*/
	/*sprintf(plotname,"Displacement%i.txt",itt); /*set name for current output file*
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Displacement writefile\n");
		}
	else{
	fprintf(outfile,"node\tx\ty\tz\trotx\troty\trotz\n"); /*column headings*
		for(n=0;n<(elemX+1)*(elemY+1);n++)
		{
			count = (6 * n);
			fprintf(outfile,"%i\t",(n));
			fprintf(outfile,"%f\t",rhs_in[count]*1000);
			fprintf(outfile,"%f\t",rhs_in[count+1]*1000);
			fprintf(outfile,"%f\t",rhs_in[count+2]*1000);
			fprintf(outfile,"%f\t",rhs_in[count+3]*1000);
			fprintf(outfile,"%f\t",rhs_in[count+4]*1000);
			fprintf(outfile,"%f\t",rhs_in[count+5]*1000);


			fprintf(outfile,"\n");
		}
	}

	fclose(outfile);*/

	/*printf("\n\nDisplacement File written");*/
	
}


