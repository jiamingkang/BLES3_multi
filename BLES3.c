/*
 *  BLES3.c
 *  Now capable of minimising 4 objectives:
 *		1) Compliance (strain energy)
 *		2) Expected Compliance
 *		3) Variance of compliance
 *		4) Weighted sum of 2) & 3)
 *
 *  Created by Peter Dunning on 06/02/2009.
 *	
 *	Version 3.0 @ 21/06/2010
 */
 
#include "ls_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "FixedGrid.h"
#include "EMatrix.h"
#include "Numbering.h"
#include "IFG.h"
#include "AFG.h"
#include "KMatrix.h"
#include "Solve.h"
#include "Strain.h"
#include "Levels.h"
#include "BCs.h"
#include "FibreSensitivity.h"

int main()
{

printf("\n\n------*---Start of BLES Version 3 Program---*------\n");
printf("Compiled without NL Velocity map & alpha/dist weighting & no smoothing at all!\nFullde-bug output!\nh_bar = 1h, holeLim = 0.02, bw=6h\n");
int i,j,k,p;	/*Incrementors*/
double ftemp;
int temp,num; /*Temperary variables*/
int count,begin; /*Other incrementors*/
FILE *outfile;	/*File varible for output files*/
char plotname[40];	/*variable to change names of plotting output files*/

/*------------------------------------------------------------------------------------------------------------
/
/		Section 1: Input Data from input file and keyboard
/
/			Inputs: * Name of Input Data File from keyboard user input
/					* Material Properties, Element Edge Length, number of X & Y elements
/						and inital hole data all from Input Data File
/					* Analysis method type, volume constraint % and lsf fix flag from keyboard
/			
/			Outputs: * Element Area and IN element stiffness matrix
/					 * Maximum Co-ordinates of struture domain conmforming to square shaped elements
/					 * Number of inital circular holes with center co-ordinate and radius data
/					 * Volume constraint in actual terms
/
/-------------------------------------------------------------------------------------------------------------*/

/*-----------Get user to enter name of input data file and open file if file exists-------------*/
/*char datafile[20] = "Femur";	/*Variable for input data file*/
char incheck[20] = "0";	/*Variable to check numerical input*/
char incheckcrit[20] = "1"; /* set objective to minimise: 1 = complience 2 = expected complience 3 = varience of complience 4 = wieghted combo of 2 & 3*/
char incheckmeth[20] = "2"; /*set method for anaylsis: 1 = BFG 2 = AFG*/
char incheckhole[20] = "0"; /*Hole Insertion? 1 = Yes 0 = No*/
char incheckvol[20] = "0.5"; /*set volume ratio between 0.1 and 0.9*/
char incheckbc[20] = "0"; /*Fix boundary nodes? 1 = yes 0 = No*/
char incheckiter[20]; /*number of iterations*/
char datafile[20], datafile2[20];
int Numlsf = 1;				/*Number of level set functions in the solver 1 or 2 depending on the model*/
double DeltaThetaMax = 0.5;		/*Maximum change in fibre angel per iteration*/

    int countt=0;
    double scale_factor=1.0;
    
FILE *infile, *infile2;		/*input file pointer*/
do
{
	printf("\nEnter Name of Data File: ");
	scanf("%s",datafile); /*get datafile name from the terminal*/
	/*infile = fopen("Femur", "r"); /*Try to open the file*/
    printf("%s", datafile);
	infile = fopen(datafile, "r"); /*Try to open the file*/
	
	/*If file does not exist then get user to re-enter*/
	if(infile == NULL){
		printf("\nCould not find data input file! - Start Again\n");
		}
	/*Otherwise the file was opened successfully*/
	else if (infile != NULL){
		printf("File opened successfully\n");
		}
}
while (infile == NULL);

int OjFlag; /*flag for objective to minimise (1 -> 4) */
do
{ 
	printf("\nEnter Objective to minimise\n1. Compliance\n2. Expected Compliance");
	printf("\n3. Variace of Compliance\n4. Weighted combo of 2. & 3.\n----->\t");
	/*scanf("%s",incheck); /*gets user input from the terminal*/
	/*check input was an integer*/
	for(i=0;incheckcrit[i] >= '0' && incheckcrit[i] <= '9'; i++);
	if(incheckcrit[i] != '\0' || i == 0)
	{
		OjFlag = 0;
	}
	else
	{
		sscanf(incheckcrit,"%i",&OjFlag); /*convert input to integer*/
		/*printf("\nMflag = %i",Mflag);*/
	}
	
	/*If Mflag out of range or invalid then get user to re-enter*/
	if((OjFlag < 1) || (OjFlag > 4)){
		OjFlag = 0;
		printf("\nInvalid type! - Start Again\n");
		}
}
while (OjFlag == 0);

int Mflag; /*flag for analysis method type*/
do
{
	printf("\nEnter Analysis Method (1=BFG, 2=AFG): ");
	/*scanf("%s",incheck); /*gets user input from the terminal*/
	/*check input was an integer*/
	for(i=0;incheckmeth[i] >= '0' && incheckmeth[i] <= '9'; i++);
	if(incheckmeth[i] != '\0' || i == 0)
	{
		Mflag = 0;
	}
	else
	{
		sscanf(incheckmeth,"%i",&Mflag); /*convert input to integer*/
		/*printf("\nMflag = %i",Mflag);*/
	}
	
	/*If Mflag out of range or invalid then get user to re-enter*/
	if((Mflag != 1) && (Mflag != 2)){
		Mflag = 0;
		printf("\nInvalid type! - Start Again\n");
		}
}
while (Mflag == 0);

int Hflag; /*flag for hole insertion*/
do
{
	printf("\nHole Insertion? (1=Yes, 0=No): ");
	/*scanf("%s",incheck); /*gets user input from the terminal*/
	/*check input was an integer*/
	for(i=0;incheckhole[i] >= '0' && incheckhole[i] <= '9'; i++);
	if(incheckhole[i] != '\0' || i == 0)
	{
		Hflag = -1;
	}
	else
	{
		sscanf(incheckhole,"%i",&Hflag); /*convert input to integer*/
		/*printf("\nMflag = %i",Mflag);*/
	}
	
	/*If Mflag out of range or invalid then get user to re-enter*/
	if((Hflag != 1) && (Hflag != 0)){
		Hflag = -1;
		printf("\nInvalid! - Start Again\n");
		}
}
while (Hflag == -1);

/*-----Read in Isotropic Material Properties from file and print to screen, e = Youngs Modulus, v = Poisons Ratio-----*/
double e1,v12,e2,v21, g12,g23,g13; 
fscanf(infile, "%lf %lf %lf %lf %lf %lf  %lf", &e1, &e2, &v12, &v21, &g12, &g23, &g13);
double e,v,e11,g33;
/*fscanf(infile, "%lf %lf",&e, &v);
printf("\nYoungs Modulus = %lf\nPoisons ratio = %lf",e,v);*/

/*Calculate material consts in [E] matrix - For Plain Stress*/
e = 1.0;
v =0.3;
e11 = e;
e11 /= 1.0 - (v * v);
/*v12 = e11 * v;*/
g33 = e11 * (1.0 - v);
g33 *= 0.5;

printf("\nPlane Stress constants:\nE11=%lf, V12=%lf, G33=%lf",e11,v12,g33);



/*Maximum values of Structures Domain & Element edge length (h): Note this program currently only cope with square shaped elements i.e. h=hx=hy*/
double maxX,maxY,h;
int elemX, elemY;

/*Read in Element Edge Length (h) and number of elements in X & Y*/
fscanf(infile,"%lf",&h);
fscanf(infile,"%i",&elemX);
fscanf(infile,"%i",&elemY);

/*Input Number of plys*/
int numply;
double plyt;
fscanf(infile,"%i",&numply);
fscanf(infile,"%lf",&plyt);

double *thetaply;
thetaply = malloc(numply*sizeof(double));
printf("\nplyt = %f", plyt);
printf("\nnumply = %i", numply);

for(i=0;i<numply;i++)
{
	fscanf(infile,"%lf",&thetaply[i]);
	printf("thetaply[%i] = %f\n",i, thetaply[i]);
}



/*Define the Tow width (for now it will be 1 so there is a garenteed one line per node) Later work out a system for tows that clear multiple elements*/
double tow = 1.0;

/*calculate maximum x and y co-ordinates*/
maxX = (double)elemX * h;
maxY = (double)elemY * h;
printf("\nNumber of Elements [%i, %i]",elemX,elemY);
printf("\nMax dimentions of domain are X=%lf, Y=%lf",maxX,maxY);

double AreaElem;
AreaElem = h * h; /*Area of elements is edge length squared*/
printf("\nElement Edge length = %lf \nElement Area = %lf\n",h,AreaElem);

/*set some parameters*/
const double hbar = 1.0 * h; /*artifical height for hole creation method (in terms of grid spacing h)*/
const double HoleLim = 0.02; /*upper limit on new hole creation volume*/
const double lBand = elemX*elemY*h; /*Narrow band width (in terms of grid spacing h)*/

/*------Calculate Number of nodes, dofs and elements in domain------*/
int NumNodes = (elemX + 1) * (elemY + 1);
int NumDof = 6 * NumNodes;
int NumElem = elemX * elemY;
double curveLenght;

printf("\nNumber: Nodes=%i, Elements=%i, Dof=%i\n",NumNodes,NumElem,NumDof);
/*--------Import Node Co-ordinates--------*/
Coord NodeCoord[NumNodes]; /*Initialize Node co-ordinate array structure*/
Coord GridLoc[NumNodes]; /*Initialize Node co-ordinate array structure*/
for(i=0;i<NumNodes;i++)
{
	fscanf(infile,"%lf",&GridLoc[i].x);
	fscanf(infile,"%lf",&GridLoc[i].y);
	GridLoc[i].z = 0.0;
	GridLoc[i].xz = 0.0;
}

for(i=0;i<NumNodes;i++)
{
	fscanf(infile,"%lf",&NodeCoord[i].x);
	fscanf(infile,"%lf",&NodeCoord[i].y);
	fscanf(infile,"%lf",&NodeCoord[i].z);
	fscanf(infile,"%lf",&NodeCoord[i].xz);
	curveLenght = NodeCoord[i].xz; 		/*Will end up as maximum*/	
}

Elem Number[elemX][elemY]; /*Initialise Numbering Array*/

for(j=0;j<elemY;j++)
{
	for(i=0;i<elemX;i++)
	{
		fscanf(infile,"%i",&Number[i][j].n);
		fscanf(infile,"%i",&Number[i][j].a);
		fscanf(infile,"%i",&Number[i][j].b);
		fscanf(infile,"%i",&Number[i][j].c);
		fscanf(infile,"%i",&Number[i][j].d);
	}
}

/*-----------Import the varible edge lenghts-------------*/
double hz[elemX];
double hxz[elemX];

printf("\n");

for(i=0;i<elemX;i++)
{
	fscanf(infile,"%lf",&hz[i]);
	printf("hz[%i] = %f\n", i, hz[i]);
}
printf("\n");

for(i=0;i<elemX;i++)
{
	fscanf(infile,"%lf",&hxz[i]);
	printf("hxz[%i] = %f\n", i, hxz[i]);
}

double maxXZ = curveLenght;
printf("\nCurveLenght = %f\n", curveLenght);

/*Get intial element area*/
double *alpha0;
alpha0 = malloc(NumElem*sizeof(double));

int q = 0;
for(j=0;j<elemY;j++)
{
	for(i=0;i<elemX;i++)
	{
		alpha0[q] = h*hxz[i];
		++q;
	}
}

/*Write Numbering info to a file*/
outfile = fopen("Numbering.txt", "w");
if(outfile == NULL){
	printf("\nFailed to open Numbering writefile\n");
	}
else{
fprintf(outfile,"n\ta\tb\tc\td\n"); /*column headings*/
	for(j=0;j<elemY;j++)
		{
		for(i=0;i<elemX;i++)
			{
			fprintf(outfile,"%i\t",Number[i][j].n);
			fprintf(outfile,"%i\t",Number[i][j].a);
			fprintf(outfile,"%i\t",Number[i][j].b);
			fprintf(outfile,"%i\t",Number[i][j].c);
			fprintf(outfile,"%i\t",Number[i][j].d);

			fprintf(outfile,"\n");
		}
	}
}

fclose(outfile);
printf("Numbering info file written");

/*--------Import the IN element stiffness matrix----------*/
double *KE; /*Initialise Element Stiffness Matrix*/
KE = malloc((576)*sizeof(double));

double *Trans;
Trans = malloc(elemX*9*sizeof(double));

for(i=0;i<(3*elemX);i++)
{
	for(j=0;j<3;j++)
	{
		temp = (3 * i) + j;
		fscanf(infile,"%lf",&Trans[temp]);
		/*printf("\nTrans[%i] = %f", temp, Trans[temp]);*/
	}
}

/*Intialise theta matrix (0 is horiziontal)*/
double *theta;
theta = malloc(elemX*elemY*numply*sizeof(double));

for(i=0;i<elemX; i++)
{
	for(j=0;j<elemY; j++)
	{
		num = Number[i][j].n-1;
		for(k=0;k<numply;k++)
		{
			theta[num*numply + k] = thetaply[k]; 
		}
	}
}
free(thetaply);

/*Print KE Matrix to file*
outfile = fopen("Element.txt", "w");
if(outfile == NULL){
	printf("\nFailed to open Element writefile\n");
	}
else{
	for(i=0;i<(24*elemX);i++)
	{
		for(j=0;j<24;j++)
		{
			temp = (24 * i) + j;
			fprintf(outfile,"%lf",KE[temp]);
			fprintf(outfile,"\t");
		}
		fprintf(outfile,"\n");
	}
}
fclose(outfile);
printf("\nElement stiffness matrix file written\n");*/

double vol_in; /*varibale for volume ratio input*/

do
{
	temp = 1;
	printf("\nEnter Volume Ratio (0.1 < Vr < 0.9): ");
	/*scanf("%s",incheck); /*gets user input from the terminal*/
	/*check input was an integer*/
	for(i=0;(incheckvol[i] >= '0' && incheckvol[i] <= '9') || incheckvol[i] == '.'; i++);
	if(incheckvol[i] != '\0' || i == 0)
	{
		vol_in = 0.0;
	}
	else
	{
		sscanf(incheckvol,"%lf",&vol_in); /*convert input to integer*/
		/*printf("\nMflag = %i",Mflag);*/
	}
	/*If volume ratio out of range then get user to re-enter*/
	if((vol_in < 0.1) || (vol_in > 0.9)){
		printf("\nvolume ratio out of range!! or invalid\n");
		temp = 0;
		}
}
while (temp == 0);

ftemp = elemY * h * curveLenght;	/*Total domain volume (area)*/
double VolRat = vol_in * ftemp;		/*set target volume for structure as the specified volume ratio * total domain volume*/
printf("Volume Constraint = %lf\n",VolRat);

/*ask user if lsf should be fixed where bc's are applied*/
int fixFlag = 0;
do
{
	temp = 1;
	printf("\nFix level set function at bc's? (1 = yes, 0 = no): ");
	/*scanf("%s",incheck); /*gets user input from the terminal*/
	/*check input was an integer*/
	for(i=0;incheckbc[i] >= '0' && incheckbc[i] <= '9'; i++);
	if(incheckbc[i] != '\0' || i == 0)
	{
		fixFlag = 2;
	}
	else
	{
		sscanf(incheckbc,"%i",&fixFlag); /*convert input to integer*/
	}

	/*If volume ratio out of range then get user to re-enter*/
	if((fixFlag != 0) && (fixFlag != 1)){
		printf("\ninput error!! try again...");
		temp = 0;
		}
}
while (temp == 0);


/*read in number of rectangular holes in domain & data*/
int NumRect;
fscanf(infile,"%i",&NumRect);
printf("\nNumber of inital rectangular holes = %i\n",NumRect);
Coord *Rect;
Rect = malloc(2 * NumRect * sizeof(Coord)); /*array to store rectangular hole data (max & min x&y coords)*/

for(i=0;i<NumRect;i++)
{
	j = 2 * i;
	fscanf(infile,"%lf",&Rect[j].x); /*min X*/
	fscanf(infile,"%lf",&Rect[j+1].x); /*max X*/
	fscanf(infile,"%lf",&Rect[j].y); /*min Y*/
	fscanf(infile,"%lf",&Rect[j+1].y); /*max Y*/

	printf("\nRect%i",i);
	printf("\n%lf",Rect[j].x); //min X
	printf("\n%lf",Rect[j+1].x); //max X
	printf("\n%lf",Rect[j].y); //min Y
	printf("\n%lf",Rect[j+1].y); //max Y
}
	
/*------------------------------------------------------------------------------------------------------------
/
/		Section 2: Number all Elements and Nodes & calculate Node Co-ordinates
/
/			Inputs: * Maximum Domain Co-ordinates, Element Edge Length from Section 1
/			
/			Outputs: * Number of elements in each direction (x, y)
/					 * Total number of Nodes, dofs and Elements in Domain
/					 * Numbering array for all Elements and Nodes in domain and data file
/					 * Second node numbering array in 2D matrix format
/					 * Co-ordinate array for all Nodes in domain and data file
/
/-------------------------------------------------------------------------------------------------------------*/


/* Write Node Co-ordinate information file */
sprintf(plotname,"%s_NCrd.txt",datafile); /*set name for current output file*/
outfile = fopen(plotname, "w");
if(outfile == NULL){
	printf("\nFailed to open Node_Coord writefile\n");
	}
else{
/*fprintf(outfile,"Node Num\tx\ty\n"); /*column headings*/
	for(i=0;i<NumNodes;i++)
	{
			fprintf(outfile,"%i\t",(i+1));
			fprintf(outfile,"%lf\t",NodeCoord[i].x);
			fprintf(outfile,"%lf\t",NodeCoord[i].y);
			fprintf(outfile,"%lf\t",NodeCoord[i].z);
			fprintf(outfile,"%lf\t",NodeCoord[i].xz);

		fprintf(outfile,"\n");
	}
}

fclose(outfile);
printf("Node Co-ordinate File written\n");

sprintf(plotname,"%sxz_NCrd.txt",datafile); /*set name for current output file*/
outfile = fopen(plotname, "w");
if(outfile == NULL){
	printf("\nFailed to open Node_Coord writefile\n");
	}
else{
/*fprintf(outfile,"Node Num\tx\ty\n"); /*column headings*/
	for(i=0;i<NumNodes;i++)
	{
			fprintf(outfile,"%i\t",(i+1));
			fprintf(outfile,"%lf\t",NodeCoord[i].xz);
			fprintf(outfile,"%lf\t",NodeCoord[i].y);

		fprintf(outfile,"\n");
	}
}

fclose(outfile);
printf("Node Co-ordinate File written\n");


int NodeX = elemX + 3; 
int NodeY = elemY + 3;	/* +1 as there is always 1 more node in a direction than elemenets, +2 for ghost nodes either side of the domain*/
int Nodes2[NodeX][NodeY]; /*Array to store ordered node numbers*/

/*Determine ordred 2D array for node numbers*/
NodeNums2(h, NodeX, NodeY, Nodes2, NumNodes, GridLoc);

/* Write Ordered Node Number information file 
outfile = fopen("Nodes2.txt", "w");
if(outfile == NULL){
	printf("\nFailed to open Nodes2 writefile\n");
	}
else{
	for(j=NodeY-1;j>=0;j--)
	{
		for(i=0;i<NodeX;i++)
		{
			fprintf(outfile,"%i\t",Nodes2[i][j]);
		}
		fprintf(outfile,"\n");
	}
}

fclose(outfile);
printf("Ordered Node Number File written");*/

/*------------------------------------------------------------------------------------------------------------
/
/		Section 3: Read in and apply Boundary conditions and Forces
/
/			Inputs: * Number of Zero Displacement Boundary conditions and co-ordinates from input file
/					* Number of Fixed Areas, max/min co-ordinates and fixed directions from input file
/					* Number of Applied Forces and co-ordinates from input file
/					* Number of Uniform forces and max/min co-ordinates from input file
/					* Total Number of Nodes and Node co-ordinates from Section 2
/			
/			Outputs: * Array storing the row and column numbers to be deleted from global stiffness matrix
/						due to Zero Displacement Boundary conditions
/					 * Array storing the force vector, i.e. right hand side (rhs) of equation
/					 * Set of nodes with fixed zero values for the signed distance function
/
/-------------------------------------------------------------------------------------------------------------*/



/*--------------------------------------First work out fixed dofs----------------------------------*/

int *fixed; /*Array to store nodes that have fixed zero lsf values thoughout the optimisation, if fixFlag == 1 */
int NumFix = 0; /*variable to count number of fixed lsf values, initally zero*/
if(fixFlag == 1)
{
	NumFix = 1; /*set to indicate that fixFlag =1*/
	fixed = calloc(NumNodes,sizeof(int)); /*Array to store nodes that have fixed zero lsf values thoughout the optimisation*/
}

Aux *Afix; /*pointer to store fixed area data (n[0],[1]=fixed direction (x,y) flag, x,y[0] = x,y min, x,y[1] = x,y max,*/
int LAfix = 0; /*Length of Afix array, initially zero*/
if( (Mflag == 1) && (fixFlag == 0) )
{
	LAfix = 50; /*Set maximum length of Afix array (if using BFG method) */
	Afix = malloc(LAfix * sizeof(Aux)); /*reserve memory for array*/
}

int TotFix; /*Total number of fixed dof*/
int *FixDofs; /*arrary to store fix dof numbers*/
FixDofs = malloc(NumDof * sizeof(int)); /*intitalise to largest possible size*/

/*function that finds all zero dispacement boundary conditions and assembles into one array */
FixDisp(h, NumNodes, GridLoc, &NumFix, fixed, &TotFix, FixDofs, &LAfix, Afix, infile);

/*Now re-allocate memory for fixed dof arrays*/
FixDofs = realloc(FixDofs, TotFix * sizeof(int)); /*fixed dof array*/

if(fixFlag == 1)
{
	fixed = realloc(fixed, NumFix * sizeof(int)); /*nodes with fixed lsf value array*/
}

else if(Mflag == 1)
{
	Afix = realloc(Afix, LAfix * sizeof(Aux)); /*fixed area data when using BFG method*/
}

/*--------------------------------Now work out required load vectors and variance matrix----------------------------*/

/*First read in loading data from the input file*/
int NumPforce;	/*number of loads applied as point forces*/
double *Pforce;	/*pointer for array that stores point force data:
					6 enties per force (x, y, mean magnitude, stand dev magnitude, mean direction, stand dev direction*/
int NumAforce;	/*number of distributed forces*/
double *Aforce;	/*pointer for array that stores point force data:
					8 enties per force (xmin, xmax, ymin, ymax, mean magnitude, stand dev magnitude, mean direction, stand dev direction*/
					
fscanf(infile, "%i", &NumPforce);	/*read in number of point loads*/
printf("\nNumber of point loads = %i",NumPforce);

if(NumPforce > 0)
{
	Pforce = malloc(NumPforce * 8 * sizeof(double));	/*allocate memory for point force data array*/
	
	for(i=0;i<NumPforce;i++)
	{
		j = i * 8; /*place indicator for Pforce array*/
		fscanf(infile, "%lf", &Pforce[j]);		/*x co-ordinate*/
		fscanf(infile, "%lf", &Pforce[j+1]);	/*y co-ordinate*/	
		fscanf(infile, "%lf", &Pforce[j+2]);	/*X magnitude of force*/
		fscanf(infile, "%lf", &Pforce[j+3]);	/*Y magnitude of force*/
		fscanf(infile, "%lf", &Pforce[j+4]);	/*Z magnitude of force*/
		fscanf(infile, "%lf", &Pforce[j+5]);	/*ROTX magnitude of force*/
		fscanf(infile, "%lf", &Pforce[j+6]);	/*ROTY magnitude of force*/
		fscanf(infile, "%lf", &Pforce[j+7]);	/*ROTZ magnitude of force*/
		/*printf("\n Pforce [%i] = (%f, %f, %f, %f, %f, %f)",j, Pforce[j], Pforce[j+1], Pforce[j+2], Pforce[j+3], Pforce[j+4], Pforce[j+5], Pforce[j+6], Pforce[j+7]);*/
	}
}

fscanf(infile, "%i", &NumAforce);	/*read in number of distributed loads*/
printf("\nNumber of distributed loads = %i",NumAforce);

if(NumAforce > 0)
{
	Aforce = malloc(NumAforce * 8 * sizeof(double));	/*allocate memory for point force data array*/
	
	for(i=0;i<NumAforce;i++)
	{
		j = i * 8; /*place indicator for Aforce array*/
		fscanf(infile, "%lf", &Aforce[j]);		/*min x co-ordinate*/
		fscanf(infile, "%lf", &Aforce[j+1]);	/*max x co-ordinate*/
		fscanf(infile, "%lf", &Aforce[j+2]);	/*min y co-ordinate*/
		fscanf(infile, "%lf", &Aforce[j+3]);	/*max y co-ordinate*/
		fscanf(infile, "%lf", &Aforce[j+4]);	/*magnitude of force*/
		fscanf(infile, "%lf", &Aforce[j+5]);	/*variance of magnitude (unused)*/
		fscanf(infile, "%lf", &Aforce[j+6]);	/*direction of load*/
		fscanf(infile, "%lf", &Aforce[j+7]);	/*variance of direction (unused)*/
	}
}
int Symetry;
fscanf(infile, "%i", &Symetry);

fclose(infile); /*Now that all the input values have been read in, close the input file*/

int NumRhs = 1; /*Number of load cases required to solve problem, initialise to one*/
int NumForces;	/*Number of dofs with applied forces, over all load cases*/
int *rhs_ind;	/*location in array pointer for mean loading conditions*/
double *rhs_val;	/*values in array pointer for mean loading conditions*/

if(OjFlag == 1) /*if objective is to minimise (deterministic) strain energy or involves variance*/
{
	/*first assign initial amount of memory*/
	rhs_ind = malloc(NumDof * sizeof(int));
	rhs_val = malloc(NumDof * sizeof(double));
	
	/*asseble mean loading array (compressed memory)*/
	LoadVector(h, NumNodes, GridLoc, &NumForces, rhs_ind, rhs_val, NumPforce, Pforce, NumAforce, Aforce);
	
	/*re-assign required memory*/
	rhs_ind = realloc(rhs_ind, NumForces * sizeof(int));
	rhs_val = realloc(rhs_val, NumForces * sizeof(double));
	
	for(i=0;i<NumForces;i++)
	{
		printf("\nrhs_ind =%i, rhs_val = %lf",rhs_ind[i],rhs_val[i]);
	}
}

/*int NumUn = 0;		/*Number of entries in all laod vectors*
int *unc_ind;		/*location in array pointer for uncertain loading conditions*
double *unc_val;	/*values in array pointer for uncertain loading conditions*
double *unc_wgt;	/*pointer array for load case weights*
unc_wgt = malloc(1000 * sizeof(double));
unc_wgt[0] = 1.0; /*weight for first load vector*/

/*if(OjFlag != 1) /*if objective involves uncertainty*
{
	/*first assign initial amount of memory*
	unc_ind = malloc(NumDof * sizeof(int));
	unc_val = malloc(NumDof * sizeof(double));
	
	/*asseble loading arrays (compressed memory)*
	UncLoads(h, NumNodes, NodeCoord, &NumUn, unc_ind, unc_val, &NumRhs, unc_wgt, NumPforce, Pforce, NumAforce, Aforce);
	
	/*re-assign required memory*
	unc_ind = realloc(unc_ind, NumUn * sizeof(int));
	unc_val = realloc(unc_val, NumUn * sizeof(double));
	
	for(i=0;i<NumUn;i++)
	{
		printf("\nunc_ind =%i, unc_val = %lf",unc_ind[i],unc_val[i]);
	}
	for(i=0;i<NumRhs;i++)
	{
		printf("\nunc_wgt =%lf",unc_wgt[i]);
	}
}*/

/*unc_wgt = realloc(unc_wgt, NumRhs * sizeof(double));

int NumVar = 0;	/*Number of nodes with co-variance matrix entries*
int *var_ind; /*pointer array to store node numbers effected by variance*
double *var_val; /*pointer array for co-variance matrix entries (3 per node) *

/*if( (OjFlag == 3) || (OjFlag == 4) ) /*if objective involves variance of strain energy*
{
	/*first assign initial amount of memory*
	var_ind = malloc(NumDof * sizeof(int));
	var_val = malloc(3 * NumDof * sizeof(double));
	
	/*compute "co-variance" matrix entries*
	VarMat(h, NumNodes, NodeCoord, &NumVar, var_ind, var_val, NumPforce, Pforce, NumAforce, Aforce);
	
	/*re-assign required memory*
	var_ind = realloc(var_ind, NumVar * sizeof(int));
	var_val = realloc(var_val, 3 * NumVar * sizeof(double));
	
	for(i=0;i<NumVar;i++)
	{
		j = 3 * i;
		printf("\nvar_ind =%i\nvar_val[1] = %lf, var_val[2] = %lf, var_val[3] = %lf",var_ind[i],var_val[j],var_val[j+1],var_val[j+2]);
	}
}*/

/*double wgt_combo = 1.0; /* weight for combined objective (no. 4), initially 1.0 for objective 3 *
				   
if(OjFlag == 4)
{
	ftemp = 0.0; /*set initial to zero*
	for(i=0;i<NumUn;i++) /*compute mean loading vector squared*
	{
		if(unc_ind[i] < NumDof)
		{
			ftemp += unc_val[i] * unc_val[i];
		}
	}
	
	wgt_combo = e / ftemp; /*weight = Youngs Modulus / f^2 *
	printf("\nNon-Dimensional Weight = %lf",wgt_combo);
	
	for(i=0;i<NumUn;i++) /*modify weights for expected compliance load cases*
	{
		unc_wgt[i] *= wgt_combo;
	}
	
	/*ask user to input weight for variance part of objective*
	do
	{
		temp = 1;
		printf("\nEnter weight for variance part of objective\n(wight for expected part = 1.0) : ");
		scanf("%s",incheck); /*gets user input from the terminal*
		/*check input was an integer*
		for(i=0;(incheck[i] >= '0' && incheck[i] <= '9') || incheck[i] == '.'; i++);
		if(incheck[i] != '\0' || i == 0)
		{
			ftemp = -1.0;
		}
		else
		{
			sscanf(incheck,"%lf",&ftemp); /*convert input to integer*
		}
		/*If weight negative then get user to re-enter*
		if(ftemp < 0.0) {
			printf("\nWeight is invalid or negative!\n");
			temp = 0;
		}
	}
	while (temp == 0);
	
	/*ftemp *= wgt_combo * wgt_combo; /*square the weight for the variance part of the objective*
	wgt_combo *= ftemp;
}*/
	
printf("\nNumber of Fixed nodes = %i",NumFix);
/*printf("\nLength of unm arrays = %i",NumUn);*/
   // for (i=0; i<NumFix; i++) {
    //    printf("fixed %d, %d", i, fixed[i] );
    //}
    
    
printf("\nNumber of RHS' = %i",NumRhs);

/*------------------------------------------------------------------------------------------------------------
/
/		Section 4: Calculate Inital Values of the signed distance function & active set of nodes
/
/			Inputs: * Maximum and minimum co-ordinates from Section 1
/					* Inital circular hole data from Section 1
/					* Node Co-ords and numbers from Section 2
/					* Set of nodes with fixed zero values for the signed distance function from Section 3
/			
/			Outputs: * Inital values for the signed distance function
/					 * Active nodes and near edge (mine) nodes within the inital narrow band
/-------------------------------------------------------------------------------------------------------------*/
double *lsf, *lsf2, *lsf2b;
lsf = malloc(NumNodes * sizeof(double)); /*Array to store all nodal values of the signed distance level set function*/
lsf2 = malloc(Numlsf*NumNodes * sizeof(double)); /*Array to store all nodal values of the signed distance level set function*/
lsf2b = malloc(NumNodes * sizeof(double)); /*Array to store all nodal values of the second signed distance level set function*/

int *active;
active = calloc(NumNodes,sizeof(int)); /*Array to track the active node set, entries = node nums*/
int *activeHoles;
activeHoles = calloc(NumNodes,sizeof(int)); /*Array to track the active node set, entries = node nums*/
int *mine;
mine = calloc(NumNodes,sizeof(int)); /*Array to track the near edge node set for narrow band method*/

/*---------First define signed distance fuction for domain without holes---------*/
double minX,minY,dist,leftX,rightX;
int dtemp;

/*For all nodes find the closest edge of the outer domain*/
for(i=0;i<NumNodes;i++)
{
	/*Is node closer to the right or left edge*/
	leftX = 0;
	rightX = 0;
	dtemp = NodeCoord[i].xz;
	for(j=0;j<dtemp;j++)
	{
		leftX += hxz[j];
	}
	for(j=dtemp;j<elemX;j++)
	{
		rightX += hxz[j];
	}

	minX = ( leftX < rightX ) ? leftX : rightX;

	/*Is node closer to the top or bottom edge*/
	ftemp = maxY - NodeCoord[i].y;
	minY = ( (NodeCoord[i].y - ftemp) < 0.000001 )  ? NodeCoord[i].y : ftemp;
	
	/*Signed distance function is then the lower of minX & minY*/
	lsf[i] = (minX < minY) ? minX : minY;
	lsf2b[i] = elemX*elemY; /*Intialise to a large number!*/
}

/*call function to update lsf if a rectangular hole exists*/
if(NumRect > 0)
{
	RectHole(NumNodes, NodeCoord, NumRect, Rect, lsf, elemX, elemY, hxz, maxXZ);
}

free(Rect); /*Inital rectangular hole data no longer required*/


/*---------Read in a lsf1 input for the second levelset-----------------*/
double *lsf_old; 
lsf_old = malloc(NumNodes * sizeof(double));

for(p=0;p<Numlsf;p++)
{
	sprintf(datafile2,"AAlsfin%i.txt",p+1);
	infile2 = fopen(datafile2, "r");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		fscanf(infile,"%lf",&lsf2b[Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fscanf(infile,"%lf",&lsf2b[Number[i][0].b-1]);
		}
		for(j=0;j<elemY;j++)
		{
			fscanf(infile,"%lf",&lsf2b[Number[0][j].d-1]);

			for(i=0;i<elemX;i++)
			{
				fscanf(infile,"%lf",&lsf2b[Number[i][j].c-1]);
			}
		}
	}
	fclose(infile2);

	ReIntlsf2(NodeX, NodeY, Nodes2, NumNodes, lsf2b, h, -1, lsf_old, fixed, NumFix,  elemX, hxz);
	ReIntlsf2(NodeX, NodeY, Nodes2, NumNodes, lsf2b, h, 1, lsf_old, fixed, NumFix, elemX, hxz);
	ArrayCopy(NumNodes, lsf_old, lsf2b);

	for(j=0;j<NumNodes;j++)
	{
		k = j + (NumNodes*p);
		//printf("\nk = %i, j = %i\n", k, j);
        
        printf("\nk = %i, lsf = %f\n", k, lsf2b[j]);
		lsf2[k] = lsf2b[j];
    
	}
    double DL;
    //DL=-0.0000001;
    //lsf2[362]=lsf2[362]+DL;  // 82, 133,224, 362, 400, 450
    //lsf2[123]=lsf2[123]-DL;

}
free(lsf_old);


	/*for(i=0;i<NumNodes;i++)
	{
		lsf2[i] += 1;
	}*/

/*---------Read in a lsf2 input-----------------*

	sprintf(datafile2,"AAlsfin.txt");
	infile2 = fopen(datafile2, "r");
	printf("\nB0\nB0");
	if(infile2 == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n\n");
		}
	else
	{
		fscanf(infile,"%lf",&lsf2[Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fscanf(infile,"%lf",&lsf2[Number[i][0].b-1]);
		}
		for(j=0;j<elemY;j++)
		{
			fscanf(infile,"%lf",&lsf2[Number[0][j].d-1]);
			for(i=0;i<elemX;i++)
			{
				fscanf(infile,"%lf",&lsf2[Number[i][j].c-1]);
			}
		}
	}
	fclose(infile2);
/*---------Read in a lsf2 input-----------------*
for(i=0;i<NumNodes;i++)
{
	lsf2[i] += 1;
}/**/

	//ThetaUpdate3(theta, elemX, elemY, Number, lsf2, h, hxz);          /* caculate the theta value*/
	
    ThetaUpdate4(theta, elemX, elemY, Number, lsf2, lsf2b, Numlsf, h, hxz);
    
    //ThetaUpdate5(theta, elemX, elemY, Number, lsf2, lsf2b, Numlsf, h, hxz);
	/*sprintf(plotname,"lsfIntial.txt"); /*set name for current topology output file*	
	/* Write Node signed distance value information file post re-initialisation*
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		for(i=0;i<NumNodes;i++){ fprintf(outfile,"lsf2[%i] = %0.3f;\n", i, lsf2[i]); }
	}
	fclose(outfile);*/

for(p=0;p<Numlsf;p++)
{
	sprintf(plotname,"PlyPlot0lsf%i.txt",p+1); /*set name for current topology output file*/	
	/* Write Node signed distance value information file post re-initialisation*/
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		k = (NumNodes*p);
		/*Bottom row of Nodes first*/
		fprintf(outfile,"%lf\t",lsf2[k+Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fprintf(outfile,"%lf\t",lsf2[k+Number[i][0].b-1]);
		}
		fprintf(outfile,"\n");
		/*Then the rest of the Nodes*/
		for(j=0;j<elemY;j++)
		{
			fprintf(outfile,"%lf\t",lsf2[k+Number[0][j].d-1]);
			for(i=0;i<elemX;i++)
			{
				fprintf(outfile,"%lf\t",lsf2[k+Number[i][j].c-1]);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);
}

	sprintf(plotname,"PlyPlot0com.txt"); /*set name for current topology output file*/	
	/* Write Node signed distance value information file post re-initialisation*/
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		/*Bottom row of Nodes first*/
		ftemp = 10000;
		for(p=0;p<Numlsf;p++)
		{
			k = (NumNodes*p);
			ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[0][0].a-1])) ? lsf2[k+Number[0][0].a-1]:ftemp;
		}
		fprintf(outfile,"%lf\t",ftemp);
		for(i=0;i<elemX;i++)
		{
			ftemp = 10000;
			for(p=0;p<Numlsf;p++)
			{
				k = (NumNodes*p);
				ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[i][0].b-1])) ? lsf2[k+Number[i][0].b-1]:ftemp;
			}
			fprintf(outfile,"%lf\t",ftemp);
		}
		fprintf(outfile,"\n");
		/*Then the rest of the Nodes*/
		for(j=0;j<elemY;j++)
		{
			ftemp = 10000;
			for(p=0;p<Numlsf;p++)
			{
				k = (NumNodes*p);
				ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[0][j].d-1])) ? lsf2[k+Number[0][j].d-1]:ftemp;
			}
			fprintf(outfile,"%lf\t",ftemp);
			for(i=0;i<elemX;i++)
			{
				ftemp = 10000;
				for(p=0;p<Numlsf;p++)
				{
					k = (NumNodes*p);
					ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[i][j].c-1])) ? lsf2[k+Number[i][j].c-1]:ftemp;
				}
				fprintf(outfile,"%lf\t",ftemp);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);



	double Cx, Cy;

	sprintf(plotname,"ThetaPlot0.txt"); 

	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
	}
	else{

		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{

				num = Number[i][j].n-1;
				Cx = NodeCoord[Number[i][j].a-1].xz + (hxz[i]/2);
				Cy = NodeCoord[Number[i][j].a-1].y + (h/2);
				fprintf(outfile,"%lf\t %lf\n",((Cx+0.5)-0.4*cos(theta[num])),((Cy+0.5)-0.4*sin(theta[num])));
				fprintf(outfile,"%lf\t %lf\n\n",((Cx+0.5)+0.4*cos(theta[num])),((Cy+0.5)+0.4*sin(theta[num])));
			}

		}
	}
	fclose(outfile);

/*---Intialise a fixed value of theta for tests---*/

	/*for(i=0;i<(elemX*elemY);i++)
	{
		theta[i] = 0.0;
	}*/

	/*for(j=0;j<elemY;j++)
	{
		for(i=0;i<elemX;i++)
		{
			num = Number[i][j].n-1;
			if(j<10){theta[num] = -0.7853981634;}
			else {theta[num] = 0.7854981634;}
		}
	}*/
/*---------------------------------------------------*/

sprintf(plotname,"Theta0.txt"); /*set name for current topology output file*/	
/* Write Theta value information file*/
outfile = fopen(plotname, "w");
if(outfile == NULL){
	printf("\nFailed to open Re-Int Implicit writefile\n");
}
else{
	/*fprintf(outfile,"Node Num\tx\ty\tphi\n"); /*column headings*/
	for(j=0;j<elemY;j++)
	{
		for(i=0;i<elemX;i++)
		{
			num = Number[i][j].n-1;
			fprintf(outfile,"%lf\t",(theta[num]*(180/3.141592654)));
		}
		fprintf(outfile,"\n");
	}
}
fclose(outfile);
/*Calculate the intial fiber continuity score*/


/* Write initial signed distance value information file *
outfile = fopen("PlotShp0.txt", "w");
if(outfile == NULL){
	printf("\nFailed to open Implicit writefile\n");
	}
else{
	/*fprintf(outfile,"Node Num\tx\ty\tphi\n"); /*column headings*
	for(i=0;i<NumNodes;i++)
	{
		fprintf(outfile,"%lf\n",lsf[i]);
	}
}

fclose(outfile);
printf("\nInital Signed Distance Info File written\n");*/

/*-------Apply the narrow band method to reduce the inital active set of nodes---------*/
int NumAct; /*Variable to count the number of active nodes*/
int NumMine; /*variable to count number of mines*/
int NumActHoles;

/*Calculate the bandwidth*/
NarBand(NumNodes, h, lBand, &NumAct, &NumMine, lsf, active, mine, NumFix, fixed, activeHoles, &NumActHoles);

/*re-size narrow band arrays*/
active = realloc(active,NumAct * sizeof(int));
mine = realloc(mine,NumMine * sizeof(int));
activeHoles = realloc(activeHoles,NumActHoles * sizeof(int));
	
/*------------------------------------------------------------------------------------------------------------
/
/		Section 5: Optimise strcuture boundary by the level set method
/
/			Inputs: * All information input / calculated so far (except inital hole data)
/			
/			Outputs: * Array for maximum, minimum boundary strain energy & lambda for each iteration
/					 * Array for maximum boundary velocity for each iteration
/					 * Array for matrix order (that is passed to MA57 solver) for each iteration
/					 * Array for total volume (area) of structure for each iteration
/					 * Array for total strain energy for each iteration
/					 * Array to store number of different element types for each iteration [IN,Pent,Quad,Tri,OUT]
/					 * Output file of the signed distance function for each iteration
/-------------------------------------------------------------------------------------------------------------*/

int itt = 0; /*Initalise number of itterations to 0*/
int numitt; /*variable to control number of iteration to run*/
int cont; /*varibale to manualy control iteration steps*/
double xCrd,yCrd; /*variables to store node x & y co-ordinates*/
double ElemU[8]; /*array for element displacement matrix*/
int ReCount = 0; /*initialize reinitialization count*/
double SensMaxOld=0;
double grad_old[NumNodes];    // this is for the conjugate gradient method.
double vnorm_old[NumNodes];
double vnorm_temp[NumNodes];
    

do
{
	temp = 1;
	printf("\n\nEnter Inital maximum number of iteration to perform (1 -> %i): ",10*NumNodes);
	scanf("%s",incheckiter); /*gets user input from the terminal*/
	/*check input was an integer*/
	for(i=0;incheckiter[i] >= '0' && incheckiter[i] <= '9'; i++);
	if(incheckiter[i] != '\0' || i == 0)
	{
		numitt = -1;
	}
	else
	{
		sscanf(incheckiter,"%i",&numitt); /*convert input to integer*/
	}

	/*If volume ratio out of range then get user to re-enter*/
	if((numitt < 1) || (numitt > 10*NumNodes)){
		printf("\nNumber out of range!!");
		temp = 0;
		}
}
while (temp == 0);

double Vmax; /*variable for maximum absolute boundary velocity*/
double lambda, lamT; /*variables for lambda values for boundary movement (lambda) and hole creation (lamT)*/
/*double gamma = 0.0; /*hard coded for now! - Convergence criteria: % of average boundary strain energy for Vnorm to deviate*/
/*double *lsfLast;
lsfLast = malloc(NumNodes * sizeof(double)); /*array to store previous lsf array (for converge test)*/

double *lsfHole;
if(Hflag == 1)
{
	lsfHole = malloc(NumNodes * sizeof(double));
	for(i=0;i<NumNodes;i++)
	{
		/*Set intial lsfHole values depending on sign of lsf*/
		lsfHole[i] = (lsf[i] > -0.000001) ? hbar : -hbar;
	}
}
/*variables for dynamincally changing the time step (if enabled)*/
double dt; /*variable for psuedo time step*/
double alt,del1,del2;	/*variables to dynamically change the time step co-efficient*/
alt = 0.1; /*initalise time step modifier*/

double *Grad; /*pointer for lsf gradient info (using upwind scheme)*/

/*Variables and arrays to store mesh data*/
short *ElemStat;	/*Element status array pointer*/
short *ElemStatb;	/*Element status array pointer*/
short *NodeStat;	/*Node status array pointer*/
short *NodeStatb;
Aux *auxNode;		/*Array to store auxillary node numbers and co-ords*/
int Ntot;	/*pointer & variable to count total nodes, inclusing auxillary ones*/
int Itotal,NIOtotal,Ototal; /*variables to count number of each element type*/
int tCount,qCount,pCount; /*Variables to count the number of the different types of B (NIO) elements*/
double *alpha; /*array to store element areas (or area ratios) */

/*Arrays to store node and element data related to the optimisation*/
double * theta_old;	/*Previous iteration's theta value, used to calculate Energy Factor*/
double *rhs_det;	/*deterministic rhs, then dispalcements*/
double *rhs_old;	/*previous iteration's dispalcements, used to calculate Energy Factor*/
double *rhs_unc;	/*uncertain rhs', then displacements*/
double *rhs_var;	/*variance rhs, then adjoint dispalcements*/
double *rhs_pnt;	/*pointer to reference different rhs' */
double *rhs_pnt2;	/* another pointer to reference different rhs' */
double *Nenergy; /*pointer for node (+ aux node) strain energy array*/
double *EF;		  /*The energy factor, caluclated from the 3 Nenergy's*/
double *Eenergy; /*pointer for array that stores element strain energy values*/
double *Vnorm;	/*pointer for node (+ aux node) normal velocity array*/
NElemDiff *DTL, *DTLb;	/*Pointer for the differential of fibre anlge with respect to the level set function at eac node*/
double *Ustrain, *Vstrain, *GXYstrain, *GXZstrain, *GYZstrain, *Ustress, *Vstress, *GXYstress, *GXZstress, *GYZstress; /*Stress and Strain varribles*/
double *ElemThetaStn, *ElemThetaStnB;	/*Element sensitivty values*/

/*Variables to calculate the lagrange multiplier such that the volume constraint is enforced*/
Aux *Lbound; /*array to store boundary length and strain energy associated with each boundary node, 
				x=bound length, y=strain enegy*/
double delVol; /*variable to calculate difference between current and target volume (area)*/

/*stuff for MA57 solver*/
int *irn, *jcn; /*pointers to Row and Column no. index arrays*/
double *A; /*pointer to Global Stiffness matrix entry index array*/
short *redun; /*pointer for array to indicate redundant dofs*/
int MatrixOrder;
int NumEntries; /*Variable for Number of index entries */

/*additional variables to fix auxillary nodes if they are in a fixed area*/
int TotFix2; /*new total fixed dofs*/
int *FixDofs2; /*new array storing dofs to be fixed*/
int IsCount, atemp,nNod;

/*arrays for pentagonal element data for the IFG method*/
/*Shape *Pshp;
Coord *Pcrds;*/

/*Set fit to vertex limit to 10% element edge length for IFG method, or -ve indicator for AFG method*/
double ftv_lim = (Mflag == 1) ? 0.1 : -0.01;

/*Arrays to store key data for each iteration*/
Astn *EmaxA;	
EmaxA = malloc(sizeof(Astn));	/*Array for maximum(x), minimum(y) boundary strain energy & lambda(e)*/
double *VolA;
VolA = malloc(sizeof(double));	/*Array for total volume (area) of structure for each iteration*/
double *CompA;
CompA = malloc(sizeof(double));	/*Array for total compilance (strain energy) for each iteration*/
double *SensA, *SensB, *FCS;
SensB = malloc(sizeof(double));	/*Array for total compilance (strain energy) for each iteration*/
SensA = malloc(sizeof(double));	/*Array for total senstitvity (strain energy differentiated wrt theta) for each iteration*/
FCS = malloc(sizeof(double)); /*fiber Continuity score stroage array*/

/*Set update calculation varribles*/
double *DeltaTheta, *Deltalsf2, *Deltalsf2b;

    
do {
	printf("\nStart iteration %i ...",itt);
    //if(itt>150){DeltaThetaMax = 17;}
    //if(itt>918){DeltaThetaMax = 5;}
	/*------------------------------------------------------------------------------------------------------------
	/
	/		Section 5.1: Analyse Structure, calculate boundary strain energy, velocities and lambda
	/
	/			Inputs: * Material properties and element edge length & maximum domain co-ordinates from Section 1
	/					* Analysis method type from Section 1
	/					* Element area & IN element stiffness matrix from Section 1
	/					* Number of: Elements, Nodes and dof from Section 2
	/					* Node co-ordinates & Mesh numbering data from Section 2
	/					* Zero Displacement Boundary conditions and force array from Section 3
	/					* Current (or inital) signed distance function from either Section 4 or 5.3/4
	/			
	/			Outputs: * Node status & Number of Auxillary nodes and co-ordinates
	/					 * Number of each element type [IN,Pent,Quad,Tri,OUT]
	/					 * Nodal displacements
	/					 * Element areas (or area ratios) and total structure area
	/					 * Pentagonal co-ordinates and shape functions (for BFG method)
	/-------------------------------------------------------------------------------------------------------------*/

	Ntot = NumNodes; /*Initalize the auxillary node number count to the current number of nodes*/
	ElemStat = calloc(Numlsf*NumElem,sizeof(short));	/*Element status array defined using dynamic memeory allocation*/
	ElemStatb = calloc(NumElem,sizeof(short));
	NodeStat = calloc(NumNodes*Numlsf,sizeof(short));	/*Node status array defined using dynamic memeory allocation*/
    
    NodeStatb = calloc(NumNodes,sizeof(short));
	/*Initialise length of auxNode array to store 2 auxillary nodes per element*/
	auxNode = calloc((2 * NumNodes),sizeof(Aux));
	
	/*---------------------Calculate node & element status, also add auxillary nodes----------------------*/
	printf("\nNtot = %i",Ntot);
    
    Aux AuxNodes[1000];  // assume there are 1000 auxlinary nodes
    Aux AuxNodesb[1000];
    
    
    
    int *NumAux=malloc( sizeof(int));
    int *NumBound=malloc( sizeof(int));
    
    int Ntot_new=NumNodes*Numlsf;
    int tempcounter=0;
    int NumofAux[Numlsf];
    //AuxNodes = malloc (NumNodes*sizeof(Aux));
    
    for (i=0; i<1000; i++) {
        AuxNodes[i].n=0;
        AuxNodes[i].x=0.0;
        AuxNodes[i].xz=0.0;
        AuxNodes[i].y=0.0;
        AuxNodes[i].z=0.0;
        
        AuxNodesb[i].n=0;
        AuxNodesb[i].x=0.0;
        AuxNodesb[i].xz=0.0;
        AuxNodesb[i].y=0.0;
        AuxNodesb[i].z=0.0;
    }

    
for(p=0;p<Numlsf;p++)
{
    // THis part is useless for my code
    Bseg *bptr=malloc(Numlsf*NumElem * sizeof(Bseg));
    int *na_conn_ind=malloc(Numlsf*(NumNodes+1) * sizeof(int));
    int *na_conn= malloc(Numlsf*4*NumNodes * sizeof(int));
    
	for(i=0;i<NumNodes;i++)
	{
		k = p*NumNodes+i;
		lsf2b[i] = lsf2[k];
        //printf("/n %d:%f,",i, lsf2b[i]); // so far , lsf value is correct
	}
    
    // come here
    
	IFG_Status(h,elemX,elemY,Number,&Ntot,NumElem,NodeX,NodeY,Nodes2,NumNodes,NodeCoord,lsf2b,NodeStatb,ElemStatb,auxNode,ftv_lim,itt,hxz,hz, AuxNodesb, bptr, na_conn_ind, na_conn, NumAux, NumBound);
    
     //Ntot_new=Ntot_new+*NumAux;   //Get the total number of nodes
    
    for(i=0; i<*NumAux; i++)
    {
  // printf("||%d:%f,%f,%d",i, AuxNodes[i].x,AuxNodes[i].y,AuxNodes[i].n);
    }
    for(i=0; i<*NumBound; i++)
    {
 //   printf("||%d: %d,%d, %d",i, bptr[i].n1, bptr[i].n2, bptr[i].e);
    }
    for(i=0; i<na_conn_ind[NumNodes]; i++)
    {
        //printf("%d:%d,",i, na_conn[i]);
    }
    
 //   printf("This is na_conn_ind %d", na_conn_ind[NumNodes]);
	for(i=0;i<NumElem;i++)
	{
		k = p*NumElem+i;
		ElemStat[k] = ElemStatb[i];
        
	}
    for(i=0;i<NumNodes;i++)
    {
        k = p*NumNodes+i;
        NodeStat[k] = NodeStatb[i];
        
    }
    
    for(i=0;i<*NumAux;i++)
    {
        k = tempcounter+i;
        AuxNodes[k].n=AuxNodesb[i].n;
        AuxNodes[k].x=AuxNodesb[i].x;
        AuxNodes[k].xz=AuxNodesb[i].xz;
        AuxNodes[k].y=AuxNodesb[i].y;
        AuxNodes[k].z=AuxNodesb[i].z;
    }
    NumofAux[p]=*NumAux;
    printf("NumofAux is %d,", NumofAux[p]);
    tempcounter+=*NumAux;
    
    for(i=0; i<*NumAux; i++)
    {
        printf("||%d:%f,%f,%d",i, AuxNodes[i].x,AuxNodes[i].y,AuxNodes[i].n);
    }
   
    
}
     Ntot_new=Ntot_new+tempcounter;

	printf("\nTotal nodes (inc. auxillary) = %i",Ntot_new);
	
	/*Algorithm to determine number of each element type*/
	elemType(NumElem, ElemStat, &Itotal, &Ototal, &NIOtotal, &tCount, &qCount, &pCount);
	
	/*sprintf(plotname,"NodeStat%i.txt",itt); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Node Strain Energy writefile\n");
	}
	else{

				fprintf(outfile,"%i\t",NodeStat[Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fprintf(outfile,"%i\t",NodeStat[Number[i][0].b-1]);
		}
		fprintf(outfile,"\n");
		for(j=0;j<elemY;j++)
		{
			fprintf(outfile,"%i\t",NodeStat[Number[0][j].d-1]);
			for(i=0;i<elemX;i++)
			{
				fprintf(outfile,"%i\t",NodeStat[Number[i][j].c-1]);
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);


	sprintf(plotname,"ElemStat%i.txt",itt);
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Element Status writefile\n");
		}
	else{

		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;
				fprintf(outfile,"%i\t",ElemStat[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	
	fclose(outfile);

	/*printf("\nElement Status File written");*/
		
	/*re-size auxillary node array to its actual size*/
	printf("\nNIOtotal = %i", NIOtotal);
	auxNode = realloc(auxNode, (2 * NIOtotal * sizeof(Aux)));
	
	/*work out if any auxilary nodes are within a fixed area - only needed for IFG*
	if((LAfix != 0) && (Mflag == 1))
	{		
		TotFix2 = TotFix-1; /*initalize to original value*
		temp = 2 * (Ntot - NumNodes); /*dofs for auxillary nodes*
		temp += TotFix; /*work out possible maximum nuber of fixed dof*
		FixDofs2 = malloc(temp * sizeof(int)); /*initally set to maximum possible*
		for(i=0;i<TotFix2;i++) /*copy the original fixed dofs array*
		{
			FixDofs2[i] = FixDofs[i];
		}
		/*work out if any auxillary nodes are with a fixed area and update data*
		AuxBound2( (LAfix / 2), Afix, NumNodes, Ntot, auxNode, NIOtotal, &TotFix2, FixDofs2);

		FixDofs2 = realloc(FixDofs2, TotFix2 * sizeof(int)); /*resize array*
		
		if(TotFix2 > TotFix)
		{
			printf("\n%i additional dofs fixed from auxillary nodes, now",TotFix2 - TotFix);
			/*for(i=0;i<(TotFix2 - 1);i++) /*copy the original fixed dofs array
			{
				printf("\nFixDofs2[%i] = %i",i,FixDofs2[i]);
			}*
		}
	}
	
	else otherwise total fixed dof unchanged*
	{*/
		TotFix2 = TotFix;	/*copy number of fixed dofs*/
		FixDofs2 = FixDofs; /*copy address across*/
	/*}
	
	/*for (i=0;i<TotFix2;i++) {
		printf("FixDofs2[%i]=%i",i,FixDofs2[i]);
	}*/
	
	/*---------------------Calcualte global stiffness matrix index arrays for the MA57 solver----------------------*/
	
	/*Calculate the maximum possible matrix order*/
	if(Mflag == 1)
	{
		MatrixOrder = 2 * Ntot; /*Order of matrix is initially the total dofs (including auxillary nodes)*/
	}
	else if(Mflag == 2)
	{
		MatrixOrder = NumDof; /*matrix order for AFG*/
	}
	
	/*add up number of OUT nodes and redefine the matrix order*/
	temp = 0;
	for(i=0;i<NumNodes;i++)
	{
		temp += (NodeStat[i] == 0) ? 1 : 0;
	}
	
	printf("\nMatrixOrder = %i",MatrixOrder);

	/*Assemble required rhs arrays*/
	if(OjFlag == 1) /*if Objective is deterministic compliance*/
	{
		/*Need to assemble rhs array from the compressed storage*/
		rhs_det = calloc((MatrixOrder), sizeof(double)); /*reallocate memory to new required length*/
		for(i=0;i<NumForces;i++)
		{
			rhs_det[rhs_ind[i]] = rhs_val[i];
			/*printf("rhs_det[%i] = %f = %f", rhs_ind[i], rhs_val[i], rhs_det[rhs_ind[i]]);*/
		}
	}
	/*printf("\nstop");
	printf("\nstop")*/;

	/*else /*if objective involves expected compliance or need the load cases for variance*
	{
		/*Need to assemble rhs array from the compressed storage*
		rhs_unc = calloc((MatrixOrder * NumRhs), sizeof(double)); /*reallocate memory to new required length*

		if(Mflag == 2)
		{
			for(i=0;i<NumUn;i++)
			{
				rhs_unc[unc_ind[i]] = unc_val[i];
			}
		}
		
		else if(Mflag == 1)
		{
			temp = MatrixOrder - NumDof; /*need to shift indexing because of auxillary nodes*
			/*printf("\ntemp = %i",temp);*
			for(i=0;i<NumUn;i++)
			{
				j = unc_ind[i] / NumDof; /*works out the rhs number*
				/*printf("\nj=%i, rhs entry %i",j,unm_ind[i] + (j * temp));*
				rhs_unc[unc_ind[i] + (j * temp)] = unc_val[i];
			}
		}
	}*/
	
	alpha = calloc(NumElem, sizeof(double)); /*Array to store element areas*/
	/*caculate the total number of entries for the MA57 solver index arrays*/
	/*if(Mflag == 1)
	{
		/*---Calculate Index Arrays for Multifrontal Solver {Also calculates stiffness matrices for Isoparametric formulation}---*
		NumEntries = (Itotal + qCount) * 36;
		NumEntries += tCount * 21;
		NumEntries += pCount * 55;
	}	
	else*/ 
	if(Mflag == 2)
	{
		/*for AFG - all non-OUT elements have an 24x24 symetric matrix*/
		NumEntries = elemX*elemY*576; /*ALL Elements are now Always in !*/  /*(Itotal + NIOtotal) * 576;*/
	}
	
	printf("\nInitial Number of Entries = %i",NumEntries);

	/*Initialize Index Arrays using dynamic memory allocation*/
	irn = calloc(NumEntries,sizeof(int));
	jcn = calloc(NumEntries,sizeof(int));
	A = malloc(NumEntries * sizeof(double));
	
	/*Initalise Pentagonal element data arrays, if some exist*/
	/*if((Mflag == 1) && (pCount != 0))
	{
		Pshp = malloc(pCount * 5 * sizeof(Shape)); /*Array to store pentagonal shape function co-efficients (20 per element)*
		Pcrds = malloc(pCount * 6 * sizeof(Coord)); /*Array to store pentagonal co-ordinates, 6 per element including center*
	}*/
	
	/*Calcualte global stiffness matrix index arrays for the MA57 solver*/
	redun = calloc(MatrixOrder, sizeof(short)); /*memory alloc for redundant dof indicator array*/
	/*if(Mflag == 1)
	{
		/*Calculate the global stiffness matrix in index format, ready for the MA57 solver - using IFG*
		IFG_Matrix(e11, v12, g33, KE, h, ElemStat, NodeStat, auxNode, irn, jcn, A, elemX, elemY, Number, NumElem, NumNodes, NodeCoord, Pshp, Pcrds, alpha);
		for(i=0;i<NumNodes;i++)
		{
			if(NodeStat[i] != 0) /*if node isn't out then dofs are not redundant*
			{
				temp = 2 * i;
				redun[temp] = 1;
				redun[temp+1] = 1;
			}
		}
		/*dofs associated with aux nodes are never redundant*
		for(i=(2 * NumNodes);i<MatrixOrder;i++)
		{
			redun[i] = 1;
		}
	}
	else*/
	/*printf("\nC");
	printf("\nC");*/


	if(Mflag == 2)
	{

	for(i=0;i<NumElem;i++)
	{
		alpha[i] = i;
	}
		/*Calculate the global stiffness matrix in index format, ready for the MA57 solver - using AFG*/
		AFG_Matrix(e1, e2, v12, v21, g12, g23, g13, KE, h, ElemStat, NodeStat, auxNode, irn, jcn, A, elemX, elemY, Number, NumElem, NumNodes, NodeCoord, alpha, redun, hxz, theta, numply, plyt, Trans);
        
        
        //printf("KEEEEEE :  %f", KE[0]);
    
    }

	/*Calculate total structure volume (area) by summing element area array*/
	VolA[itt] = 0;
	for(i=0;i<elemX;i++)
	{
		for(j=0;j<elemY;j++)
		{
			VolA[itt] += h*hxz[i]; /*PolyArea2(lsf,elemX,elemY,Number,NumNodes,NodeCoord,NodeStat, ElemStat, alpha0, hxz, h, alpha);*/ /*All elements are in for now*/
		}
	}

	/*---------------------Solve the FE equation using the MA57 solver-----------------*/
	/*if(OjFlag == 1) /*if Objective is deterministic compliance -> solve deterministic f=Ku*
	{*/
    
    
		IFG_Solve(irn, jcn, A, TotFix2, FixDofs2, rhs_det, NumEntries, MatrixOrder, redun, 1, itt);
    
    
	/*}*/
	
	/*else /*Otherwise find displacements for multiple load cases*
	{
		IFG_Solve(irn, jcn, A, TotFix2, FixDofs2, rhs_unc, NumEntries, MatrixOrder, redun, NumRhs);
	}*/
	
	minX = 0.0; /*Use for variance part of objective*/
	/*if( (OjFlag == 3) || (OjFlag == 4) ) /*if objective involves variance of compliance*
	{
		rhs_var = calloc((MatrixOrder * NumRhs), sizeof(double)); /*allocate memory for adjoint force vector*
		/*construct adjoint load vectors using "co-variance" matrix entries and multiple load case displacements*
		for(j=0;j<NumRhs;j++)
		{
			ftemp = (j == 0) ? 4.0 : 2.0; /*factor to scale the load cases for variance computation*
			for (i=0;i<NumVar;i++)
			{
				temp = (2 * var_ind[i]) + (j * MatrixOrder);	/* X dof for j'th rhs*
				k = 3 * i;	/*place indicator for var_val array*
				/*printf("\nDisplacements for force %i are: Ux=%lf, Uy=%lf",i+1,rhs_det[temp],rhs_det[temp+1]);*
				rhs_var[temp] = ftemp * ( (rhs_unc[temp] * var_val[k]) + (rhs_unc[temp+1] * var_val[k+2]) );
				rhs_var[temp+1] = ftemp * ( (rhs_unc[temp+1] * var_val[k+1]) + (rhs_unc[temp] * var_val[k+2]) );
			}
		}
		
		for(i=0;i<NumRhs;i++)
		{
			ftemp = unc_wgt[i]; /* weighting for load case *
			for(j=0;j<MatrixOrder;j++)
			{
				k = j + (MatrixOrder * i); /*rhs array indicator*
				if(fabs(rhs_var[k]) > 0.000001)
				{
					/*printf("\nAdjoint Force at dof %i for load case %i = %lf",j+1,i+1,rhs_var[k]);*
					minX += ftemp * rhs_var[k] * rhs_unc[k]; /*compute analytical variance*
				}
			}
		}
		minX *= 0.5 * wgt_combo; /*need to divide by 2 because sensitivity twice objective*
		printf("\nVariance (normalised for combined objective) = %lf",minX);
		
		/*Now solve the adjoint system*
		IFG_Solve(irn, jcn, A, TotFix2, FixDofs2, rhs_var, NumEntries, MatrixOrder, redun, NumRhs);
	}*/

	/*sprintf(plotname,"alpha%i.txt",itt); /*set name for current output file*
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Matlab Area ratio writefile\n");
		}
	else{

		for(i=0;i<elemX;i++)
		{
			for(j=0;j<elemY;j++)
			{
				num = Number[i][j].n-1;
				fprintf(outfile,"%lf\t",alpha[num]/**alpha0[num]*);
			}
			fprintf(outfile,"\n");
		}
	}
	
	fclose(outfile);
	printf("\nMatlab Area Ratio File written");*/


	free(irn);
	free(jcn);
	free(A);
	free(redun);
	
	
	/* Write Element Area Ratio information file
	sprintf(plotname,"AreaRatio%i.txt",itt); /*set name for current output file
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Area ratio plotting writefile\n");
		}
	else{

		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n - 1;	/*Element number, -1 for C array storage
				fprintf(outfile,"%lf\t",alpha[num]);
			}
			fprintf(outfile,"\n");
		}
	}
	
	fclose(outfile);

	printf("\nArea Ratio plotting File written");*/
	
	/* Write Element Area Ratio information file for matlab*/
	
	/*------------------------------------------------------------------------------------------------------------
	/
	/		Section 5.2: Calculate boundary strain energy, velocities and lambda
	/
	/			Inputs: * Material properties and element edge length & maximum domain co-ordinates from Section 1
	/					* Analysis method type and Volume constraint from Section 1
	/					* Element edge length from Section 1
	/					* Number of: Elements and Nodes from Section 2
	/					* Node co-ordinates & Mesh numbering data from Section 2
	/					* Nodes with fixed zero values for the signed distance function from Section 3
	/					* Current signed distance function and narrow band info from either Section 4 or 5.3/4
	/					* Node status & Number of Auxillary nodes and co-ordinates from Section 5.1
	/					* Number of each element type [IN,Pent,Quad,Tri,OUT] from Section 5.1
	/					* Nodal displacements from Section 5.1
	/					* Element areas (or area ratios) and total structure area from Section 5.1
	/					* Pentagonal co-ordinates and shape functions (for BFG method) from Section 5.1
	/					* Current (or inital) lsfHole values from Section 5.3/4
	/			
	/			Outputs: * Boundary normal velocities
	/					 * Total strain energy
	/					 * Maximum boundary velocity, time step, lambda and lamT (for hole creation method)
	/					 * Updated lsfHole values (for hole insertion method)
	/-------------------------------------------------------------------------------------------------------------*/
	
	/*--------------------------------------Calculate boundary node strain energy------------------------------------*/
	Eenergy = calloc(NumElem,sizeof(double));	/*Dynamic allocation of element strain energy array*/
	Nenergy = calloc(Ntot,sizeof(double));	/*Dynamic allocation of node strain energy array (including auxillary nodes)*/
	
	/*if(Mflag == 1)
	{		
		/*work out number of gauss / int points for strain energy calc*
		temp = (Itotal + qCount) * 4;
		temp += pCount * 5;
		temp += tCount;
		printf("\nTotal Number of Sensitivity evaluations =%i", (temp * NumRhs));	
		dist = 0.0;
		
		/*calcualte nodal strain energy using 1st order least squares of intergration point values*
		if(OjFlag == 1) /*deterministic objective*
		{
			IFG_StrainLS(elemX, elemY, Number, ElemStat, NodeStat, alpha, e11, v12, g33, h, auxNode, Nenergy, Pshp, Pcrds, rhs_det, rhs_det, NumNodes, NodeCoord, temp, Ntot,NIOtotal,Eenergy, Hflag, 1.0);
			
			/*Calculate total compliance by summing element value array*
			for(j=0;j<NumElem;j++)
			{
				dist += Eenergy[j];
			}
			printf("\nDeterministic Compliance = %lf",dist);
		}
		
		/*if( (OjFlag == 2) || (OjFlag == 4) ) /*objective involves expected compliance*
		{
			printf("\nExpected Compliance (normalised for combined objective):");
			for(i=0;i<NumRhs;i++)
			{
				j = MatrixOrder * i;
				ftemp = unc_wgt[i]; /* weighting for laod case *
				rhs_pnt = &rhs_unc[j]; /*set pointer to start of next rhs*
			
				IFG_StrainLS(elemX, elemY, Number, ElemStat, NodeStat, alpha, e11, v12, g33, h, auxNode, Nenergy, Pshp, Pcrds, rhs_pnt, rhs_pnt, NumNodes, NodeCoord, temp, Ntot,NIOtotal,Eenergy, Hflag, ftemp);
				
				/*Calculate total compliance by summing element value array*
				ftemp = 0.0;
				for(j=0;j<NumElem;j++)
				{
					ftemp += Eenergy[j];
				}
				printf("\nAfter %i load cases = %lf",i+1,ftemp);
			}
			dist = ftemp;
		}
		
		if( (OjFlag == 3) || (OjFlag == 4) ) /*objective involves compliance variance*
		{
			for(i=0;i<NumRhs;i++)
			{
				j = MatrixOrder * i;
				ftemp = unc_wgt[i] * wgt_combo; /* weighting for load case *
				rhs_pnt = &rhs_unc[j]; 
				rhs_pnt2 = &rhs_var[j]; /*set pointers to start of next rhs*
			
				IFG_StrainLS(elemX, elemY, Number, ElemStat, NodeStat, alpha, e11, v12, g33, h, auxNode, Nenergy, Pshp, Pcrds, rhs_pnt, rhs_pnt2, NumNodes, NodeCoord, temp, Ntot,NIOtotal,Eenergy, Hflag, ftemp);
			}
			
			/*Work out current Objective value*
			dist += minX; /*Add variance to expected value (or just copy variance for OjFlag == 3) *
			if(OjFlag == 4)
			{
				printf("\nNormalised Combined Objective = %lf",dist);
			}
		}
		
		if(pCount != 0)
		{
			free(Pshp);
			free(Pcrds);
		}
	}
	
	else*/
	 if(Mflag == 2)
	{	
		temp = (NumElem) * 4; /*number of gauss points*/
		printf("\nTotal Number of Sensitivity evaluations =%i", (temp * NumRhs));	
		dist = 0.0;
		/*calcualte nodal strain energy using 1st order least squares of intergration point values*/
		if(OjFlag == 1) /*deterministic objective*/
		{
			AFG_Strain_LS(elemX,elemY,Number,ElemStat,NodeStat,alpha,e1, e2,v12,v21,g12,g23,g13,h,auxNode,Nenergy,rhs_det, rhs_det,NumNodes, NodeCoord, temp, Ntot, NIOtotal, Eenergy, Hflag, 1.0, hxz, hz, theta, Trans, numply, plyt);
            
            
            
            
			
			/*Calculate total compliance by summing element value array*/
			for(j=0;j<NumElem;j++)
			{
				dist += Eenergy[j];
                
                
                
                if (j==80) {
                    double vvvvv;
                    vvvvv=Eenergy[j];
                    printf("\nDeterministic Compliance of number 40= %.16lf",vvvvv);
                }
			}
			printf("\nDeterministic Compliance = %.16lf",dist);
			}
        
       // double dist_old;
        
         //   dist_old=dist;
        

        
        
        
        
        
        
        
			sprintf(plotname,"ElemStn%i.txt",itt);
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
						fprintf(outfile,"%0.9f\t",0.5 * Eenergy[num]);
					}
					fprintf(outfile,"\n");
				}
			}
			
			fclose(outfile);
	
			printf("\nElem Energy plotting File written");/**/

		
		/*if( (OjFlag == 2) || (OjFlag == 4) ) /*objective involves expected compliance*
		{
			printf("\nExpected Compliance (normalised for combined objective):");
			for(i=0;i<NumRhs;i++)
			{
				j = MatrixOrder * i;
				ftemp = unc_wgt[i]; /* weighting for laod case *
				/*printf("\nweight = %lf",ftemp);*
				rhs_pnt = &rhs_unc[j]; /*set pointer to start of next rhs*
				
				AFG_Strain_LS(elemX, elemY, Number, ElemStat, NodeStat, alpha, e11, v12, g33, h, auxNode, Nenergy, rhs_pnt, rhs_pnt, NumNodes, NodeCoord, temp, Ntot, NIOtotal, Eenergy, Hflag, ftemp);
				
				/*Calculate total compliance by summing element value array*
				ftemp = 0.0;
				for(j=0;j<NumElem;j++)
				{
					ftemp += Eenergy[j];
				}
				printf("\nAfter %i load cases = %lf",i+1,ftemp);
			}
			dist = ftemp;
		}*/
		
		/*if( (OjFlag == 3) || (OjFlag == 4) ) /*objective involves compliance variance*
		{
			for(i=0;i<NumRhs;i++)
			{
				j = MatrixOrder * i;
				ftemp = unc_wgt[i] * wgt_combo; /* weighting for laod case *
				/*printf("\nweight for laod case %i = %lf",i+1,ftemp);*
				rhs_pnt = &rhs_unc[j]; 
				rhs_pnt2 = &rhs_var[j]; /*set pointers to start of next rhs*

				AFG_Strain_LS(elemX, elemY, Number, ElemStat, NodeStat, alpha, e11, v12, g33, h, auxNode, Nenergy, rhs_pnt, rhs_pnt2, NumNodes, NodeCoord, temp, Ntot, NIOtotal, Eenergy, Hflag, ftemp);
				
				/*Calculate variance by summing element value array*
				ftemp = 0.0;
				for(j=0;j<NumElem;j++)
				{
					ftemp += Eenergy[j];
				}
				printf("\nAfter %i load cases = %lf",i+1,ftemp);
			}
			
			/*Work out current Objective value*
			dist += minX; /*Add variance to expected value (or just copy variance for OjFlag == 3) *
			if(OjFlag == 4)
			{
				printf("\nNormalised Combined Objective = %lf",dist);
			}
			
			/*Check variance
			ftemp = 0.0;
			for(j=0;j<NumElem;j++)
			{
				ftemp += Eenergy[j];
			}
			printf("\nVariance check = %lf",0.5*ftemp);*
		}*/
	}
	
	CompA[itt] = dist; /* Store total objective value*/

	/*----------------------------------------------------------------------------------------------------------------------------------*/
	/*----------------------------------------------------------------------------------------------------------------------------------*/
	/*			IN THE FOLLOWING SECTION THE SENSITVITY OF THE COMPLIIANCE TO THE FIBRE ANGEL (AS DEFINED BY THE LSF FUNCTION)			*/
	/*----------------------------------------------------------------------------------------------------------------------------------*/

	/*DiffThetaLsf = malloc(NumNodes*sizeof(double));*/

	EF = malloc(NumNodes*sizeof(double));

	/*Differentiate theta wrt to the level set function*/
	/*DthetaDlsfCalcOLDVERSION(DiffThetaLsf, elemX, elemY, Number, NumNodes, NodeCoord, lsf2, h, hxz, NodeX, NodeY, Nodes2);*/

	DTL = malloc(Numlsf*NumNodes*sizeof(NElemDiff));
	DTLb = malloc(NumNodes*sizeof(NElemDiff));
    double *NodeDTL;
    NodeDTL= malloc(Ntot_new*sizeof(double));
    
    
    double *Nodesenpp;
    double *NodeThetaStn;
    NodeThetaStn = malloc(Ntot_new*sizeof(double));
    Nodesenpp = malloc(Ntot_new*sizeof(double));
    Gstn *NodeDTLg;
    Gstn *Nodeseng;
    NodeDTLg = malloc(4*elemX*elemY*sizeof(Gstn));
    Nodeseng = malloc(4*elemX*elemY*sizeof(Gstn));
    
for(p=0;p<Numlsf;p++)
{
	for(i=0;i<NumNodes;i++)	{k = p*NumNodes+i;	lsf2b[i] = lsf2[k];
        printf("\nlsf value is %f", lsf2b[i]);
    }


    
	DthetaDlsfCalc(DTLb, elemX, elemY, Number, lsf2b, h, hxz, NumNodes, NodeCoord,Ntot_new, AuxNodes, NodeDTL, NodeDTLg);
    
	for(i=0;i<NumNodes;i++)	
	{
		k = p*NumNodes+i;
		DTL[k].ai = DTLb[i].ai;	DTL[k].an = DTLb[i].an;  DTL[k].ad = DTLb[i].ad;
		DTL[k].bi = DTLb[i].bi;	DTL[k].bn = DTLb[i].bn;  DTL[k].bd = DTLb[i].bd;
		DTL[k].ci = DTLb[i].ci;	DTL[k].cn = DTLb[i].cn;  DTL[k].cd = DTLb[i].cd;
		DTL[k].di = DTLb[i].di;	DTL[k].dn = DTLb[i].dn;  DTL[k].dd = DTLb[i].dd;
        
        printf("\nDTL value is %f", DTL[k].ad);
	}
	
}

	/*Now calculate the Energy factor.*/
	EnergyFactor(elemX, elemY, Number, NodeX, NodeY, Nodes2, h, hxz, rhs_det, rhs_old, theta, theta_old, e1, e2, v12, v21, g12, g23, g13, EF, Trans, numply, plyt, itt);
	
	if(itt>0)
	{
		free(rhs_old);
		free(theta_old);
	}
	ElemThetaStn = malloc(elemX*elemY*sizeof(double));					/*Don't forget to free this value*/
	ElemThetaStnB = malloc(elemX*elemY*sizeof(double));
	/*Now calculate the element sensitvities*/

    
	SensitivityCalculation(elemX, elemY, Number, NumNodes, NodeCoord, h, hxz, NodeX, NodeY, Nodes2, rhs_det, Trans, numply, plyt, theta, e1, e2, v12, v21, g12, g23, g13, ElemThetaStn, ElemThetaStnB, EF, itt, Ntot_new, AuxNodes, NodeThetaStn, Nodeseng);
    
    double *NodesenTemp;
    NodesenTemp = malloc(Ntot_new*sizeof(double));
    Gstn Sensofc[elemX*elemY*4];  // store sensitivity at each element.
    
    // got the sensitivity at each node point of each element, not in Gausspoint
    for (i=0; i<elemX*elemY; i++) {
        Sensofc[i*4].u  =ElemThetaStn[i]*NodeDTLg[i*4].u;
        Sensofc[i*4+1].u=ElemThetaStn[i]*NodeDTLg[i*4].u;
        Sensofc[i*4+2].u=ElemThetaStn[i]*NodeDTLg[i*4].u;
        Sensofc[i*4+3].u=ElemThetaStn[i]*NodeDTLg[i*4].u;
        
    }
    
    
    
    
    Gstn *GPsens;
    
    GPsens = malloc(4*elemX*elemY*sizeof(Gstn));
    
    int kk;
    // Get the ture sensitivity at each single node point, consist of de/dtheta dtheta/dfi.....
    // so, we should consider multi fiber path in here now.
    //de/dfi is different in different set of mesh, so we now numlsf set of de/dfi which is nodesesntemp
    int nodecounter;
    int auxnodecounter;
    auxnodecounter=0;
    nodecounter=0;
for(p=0;p<Numlsf;p++)
    {
       
        for(i=0;i<NumNodes;i++)
        {
            k = p*NumNodes+i;
            DTLb[i].ai = DTL[k].ai;	DTLb[i].an = DTL[k].an;  DTLb[i].ad = DTL[k].ad;
            DTLb[i].bi = DTL[k].bi;	DTLb[i].bn = DTL[k].bn;  DTLb[i].bd = DTL[k].bd;
            DTLb[i].ci = DTL[k].ci;	DTLb[i].cn = DTL[k].cn;  DTLb[i].cd = DTL[k].cd;
            DTLb[i].di = DTL[k].di;	DTLb[i].dn = DTL[k].dn;  DTLb[i].dd = DTL[k].dd;
        }
        for(i=0;i<NumElem;i++)	{k = p*NumElem+i;	ElemStatb[i] = ElemStat[k];}
        
        NodeSen(NumNodes, DTLb, ElemStatb,  NodesenTemp, ElemThetaStn);

    for (kk=0; kk<(elemX+1)*(elemY+1); kk++) {
        //printf("\n Node sens %d: %f", kk, NodesenTemp[kk]);
    }
        
    for (i=0; i<NumofAux[p]; i++) {
            k=auxnodecounter+i;
            AuxNodesb[i].n=AuxNodes[k].n;
            AuxNodesb[i].x=AuxNodes[k].x;
            AuxNodesb[i].xz=AuxNodes[k].xz;
            AuxNodesb[i].y=AuxNodes[k].y;
            AuxNodesb[i].z=AuxNodes[k].z;
    }
        auxnodecounter=auxnodecounter+NumofAux[p];
    
    //Nodesensitivity( elemX,  elemY, Number,  NumNodes, NodeCoord,h, hxz, Ntot_new, AuxNodes, NodeDTL, NodesenTemp, NodeDTLg);
   
    //NodesenP(elemX, elemY, Number,NumNodes,NodeCoord,h, hxz,  Ntot_new, AuxNodes, NodeDTL, NodesenTemp, NodeDTLg, Nodeseng, Nodesenpp);

    // use the shape function to get the sensitivity at Gauss point
    //store the sensitivity at GPsens
    
    GPsensitivity( elemX, elemY, Number,  NumNodes, NodeCoord, h, hxz,  NumNodes+NumofAux[p], AuxNodesb, NodeDTL, NodesenTemp, GPsens);
    
    // Get the sensitivity at Auxulary node and node point as well, the later one can used to compare with the sensitivity we got earlier to check if it is correct.
    // store the sensitivit at  Nodesenpp
    double *Nodesenpptemp = malloc((NumNodes+NumofAux[p])*sizeof(double));
    
        
    
        
    NodeAuxsens( elemX, elemY, Number, NumNodes,  NodeCoord, h,  hxz, NumNodes+NumofAux[p], AuxNodesb,GPsens,Nodesenpptemp);
    
    
    //NodeAuxsens2( elemX, elemY, Number, NumNodes,  NodeCoord, h,  hxz, Ntot_new, AuxNodes,Nodesenpp,NodesenTemp);

    
    
    
    //Get the sensitivity at intersection point, simply add two sensitivity together
    
    //NodeAuxsens3(elemX, elemY,  Number,  NumNodes, NodeCoord,  h,  Ntot_new, AuxNodes, Nodesenpp,NodesenTemp);
    
    // GEt the sens at intersection point, use linear interploation
    // use sens only if the element was intersected, and no Gauss sens
    
    //NodeAuxsens4(elemX, elemY,  Number,  NumNodes, NodeCoord,  h,  Ntot_new, AuxNodes, Nodesenpp,NodesenTemp);
        
        for(i=0;i<NumNodes+NumofAux[p];i++)
        {
            
            k = nodecounter+i;
            Nodesenpp[k] = Nodesenpptemp[i];
        }
        nodecounter=nodecounter+NumNodes+NumofAux[p];
}

    
    
   
    free(NodeDTLg);
    free(Nodeseng);
    free(GPsens);
    
    
   

    for (kk=0; kk<Ntot_new; kk++) {
               // printf("Sensitivity%d: %f and %f", kk,NodeThetaStn[kk], NodeDTL[kk]);


    }
    
	/*Calculate the total sensitvity*/
	SensA[itt] = 0;
	SensB[itt] = 0;
	for(i=0;i<(elemX*elemY);i++)
	{
		SensA[itt] += ElemThetaStn[i];
		/*SensB[itt] += ElemThetaStnB[i];*/
		/*printf("\nabs(ElemThetaStn[%i] = %f \n SensA[%i] = %f", i, fabs(ElemThetaStn[i]), itt, SensA[itt]);*/
	}

	/*sprintf(plotname,"Energyfactor%i.txt",itt); /*set name for current topology output file*
	/* Write Node signed distance value information file post re-initialisation*
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		/*Bottom row of Nodes first*
		fprintf(outfile,"%lf\t",EF[Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fprintf(outfile,"%lf\t",EF[Number[i][0].b-1]);
		}
		fprintf(outfile,"\n");
		/*Then the rest of the Nodes*
		for(j=0;j<elemY;j++)
		{
			fprintf(outfile,"%lf\t",EF[Number[0][j].d-1]);
			for(i=0;i<elemX;i++)
			{
				fprintf(outfile,"%lf\t",EF[Number[i][j].c-1]);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);*/

	/*sprintf(plotname,"DTL%i.txt",itt); /*set name for current topology output file*
	/* Write Node signed distance value information file post re-initialisation*
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		fprintf(outfile,"Node\tElem1Index\tNum\tDTheta\tElem2Index\tNum\tDTheta\tElem3Index\tNum\tDTheta\tElem4Index\tNum\tDTheta");
		/*These have to be out put as a list*
		for(i=0;i<NumNodes;i++)
		{
			fprintf(outfile,"\n%i\t%i\t%i\t%f\t%i\t%i\t%f\t%i\t%i\t%f\t%i\t%i\t%f", i, DTL[i].ai, DTL[i].an, DTL[i].ad, DTL[i].bi, DTL[i].bn, DTL[i].bd, DTL[i].ci, DTL[i].cn, DTL[i].cd, DTL[i].di, DTL[i].dn, DTL[i].dd);
		}
	}
	fclose(outfile);*/

	/*----------------------------------------------------------------------------------------------------------------------------------*/
	/*-------------------------Here we use the sensitivity values to update the lsf functions-------------------------------------------*/

	DeltaTheta = malloc(elemX*elemY*sizeof(double));
	Deltalsf2 = malloc(Numlsf*NumNodes*sizeof(double));
	Deltalsf2b = malloc(NumNodes*sizeof(double));
    
    
    double *DeltaT;
    DeltaT = malloc(1*sizeof(double));
	/*First get the change in fibre angle theta desired in each element, scaled from the sensitivities and the user set deltathetamax*/
	//ScaleDeltaTheta(elemX, elemY, ElemThetaStn, DeltaTheta, DeltaThetaMax, itt, CompA, Number,ElemStat, ElemStatb, Numlsf);

	/*Now multiply each of the nodal sensitvites of lsf to a change in fibre angel by the desired change in fibre angle to get the change in the lsf value.*/

    nodecounter=0;
    auxnodecounter=0;
for(p=0;p<Numlsf;p++)
{
	for(i=0;i<NumNodes;i++)	
	{
		k = p*NumNodes+i;
		DTLb[i].ai = DTL[k].ai;	DTLb[i].an = DTL[k].an;  DTLb[i].ad = DTL[k].ad;
		DTLb[i].bi = DTL[k].bi;	DTLb[i].bn = DTL[k].bn;  DTLb[i].bd = DTL[k].bd;
		DTLb[i].ci = DTL[k].ci;	DTLb[i].cn = DTL[k].cn;  DTLb[i].cd = DTL[k].cd;
		DTLb[i].di = DTL[k].di;	DTLb[i].dn = DTL[k].dn;  DTLb[i].dd = DTL[k].dd;
	}
	for(i=0;i<NumElem;i++)	{k = p*NumElem+i;	ElemStatb[i] = ElemStat[k];}
    for(i=0;i<NumNodes;i++)	{k = p*NumNodes+i;	NodeStatb[i] = NodeStat[k];}
	/*First get the change in fibre angle theta desired in each element, scaled from the sensitivities and the user set deltathetamax*/
	//ScaleDeltaTheta(elemX, elemY, ElemThetaStn, DeltaTheta, DeltaThetaMax, itt, CompA, Number, ElemStat, ElemStatb, Numlsf, NodeThetaStn, NodeDTL, DeltaT, Ntot_new );

	/*Now multiply each of the nodal sensitvites of lsf to a change in fibre angel by the desired change in fibre angle to get the change in the lsf value.*/
    
    
    Vnorm = malloc((NumNodes+NumofAux[p])*sizeof(double));
    for(i=0;i<NumNodes;i++)	{k = p*NumNodes+i;	lsf2b[i] = lsf2[k];}
    for(i=0; i<NumNodes+NumofAux[p]; i++)
    {
        //Vnorm[i]=NodeThetaStn[i]*NodeDTL[i];
        k=nodecounter+i;
        Vnorm[i]=Nodesenpp[k];
        //Vnorm[i]=NodesenTemp[i];
        printf("\n Vnorm_before %d: %f ", i,  Vnorm[i]);
        //printf("\n %d: %f, %f, %f", i, NodeThetaStn[i], NodeDTL[i], Vnorm[i]);  // print out sensitivityies
       
    }
    for (i=0; i<NumofAux[p]; i++) {
        k=auxnodecounter+i;
        AuxNodesb[i].n=AuxNodes[k].n;
        AuxNodesb[i].x=AuxNodes[k].x;
        AuxNodesb[i].xz=AuxNodes[k].xz;
        AuxNodesb[i].y=AuxNodes[k].y;
        AuxNodesb[i].z=AuxNodes[k].z;
    }
    
    auxnodecounter=auxnodecounter+NumofAux[p];
    nodecounter=nodecounter+NumNodes+NumofAux[p];

    
  
    // put the Vnorm value in mesh nodes which is not the boundary node, equal to 0
    for (i=0; i<NumNodes; i++) {
        
        if (NodeStatb[i]!=2) {
            Vnorm[i]=0.0;
        }
    
    }
    
   // ScaleDeltaTheta(elemX, elemY, ElemThetaStn, DeltaTheta, DeltaThetaMax, itt, CompA, Number, ElemStat, ElemStatb, Numlsf, NodeThetaStn, NodeDTL, DeltaT, Ntot_new, Vnorm );
    
    Vext(NodeX, NodeY,Nodes2,NodeStatb,Vnorm,lsf2b,active,NumAct,NumNodes, NumNodes+NumofAux[p],NodeCoord, h, AuxNodesb, -1 , dt,elemX,hxz);
    Vext(NodeX, NodeY,Nodes2,NodeStatb,Vnorm,lsf2b,active,NumAct,NumNodes, NumNodes+NumofAux[p],NodeCoord, h, AuxNodesb, 1 , dt,elemX,hxz);
    
    // numNodes+NumofAux[p], this was Ntot_new
   
    //ScaleDeltaTheta(elemX, elemY, ElemThetaStn, DeltaTheta, DeltaThetaMax, itt, CompA, Number, ElemStat, ElemStatb, Numlsf, NodeThetaStn, NodeDTL, DeltaT, Ntot_new, NodesenTemp, &SensMaxOld);
    *DeltaT=0.0;
    
    
    ScaleVnorm(elemX, elemY, ElemThetaStn, DeltaTheta, DeltaThetaMax, itt, CompA, Number, ElemStat, ElemStatb, Numlsf, NodeThetaStn, NodeDTL, DeltaT, NumNodes+NumofAux[p], Vnorm, &SensMaxOld);

    
//I want to noemalize the sensiticity to make sure the change of level set value not too werid
    //printf("SensMaxlalal=%f,",SensMaxOld);
    
    for (i=0; i<NumNodes; i++) {
        
        
            Vnorm[i]=Vnorm[i]*(*DeltaT);
        
        
    }
    for (i=0; i<NumNodes; i++) {
        
        //printf("\nlsf: %f, %f", lsf2[i], lsf[i] );
        
    }
    
    Grad = malloc(NumNodes*sizeof(double));// array for gradinet info
    int i2,j2;
    // Calculate lsf gradients & update
    j2=NodeY-1;
    i2=NodeX-1;
    for(j=1;j<j2;j++)
    {
        for(i=1;i<i2;i++)
        {
            k = Nodes2[i][j]; // read in node number
           // printf("kkk: %d", k);
            //if(active[k]) // If node is active (not fixed and within narrow band)
           // {
            Grad[k]=0.0;  // initialise the grad value
                if(fabs(Vnorm[k]) > 1.0e-6) // If normal velocity not zero
                {
                    // work out which sign to use for gradient calculation - upwind scheme
                    temp = (Vnorm[k] < -0.0) ? -1 : 1;
                    //double GradWENO(int Xi, int Yj, int num, double *lsf, int NodeX, int NodeY, int Nodes2[NodeX][NodeY], int sign,double Vn, double dt, double h, int elemX, double hxz[elemX])
                    //Grad[k]=GradWENO(i, j, k, lsf2, NodeX, NodeY, Nodes2, temp,Vnorm[k], dt, h,elemX, hxz);
                    // I don't need this part now...so just comment it out
                                   }
          //  }
        }
    }
    
  /*
    
    for(i=0; i<(elemX+1)*(elemY+1); i++)
    {
        if (i==82||i==123) {
            Vnorm[i]=Vnorm[i];
        }
        else{
            
            Vnorm[i]=0.0;
        
        }
        
        
        // printf("\n %d", active[i]);
         printf("%d, %lf \n", i,Vnorm[i]);
    }
    
*/
    
    
    
    // prepare for the conjugate method
    
    
    if (itt>1)
    {
        
        double betak;
        double tempb,tempb1;
        tempb=0.0;
        tempb1=0.0;
        for (i=0; i<NumNodes; i++) {
            
            //tempb+=grad_old[i]*grad_old[i];
            //tempb1+=Grad[i]*Grad[i];
            
            tempb+=vnorm_old[i]*vnorm_old[i];
            tempb1+=Vnorm[i]*Vnorm[i];
        }
        tempb=sqrt(tempb);
        tempb1=sqrt(tempb1);
        betak=(tempb1/tempb)*(tempb1/tempb);
        
        
        for (i=0; i<NumNodes; i++) {
            // if (fabs(Vnorm[i])>0.000000001) {
            //  Grad[i]=Vnorm[i]/fabs(Vnorm[i]);   // get the sign of this number.
            //}
            
            //Grad[i]=-Grad[i]+betak*(-grad_old[i]);
            //vnorm_temp[i]=-Vnorm[i];
            vnorm_temp[i]=-Vnorm[i]+betak*(-vnorm_old[i]);
            printf("\vnormtemp %d lsf: %f, grad:", i,vnorm_temp[i]);
            
        }
        
        
    }
    else
    {
        
        for (i=0; i<NumNodes; i++) {
            // if (fabs(Vnorm[i])>0.000000001) {
            //  Grad[i]=Vnorm[i]/fabs(Vnorm[i]);   // get the sign of this number.
            //}
            
            //Grad[i]=-Grad[i]+betak*(-grad_old[i]);
            
            vnorm_temp[i]=-Vnorm[i];
            printf("\n vnormtemp %d lsf: %f, grad:", i,vnorm_temp[i]);
            
        }
    
    }
    
    for (i=0; i<NumNodes; i++) {
        //grad_old[i]=Grad[i];
        //if (fabs(Vnorm[i])>0.000000001) {
        vnorm_old[i]=Vnorm[i];
        //vnorm_old[i]=Vnorm[i]/fabs(Vnorm[i]);
        // grad_old[i]=vnorm_old[i]*Grad[i];
        // what I want in this line is the sign of the velocity;  vnorm_old is just the -1 or +1.
        
        // }
        
    }
    
   
    
    DeltaLsfCalc(NumNodes, Deltalsf2b, DeltaTheta, DTLb, ElemStatb, vnorm_temp);
    //DeltaLsfCalc(NumNodes, Deltalsf2b, DeltaTheta, DTLb, ElemStatb,Vnorm);
    
  
    
    for(i=0;i<NumNodes;i++)	{k = p*NumNodes+i;	Deltalsf2[k] = Deltalsf2b[i];}
    free(Vnorm);

    
}
    
    
    

    
    
    
    
	/*Finally update the lsf function values*/
	for(i=0;i<(Numlsf*NumNodes);i++)
	{
		
    /*    if (Deltalsf2[i]>0)
        {
            Deltalsf2[i] = atan(Deltalsf2[i]);
        }
        else if (Deltalsf2[i]<0) {
            
            Deltalsf2[i] = atan(Deltalsf2[i]);
        }
        
     
     */
         //lsf2[i]-= *DeltaT*Deltalsf2[i]*Grad[i];
       // lsf2[i]-= *DeltaT*Deltalsf2[i];
        
        
        //lsf2[i]-= *DeltaT*Deltalsf2[i]*Grad[i];
        //printf("\nVornm: %d : %f, %f, %f, %f", i, *DeltaT*Deltalsf2[i]*Grad[i], *DeltaT, Deltalsf2[i], Grad[i]);
        
        //lsf2[i]-= *DeltaT*NodesenTemp[i]*scale_factor;
        
        //lsf2[i]-= *DeltaT*Deltalsf2[i]*Grad[i];
       // lsf2[i]-= *DeltaT*fabs(Deltalsf2[i])*Grad[i];
        //lsf2[i]+= *DeltaT*Deltalsf2[i];
        /*
        if (0.15*Deltalsf2[i]>0.8)
        {
            lsf2[i]+= 0.8;
            printf("lsf value just too large!!!!!!");
        }
        else if (0.15*Deltalsf2[i]<-0.8) {
            
            lsf2[i]+= -0.8;
            printf("lsf value just too large!!!!!!");
        }
        else{
        
            lsf2[i]+= 0.15*Deltalsf2[i];
        }
        */
        lsf2[i]+= 0.15*Deltalsf2[i];
        
        //lsf2[i]+= *DeltaT*fabs(Deltalsf2[i])*Deltalsf2[i];
        double Templsf;
        Templsf=0.0;
        //Templsf=*DeltaT*Deltalsf2[i]*Grad[i];
        if (fabs(Templsf)-0.1>0.0) {
            //Templsf=Templsf*0.25;
        }
        //lsf2[i]+=Templsf;
        //lsf2[i]+= (*DeltaT)*Deltalsf2[i]*Grad[i];
        
        printf("\ndelta %d lsf: %f, grad:", i,Deltalsf2[i]*0.1);
        
        //printf("\ndelta %d lsf: %f, grad: %f", i,(*DeltaT)*Deltalsf2[i]*Grad[i],Grad[i]);
	}
    
   
    
    
    
    
	for(i=0;i<(elemX*elemY);i++)
	{
	//	SensB[itt] += DeltaTheta[i];
	}

	/*Reintiailise leveset function to update all the other values*/
	lsf_old = malloc(NumNodes * sizeof(double));
    
for(p=0;p<Numlsf;p++)
{
	for(i=0;i<NumNodes;i++)	{k = p*NumNodes+i;	lsf2b[i] = lsf2[k];}

	ReIntlsf2(NodeX, NodeY, Nodes2, NumNodes, lsf2b, h, -1, lsf_old, fixed, NumFix,  elemX, hxz);
	ReIntlsf2(NodeX, NodeY, Nodes2, NumNodes, lsf2b, h, 1, lsf_old, fixed, NumFix,  elemX, hxz);
	ArrayCopy(NumNodes, lsf_old, lsf2b);

	for(i=0;i<NumNodes;i++)	{k = p*NumNodes+i;	lsf2[k] = lsf2b[i];}
}
	free(lsf_old);


    

    
    
    
    
    

    
    
    
    
    
    
    
for(p=0;p<Numlsf;p++)
{
	sprintf(plotname,"PlyPlot%ilsf%i.txt",itt+1,p+1);
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		k = (NumNodes*p);

		fprintf(outfile,"%lf\t",lsf2[k+Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fprintf(outfile,"%lf\t",lsf2[k+Number[i][0].b-1]);
		}
		fprintf(outfile,"\n");

		for(j=0;j<elemY;j++)
		{
			fprintf(outfile,"%lf\t",lsf2[k+Number[0][j].d-1]);
			for(i=0;i<elemX;i++)
			{
				fprintf(outfile,"%lf\t",lsf2[k+Number[i][j].c-1]);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);
}

	sprintf(plotname,"PlyPlot%icom.txt",itt+1); /*set name for current topology output file*/	
	/* Write Node signed distance value information file post re-initialisation*/
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		/*Bottom row of Nodes first*/
		ftemp = 10000;
		for(p=0;p<Numlsf;p++)
		{
			k = (NumNodes*p);
			ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[0][0].a-1])) ? lsf2[k+Number[0][0].a-1]:ftemp;
		}
		fprintf(outfile,"%lf\t",ftemp);
		for(i=0;i<elemX;i++)
		{
			ftemp = 10000;
			for(p=0;p<Numlsf;p++)
			{
				k = (NumNodes*p);
				ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[i][0].b-1])) ? lsf2[k+Number[i][0].b-1]:ftemp;
			}

			fprintf(outfile,"%lf\t",ftemp);
		}
		fprintf(outfile,"\n");
		/*Then the rest of the Nodes*/
		for(j=0;j<elemY;j++)
		{
			ftemp = 10000;
			for(p=0;p<Numlsf;p++)
			{
				k = (NumNodes*p);
				ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[0][j].d-1])) ? lsf2[k+Number[0][j].d-1]:ftemp;
			}
			fprintf(outfile,"%lf\t",ftemp);
			for(i=0;i<elemX;i++)
			{
				ftemp = 10000;
				for(p=0;p<Numlsf;p++)
				{
					k = (NumNodes*p);
					ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[i][j].c-1])) ? lsf2[k+Number[i][j].c-1]:ftemp;
				}
				fprintf(outfile,"%lf\t",ftemp);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);


	/*sprintf(plotname,"PlyDeltaPlot%i.txt",itt);
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		fprintf(outfile,"%lf\t",Deltalsf2[Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fprintf(outfile,"%lf\t",Deltalsf2[Number[i][0].b-1]);
		}
		fprintf(outfile,"\n");
		for(j=0;j<elemY;j++)
		{
			fprintf(outfile,"%lf\t",Deltalsf2[Number[0][j].d-1]);
			for(i=0;i<elemX;i++)
			{
				fprintf(outfile,"%lf\t",Deltalsf2[Number[i][j].c-1]);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);

	/*----------------------------------------------------------------------------------------------------------------------------------*/
	
	/*---------Store theta and displacement from this iteration for the energy factor calcualtion------------*/
	rhs_old = malloc(MatrixOrder*sizeof(double));
	theta_old = malloc(elemX*elemY*sizeof(double));
	
	for(i=0;i<MatrixOrder;i++)
	{
		rhs_old[i] = rhs_det[i];
	}

	for(i=0;i<elemX*elemY;i++)
	{
		theta_old[i] = theta[i];
	}

	free(DTL);
	free(DTLb);
	free(EF);
	free(ElemThetaStn);
	free(ElemThetaStnB);
	free(DeltaTheta);
	free(Deltalsf2);
    
	free(Deltalsf2b);
    free(NodesenTemp);
    free(Grad);
    
	free(ElemStat);
	free(ElemStatb);
	free(auxNode);
	free(rhs_det);
	free(alpha);
	free(NodeStat);
    free(NodeStatb);
	free(Eenergy);
	free(Nenergy);
    
   // free(AuxNodes);
    
    //free(bptr);     // Just simply commment out this line, I don't know if it will affect the programm,
                        // July 2, 17:56
    
    
    //free(na_conn_ind);
    free(NumBound);
    free(NumAux);
    free(NodeThetaStn);
    free(NodeDTL);
    free(DeltaT);
    free(Nodesenpp);
  
    
    //free(na_conn);

	FCS[itt] = FibreContScore(NumNodes, NodeCoord,elemX, elemY, Number, h, hxz, 10.0, theta, itt);

	//ThetaUpdate3(theta, elemX, elemY, Number, lsf2, h, hxz);
    
	ThetaUpdate4(theta, elemX, elemY, Number, lsf2, lsf2b, Numlsf, h, hxz);
	printf("\nFCS = %f", FCS[itt]);

	sprintf(plotname,"Theta%i.txt",itt+1); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
	}
	else{
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n-1;
				fprintf(outfile,"%lf\t",(theta[num]*(180/3.141592654)));
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);

	sprintf(plotname,"ThetaPlot%i.txt",itt+1); 
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
	}
	else{
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n-1;
				Cx = NodeCoord[Number[i][j].a-1].xz + (hxz[i]/2);
				Cy = NodeCoord[Number[i][j].a-1].y + (h/2);
				fprintf(outfile,"%lf\t %lf\n",((Cx+0.5)-0.4*cos(theta[num])),((Cy+0.5)-0.4*sin(theta[num])));
				fprintf(outfile,"%lf\t %lf\n\n",((Cx+0.5)+0.4*cos(theta[num])),((Cy+0.5)+0.4*sin(theta[num])));
			}
		}
	}
	fclose(outfile);

	/*for(i=0;i<(elemX*elemY);i++)
	{
		theta[i] += (3.141592654/180);
	}*/


	/*--------------Small section for interactive control of number of iterations to perform--------------*/
	if(itt >= numitt)
	{
		printf("\nDo you wish to continue (1=yes): ");
		/*scanf("%s",incheck); /*gets user input from the terminal*/
		/*check input was an integer*/
		for(i=0;incheck[i] >= '0' && incheck[i] <= '9'; i++);

		if(incheck[i] != '\0' || i == 0)
		{
			cont = -1;
		}
		else	{
			sscanf(incheck,"%i",&cont); /*convert input to integer*/
		}

		if(cont == 1)
		{
			do
			{
				temp = 1;
				printf("\nFor how many extra iterations (1 -> %i): ", NumNodes);
				scanf("%s",incheck); /*gets user input from the terminal*/
				for(i=0;incheck[i] >= '0' && incheck[i] <= '9'; i++);
				if(incheck[i] != '\0' || i == 0)
				{
					j = -1;
				}
				else
				{
					sscanf(incheck,"%i",&j); /*convert input to integer*/
				}

				/*If volume ratio out of range then get user to re-enter*/
				if((j < 1) || (j > NumNodes)){
					printf("\nNumber out of range!!");
					temp = 0;
					}
				else {
					numitt += j;
				}
			}
			while (temp == 0);
		}
	}
	
	else	{
		cont = 1;
		}
	
	if(cont == 1) /*if continuing, increase size of data storage arrays*/
	{
		itt++; /*Update itterartion number*/
		temp = itt+1;
		EmaxA = realloc(EmaxA, temp * sizeof(Astn));	/*Array for maximum, minimum boundary strain energy & lambda*/
		VolA = realloc(VolA, temp * sizeof(double));	/*Array for total volume (area) of structure for each iteration*/
		CompA = realloc(CompA, temp * sizeof(double));	/*Array for objective for each iteration*/
		SensA = realloc(SensA, temp * sizeof(double));
		SensB = realloc(SensB, temp * sizeof(double));
		FCS = realloc(FCS, temp * sizeof(double));
	}
	/*double X,Y;
	X = 7/0;*/
} while (cont == 1); /*If solution hasn't converged after a number of itteration equal to the number of elements, then abort*/
printf("\nSolution stopped after %i iterations\nTotal Objective Value = %.16lf\n",itt,CompA[itt]);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// end of loop
    

/* Write Node signed distance value information file post re-initialisation*/
sprintf(plotname,"%s_PlotShpFinal.txt",datafile); /*set name for output file*/	
outfile = fopen(plotname, "w");
if(outfile == NULL){
	printf("\nFailed to open Final Implicit writefile\n");
	}
else{
	for(i=0;i<NumNodes;i++)
	{
		fprintf(outfile,"%lf\n",lsf2[i]);
	}
}

fclose(outfile);
printf("\nFinal Signed Distance Info File written\n");

	sprintf(plotname,"XThetaPlotFinal.txt"); 

	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
	}
	else{

		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{

				num = Number[i][j].n-1;
				Cx = NodeCoord[Number[i][j].a-1].xz + (hxz[i]/2);
				Cy = NodeCoord[Number[i][j].a-1].y + (h/2);
				fprintf(outfile,"%lf\t %lf\n",((Cx+0.5)-0.4*cos(theta[num])),((Cy+0.5)-0.4*sin(theta[num])));
				fprintf(outfile,"%lf\t %lf\n\n",((Cx+0.5)+0.4*cos(theta[num])),((Cy+0.5)+0.4*sin(theta[num])));
			}

		}
	}
	fclose(outfile);

for(p=0;p<Numlsf;p++)
{
	sprintf(plotname,"XPlyPlotFlsf%i.txt",p+1); /*set name for current topology output file*/	
	/* Write Node signed distance value information file post re-initialisation*/
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		k = (NumNodes*p);
		/*Bottom row of Nodes first*/
		fprintf(outfile,"%lf\t",lsf2[k+Number[0][0].a-1]);
		for(i=0;i<elemX;i++)
		{
			fprintf(outfile,"%lf\t",lsf2[k+Number[i][0].b-1]);
		}
		fprintf(outfile,"\n");
		/*Then the rest of the Nodes*/
		for(j=0;j<elemY;j++)
		{
			fprintf(outfile,"%lf\t",lsf2[k+Number[0][j].d-1]);
			for(i=0;i<elemX;i++)
			{
				fprintf(outfile,"%lf\t",lsf2[k+Number[i][j].c-1]);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);
}

	sprintf(plotname,"XPlyPlotFcom.txt"); /*set name for current topology output file*/	
	/* Write Node signed distance value information file post re-initialisation*/
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
		}
	else
	{
		/*Bottom row of Nodes first*/
		ftemp = 10000;
		for(p=0;p<Numlsf;p++)
		{
			k = (NumNodes*p);
			ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[0][0].a-1])) ? lsf2[k+Number[0][0].a-1]:ftemp;
		}
		fprintf(outfile,"%lf\t",ftemp);
		for(i=0;i<elemX;i++)
		{
			ftemp = 10000;
			for(p=0;p<Numlsf;p++)
			{
				k = (NumNodes*p);
				ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[i][0].b-1])) ? lsf2[k+Number[i][0].b-1]:ftemp;
			}

			fprintf(outfile,"%lf\t",ftemp);
		}
		fprintf(outfile,"\n");
		/*Then the rest of the Nodes*/
		for(j=0;j<elemY;j++)
		{
			ftemp = 10000;
			for(p=0;p<Numlsf;p++)
			{
				k = (NumNodes*p);
				ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[0][j].d-1])) ? lsf2[k+Number[0][j].d-1]:ftemp;
			}
			fprintf(outfile,"%lf\t",ftemp);
			for(i=0;i<elemX;i++)
			{
				ftemp = 10000;
				for(p=0;p<Numlsf;p++)
				{
					k = (NumNodes*p);
					ftemp = (fabs(ftemp)>fabs(lsf2[k+Number[i][j].c-1])) ? lsf2[k+Number[i][j].c-1]:ftemp;
				}
				fprintf(outfile,"%lf\t",ftemp);
			}
			fprintf(outfile,"\n");
		}

	}
	fclose(outfile);


	sprintf(plotname,"XThetaFINAL.txt",itt+1); /*set name for current topology output file*/
	/* Write Theta value information file*/
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Re-Int Implicit writefile\n");
	}
	else{
		/*fprintf(outfile,"Node Num\tx\ty\tphi\n"); /*column headings*/
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				num = Number[i][j].n-1;
				fprintf(outfile,"%lf\t",(theta[num]*(180/3.141592654)));
			}
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);

/* Write Objetcive and Constraint convergence, per iteration, info file*/
sprintf(plotname,"%s_Convergence.txt",datafile); /*set name for output file*/	
outfile = fopen(plotname, "w");
if(outfile == NULL){
	printf("\nFailed to open TotStn writefile\n");
	}
else{
	for(i=0;i<=itt;i++)
	{
		fprintf(outfile,"%i\t",i); /*iteration*/
		fprintf(outfile,"%lf\t",CompA[i]); /*Objective*/
		fprintf(outfile,"%0.2f\t",FCS[i]); /*Fiber Continuity score*/
		fprintf(outfile,"%lf\t",SensB[i]); /*Compliance*/
		fprintf(outfile,"%lf\n",SensA[i]); /*Sensitivity*/
	}
}

fclose(outfile);
printf("\nConvergence History file written\n");

printf("\n\n------*---End of BLES Version 3 Program---*------\n\n");
return (0);
}
