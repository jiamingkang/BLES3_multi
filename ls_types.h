/*
 *  ls_types.h
 *  
 *
 *  Created by Peter Dunning on 06/02/2009.
 *
 *	Based on fg_types.h
 *	Version beta 1.0 @ 06/02/2009
 */

#ifndef __LS_TYPES_H
#define __LS_TYPES_H

/*Element and node Numbering information structure*/
typedef struct 
				{
				int n;	/* element number */
				int a,b,c,d; /* node numbers */
				} Elem;

/*Co-ordinate inforamtion structure*/
typedef struct
				{
				double x,y,z,xz; /*Co-ordinates*/
				} Coord;
				
/*Structure for inital circular hole data*/
typedef struct
				{
				double x,y; /*center co-ordinates*/
				double r;	/*hole radius*/
				} CirH;

/*Structure for line selection set data*/
typedef struct
				{
				double x,y; /*center co-ordinates*/
				int r;	/*hole radius*/
				} SET;

/*Structure to hold four doubleing point data, i.e. node or gauss point data per element*/				
typedef struct 
				{
				double a,b,c,d; /* anti-clockwise numbers */
				} FourF;
				
/*Structure to hold Pentagonal shape function co-efficient data*/				
typedef struct 
				{
				double a,b,c,d; /* Co-efficients*/
				} Shape;
				
/*Structure for position of nodes */				
typedef struct {
				short x,y;
				} Pos;

/*Structure for strain energy at a co-ordinate*/
typedef struct
				{
				double x,y; /*Co-ordinates*/
				double e;	/*Strain Energy*/
				} Astn;

/*Structure for auxillary node data*/
typedef struct
				{
				double x,y,z,xz;
				int n;
				} Aux;

/*Structure for auxillary node data*/
typedef struct
                {
				double x,y,z,xz;
				int n;
                } Auxreal;
// Structure to hold data for a boundary segment
typedef struct
                {
                int n1, n2; // node numbers for end points
                int e;  // associated element
                }	Bseg;

/*Structure for strain energy at a co-ordinate - with element area data*/
typedef struct
				{
				double x,y;	/*Co-ordinates*/
				double u;	/*Strain Energy*/
				double a;	/*assosiated element area*/
				} Gstn;
typedef struct
                {
				int an,bn,cn,dn;	/*The Neighbouring element numbers*/
				int ai,bi,ci,di;	/*Idex stating weather or not each element exists*/
				double ad,bd,cd,dd;			/*assosiated element differnetial of theta wrt to the node lsf value*/
                } NElemDiff;

/*Structure to store the diffential of theta wrt lsf in the neighbouring elements to the node*/

typedef struct
                {
				int an, bn, cn, dn;	/*The Neighbouring element numbers*/
				int ai, bi, ci, di;	/*Idex stating weather or not each element exists*/
				double ad, bd, cd, dd;			/*assosiated element differnetial of theta wrt to the node lsf value*/
				} NElemDiffFomer;


/*structure to store 3 short interger numbers*/
typedef struct {
				short a,b,c;
				} Three;

/*structure to store 4 short interger numbers*/
typedef struct {
				short a,b,c,d;
				} Four;
				
#endif
