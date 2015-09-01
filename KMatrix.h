/*
 *  KMatrix.h
 *  
 *  Created by Peter Dunning on 26/02/2007.
 *
 */
#include "ls_types.h"
#include <stdio.h>
#include <math.h>

/*Function that assemebles the element stiffness matricies into the Index arrays used by the MA57 solver*/
void Assemble2(int begin, int nNod, int *tnodes, double *KE, int *irn, int *jrn, double *A);
