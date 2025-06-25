

#ifndef _GIBBS_
#define _GIBBS_

#include <stdio.h>
#include "matrix.h"
#include "unit.h"
#include "angl.h"
#include "SAT_const.h"

void gibbs(Matrix r1,Matrix r2,Matrix r3, Matrix & v2, double & theta, double & theta1, double & cop, string & error);

#endif