

#ifndef _HGIBBS_
#define _HGIBBS_

#include "matrix.h"
#include <stdio.h>
#include "SAT_const.h"
#include "unit.h"
#include "angl.h"


void hgibbs(Matrix r1, Matrix r2, Matrix r3, double Mjd1, double Mjd2, double Mjd3, Matrix & v2,
 double & theta, double &theta1, double & cop, string & error);

 #endif