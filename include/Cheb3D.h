#ifndef _CHEB3D_
#define _CHEB3D_


#include <math.h>
#include "matrix.h"

Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz);

#endif