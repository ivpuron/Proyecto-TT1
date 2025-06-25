#ifndef _IERS_
#define _IERS_

#include <math.h>
#include "matrix.h"

/*Método para quedarse con la parte decimal de un número.
Por ejemplo, del número 2.5 devuelve 0.5.*/

void IERS(Matrix eop,double  Mjd_UTC,char interp, double &x_pole, double &y_pole, double &UT1_UTC, double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC);
void IERS(Matrix eop,double  Mjd_UTC, double &x_pole, double &y_pole, double &UT1_UTC, double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC);



#endif 