

#ifndef _JPL_
#define _JPL_


#include <math.h>
#include "globals.h"
#include "Cheb3D.h"



void JPL_Eph_DE430(Matrix & r_Mercury,Matrix & r_Venus,Matrix & r_Earth,Matrix & r_Mars,Matrix & r_Jupiter,Matrix & r_Saturn,Matrix & r_Uranus, 
          Matrix & r_Neptune,Matrix & r_Pluto,Matrix & r_Moon,Matrix & r_Sun, double Mjd_TDB);


#endif