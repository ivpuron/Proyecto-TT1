#ifndef PROYECTO_GLOBALS_H
#define PROYECTO_GLOBALS_H

#include "matrix.h"
#include "Mjday.h"
#include "SAT_Const.h"
#include <cstdio>
#include <cstdlib>
#include <string.h>

typedef struct{
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

class Global{
    public:
     static Matrix *eopdata;
     static Param auxparam;
     static Matrix *Cnm;
     static Matrix *Snm;
     static Matrix *PC;
     static Matrix *obs;


     static void eop19620101(int c);
     static void GGM03S();
     static void DE430Coeff(int f, int c);
     static void GEOS3(int f);
};

#endif //PROYECTO_GLOBALS_H