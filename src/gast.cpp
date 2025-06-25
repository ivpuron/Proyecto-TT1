


#include "../include/gast.h"

double gast (double Mjd_UT1){

    return custom_mod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), Const::pi2 );
}