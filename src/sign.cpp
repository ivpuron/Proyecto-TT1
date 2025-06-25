#include "../include/sign.h"

double sign(double a, double b){
    if (b>=0)
        return fabs(a);
    else
        return - fabs(a);
}
