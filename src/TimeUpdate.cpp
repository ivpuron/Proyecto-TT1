


#include "../include/TimeUpdate.h"

void TimeUpdate(Matrix & P, Matrix Phi, double Qdt){

    P = Phi*P*transpose(Phi) + Qdt;
}

void TimeUpdate(Matrix & P, Matrix Phi){

    P = Phi*P*transpose(Phi);
}