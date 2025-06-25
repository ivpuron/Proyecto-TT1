#include "../include/EccAnom.h"


double EccAnom (double M, double e)
{
    int maxit = 15;
    int i = 1;

    // Starting value
    M = fmod(M, Const::pi2);
    double E;
    if (e<0.8){
        E = M; 
    }
    else{
        E = Const::pi2;
    }

    double f = E - e*sin(E) - M;
    E = E - f / ( 1.0 - e*cos(E) );

    // Iteration
    while (fabs(f) > 1e2*pow(2,-52)){    
        f = E - e*sin(E) - M;
        E = E - f / ( 1.0 - e*cos(E) );
        i = i+1;
        if (i==maxit){
            cout << "convergence problems in EccAnom";
            exit(EXIT_FAILURE);
        }  
    }
    return E;
}