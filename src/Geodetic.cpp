

#include "../include/Geodetic.h"



void Geodetic(Matrix r, double & lon, double & lat, double & h){
    double R_equ = Const::R_Earth;
    double f = Const::f_Earth;

    double epsRequ = 2.22044604925031e-16*R_equ;        // Convergence criterion
    double e2      = f*(2.0-f);        // Square of eccentricity

    double X = r(1);                   // Cartesian coordinates
    double Y = r(2);
    double Z = r(3);
    double rho2 = X*X + Y*Y;           // Square of distance from z-axis

    // Check validity of input data
    if (r.norm()==0.0){
        cout << " invalid input in Geodetic constructor\n";
        lon = 0.0;
        lat = 0.0;
        h   = -R_equ;
    }

    // Iteration 
    double dZ = e2*Z;
    double ZdZ, Nh, SinPhi, N, dZ_new;

    while(true){
        ZdZ    =  Z + dZ;
        Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
        SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
        N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        dZ_new =  N*e2*SinPhi;
        if ( fabs(dZ-dZ_new) < epsRequ ){
            break;
        }
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;
}