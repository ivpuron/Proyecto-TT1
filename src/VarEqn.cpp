


#include "../include/VarEqn.h"


Matrix VarEqn(double x, Matrix & yPhi){

    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(*Global::eopdata, Global::auxparam.Mjd_UTC,'l',x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC );
    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    double Mjd_UT1 = Global::auxparam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    // Transformation matrix
    Matrix P = PrecMatrix(Const::MJD_J2000,Global::auxparam.Mjd_TT + x/86400);
    Matrix N = NutMatrix(Global::auxparam.Mjd_TT + x/86400.0);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    Matrix r(3);
    r(1)=yPhi(1);
    r(2)=yPhi(2);
    r(3)=yPhi(3);    
    Matrix v(4);                                       // Position
    v(1) = yPhi(4);
    v(2)= yPhi(5);
    v(3) = yPhi(6);   
    Matrix Phi = zeros(6,6);

    // State transition matrix
    for(int j=1 ; j<=6;j++){
        for(int i=1; i<=6;i++){
            Phi(i,j) = yPhi(6*j+i);
        }
    }

    // Acceleration and gradient
    Matrix a = AccelHarmonic ( r, E, Global::auxparam.n, Global::auxparam.m );
    Matrix G = G_AccelHarmonic ( r, E, Global::auxparam.n, Global::auxparam.m );

    // Time derivative of state transition matrix
    Matrix yPhip = zeros(42,1);
    Matrix dfdy = zeros(6,6);

    for (int i=1; i<=3;i++){
        for (int j=1; j<=3; j++){
            dfdy(i,j) = 0.0;                 // dv/dr(i,j)
            dfdy(i+3,j) = G(i,j);            // da/dr(i,j)
            if ( i==j ){
                dfdy(i,j+3) = 1;
            }
            else{
                dfdy(i,j+3) = 0;             // dv/dv(i,j)
            }
            dfdy(i+3,j+3) = 0.0;             // da/dv(i,j)
        }
    }

    Matrix Phip (6,6);
    Phip = dfdy*Phi;

    // Derivative of combined state vector and state transition matrix
    for (int i=1;i<=3;i++){
        yPhip(i)   = v(i);                 // dr/dt(i)
        yPhip(i+3) = a(i);                 // dv/dt(i)
    }

    for (int i=1; i<=6; i++){
        for (int j=1; j<=6;j++){
            yPhip(6*j+i) = Phip(i,j);     // dPhi/dt(i,j)
        }
    }
    return yPhip;
}