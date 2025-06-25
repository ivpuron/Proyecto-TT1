#include <iostream>
#include <string>
#include <stdio.h>

#include "../include/globals.h"
#include "../include/Matrix.h"
#include "../include/SAT_const.h"
#include "../include/Position.h"
#include "../include/Accel.h"
#include "../include/DEInteg.hpp"
#include "../include/LTC.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/VarEqn.h"
#include "../include/gmst.h"
#include "../include/TimeUpdate.h"
#include "../include/AzElPa.h"
#include "../include/MeasUpdate.h"
#include "../include/anglesg.h"

using namespace std;

int main(){
    Global::eop19620101(6);
    Global::eop19620101(21413);
    Global::DE430Coeff(2285, 1020);
    Global::GEOS3(46);
    Global::GGM03S();

    int nobs=46;

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Const::Rad; // [rad]
    double sigma_el = 0.0139*Const::Rad; // [rad]

    // Kaena Point station
    double lat = Const::Rad*21.5748;     // [rad]
    double lon = Const::Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix Rs = Position(lon, lat, alt);

    double Mjd1 = (*Global::obs)(1,1);
    double Mjd2 = (*Global::obs)(9,1);
    double Mjd3 = (*Global::obs)(18,1);

    Matrix r2(3), v2(3);

    anglesg((*Global::obs)(1,2),(*Global::obs)(9,2),
            (*Global::obs)(18,2),(*Global::obs)(1,3),
            (*Global::obs)(9,3),(*Global::obs)(18,3),
            Mjd1,Mjd2,Mjd3,Rs,Rs,Rs,r2,v2);
    // [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
    //                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

    Matrix Y0_apr(6);
    for(int i = 1; i <= 3; i++){
        Y0_apr(i) = r2(i);
        Y0_apr(i+3) = v2(i);
    }

    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = (*Global::obs)(9,1);

    Global::auxparam.Mjd_UTC = Mjd_UTC;
    Global::auxparam.n      = 20;
    Global::auxparam.m      = 20;
    Global::auxparam.sun     = 1;
    Global::auxparam.moon    = 1;
    Global::auxparam.planets = 1;

    int n_eqn  = 6;

    Matrix Y = DEInteg(Accel, n_eqn, Y0_apr, 0, -((*Global::obs)(9,1)-Mjd0)*86400.0, 1e-13, 1e-6);

    Matrix P = zeros(6,6);
    
    for(int i=1; i<=3;i++){
        P(i,i)=1e8;
    }
    for(int i=4;i<=6;i++){
        P(i,i)=1e3;
    }

    Matrix LT = LTC(lon,lat);

    Matrix yPhi = zeros(42,1);
    Matrix Phi  = zeros(6,6);

    for(int i = 1; i <= 3; i++){
        P(i,i)=1e8;
        P(i+3,i+3)=1e3;
    }
    


    // Measurement loop
    double t = 0, t_old;
    Matrix Y_old(6);
    Matrix eop = *Global::eopdata;
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;

    double Mjd_TT, Mjd_UT1, Azim, theta, Elev, Dist;

    Matrix r(3), U(3), s(3), dAds(3), dEds(3), aux1(3), aux2(3), aux3(3), aux4(3), dAdY(6),
    dEdY(6), K(6), dDds(3), dDdY(6);

    for (int i=1;i<=nobs;i++){    
        // Previous step
        t_old = t;
        Y_old = Y;
        
        // Time increment and propagation
        Mjd_UTC = (*Global::obs)(i,1);                       // Modified Julian Date
        t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
        
        IERS(*Global::eopdata, Global::auxparam.Mjd_UTC,'l',x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC );
        timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        Global::auxparam.Mjd_UTC = Mjd_UTC;
        Global::auxparam.Mjd_TT = Mjd_TT;
            
        for(int ii=1;ii<=6;ii++){
            yPhi(ii) = Y_old(ii);
            for (int j=1;j<=6;j++){  
                if (ii==j){ 
                    yPhi(6*j+ii) = 1; 
                }
                else{
                    yPhi(6*j+ii) = 0;
                }
            }
        }
        
        yPhi = DEInteg (VarEqn, 42, yPhi, 0, t-t_old, 1e-13, 1e-6);
        
        // Extract state transition matrices
        for (int j=1;j<=6;j++){
            for(int k = 1; k <=6; k++){
                Phi(k,j) = yPhi(6*j+k);
            }
        }
        
        Y = DEInteg (Accel, 6, Y_old, 0, t-t_old, 1e-13, 1e-6);
        
        // Topocentric coordinates
        theta = gmst(Mjd_UT1);                    // Earth rotation
        U = R_z(theta);
        r(1)=Y(1); r(2)=Y(2); r(3)=Y(3);
        s = transpose(LT*(U*transpose(r)-transpose(Rs)));                          // Topocentric position [m]
        
        // Time update
        TimeUpdate(P, Phi);
            
        // Azimuth and partials
        AzElPa(s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation
        dAdY = zeros(1,6);
        Matrix aux=dAds*LT*U;
        dAdY(1)=aux(1); dAdY(2)=aux(2); dAdY(3)=aux(3);
        
        // Measurement update
        MeasUpdate(Y, (*Global::obs)(i,2), Azim,  dAdY, sigma_az,P, 6,K);
        
        // Elevation and partials
        r(1)=Y(1); r(2)=Y(2); r(3)=Y(3);
        s = transpose(LT*(U*transpose(r)-transpose(Rs)));                          // Topocentric position [m]
        AzElPa(s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation
        dEdY = zeros(1,6);
        aux=dEds*LT*U;
        for (int k = 1; k <= 3; k++){
            dEdY(k) = aux(k);
        }
        
        // Measurement update
        MeasUpdate ( Y, (*Global::obs)(i,3), Elev, dEdY,sigma_el,  P, 6,K );

        // Range and partials
        r(1)=Y(1); r(2)=Y(2); r(3)=Y(3);
        s = transpose(LT*(U*transpose(r)-transpose(Rs)));                         // Topocentric position [m]
        Dist = s.norm(); dDds = (s/Dist);         // Range
        dDdY = zeros(1,6);
        aux = dDdY*LT*U;
        for (int k = 1; k <= 3; k++){
            dDdY(k) = aux(k);
        }
        
        
        // Measurement update
        MeasUpdate ( Y, (*Global::obs)(i,4), Dist, dDdY,sigma_range,  P, 6,K );
    }

    IERS(eop,(*Global::obs)(46,1),'l',x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
        timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    Global::auxparam.Mjd_UTC = Mjd_UTC;
    Global::auxparam.Mjd_TT = Mjd_TT;

    Matrix Y0 = DEInteg (Accel, n_eqn, Y, 0, -((*Global::obs)(46,1)-(*Global::obs)(1,1))*86400.0, 1e-13, 1e-6);

    Matrix Y_true(6);Y(1)=5753.173e3; Y(2)=2673.361e3; Y(3)=3440.304e3, Y(4)=4.324207e3, Y(5)=-1.924299e3, Y(6)=-5.728216e3;

    printf("\nError of Position Estimation\n");
    printf("dX      %10.1f [m]\n",Y0(1)-Y_true(1));
    printf("dY      %10.1f [m]\n",Y0(2)-Y_true(2));
    printf("dZ      %10.1f [m]\n",Y0(3)-Y_true(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx     %8.1f [m/s]\n",Y0(4)-Y_true(4));
    printf("dVy     %8.1f [m/s]\n",Y0(5)-Y_true(5));
    printf("dVz     %8.1f [m/s]\n",Y0(6)-Y_true(6));

    return 0;
}

