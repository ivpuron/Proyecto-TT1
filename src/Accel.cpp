


#include "../include/Accel.h"

Matrix Accel(double x, Matrix & Y){

    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(*Global::eopdata, Global::auxparam.Mjd_UTC+ x/86400,'l',x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC );
    double UT1_TAI,  UTC_GPS,  UT1_GPS,  TT_UTC,  GPS_UTC;
    timediff( UT1_UTC,  TAI_UTC,  UT1_TAI,  UTC_GPS,  UT1_GPS,  TT_UTC,  GPS_UTC);
    double Mjd_UT1 =  Global::auxparam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    double Mjd_TT =  Global::auxparam.Mjd_UTC + x/86400 + TT_UTC/86400;

    Matrix P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    Matrix r_Mercury(3),r_Venus (3),r_Earth (3),r_Mars (3),r_Jupiter (3),r_Saturn (3),r_Uranus (3), 
    r_Neptune (3),r_Pluto (3),r_Moon (3),r_Sun (3);
     JPL_Eph_DE430(r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, 
    r_Neptune,r_Pluto,r_Moon,r_Sun, MJD_TDB);

    // Acceleration due to harmonic gravity field
    Matrix aux(3);
    aux(1)=Y(1); aux(2)=Y(2); aux(3)=Y(3);
    Matrix a = AccelHarmonic(aux, E,  Global::auxparam.n,  Global::auxparam.m);

    // Luni-solar perturbations
    if ( Global::auxparam.sun){
        a = a + AccelPointMass(aux,r_Sun,Const::GM_Sun);
    }

    if ( Global::auxparam.moon){
        a = a + AccelPointMass(aux,r_Moon,Const::GM_Moon);
    }

    // Planetary perturbations
    if ( Global::auxparam.planets){
        a = a + AccelPointMass(aux,r_Mercury,Const::GM_Mercury);
        a = a + AccelPointMass(aux,r_Venus,Const::GM_Venus);
        a = a + AccelPointMass(aux,r_Mars,Const::GM_Mars);
        a = a + AccelPointMass(aux,r_Jupiter,Const::GM_Jupiter);
        a = a + AccelPointMass(aux,r_Saturn,Const::GM_Saturn);
        a = a + AccelPointMass(aux,r_Uranus,Const::GM_Uranus);    
        a = a + AccelPointMass(aux,r_Neptune,Const::GM_Neptune);
        a = a + AccelPointMass(aux,r_Pluto,Const::GM_Pluto);
    }
    Matrix dY(6);
    dY(1)=Y(4); dY(2)=Y(5); dY(3)=Y(6);
    dY(4)=a(1); dY(5) = a(2); dY(6)=a(3);
    return dY;
}