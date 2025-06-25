#include "..\include\R_x.h"
#include "..\include\R_y.h"
#include "..\include\R_z.h"
#include "..\include\sign.h"
#include "..\include\Frac.h"
#include "..\include\Legendre.h"
#include "..\include\timediff.h"
#include "..\include\unit.h"
#include "..\include\globals.h"
#include "..\include\accelPointMass.h"
#include "..\include\AzElPa.h"
#include "..\include\Cheb3D.h"
#include "..\include\EccAnom.h"
#include "..\include\Geodetic.h"
#include "..\include\IERS.h"
#include "..\include\MeanObliquity.h"
#include "..\include\Mjday.h"
#include "..\include\Mjday_TDB.h"
#include "..\include\Position.h"
#include "..\include\NutAngles.h"
#include "..\include\EqnEquinox.h"
#include "..\include\NutMatrix.h"
#include "..\include\PoleMatrix.h"
#include "..\include\PrecMatrix.h"
#include "..\include\TimeUpdate.h"
#include "..\include\angl.h"
#include "..\include\gmst.h"
#include "..\include\LTC.h"
#include "..\include\MeasUpdate.h"
#include "..\include\elements.h"
#include "..\include\AccelHarmonic.h"
#include "..\include\gast.h"
#include "..\include\GHAMatrix.h"
#include "..\include\G_AccelHarmonic.h"
#include "..\include\JPL_Eph_DE430.h"
#include "..\include\Accel.h"
#include "..\include\gibbs.h"
#include "..\include\hgibbs.h"
#include "..\include\VarEqn.h"
#include "..\include\DEInteg.hpp"
#include "..\include\anglesg.h"
#include "..\include\Matrix.h"

#include <stdio.h>
#include <cmath>
#include <string>


int tests_run=0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__,__LINE__)
#define _assert(test) do {if(!(test)){FAIL(); return -1;}} while(0)
#define _verify(test) do {int r=test(); tests_run++; if(r) return r;} while(0)

int m_equals(Matrix A, Matrix B, int f, int c, double p){
    int equal=1;
    for(int i=1;i<=f;i++){
        for(int j=1;j<=c;j++){
            if(fabs(A(i,j)-B(i,j))>p){
                printf("%2.20lf %2.20lf\n", A(i,j),B(i,j));
                equal=0;
            }
        }
    }
    return equal;
}

int Frac1(){
    double a = -2.1;
    double res = Frac(a);
    _assert(res==-0.1);

    return 0;
}

int Frac2(){
    double a = 3.254;
    double res = Frac(a);
    _assert(res==0.254);

    return 0;
}

int R_x1(){

    Matrix rx = R_x(1);

    Matrix A(3,3);

    A(1,1) = 1.0;  A(1,2) = 0;        A(1,3) = 0;
    A(2,1) = 0;  A(2,2) = 0.8776;   A(2,3) = 0.4794;
    A(3,1) = 0;  A(3,2) = -0.4794;  A(3,3) = 0.8776;

    _assert(m_equals(A, rx, 3, 3, 1e-4));

    return 0;
}

int R_x2(){

    Matrix rx = R_x(0.48);

    Matrix A(3,3);

    A(1,1) = 1.0;  A(1,2) =  0;       A(1,3) = 0;
    A(2,1) = 0;  A(2,2) = 0.5403;   A(2,3) = 0.8415;
    A(3,1) = 0;  A(3,2) = -0.8415;  A(3,3) = 0.5403;

    _assert(m_equals(A, rx, 3, 3, 1e-4));

    return 0;
}

int R_y1(){

    Matrix ry = R_y(0.5);

    Matrix A(3,3);

    A(1,1) = 0.8776;  A(1,2) = 0;        A(1,3) = -0.4794;
    A(2,1) = 0;  A(2,2) = 0.8776;   A(2,3) = 0;
    A(3,1) = 0.4794;  A(3,2) = 0;  A(3,3) = 0.8776;

    _assert(m_equals(A, ry, 3, 3, 1e-4));

    return 0;
}

int R_y2(){

    Matrix ry = R_y(0.75);

    Matrix A(3,3);

    A(1,1) = 0.7317;  A(1,2) =  0;       A(1,3) = -0.6816;
    A(2,1) = 0;  A(2,2) = 1.0;   A(2,3) = 0;
    A(3,1) = 0.6816;  A(3,2) = 0;  A(3,3) = 0.7317;

    _assert(m_equals(A, ry, 3, 3, 1e-4));

    return 0;
}

int R_z1(){

    Matrix rz = R_z(0.65);

    Matrix A(3,3);

    A(1,1) = 0.7961;  A(1,2) = 0;        A(1,3) = 0.6052;
    A(2,1) = -0.6052;  A(2,2) = 0.7961;   A(2,3) = 0;
    A(3,1) = 0;  A(3,2) = 0;  A(3,3) = 1.0;

    _assert(m_equals(A, rz, 3, 3, 1e-4));

    return 0;
}

int R_z2(){

    Matrix rz = R_z(0.48);

    Matrix A(3,3);

    A(1,1) = 0.9737;  A(1,2) = 0;        A(1,3) = 0.2280;
    A(2,1) = -0.2280;  A(2,2) = 0.9737;   A(2,3) = 0;
    A(3,1) = 0;  A(3,2) = 0;  A(3,3) = 1.0;

    _assert(m_equals(A, rz, 3, 3, 1e-4));

    return 0;
}

int sign1(){
    double a = - 2.5;
    double b = 0.74;
    double res = sign(a,b);
    _assert(res==2.5);

    return 0;
}

int sign2(){
    double a = - 2.5;
    double b = - 0.23;
    double res = sign(a,b);
    _assert(res==-2.5);

    return 0;
}

int legendre1(){
    Matrix aux1(4,4);
    Matrix aux2(4,4);

    Matrix A(4,4);
    Matrix B(4,4);

    A(1,1) = 1;         A(1,2) = 0;         A(1,3) = 0;         A(1,4) = 0;
    A(2,1) = 1.0891;    A(2,2) = 1.3468;    A(2,3) = 0;         A(2,4) = 0;
    A(3,1) = 0.2081;    A(3,2) = 1.8936;    A(3,3) = 1.1708;    A(3,4) = 0;
    A(4,1) = -0.8510;    A(4,2) = 1.2307;    A(4,3) = 1.9478;    A(4,4) = 0.9834;

    B(1,1) = 0;         B(1,2) = 0;         B(1,3) = 0;         B(1,4) = 0;
    B(2,1) = 1.3468;    B(2,2) = -1.0891;   B(2,3) = 0;         B(2,4) = 0;
    B(3,1) = 3.2799;    B(3,2) = 0.8104;   B(3,3) = -1.8936;   B(3,4) = 0;
    B(4,1) = 3.0146;    B(4,2) = 5.1644;   B(4,3) = -0.7416;   B(4,4) = -2.3856;

    Legendre(3, 3, 0.68, aux1, aux2);
    
    _assert(m_equals(A, aux1, 4, 4, 1e-4));
    _assert(m_equals(B, aux2, 4, 4, 1e-4));

    return 0;
}


int legendre2(){
    Matrix aux1(3,4);
    Matrix aux2(3,4);

    Matrix A(3,4);
    Matrix B(3,4);

    A(1,1) = 1;         A(1,2) = 0;         A(1,3) = 0;         A(1,4)=0;
    A(2,1) = 0.5939;    A(2,2) = 1.6270;    A(2,3) = 0;         A(2,4)=0;
    A(3,1) = -0.7237;    A(3,2) = 1.2475;    A(3,3) = 0.7088;   A(3,4)=0;

    B(1,1) = 0;         B(1,2) = 0;         B(1,3) = 0;    B(1,4)=0;    
    B(2,1) = 1.6270;    B(2,2) = -0.5939;   B(2,3) = 0;        B(2,4)=0; 
    B(3,1) = 2.1608;    B(3,2) = 2.9622;   B(3,3) = -1.2475;   B(3,4)=0;

    Legendre(2, 3, 0.35, aux1, aux2);

    _assert(m_equals(A, aux1, 3, 2, 1e-4));
    _assert(m_equals(B, aux2, 3, 2, 1e-4));

    return 0;
}

int timediff1(){
    double UT1_UTC = 11.3;
    double TAI_UTC = 5.4;
    
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS; 
    double TT_UTC;
    double GPS_UTC;

    double UT1_TAI_m = 5.9;
    double UTC_GPS_m = 13.6;
    double UT1_GPS_m = 24.9; 
    double TT_UTC_m = 37.5840;
    double GPS_UTC_m = -13.6;

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double eps = 1e-4;

    _assert(fabs(UT1_TAI-UT1_TAI_m)<eps && fabs(UTC_GPS-UTC_GPS_m)<eps && fabs(UT1_GPS-UT1_GPS_m)<eps
            && fabs(TT_UTC-TT_UTC_m)<eps && fabs(GPS_UTC-GPS_UTC_m)<eps);

    return 0;
}

int unit1(){

    Matrix vec(3);
    Matrix outvec(3);
    Matrix result(3);

    vec(1) = 1;             vec(2) = 2;             vec(3) = 3;
    result(1) = 0.2673;     result(2) = 0.5345;     result(3) = 0.8018;
    outvec = unit(vec);

    _assert(m_equals(result, outvec, 1, 3, 1e-4));

    return 0;
}

int unit2(){
    Matrix vec(3);
    Matrix outvec(3);
    Matrix result(3);
    
    vec(1) = 1.3;           vec(2) = -0.64;          vec(3) = 6.87;
    result(1) = 0.1852;     result(2) = -0.0912;     result(3) = 0.9785;

    outvec = unit(vec);
    
    _assert(m_equals(result, outvec, 1, 3, 1e-4));

    return 0;
}

int accelPointMass1(){
    Matrix result(3);
    Matrix r(3);        r(1) = 1;       r(2) = 2;       r(3) = 3;
    Matrix s(3);        s(1) = 4;       s(2) = 5;       s(3) = 6;

    Matrix compare(3);       compare(1) = 0.0077;    compare(2) = 0.0070; compare(3) = 0.0063;

    result = AccelPointMass(r, s, 0.5);
    _assert(m_equals(compare, result, 1, 3, 1e-4));

    return 0;
}

int accelPointMass2(){
    Matrix result(3);
    Matrix r(3);        r(1) = 2.3;       r(2) = 4.5;       r(3) = -1.6;
    Matrix s(3);        s(1) = 8.54;      s(2) = -6.7;        s(3) = 3.2;

    Matrix compare(3);              compare(1) = -0.0022;
    compare(2) = -0.0002; compare(3) = -0.0002;

    result = AccelPointMass(r, s, 0.64);
    _assert(m_equals(compare, result, 1, 3, 1e-4));

    return 0;
}

int AzElPa1(){
    Matrix dAds(3);
    Matrix dEds(3);  
    double Az;  double El;    
    Matrix s(3);        s(1) = 1;      s(2) = 2;        s(3) = 3;

    Matrix dAds_comp(3);  dAds_comp(1) = 0.4;   dAds_comp(2) = -0.2;   dAds_comp(3) = 0.0;
    Matrix dEds_comp(3);  dEds_comp(1) = -0.0958;   dEds_comp(2) = -0.1917;   dEds_comp(3) = 0.1597;
    double Az_comp = 0.4636;  double El_comp = 0.9303;
    AzElPa(s, Az, El, dAds, dEds);
    
    _assert(m_equals(dAds_comp , dAds , 1, 3, 1e-4));
    _assert(m_equals(dEds_comp , dEds , 1, 3, 1e-4));
    _assert(fabs(Az-Az_comp)<1e-5);
    _assert(fabs(El-El_comp)<1e-5);

    return 0;
}


int AzElPa2(){
    Matrix dAds(3);
    Matrix dEds(3);  
    double Az;  double El;    
    Matrix s(3);        s(1) = 5.34;      s(2) = -9.65;        s(3) = 3.21;

    Matrix dAds_comp(3);  dAds_comp(1) = -0.0793;   dAds_comp(2) = -0.0439;    dAds_comp(3) = 0.0;
    Matrix dEds_comp(3);  dEds_comp(1) = -0.0118;    dEds_comp(2) = 0.0213;    dEds_comp(3) = 0.0836;
    double Az_comp = 2.6362;  double El_comp = 0.2832;

    AzElPa(s, Az, El, dAds, dEds);

    _assert(m_equals(dAds_comp , dAds , 1, 3, 1e-4));
    _assert(m_equals(dEds_comp , dEds , 1, 3, 1e-4));
    _assert(fabs(Az-Az_comp)<1e-5);
    _assert(fabs(El-El_comp)<1e-5);

    return 0;
}

int Cheb3D1(){

    double t = 2.5; 
    int N = 3; 
    double Ta = 1.3;
    double Tb = 4.5;
    Matrix Cx(3); Matrix Cy(3); Matrix Cz(3);
    Matrix result(3);
    Matrix compare(3);  
    
    Cx(1) = 1;      Cx(2) = 2;        Cx(3) = 3;
    Cy(1) = 4;      Cy(2) = 5;        Cy(3) = 6;
    Cz(1) = 7;      Cz(2) = 8;        Cz(3) = 9;

    compare(1) = -2.1250;   compare(2) = -2.5;    compare(3) = -2.8750;

    result = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    _assert(m_equals(compare , result , 1, 3, 1e-4));

    return 0;
}

int EccAnom1(){
    double M = 5.0; 
    double e = 0.7; 
    double compare = 4.3464;
    double result = EccAnom(M, e);
    _assert(fabs(compare - result)<1e-4);

    return 0;
}

int EccAnom2(){
    double M = -2.8; 
    double e = 1.0; 
    double compare = 3.3128;
    double result = EccAnom(M, e);
    _assert(fabs(compare - result)<1e-4);

    return 0;
}




int all_test(){
    _verify(Frac1);
    _verify(Frac2);
    _verify(R_x1);
    _verify(R_x2);
    _verify(R_y1);
    _verify(R_y2);
    _verify(R_z1);
    _verify(R_z2);
    _verify(sign1);
    _verify(sign2);
    _verify(legendre1);
    _verify(legendre2);
    _verify(timediff1);
    _verify(unit1);
    _verify(unit2);
    _verify(accelPointMass1);
    _verify(accelPointMass2);
    _verify(AzElPa1);
    _verify(AzElPa2);
    _verify(Cheb3D1);
    _verify(EccAnom1);
    _verify(EccAnom2);

    return 0;
}

int main(){
    return all_test();
}
