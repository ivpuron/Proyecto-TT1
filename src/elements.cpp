

#include "../include/elements.h"


 void elements (double & p, double & a, double & e, double & i, double & Omega, double & omega, double & M, Matrix y){
    Matrix r(3);
    r(1)=y(1);
    r(2)=y(2);
    r(3)=y(3);    
    Matrix v(4);                                       // Position
    v(1) = y(4);
    v(2)= y(5);
    v(3) = y(6);                                        // Velocity
    
    Matrix h = cross(r,v);                                    // Areal velocity
    double magh = h.norm();
    p = magh*magh/Const::GM_Earth;
    double H = h.norm();
    
    Omega = atan2 ( h(1), -h(2) );                     // Long. ascend. node 
    Omega = custom_mod(Omega,Const::pi2);
    i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) ); // Inclination        
    double u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );    // Arg. of latitude   
    
    double R  = r.norm();                                      // Distance           
    
    a = 1/(2/R-dot(v,v)/Const::GM_Earth);               // Semi-major axis    
    
    double eCosE = 1-R/a;                                     // e*cos(E)           
    double eSinE = dot(r,v)/sqrt(Const::GM_Earth*a);           // e*sin(E)           
    
    double e2 = eCosE*eCosE +eSinE*eSinE;
    e  = sqrt(e2);                                     // Eccentricity 
    double E  = atan2(eSinE,eCosE);                           // Eccentric anomaly  
    
    M  = custom_mod(E-eSinE,Const::pi2);                             // Mean anomaly
    
    double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          // True anomaly
    
    omega = custom_mod(u-nu,Const::pi2);                             // Arg. of perihelion 
 }