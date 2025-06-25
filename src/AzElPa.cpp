#include "../include/AzElPa.h"


void AzElPa(Matrix s,double & Az, double & El, Matrix &dAds, Matrix &dEds){


    double rho = sqrt(s(1)*s(1)+s(2)*s(2));

    // Angles
    Az = atan2(s(1,1),s(1,2));

    if (Az<0.0) 
        Az = Az+Const::pi2;

    El = atan ( s(3) / rho );

    // Partials
    dAds(1) =  s(2)/(rho*rho);
    dAds(2) = -s(1)/(rho*rho);
    dAds(3) =  0.0;
    dEds(1) =  -s(1)*s(3)/rho;
    dEds(2) = -s(2)*s(3)/rho;
    
    dEds(3) = rho;
    dEds=dEds/dot(s,s); 
}