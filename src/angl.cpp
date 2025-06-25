


#include "../include/angl.h"


double angl ( Matrix vec1,Matrix vec2 )
{
    double small     = 0.00000001;
    double undefined = 999999.1;

    double magv1 = vec1.norm();
    double magv2 = vec2.norm();

    double theta;

    if (magv1*magv2 > small*small){
        double temp= dot(vec1,vec2) / (magv1*magv2);
        if (fabs( temp ) > 1.0){
            temp= sign(temp,temp) * 1.0;
        }
        theta= acos( temp );
    }
    else{
        theta= undefined;
    }
    return theta;
}