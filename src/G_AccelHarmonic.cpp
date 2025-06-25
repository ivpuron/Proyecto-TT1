


#include "../include/G_AccelHarmonic.h"

Matrix G_AccelHarmonic( Matrix r, Matrix U, int n_max, int m_max ){

    double d = 1.0;   // Position increment [m]

    Matrix G = zeros(3,3);
    Matrix dr = zeros(1,3);

    Matrix da (3);

    // Gradient
    for (int i=1; i<=3 ; i++){
        // Set offset in i-th component of the position vector
        dr = zeros(1,3);
        dr(i) = d;
        // Acceleration difference
        da = AccelHarmonic ( r+dr/2,U, n_max, m_max ) - 
            AccelHarmonic ( r-dr/2,U, n_max, m_max );
        // Derivative with respect to i-th axis
        G(1,i) = da(1)/d;   
        G(2,i) = da(2)/d; 
        G(3,i) = da(3)/d;  
    }
    return G;
}