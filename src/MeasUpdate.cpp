

#include "../include/MeasUpdate.h"

void MeasUpdate(Matrix & x, Matrix z, double g, Matrix s, Matrix G, Matrix & P, double n, Matrix & K){

    int m = z.n_cols();
    Matrix Inv_W = zeros(m,m);

    for (int i=1; i<=m; i++){
        Inv_W(i,i) = s(i)*s(i);    // Inverse weight (measurement covariance)
    }

    // Kalman gain
    K = P*transpose(G)*inv(Inv_W+G*P*transpose(G));

    // State update
    x = x + K*(z-g);

    // Covariance update
    P = (eye(n)-transpose(K)*G)*P;
}