#include "../include/AccelPointMass.h"
/*%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%				  mass
%
% Inputs:
%   r           Satellite position vector 
%   s           Point mass position vector
%   GM          Gravitational coefficient of point mass
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------*/
Matrix AccelPointMass(Matrix r, Matrix s, double GM){

// Relative position vector of satellite w.r.t. point mass 
Matrix d = r - s;

// Acceleration 
return (d/(d.norm()*d.norm()*d.norm()) + s/(s.norm()*s.norm()*s.norm()) )*(-GM);
}