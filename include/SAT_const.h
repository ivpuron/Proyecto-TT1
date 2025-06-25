#ifndef _CONST_
#define _CONST_

#include "matrix.h"
#include <math.h>

class Const{
public:
    // Mathematical constants
    static constexpr double pi = 3.14159265358979323846f;
    static constexpr double pi2 = 2.0*pi;
    static constexpr double Rad = pi/180.0;
    static constexpr double Deg = 180.0/pi;
    static constexpr double Arcs = 3600*180.0/pi;

    // General
    static constexpr double MJD_J2000 = 51544.5;
    static constexpr double T_B1950 = -0.500002108;
    static constexpr double c_light = 299792458.000000000;
    static constexpr double AU = 149597870700.000000;

    // Physical parameters of the Earth, Sun and Moon

    // Equatorial radius anf flattening
    static constexpr double R_Earth = 6378.1363e3;
    static constexpr double f_Earth = 1/298.257223563;
    static constexpr double R_Sun = 696000e3;
    static constexpr double R_Moon = 1738e3;

    //FALTAN
    
    // Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
    static constexpr double omega_Earth = 15.04106717866910/3600*Rad;   // [rad/s]; WGS-84

    // Gravitational coefficients
    static constexpr double GM_Earth    = 398600.435436e9;                  // [m^3/s^2]; DE430
    static constexpr double GM_Sun      = 132712440041.939400e9;            // [m^3/s^2]; DE430
    static constexpr double GM_Moon     = GM_Earth/81.30056907419062;       // [m^3/s^2]; DE430
    static constexpr double GM_Mercury  = 22031.780000e9;                   // [m^3/s^2]; DE430
    static constexpr double GM_Venus    = 324858.592000e9;                  // [m^3/s^2]; DE430
    static constexpr double GM_Mars     = 42828.375214e9;                   // [m^3/s^2]; DE430
    static constexpr double GM_Jupiter  = 126712764.800000e9;               // [m^3/s^2]; DE430
    static constexpr double GM_Saturn   = 37940585.200000e9;                // [m^3/s^2]; DE430
    static constexpr double GM_Uranus   = 5794548.600000e9;                 // [m^3/s^2]; DE430
    static constexpr double GM_Neptune  = 6836527.100580e9;                 // [m^3/s^2]; DE430
    static constexpr double GM_Pluto    = 977.0000000000009e9;              // [m^3/s^2]; DE430

    // Solar radiation pressure at 1 AU 
    static constexpr double P_Sol       = 1367/c_light; // [N/m^2] (~1367 W/m^2); IERS 96


};

#endif