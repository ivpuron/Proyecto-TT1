//--------------------------------------------------------------------------
//
// IERS: Management of IERS time and polar motion data
//  
// Last modified:   2018/02/01   M. Mahooti
// 
//--------------------------------------------------------------------------
#include <math.h>
#include "../include/IERS.h"
#include "../include/SAT_const.h"

void IERS(Matrix eop, double  Mjd_UTC, double & x_pole, double & y_pole, double & UT1_UTC, double & LOD, double & dpsi, double & deps, double & dx_pole, double & dy_pole, double & TAI_UTC){    
        double mjd = (floor(Mjd_UTC));
        //Matrix aux = eop.fila(4);
        int n = eop.n_cols();
        int i=0;
        for(int j=1; j<=n;j++){
            if(mjd==eop(4,j)){
                i=j;
                break;
            }
        }

        Matrix aux = eop.columna(i);
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = aux(5)/Const::Arcs;  // Pole coordinate [rad]
        y_pole  = aux(6)/Const::Arcs;  // Pole coordinate [rad]
        UT1_UTC = aux(7);             // UT1-UTC time difference [s]
        LOD     = aux(8);             // Length of day [s]
        dpsi    = aux(9)/Const::Arcs;
        deps    = aux(10)/Const::Arcs;
        dx_pole = aux(11)/Const::Arcs; // Pole coordinate [rad]
        dy_pole = aux(12)/Const::Arcs; // Pole coordinate [rad]
        TAI_UTC = aux(13);            // TAI-UTC time difference [s]
}

void IERS(Matrix eop,double Mjd_UTC,char interp, double &x_pole, double &y_pole, double &UT1_UTC, double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC){
    
    if(interp =='l'){
        // linear interpolation
        double mjd = (floor(Mjd_UTC));
        int i = 0;
        for(int j=1;j<=eop.n_cols();j++){
            if(mjd==eop(4,j)){
                i=j;
                break;
            }
        }
        Matrix preeop = eop.columna(i);
        Matrix nexteop = eop.columna(i+1);
        double mfme = 1440.0*(Mjd_UTC-floor(Mjd_UTC));
        double fixf = mfme/1440.0;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole = preeop(5) + (nexteop(5) - preeop(5)) * fixf;
        y_pole = preeop(6) + (nexteop(6) - preeop(6)) * fixf;
        UT1_UTC = preeop(7) + (nexteop(7) - preeop(7)) * fixf;
        LOD = preeop(8) + (nexteop(8) - preeop(8)) * fixf;
        dpsi = preeop(9) + (nexteop(9) - preeop(9)) * fixf;
        deps = preeop(10) + (nexteop(10) - preeop(10)) * fixf;
        dx_pole = preeop(11) + (nexteop(11) - preeop(11)) * fixf;
        dy_pole = preeop(12) + (nexteop(12) - preeop(12)) * fixf;
        TAI_UTC = preeop(13);
        x_pole  = x_pole/Const::Arcs;  //% Pole coordinate [rad]
        y_pole  = y_pole/Const::Arcs;  //% Pole coordinate [rad]
        dpsi    = dpsi/Const::Arcs;
        deps    = deps/Const::Arcs;
        dx_pole = dx_pole/Const::Arcs; //% Pole coordinate [rad]
        dy_pole = dy_pole/Const::Arcs; //% Pole coordinate [rad]
    }

        else if (interp =='n'){
        double mjd = (floor(Mjd_UTC));
        Matrix aux = eop.fila(4);
        int n = eop.n_cols();
        int i=0;
        for(int j=1; j<=n;j++){
            if(mjd==aux(j)){
                i=j;
                break;
            }
        }

        aux = eop.columna(i);
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = aux(5)/Const::Arcs;  // Pole coordinate [rad]
        y_pole  = aux(6)/Const::Arcs;  // Pole coordinate [rad]
        UT1_UTC = aux(7);             // UT1-UTC time difference [s]
        LOD     = aux(8);             // Length of day [s]
        dpsi    = aux(9)/Const::Arcs;
        deps    = aux(10)/Const::Arcs;
        dx_pole = aux(11)/Const::Arcs; // Pole coordinate [rad]
        dy_pole = aux(12)/Const::Arcs; // Pole coordinate [rad]
        TAI_UTC = aux(13);            // TAI-UTC time difference [s]
    }
}
