

#include "..\include\Mjday.h"

double Mjday(int yr, int mon, int day, int hr, int min, int sec){

    double jd = 367.0 * (double)yr
            - floor( (7.0 * ((double)yr + floor( ((double)mon + 9.0) / 12.0) ) ) * 0.25 )
            + floor( 275.0 * (double)mon / 9.0 )
            + (double)day + 1721013.5
            + ( ((double)sec/60.0 + (double)min ) / 60.0 + (double)hr ) / 24.0;

    return jd-2400000.5;
}