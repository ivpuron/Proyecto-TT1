

#ifndef _ANGLESG_
#define _ANGLESG_

#include <math.h>
#include "matrix.h"
#include "globals.h"
#include "Geodetic.h"
#include "LTC.h"
#include "IERS.h"
#include "timediff.h"
#include "PoleMatrix.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "GHAMatrix.h"
#include "utils.h"
#include "gibbs.h"
#include "hgibbs.h"
#include"elements.h"


void anglesg (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3,
 Matrix Rs1, Matrix Rs2, Matrix Rs3, Matrix & r2, Matrix & v2);

#endif