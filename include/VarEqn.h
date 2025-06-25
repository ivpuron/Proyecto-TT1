


#ifndef _VAREQN_
#define _VAREQN_

#include "matrix.h"
#include "SAT_const.h"
#include "IERS.h"
#include "globals.h"
#include "timediff.h"
#include "PoleMatrix.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h" 
#include "G_AccelHarmonic.h" 

Matrix VarEqn(double x, Matrix & yPhi);

#endif