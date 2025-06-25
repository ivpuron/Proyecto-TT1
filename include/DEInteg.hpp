
// $Header$
//----------------------------------------------------------------------------------------
//                          DEInteg
//----------------------------------------------------------------------------------------
// Under the MIT License 2020
//
// Created: 2024/04/26
//
/*
 * @file DEInteg.hpp
 * @brief Header file for the ODE solver functions
 *
 * @details This header file contains the declarations for solving ordinary differential equations (ODEs) stepsize multistep method of Shampine & Gordon.
 *
 * @author Lorena Remacha Bordallo
*/

# include "matrix.h"


void de ( Matrix (*f) (double, Matrix &), int neqn, double y[],
          double &t, double tout, double relerr, double abserr, int &iflag, double yy[],
          double wt[], double p[], double yp[], double ypout[], double phi[],
          double alpha[], double beta[], double sig[], double v[], double w[],
          double g[], bool &phase1, double psi[], double &x, double &h, double &hold,
          bool &start, double &told, double &delsgn, int &ns, bool &nornd, int &k, int &kold,
          int &isnold );

int i4_sign ( int i );


void intrp ( double x, double y[], double xout, double yout[], double ypout[],
             int neqn, int kold, double phi[], double psi[] );

void ode ( Matrix (*f) (double, Matrix &), int neqn, double y[],
           double &t, double tout, double relerr, double abserr, int &iflag,
           double work[], int iwork[] );


double r8_sign ( double x );

void step ( double &x, double y[], Matrix (*f) (double, Matrix &),
            int neqn, double &h, double &eps, double wt[], bool &start, double &hold,
            int &k, int &kold, bool &crash, double phi[], double p[], double yp[],
            double psi[], double alpha[], double beta[], double sig[], double v[],
            double w[], double g[], bool &phase1, int &ns, bool &nornd );


void timestamp ( );


Matrix DEInteg ( Matrix (*f) (double, Matrix &), int neqn, Matrix y,
           double t, double tout, double relerr, double abserr);
