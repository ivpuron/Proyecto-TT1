

#ifndef _UTILS_
#define _UTILS_

#include <math.h>

#include "matrix.h"
int rpoly(double *op, int degree, double *zeror, double *zeroi);
void quad(double a,double b1,double c,double *sr,double *si,
          double *lr,double *li);

void fxshfr(int l2, int *nz);

void quadit(double *uu,double *vv,int *nz);

void realit(double sss, int *nz, int *iflag);

void calcsc(int *type);

void nextk(int *type);

void newest(int type,double *uu,double *vv);

void quadsd(int n,double *u,double *v,double *p,double *q,
            double *a,double *b);

void roots(Matrix coef, Matrix & real, Matrix & im);
double custom_mod(double x, double y);

#endif