/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

/*
 * Prototypes ATLAS Level 1 functions not defined in atlas_aux.h
 */
#ifndef ATLAS_LEVEL1_H
#define ATLAS_LEVEL1_H

/*
 * Many level one blas routines actually taken care of by atlas auxiliary
 */
#include "atlas_aux.h"

float ATL_sdsdot(const int N, const float alpha, const float *X,
                 const int incX, const float *Y, const int incY);
double ATL_dsdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY);
/*
 * Routines with all four types
 */
void ATL_sswap(const int N, float *X, const int incX,
               float *Y, const int incY);
int ATL_isamax(const int N, const float *X, const int incX);

void ATL_dswap(const int N, double *X, const int incX,
               double *Y, const int incY);
int ATL_idamax(const int N, const double *X, const int incX);

void ATL_cswap(const int N, float *X, const int incX,
               float *Y, const int incY);
int ATL_icamax(const int N, const float *X, const int incX);

void ATL_zswap(const int N, double *X, const int incX,
               double *Y, const int incY);
int ATL_izamax(const int N, const double *X, const int incX);

/*
 * Routines with real types
 */
void ATL_srotg(float *a, float *b, float *c, float *s);
void ATL_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void ATL_srot(const int N, float *X, const int incX,
              float *Y, const int incY, const float c, const float s);
void ATL_srotm(const int N, float *X, const int incX,
               float *Y, const int incY, const float *P);
float ATL_sdot(const int N, const float *X, const int incX,
                     const float *Y, const int incY);
void ATL_sssq(const int N, const float *X, const int incX,
              float *scal0, float *ssq0);
float ATL_snrm2(const int N, const float *X, const int incX);
float ATL_sasum(const int N, const float *X, const int incX);

void ATL_drotg(double *a, double *b, double *c, double *s);
void ATL_drotmg(double *d1, double *d2, double *b1, const double b2, double *P);
void ATL_drot(const int N, double *X, const int incX,
              double *Y, const int incY, const double c, const double s);
void ATL_drotm(const int N, double *X, const int incX,
               double *Y, const int incY, const double *P);
double ATL_ddot(const int N, const double *X, const int incX,
                     const double *Y, const int incY);
void ATL_dssq(const int N, const double *X, const int incX,
              double *scal0, double *ssq0);
double ATL_dnrm2(const int N, const double *X, const int incX);
double ATL_dasum(const int N, const double *X, const int incX);

/*
 * Routines with complex types
 */
void ATL_csrot(const int N, float *X, const int incX,
               float *Y, const int incY, const float c, const float s);
void ATL_crotg(float *a, const float *b, float *c, float *s);
void ATL_cdotu_sub(const int N, const float *X, const int incX,
                   const float *Y, const int incY, float *dot);
void ATL_cdotc_sub(const int N, const float *X, const int incX,
                   const float *Y, const int incY, float *dot);
void ATL_cssq(const int N, const float *X, const int incX,
              float *scal0, float *ssq0);
float ATL_scnrm2(const int N, const float *X, const int incX);
float ATL_scasum(const int N, const float *X, const int incX);

void ATL_zdrot(const int N, double *X, const int incX,
               double *Y, const int incY, const double c, const double s);
void ATL_zrotg(double *a, const double *b, double *c, double *s);
void ATL_zdotu_sub(const int N, const double *X, const int incX,
                   const double *Y, const int incY, double *dot);
void ATL_zdotc_sub(const int N, const double *X, const int incX,
                   const double *Y, const int incY, double *dot);
void ATL_zssq(const int N, const double *X, const int incX,
              double *scal0, double *ssq0);
double ATL_dznrm2(const int N, const double *X, const int incX);
double ATL_dzasum(const int N, const double *X, const int incX);


#define ATL_casum ATL_scasum
#define ATL_zasum ATL_dzasum
#define ATL_cnrm2 ATL_scnrm2
#define ATL_znrm2 ATL_dznrm2

#endif
