/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define DREAL
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias1.h"
#endif
#include "atlas_level1.h"
#include "cblas.h"

void cblas_drotmg(double *d1, double *d2, double *b1, const double b2,
                  double *P)
{
   Mjoin(PATL,rotmg)(d1, d2, b1, b2, P);
}

