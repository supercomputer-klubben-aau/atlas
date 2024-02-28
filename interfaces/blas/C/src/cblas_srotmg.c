/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define SREAL
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias1.h"
#endif
#include "atlas_level1.h"
#include "cblas.h"

void cblas_srotmg(float *d1, float *d2, float *b1, const float b2,
                  float *P)
{
   Mjoin(PATL,rotmg)(d1, d2, b1, b2, P);
}

