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

void catlas_sset(const int N, const float alpha, float *X, const int incX)
{
   ATL_sset( N, alpha, X, (incX >= 0 ? incX : -incX) );
}
