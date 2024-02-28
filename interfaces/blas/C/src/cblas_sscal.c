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

void cblas_sscal(const int N, const float alpha, float *X, const int incX)
{
   if (N > 0 && incX > 0)
      ATL_sscal(N, alpha, X, incX);
}
