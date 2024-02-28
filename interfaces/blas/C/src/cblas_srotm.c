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

void cblas_srotm(const int N, float *X, const int incX,
                 float *Y, const int incY, const float *P)
{
   if (N > 0)
   {
      if (incX < 0)
      {
         if (incY < 0) ATL_srotm(N, X, -incX, Y, -incY, P);
         else ATL_srotm(N, X+(1-N)*incX, incX, Y, incY, P);
      }
      else if (incY < 0) ATL_srotm(N, X+(N-1)*incX, -incX, Y, -incY, P);
      else ATL_srotm(N, X, incX, Y, incY, P);
   }
}
