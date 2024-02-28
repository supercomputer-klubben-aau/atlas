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

void cblas_drotm(const int N, double *X, const int incX,
                 double *Y, const int incY, const double *P)
{
   if (N > 0)
   {
      if (incX < 0)
      {
         if (incY < 0) ATL_drotm(N, X, -incX, Y, -incY, P);
         else ATL_drotm(N, X+(1-N)*incX, incX, Y, incY, P);
      }
      else if (incY < 0) ATL_drotm(N, X+(N-1)*incX, -incX, Y, -incY, P);
      else ATL_drotm(N, X, incX, Y, incY, P);
   }
}
