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

void catlas_saxpby(const int N, const float alpha, const float *X,
                   const int incX, const float beta, float *Y, const int incY)
{
   int incx=incX, incy=incY;

   if (N > 0)
   {
      if (incX >= 0 && incY >= 0) goto L1;
      else if (incY < 0)
      {
         if (incX < 0) { incx = -incX; incy = -incY; }
         else Y -= ((N-1))*incY;
      }
      else if (incX < 0) X -= ((N-1))*incX;
L1:
      ATL_saxpby(N, alpha, X, incx, beta, Y, incy);
   }
}
