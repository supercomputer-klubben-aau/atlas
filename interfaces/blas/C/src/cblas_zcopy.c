/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define DCPLX
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias1.h"
#endif
#include "atlas_level1.h"
#include "cblas.h"

void cblas_zcopy(const int N, const void *X, const int incX,
                      void *Y, const int incY)
{
   const double *x = X;
   double *y = Y;
   int incx = incX, incy = incY;

   if (N > 0)
   {
      if (incX < 0)
      {
         if (incY < 0) { incx = -incx; incy = -incY; }
         else x += -incX * ((N-1)<<1);
      }
      else if (incY < 0)
      {
         incy = -incy;
         incx = -incx;
         x += (N-1)*(incX<<1);
      }
      ATL_zcopy(N, x, incx, y, incy);
   }
}
