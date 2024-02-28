/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define SCPLX
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias1.h"
#endif
#include "atlas_level1.h"
#include "cblas.h"

void cblas_cdotc_sub(const int N, const void *X, const int incX,
                      const void *Y, const int incY, void * dot)
{
   const float *x = X, *y = Y;
   float *tdot=dot;
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
      ATL_cdotc_sub(N, x, incx, y, incy, dot);
   }
   else *tdot = tdot[1] = 0.0f;
}
