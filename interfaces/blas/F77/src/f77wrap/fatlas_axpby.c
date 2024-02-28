/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_f77blas.h"
#include "atlas_f77wrap.h"
void F77axpby(const F77_INTEGER *N, const TYPE *alpha, const TYPE *X,
              const F77_INTEGER *incX, const TYPE *beta,
              TYPE *Y, const F77_INTEGER *incY)
{
   int incx=(*incX), incy=(*incY);

   if (*N > 0)
   {
      if (incx >= 0 && incy >= 0) goto L1;
      else if (incy < 0)
      {
         if (incx < 0) { incx = -incx; incy = -incy; }
         else Y -= ((*N-1)SHIFT)*incy;
      }
      else if (incx < 0) X -= ((*N-1)SHIFT)*incx;
L1:
      #ifdef TREAL
         Mjoin(PATL,axpby)(*N, *alpha, X, incx, *beta, Y, incy);
      #else
         Mjoin(PATL,axpby)(*N, alpha, X, incx, beta, Y, incy);
      #endif
   }
}
