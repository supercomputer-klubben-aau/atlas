/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,cplxinvert)(ATL_CINT N, TYPE *X, ATL_CINT incx,
                            TYPE *Y, ATL_CINT incy)
/*
 * Y(:) = 1 / X(:)
 * Invert N complex scalars held in X, and write answer to Y.
 * X & Y can be same space
 */
{
   int i;
   const TYPE one=1.0, none=(-1.0);
   ATL_CINT incX=incx+incx, incY=incy+incy;
   register TYPE rtmp, itmp, t0;

   for (i=N; i; i--, X += incX, Y += incY)
   {
      rtmp = *X;
      itmp = X[1];
      if (Mabs(itmp) <= Mabs(rtmp))
      {
         t0 = itmp / rtmp;
         *Y = rtmp = one / (rtmp + itmp*t0);
         Y[1] = -rtmp * t0;
      }
      else
      {
         t0 = rtmp / itmp;
         Y[1] = rtmp = none / (itmp + rtmp*t0);
         *Y = -t0 * rtmp;
      }
   }
}
