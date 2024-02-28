/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,cplxdivide)
   (ATL_CINT N, const TYPE *beta, TYPE *X, ATL_CINT incx, TYPE *Y, ATL_CINT incy)
/*
 * Y(:) = X(:)/beta, wt division done with safe complex arithmetic.
 * It is OK for Y & X to be the same pointer
 * This code is straight adaptation of LAPACK's DLADIV, which comes from
 * the algorithm developed by Robert L. Smith (Art of Comp Prog, Vol.2 p.195)
 */
{
   ATL_CINT incY=incy+incy, incX = incx+incx;
   ATL_INT i;
   const register TYPE rb = beta[0], ib = beta[1];
   register TYPE rx, ix, e, f;

   if (Mabs(ib) < Mabs(rb))
   {
      e = ib / rb;
      f = rb + ib*e;
      for (i=N; i; i--, X += incX, Y += incY)
      {
         rx = *X; ix = X[1];
         Y[0] = (rx + ix*e) / f;
         Y[1] = (ix - rx*e) / f;
      }
   }
   else
   {
      e = rb / ib;
      f = ib + rb*e;
      for (i=N; i; i--, X += incX, Y += incY)
      {
         rx = *X; ix = X[1];
         Y[0] = (ix + rx*e) / f;
         Y[1] = (ix*e - rx) / f;
      }
   }
}
