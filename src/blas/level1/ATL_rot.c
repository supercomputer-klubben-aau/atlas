/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#include "atlas_misc.h"
#include "atlas_level1.h"

#ifdef TREAL

void Mjoin(PATL,rot)(const int N, TYPE *X, const int incX,
                     TYPE *Y, const int incY, const TYPE c, const TYPE s)
{
   int i;
   TYPE tmp;

   if (c != ATL_rone || s != ATL_rzero)
   {
      if (incX == 1 && incY == 1)
      {
         for (i=0; i != N; i++)
         {
            tmp = c * X[i] + s * Y[i];
            Y[i] = c * Y[i] - s * X[i];
            X[i] = tmp;
         }
      }
      else
      {
         for (i=N; i; i--, Y += incY, X += incX)
         {
            tmp = c * *X + s * *Y;
            *Y = c * *Y - s * *X;
            *X = tmp;
         }
      }
   }
}

#else /* complex routine */
void Mjoin(Mjoin(PATL,UPR),rot)
   (const int N, TYPE *X, const int incx, TYPE *Y, const int incy,
    const TYPE c0, const TYPE s0)
{
   const register TYPE c = c0, s = s0;
   register TYPE rx, ix, ry, iy;
   const int incX = incx<<1, incY = incy<<1;
   int i;

   if (N > 0)
   {
      if (incx == 1 && incy == 1)
      {
         for (i=N; i; i--, X += 2, Y += 2)  /* maybe compiler unrolls */
         {
            rx = *X;  ix = X[1];
            ry = *Y;  iy = Y[1];
            *X   = c * rx + s * ry;
            X[1] = c * ix + s * iy;
            *Y   = c * ry - s * rx;
            Y[1] = c * iy - s * ix;
         }
      }
      else
      {
         for (i=N; i; i--, X += incX, Y += incY)  /* maybe compiler unrolls */
         {
            rx = *X;  ix = X[1];
            ry = *Y;  iy = Y[1];
            *X   = c * rx + s * ry;
            X[1] = c * ix + s * iy;
            *Y   = c * ry - s * rx;
            Y[1] = c * iy - s * ix;
         }
      }
   }
}
#endif
