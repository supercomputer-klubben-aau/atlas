/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"

void Mjoin(PATL,vdiff)(const int N, const TYPE *X, const int incX,
                       const TYPE *Y, const int incY, TYPE *Z, const int incZ)
/*
 * Z <- X - Y
 */
{
   int i;
   #ifdef TREAL
      for (i=N; i; i--, X += incX, Y += incY, Z += incZ) *Z = *X - *Y;
   #else
      const int incx = incX<<1, incy=incY<<1, incz = incZ<<1;
      for (i=N; i; i--, X += incx, Y += incy, Z += incz)
      {
         *Z = *X - *Y;
         Z[1] = X[1] - Y[1];
      }
   #endif
}
