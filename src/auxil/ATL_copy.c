/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

#ifdef Conj_
void Mjoin(PATL,copyConj)
/*
 * Y <- conj(X)
 */
#else
void Mjoin(PATL,copy)
/*
 * Y <- X
 */
#endif
   (const int N, const TYPE *X, const int incX, TYPE *Y, const int incY)
{
   #ifdef TREAL
      #define M N
   #else
      const int M = N<<1, incx = incX<<1, incy = incY<<1;
   #endif
   int i;

   #ifndef Conj_
      if (incY == 1 && incX == 1) for (i=0; i != M; i++) Y[i] = X[i];
      else
   #endif
   #ifdef TREAL
         for (i=N; i; i--, X += incX, Y += incY) *Y = *X;
   #else
      {
         for (i=N; i; i--, X += incx, Y += incy)
         {
            *Y = *X;
            #ifdef Conj_
               Y[1] = -X[1];
            #else
               Y[1] = X[1];
            #endif
         }
      }
   #endif
}
