/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

#ifdef Conj_
void Mjoin(PATL,scalConj)
/*
 * X <- alpha * conj(X)
 */
#else
void Mjoin(PATL,scal)
/*
 * X <- alpha * X
 */
#endif
   (const int N, const SCALAR alpha, TYPE *X, const int incX)
{
   int i;
   #ifdef TREAL
      if (alpha != ATL_rzero)
      {
         if (incX == 1) for (i=0; i != N; i++) X[i] *= alpha;
         else for (i=N; i; i--, X += incX) *X *= alpha;
      }
      else Mjoin(PATL,zero)(N, X, incX);
   #else
      const register TYPE ralpha = *alpha, ialpha = alpha[1];
      #ifdef Conj_
         const register TYPE conjal = -ralpha;
      #else
         #define conjal ralpha
      #endif
      register TYPE rtmp, itmp;
      int incx = incX<<1;

      if (ialpha == ATL_rzero)
      {
         if (ralpha != ATL_rzero)
         {
      #ifndef Conj_
            if (incX == 1) Mjoin(PATLU,scal)(N<<1, ralpha, X, incX);
            else
      #endif
            {
               for (i=N; i; i--, X += incx)
               {
                  *X *= ralpha;
                  X[1] *= conjal;
               }
            }
         }
         else Mjoin(PATL,zero)(N, X, incX);
      }
      else
      {
         if (incX == 1)
         {
            for (i=N; i; i--, X += 2)
            {
               rtmp = X[0];
               itmp = X[1];
               #ifdef Conj_
                  X[0] = rtmp * ralpha + itmp * ialpha;
                  X[1] = rtmp * ialpha - itmp * ralpha;
               #else
                  X[0] = rtmp * ralpha - itmp * ialpha;
                  X[1] = rtmp * ialpha + itmp * ralpha;
               #endif
            }
         }
         else
         {
            for (i=N; i; i--, X += incx)
            {
               rtmp = X[0];
               itmp = X[1];
               #ifdef Conj_
                  X[0] = rtmp * ralpha + itmp * ialpha;
                  X[1] = rtmp * ialpha - itmp * ralpha;
               #else
                  X[0] = rtmp * ralpha - itmp * ialpha;
                  X[1] = rtmp * ialpha + itmp * ralpha;
               #endif
            }
         }
      }
   #endif
}
