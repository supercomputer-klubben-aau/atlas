/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,infnrm)(const int N, const TYPE *X, const int incX)
{
   register TYPE max=0.0;
   register int i;
   if (N > 0)
   {
   #ifdef TCPLX
      const int incX2=incX+incX;
      for (i=0; i < N; i++, X += incX2)
      {
         register TYPE t0;
         t0 = Mabs(*X) + Mabs(X[1]);
         if (t0 != t0)
            return(t0);
         max = (max >= t0) ? max : t0;
      }
   #else
      for (i=0; i < N; i++, X += incX)
      {
         register TYPE t0;
         t0 = Mabs(*X);
         if (t0 != t0)
            return(t0);
         max = (max >= t0) ? max : t0;
      }
   #endif
   }
   return(max);
}
