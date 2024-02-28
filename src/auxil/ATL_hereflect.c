/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,hereflect)
   (const enum ATLAS_UPLO Uplo, const int N, TYPE *C, const int ldc)
/*
 * If Uplo == Lower, reflect lower triangle into upper,
 * If Uplo == Upper, reflect upper triangle into lower.
 */
{
   int j;
   #ifdef TCPLX
      const size_t ldc2 = ldc+ldc, incC = ldc2+2;
   #else
      const size_t incC = ldc+1, ldc2=ldc;
   #endif
   TYPE *pC;
   if (Uplo == AtlasLower)
   {
      for (j=0; j < N-1; j++, C += incC)
         Mjoin(PATL,copyConj)(N-j-1, C+(1 SHIFT), 1, C+ldc2, ldc);
   }
   else
   {
      pC = C + (((size_t)(N-1))SHIFT);
      C += ldc2*(N-1);
      for (j=0; j < N-1; j++, C -= ldc2, pC -= (1 SHIFT))
         Mjoin(PATL,copyConj)(N-j-1, C, 1, pC, ldc);
   }
}
#ifdef ldc2
   #undef ldc2
#endif
