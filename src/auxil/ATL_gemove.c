/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,gemove)
   (const int M, const int N, const SCALAR alpha, const TYPE *A, const int lda,
    TYPE *C, const int ldc)
/*
 * C <- alpha * A
 */
{
#ifdef TREAL
   if (alpha == ATL_rone) Mjoin(PATL,gecopy)(M, N, A, lda, C, ldc);
   else if (alpha == ATL_rzero) Mjoin(PATL,gezero)(M, N, C, ldc);
   else Mjoin(PATL,gemove_aX)(M, N, alpha, A, lda, C, ldc);
#else
   TYPE ralpha = *alpha;

   if (alpha[1] == ATL_rzero)
   {
      if (ralpha == ATL_rone) Mjoin(PATL,gecopy)(M, N, A, lda, C, ldc);
      else if (ralpha == ATL_rzero) Mjoin(PATL,gezero)(M, N, C, ldc);
      else Mjoin(PATL,gemove_aXi0)(M, N, alpha, A, lda, C, ldc);
   }
   else Mjoin(PATL,gemove_aX)(M, N, alpha, A, lda, C, ldc);
#endif
}
