/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,geadd)
   (const int M, const int N, const SCALAR alpha, const TYPE *A, const int lda,
    const SCALAR beta, TYPE *C, const int ldc)
/*
 * C <- alpha*A + beta*C
 */
{
#ifdef TREAL
   if (beta == ATL_rzero) Mjoin(PATL,gemove)(M, N, alpha, A, lda, C, ldc);
   else if (alpha == ATL_rzero) Mjoin(PATL,gescal)(M, N, beta, C, ldc);
   else if (alpha == ATL_rone)
   {
      if (beta == ATL_rone)
         Mjoin(PATL,geadd_a1_b1)(M, N, alpha, A, lda, beta, C, ldc);
      else if (beta == ATL_rzero)
         Mjoin(PATL,geadd_a1_b0)(M, N, alpha, A, lda, beta, C, ldc);
      else Mjoin(PATL,geadd_a1_bX)(M, N, alpha, A, lda, beta, C, ldc);
   }
   else
   {
      if (beta == ATL_rone)
         Mjoin(PATL,geadd_aX_b1)(M, N, alpha, A, lda, beta, C, ldc);
      else if (beta == ATL_rzero)
         Mjoin(PATL,geadd_aX_b0)(M, N, alpha, A, lda, beta, C, ldc);
      else Mjoin(PATL,geadd_aX_bX)(M, N, alpha, A, lda, beta, C, ldc);
   }
#else
   const int AlphaIsReal = (alpha[1] == ATL_rzero);
   const int AlphaIsOne = (AlphaIsReal && *alpha == ATL_rone);
   const int AlphaIsZero = (AlphaIsReal && *alpha == ATL_rzero);
   const int BetaIsReal = (beta[1] == ATL_rzero);
   const int BetaIsOne = (BetaIsReal && *beta == ATL_rone);
   const int BetaIsZero = (BetaIsReal && *beta == ATL_rzero);

   if (BetaIsZero) Mjoin(PATL,gemove)(M, N, alpha, A, lda, C, ldc);
   else if (AlphaIsZero) Mjoin(PATL,gescal)(M, N, beta, C, ldc);
   else if (AlphaIsOne)
   {
      if (BetaIsOne)
         Mjoin(PATL,geadd_a1_b1)(M, N, alpha, A, lda, beta, C, ldc);
      else if (BetaIsZero)
         Mjoin(PATL,geadd_a1_b0)(M, N, alpha, A, lda, beta, C, ldc);
      else if (BetaIsReal)
         Mjoin(PATL,geadd_a1_bXi0)(M, N, alpha, A, lda, beta, C, ldc);
      else Mjoin(PATL,geadd_a1_bX)(M, N, alpha, A, lda, beta, C, ldc);
   }
   else if (AlphaIsReal)
   {
      if (BetaIsOne)
         Mjoin(PATL,geadd_aXi0_b1)(M, N, alpha, A, lda, beta, C, ldc);
      else if (BetaIsZero)
         Mjoin(PATL,geadd_aXi0_b0)(M, N, alpha, A, lda, beta, C, ldc);
      else if (BetaIsReal)
         Mjoin(PATL,geadd_aXi0_bXi0)(M, N, alpha, A, lda, beta, C, ldc);
      else Mjoin(PATL,geadd_aXi0_bX)(M, N, alpha, A, lda, beta, C, ldc);
   }
   else
   {
      if (BetaIsOne)
         Mjoin(PATL,geadd_aX_b1)(M, N, alpha, A, lda, beta, C, ldc);
      else if (BetaIsZero)
         Mjoin(PATL,geadd_aX_b0)(M, N, alpha, A, lda, beta, C, ldc);
      else if (BetaIsReal)
         Mjoin(PATL,geadd_aX_bXi0)(M, N, alpha, A, lda, beta, C, ldc);
      else Mjoin(PATL,geadd_aX_bX)(M, N, alpha, A, lda, beta, C, ldc);
   }
#endif
}
