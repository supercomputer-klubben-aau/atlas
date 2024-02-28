/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lvl2.h"
#include "atlas_lvl3.h"
#include "atlas_cacheedge.h"
#include "atlas_cache.h"
/*
 * If I don't believe CacheEdge setting (or not set), set L2 size to 4*L1
 */
#ifdef CacheEdge
   #if CacheEdge > 4194304 || CacheEdge == 0
      #define MY_CE (4*ATL_MulBySize(L1C_ELTS))
   #else
      #define MY_CE CacheEdge
   #endif
#else
   #define MY_CE (4*ATL_MulBySize(L1C_ELTS))
#endif

void Mjoin(PATL,gemv)
   (const enum ATLAS_TRANS TA, ATL_CINT M, ATL_CINT N, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *X, ATL_CINT incX,
    const SCALAR beta, TYPE *Y, ATL_CINT incY)
/*
 * ATL_gemv is a wrapper that chooses which context-tuned GEMV to call.
 * Supported contexts are tuned for in-L1 performance (_L1), tuned for
 * in-L2 performance (_L2), and tuned for out-of-cache (no suffix).
 * Right now, we do this based sheerly on operand size, since this matches
 * most common (we think!) LAPACK usage, like in factorizations.  If we
 * had access to the hardware counters, a better way to would be to access
 * say 10/100 elements of A, and call one of these contexts based on how
 * many cache misses were generated
 */
{
   const size_t opsize = (M*N+M+N)*sizeof(TYPE)SHIFT, Abeg, Aend, Obeg, Oend;
   int OL;

   if (TA == AtlasNoTrans)
   {
      if (opsize > MY_CE)
        Mjoin(PATL,gemvN)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
      else if (opsize > ATL_MulBySize(L1C_ELTS))
        Mjoin(PATL,gemvN_L2)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
      else
        Mjoin(PATL,gemvN_L1)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
   }
   else
   {
   #ifdef TCPLX
      if (TA == AtlasTrans)
      {
   #endif
      if (opsize > MY_CE)
         Mjoin(PATL,gemvT)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
      else if (opsize > ATL_MulBySize(L1C_ELTS))
         Mjoin(PATL,gemvT_L2)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
      else
         Mjoin(PATL,gemvT_L1)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
   #ifdef TCPLX
      }
      else if (TA == AtlasConjTrans)
      {
         if (opsize > MY_CE)
            Mjoin(PATL,gemvCT)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
         else if (opsize > ATL_MulBySize(L1C_ELTS))
            Mjoin(PATL,gemvCT_L2)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
         else
            Mjoin(PATL,gemvCT_L1)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
      }
      else /* if (TA == AtlasConj) */
      {
         if (opsize > MY_CE)
            Mjoin(PATL,gemvCN)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
         else if (opsize > ATL_MulBySize(L1C_ELTS))
            Mjoin(PATL,gemvCN_L2)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
         else
            Mjoin(PATL,gemvCN_L1)(M, N, alpha, A, lda, X, incX, beta, Y, incY);
      }
   #endif
   }
}
