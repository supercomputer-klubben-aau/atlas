/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#ifdef TREAL
   void Mjoin(Mjoin(Mjoin(PATL,sycopy),UploNM),_a1)
      (ATL_CSZT  N, const SCALAR alpha0, const TYPE *A, ATL_CSZT lda, TYPE *C);
   void Mjoin(Mjoin(Mjoin(PATL,sycopy),UploNM),_aX)
      (ATL_CSZT  N, const SCALAR alpha0, const TYPE *A, ATL_CSZT lda, TYPE *C);
#else
   void Mjoin(Mjoin(PATL,sycopy),UploNM)
      (ATL_CSZT  N, const TYPE *A, ATL_CSZT lda, TYPE *C);
#endif

void Mjoin(Mjoin(PATL,symmL),UploNM)
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc)
{
   #ifdef TREAL
      const SCALAR alpha=*( (const SCALAR *)valpha );
      const SCALAR beta =*( (const SCALAR *)vbeta  );
      const SCALAR one=1.0;
   #else
      #define alpha valpha
      #define beta vbeta
   #endif
   TYPE *a;
   void *va;

   if (N > SYMM_Xover)
   {
      va = malloc(ATL_Cachelen + (ATL_MulBySize(M)*M));
      ATL_assert(va);
      a = ATL_AlignPtr(va);
      #ifdef TREAL
         if ( SCALAR_IS_ONE(alpha) )
            Mjoin(Mjoin(Mjoin(PATL,sycopy),UploNM),_a1)(M, alpha, A, lda, a);
         else Mjoin(Mjoin(Mjoin(PATL,sycopy),UploNM),_aX)(M, alpha, A, lda, a);
         ATL_gemm(AtlasTrans, AtlasNoTrans, M, N, M, one, a, M, B, ldb,
                  beta, C, ldc);
      #else
         Mjoin(Mjoin(PATL,sycopy),UploNM)(M, A, lda, a);
         ATL_gemm(AtlasTrans, AtlasNoTrans, M, N, M, valpha, a, M, B, ldb,
                  vbeta, C, ldc);
      #endif
      free(va);
   }
   else Mjoin(PATL,refsymm)(AtlasLeft, Uplo_, M, N, alpha, A, lda, B, ldb,
                            beta, C, ldc);
}
