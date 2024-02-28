/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

void Mjoin(Mjoin(PATL,hecopy),UploNM)
   (ATL_CSZT  N, const TYPE *A, ATL_CSZT lda, TYPE *C);

void Mjoin(Mjoin(PATL,hemmL),UploNM)
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc)
{
   TYPE *a;
   void *va;

   if (N > HEMM_Xover)
   {
      va = malloc(ATL_Cachelen + (ATL_MulBySize(M)*M));
      ATL_assert(va);
      a = ATL_AlignPtr(va);
      Mjoin(Mjoin(PATL,hecopy),UploNM)(M, A, lda, a);
      ATL_gemm(AtlasNoTrans, AtlasNoTrans, M, N, M, alpha, a, M, B, ldb,
               beta, C, ldc);
      free(va);
   }
   else Mjoin(PATL,refhemm)(AtlasLeft, Uplo_, M, N, alpha, A, lda, B, ldb,
                            beta, C, ldc);
}
