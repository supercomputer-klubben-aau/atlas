/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#ifdef Upper_
   #define her2k_put Mjoin(PATL,her2k_putU)
   #ifdef Transpose_
      int Mjoin(PATL,her2kUC)
   #else
      int Mjoin(PATL,her2kUN)
   #endif
#else
   #define her2k_put Mjoin(PATL,her2k_putL)
   #ifdef Transpose_
      int Mjoin(PATL,her2kLC)
   #else
      int Mjoin(PATL,her2kLN)
   #endif
#endif
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc)
{
   int i;
   void *vc=NULL;
   TYPE *c;
   const TYPE beta =*( (const TYPE *)vbeta  );
   const TYPE zero[2]={0.0, 0.0};

   i = ATL_MulBySize(N)*N;
   if (i <= ATL_MaxMalloc) vc = malloc(ATL_Cachelen+i);
   if (vc == NULL) return(1);
   c = ATL_AlignPtr(vc);
   #ifdef Transpose_
      ATL_gemm(AtlasConjTrans, AtlasNoTrans, N, N, K, valpha, A, lda, B, ldb,
   #else
      ATL_gemm(AtlasNoTrans, AtlasConjTrans, N, N, K, valpha, A, lda, B, ldb,
   #endif
               zero, c, N);
   if ( beta == 1.0 ) Mjoin(her2k_put,_b1)(N, c, vbeta, C, ldc);
   else if ( beta == 0.0 ) Mjoin(her2k_put,_b0)(N, c, vbeta, C, ldc);
   else Mjoin(her2k_put,_bXi0)(N, c, vbeta, C, ldc);
   free(vc);
   return(0);
}
