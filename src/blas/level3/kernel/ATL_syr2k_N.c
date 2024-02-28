/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#ifdef Upper_
   #define syr2k_put Mjoin(PATL,syr2k_putU)
   int Mjoin(PATL,syr2kUN)
#else
   #define syr2k_put Mjoin(PATL,syr2k_putL)
   int Mjoin(PATL,syr2kLN)
#endif
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc)
{
   int i;
   void *vc=NULL;
   TYPE *c;
   #ifdef TREAL
      const SCALAR alpha=*( (const SCALAR *)valpha );
      const SCALAR beta =*( (const SCALAR *)vbeta  );
      const SCALAR one=1.0, zero=0.0;
   #else
      #define alpha valpha
      const TYPE *beta=vbeta;
      const TYPE one[2]={1.0,0.0}, zero[2]={0.0,0.0};
   #endif

   i = ATL_MulBySize(N)*N;
   if (i <= ATL_MaxMalloc) vc = malloc(ATL_Cachelen+i);
   if (vc == NULL) return(1);
   c = ATL_AlignPtr(vc);
   ATL_gemm(AtlasNoTrans, AtlasTrans, N, N, K, alpha, A, lda, B, ldb,
            zero, c, N);
   if ( SCALAR_IS_ONE(beta) ) Mjoin(syr2k_put,_b1)(N, c, beta, C, ldc);
   else if ( SCALAR_IS_ZERO(beta) ) Mjoin(syr2k_put,_b0)(N, c, beta, C, ldc);
   #ifdef TCPLX
      else if (SCALAR_IS_NONE(beta)) Mjoin(syr2k_put,_bn1)(N, c, beta, C, ldc);
      else if (beta[1] == *zero) Mjoin(syr2k_put,_bXi0)(N, c, beta, C, ldc);
   #endif
   else Mjoin(syr2k_put,_bX)(N, c, beta, C, ldc);
   free(vc);
   return(0);
}
