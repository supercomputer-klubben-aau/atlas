/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#define syr_put Mjoin(Mjoin(PATL,trput),UploNM)
int Mjoin(Mjoin(Mjoin(PATL,syrk),UploNM),N)
      (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A,
       ATL_CSZT lda, const void *vbeta, void *C, ATL_CSZT ldc)
{
#ifdef USE_AMM
   #ifdef TREAL
      Mjoin(PATL,syrk_IP)(Uplo_, AtlasNoTrans, N, K, *((TYPE*)valpha),
                          A, lda, *((TYPE*)vbeta), C, ldc);
   #else
      Mjoin(PATL,syrk_IP)(Uplo_, AtlasNoTrans, N, K, valpha, A, lda,
                          vbeta, C, ldc);
   #endif
#else
   void *vc;
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

   if (K > SYRK_Xover)
   {
#if defined(USE_AMM)
      int Mjoin(PATL,ammm_syrk)
         (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT N, ATL_CSZT K,
          const SCALAR alpha, const TYPE *A, ATL_CSZT lda, const SCALAR beta,
          TYPE *C, ATL_CSZT ldc);
      return(Mjoin(PATL,ammm_syrk)(Uplo_, AtlasNoTrans, N, K, alpha,
                                   A, lda, beta, C, ldc));
#else
      vc = malloc(ATL_Cachelen+ATL_MulBySize(N)*N);
      ATL_assert(vc);
      c = ATL_AlignPtr(vc);
      ATL_gemm(AtlasNoTrans, AtlasTrans, N, N, K, alpha, A, lda, A, lda,
               zero, c, N);
      if ( SCALAR_IS_ONE(beta) ) Mjoin(syr_put,_b1)(N, c, beta, C, ldc);
      else if ( SCALAR_IS_ZERO(beta) ) Mjoin(syr_put,_b0)(N, c, beta, C, ldc);
      #ifdef TCPLX
         else if ( SCALAR_IS_NONE(beta) )
              Mjoin(syr_put,_bn1)(N, c, beta, C, ldc);
         else if (beta[1] == *zero) Mjoin(syr_put,_bXi0)(N, c, beta, C, ldc);
      #endif
      else Mjoin(syr_put,_bX)(N, c, beta, C, ldc);
      free(vc);
#endif
   }
   else Mjoin(PATL,refsyrk)(Uplo_, AtlasNoTrans, N, K, alpha, A, lda,
                            beta, C, ldc);
#endif
   return(0);
}
