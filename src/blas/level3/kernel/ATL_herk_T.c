/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#define her_put Mjoin(Mjoin(PATL,heput),UploNM)
int Mjoin(Mjoin(Mjoin(PATL,herk),UploNM),C)
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc)
{
#if 0
   Mjoin(PATL,herk_IP)(Uplo_, AtlasConjTrans, N, K, *((TYPE*)valpha),
                       A, lda, *((TYPE*)vbeta), C, ldc);
#else
   void *vc;
   TYPE *c;
   TYPE alpha[2];
   const TYPE beta =*( (const TYPE *)vbeta  );
   const TYPE zero[2]={0.0, 0.0};

   alpha[0] = *( (const TYPE *)valpha );
   if (K > HERK_Xover)
   {
#if defined(USE_AMM)
      TYPE BE[2] = {beta, ATL_rzero};
      int Mjoin(PATL,ammm_syrk)
         (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT N, ATL_CSZT K,
          const SCALAR alpha, const TYPE *A, ATL_CSZT lda, const SCALAR beta,
          TYPE *C, ATL_CSZT ldc);
      alpha[1] = 0.0;
      return(Mjoin(PATL,ammm_syrk)(Uplo_, AtlasConjTrans, N, K, alpha,
                                   A, lda, BE, C, ldc));
#else
      alpha[1] = 0.0;
      vc = malloc(ATL_Cachelen+ATL_MulBySize(N)*N);
      ATL_assert(vc);
      c = ATL_AlignPtr(vc);
      ATL_gemm(AtlasConjTrans, AtlasNoTrans, N, N, K, alpha, A, lda, A, lda,
               zero, c, N);
      if (beta == 1.0) Mjoin(her_put,_b1)(N, c, vbeta, C, ldc);
      else if (beta == 0.0) Mjoin(her_put,_b0)(N, c, vbeta, C, ldc);
      else Mjoin(her_put,_bXi0)(N, c, vbeta, C, ldc);
      free(vc);
#endif
   }
   else Mjoin(PATL,refherk)(Uplo_, AtlasTrans, N, K, *alpha, A, lda,
                            beta, C, ldc);
#endif
   return(0);
}
