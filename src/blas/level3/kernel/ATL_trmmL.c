/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#ifdef Upper_
   #ifdef Transpose_
      #ifdef UnitDiag_
         #define ATLP UTU
         #define ATL_trcopy Mjoin(PATL,trcopyU2U_U)
      #else
         #define ATLP UTN
         #define ATL_trcopy Mjoin(PATL,trcopyU2U_N)
      #endif
   #elif defined(ConjTrans_)
      #ifdef UnitDiag_
         #define ATLP UCU
         #define ATL_trcopy Mjoin(PATL,trcopyU2Uc_U)
      #else
         #define ATLP UCN
         #define ATL_trcopy Mjoin(PATL,trcopyU2Uc_N)
      #endif
   #else
      #ifdef UnitDiag_
         #define ATL_trcopy Mjoin(PATL,trcopyU2L_U)
         #define ATLP UNU
      #else
         #define ATL_trcopy Mjoin(PATL,trcopyU2L_N)
         #define ATLP UNN
      #endif
   #endif
#else
   #ifdef Transpose_
      #ifdef UnitDiag_
         #define ATL_trcopy Mjoin(PATL,trcopyL2L_U)
         #define ATLP LTU
      #else
         #define ATL_trcopy Mjoin(PATL,trcopyL2L_N)
         #define ATLP LTN
      #endif
   #elif defined(ConjTrans_)
      #ifdef UnitDiag_
         #define ATL_trcopy Mjoin(PATL,trcopyL2Lc_U)
         #define ATLP LCU
      #else
         #define ATL_trcopy Mjoin(PATL,trcopyL2Lc_N)
         #define ATLP LCN
      #endif
   #else
      #ifdef UnitDiag_
         #define ATL_trcopy Mjoin(PATL,trcopyL2U_U)
         #define ATLP LNU
      #else
         #define ATL_trcopy Mjoin(PATL,trcopyL2U_N)
         #define ATLP LNN
      #endif
   #endif
#endif

void Mjoin(Mjoin(PATL,trmmL),ATLP)
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc)
{
   #ifdef TREAL
      const SCALAR alpha=*( (const SCALAR *)valpha );
      const SCALAR one=1.0, zero=0.0;
   #else
      const TYPE zero[2]={0.0,0.0};
      #define alpha valpha
   #endif
   void *va;
   TYPE *a;

   if (N > TRMM_Xover)
   {
      va = malloc(ATL_Cachelen + ATL_MulBySize(M)*M);
      ATL_assert(va);
      a = ATL_AlignPtr(va);
      #ifdef TREAL
         if ( SCALAR_IS_ONE(alpha) ) Mjoin(ATL_trcopy,_a1)(M, alpha, A, lda, a);
         else Mjoin(ATL_trcopy,_aX)(M, alpha, A, lda, a);
         ATL_assert(M <= ATL_VWopgen_MAX_KB); /* aliasing chk */
         ATL_almm(AtlasTrans, AtlasNoTrans, M, N, M, one, a, M, C, ldc,
                  zero, C, ldc);
      #else
         ATL_trcopy(M, A, lda, a);
         ATL_assert(M <= ATL_VWopgen_MAX_KB); /* aliasing chk */
         ATL_almm(AtlasTrans, AtlasNoTrans, M, N, M, valpha, a, M, C, ldc,
                  zero, C, ldc);
      #endif
      free(va);
   }
   else Mjoin(PATL,reftrmm)(AtlasLeft, Uplo_, Trans_, Unit_, M, N, alpha,
                            A, lda, C, ldc);
}
