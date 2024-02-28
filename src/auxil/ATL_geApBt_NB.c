/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2008 R. Clint Whaley
 */
#include "atlas_misc.h"
#ifdef Conj_
   #define axpbyT Mjoin(PATL,axpbyConj)
void Mjoin(PATL,geApBc_NB)
#else
   #define axpbyT Mjoin(PATL,axpby)
void Mjoin(PATL,geApBt_NB)
#endif
   (ATL_CINT M, ATL_CINT N, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc)
/*
 * C <- beta*C + A + B'; this routine needs small N, or you have TLB issues
 */
{
   #ifdef TREAL
      #define lda2 lda
      #define ldc2 ldc
      #define ONE ATL_rone
   #else
      ATL_CINT lda2 = lda+lda, ldc2 = ldc+ldc;
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
   #endif
   ATL_INT j;

   for (j=0; j < N; j++)
   {
      Mjoin(PATL,axpby)(M, ONE, A, 1, beta, C, 1);
      axpbyT(M, ONE, B, ldb, ONE, C, 1);
      C += ldc2;
      A += lda2;
      B += 1 SHIFT;
   }
}
#ifdef TREAL
   #undef lda2
   #undef ONE
#endif
