/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2008 R. Clint Whaley
 */
#include "atlas_misc.h"
#define NB 40  /* all three matrices fit in 128 entry TLB */
#ifdef Conj_
   #define ApBt_NB Mjoin(PATL,geApBc_NB)
   #define ApAt_NB Mjoin(PATL,heApAc_NB)
void Mjoin(PATL,heApAc)
#else
   #define ApBt_NB Mjoin(PATL,geApBt_NB)
   #define ApAt_NB Mjoin(PATL,syApAt_NB)
void Mjoin(PATL,syApAt)
#endif
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const TYPE *A, ATL_CINT lda,
    const SCALAR beta, TYPE *C, ATL_CINT ldc)
/*
 * C <- A + A', C is Upper or Lower symmetric
 * Does by blocks in order to prevent massive TLB problems.
 */
{
   ATL_INT i, j, m, n, istart, iend;

   for (j=0; j < N; j += NB)
   {
      n = N - j;
      n = Mmin(NB, n);
      if (Uplo == AtlasLower) { istart = j; iend = N; }
      else { istart = 0; iend = j+NB; }

      for (i=istart; i < iend; i += NB)
      {
         m = N - i;
         m = Mmin(NB, m);
         if (i != j)
            ApBt_NB(m, n, A+((i+j*lda)SHIFT), lda, A+((j+i*lda)SHIFT),
                    lda, beta, C+((i+j*ldc)SHIFT), ldc);
         else
            ApAt_NB(Uplo, n, A+((i+j*lda)SHIFT), lda, beta,
                    C+((i+j*ldc)SHIFT), ldc);
      }
   }
}
