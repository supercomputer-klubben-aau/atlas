/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2008 R. Clint Whaley
 */
#include "atlas_misc.h"
#ifdef Conj_
   #define axpbyT Mjoin(PATL,axpbyConj)
void Mjoin(PATL,heApAc_NB)
#else
   #define axpbyT Mjoin(PATL,axpby)
void Mjoin(PATL,syApAt_NB)
#endif
   (const enum ATLAS_UPLO uplo, ATL_CINT N, const TYPE *A, ATL_CINT lda,
    const SCALAR beta, TYPE *C, ATL_CINT ldc)
/*
 * C <- beta*C + A + A', C is Upper or Lower symmetric
 */
{
   const TYPE *Ac, *Ar;
   #ifdef Conj_
      TYPE *C0 = C;
   #endif
   #ifdef TREAL
      #define ldc2 ldc
      #define lda2 lda
      #define ONE ATL_rone
      const int ldap1 = lda+1, ldcp1 = ldc+1;
   #else
      const int lda2 = lda+lda, ldc2 = ldc+ldc, ldap1 = lda2+2, ldcp1 = ldc2+2;
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
   #endif
   int j;

   if (uplo == AtlasUpper)
   {
      Ac = Ar = A;
      for (j=0; j < N; j++)
      {
         Mjoin(PATL,axpby)(j+1, ONE, Ac, 1, beta, C, 1);
         axpbyT(j+1, ONE, Ar, lda, ONE, C, 1);
         C += ldc2;
         Ar += 1 SHIFT;
         Ac += lda2;
      }
   }
   else
   {
      for (j=0; j < N; j++)
      {
         Mjoin(PATL,axpby)(N-j, ONE, A, 1, beta, C, 1);
         axpbyT(N-j, ONE, A, lda, ONE, C, 1);
         C += ldcp1;
         A += ldap1;
      }
   }
#ifdef Conj_
   Mjoin(PATLU,zero)(N, C0+1, ldcp1);  /* zero imag part of diagonal */
#endif
}
#ifdef TREAL
   #undef lda2
   #undef ldb2
   #undef ONE
#endif

