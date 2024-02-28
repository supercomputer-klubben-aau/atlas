/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,tradd)(const enum ATLAS_UPLO Uplo, ATL_CINT N, const TYPE *A,
                       ATL_CINT lda, const SCALAR beta, TYPE *C, ATL_CINT ldc)
/*
 * C <- beta*C + A, A & C are triangular matrices of order N
 */
{
   int j;

   #ifdef TCPLX
      ATL_CINT lda2 = lda+lda, ldc2 = ldc+ldc;
      ATL_CINT ldap1 = lda2 + 2, ldcp1 = ldc2+2;
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
   #else
      TYPE ONE = ATL_rone;
      ATL_CINT ldap1 = lda+1, ldcp1 = ldc+1;
      #define lda2 lda
      #define ldc2 ldc
   #endif
   if (Uplo == AtlasLower)
   {
      for (j=0; j < N; j++, A += ldap1, C += ldcp1)
         Mjoin(PATL,axpby)(N-j, ONE, A, 1, beta, C, 1);
   }
   else /* Upper */
   {
      for (j=0; j < N; j++, A += lda2, C += ldc2)
         Mjoin(PATL,axpby)(j+1, ONE, A, 1, beta, C, 1);
   }
}
