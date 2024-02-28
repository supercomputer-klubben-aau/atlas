/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,trsetL)
   (ATL_CINT M, ATL_CINT N, const SCALAR alpha, const SCALAR beta,
    TYPE *A, ATL_CINT lda)
/*
 * Sets main diagonal to beta, rest of lower triangle to alpha, does not
 * touch upper triangle
 */
{

   ATL_INT j;
   ATL_CINT ldap1 = (lda+1)SHIFT;

   for (j=0; j < N; j++, A += ldap1)
   {
      #ifdef TCPLX
         *A = *beta;
         A[1] = beta[1];
      #else
         *A = beta;
      #endif
      if (N-j-1)
         Mjoin(PATL,set)(N-j-1, alpha, A+(1 SHIFT), 1);
   }
}
