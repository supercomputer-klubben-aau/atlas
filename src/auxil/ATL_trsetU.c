/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,trsetU)
   (ATL_CINT M, ATL_CINT N, const SCALAR alpha, const SCALAR beta,
    TYPE *A, ATL_CINT lda)
/*
 * Sets main diagonal to beta, rest of upper triangle to alpha, does not
 * touch lower triangle
 */
{
   ATL_INT j;
   #ifdef TCPLX
      ATL_CINT lda2 = lda+lda;
   #else
      #define lda2 lda
   #endif
   for (j=0; j < N; j++, A += lda2)
   {
      if (j)
         Mjoin(PATL,set)(j, alpha, A, 1);
      #ifdef TCPLX
         A[j+j] = *beta;
         A[j+j+1] = beta[1];
      #else
         A[j] = beta;
      #endif
   }
}
#ifndef TCPLX
   #undef lda2
#endif
