/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2011, 2010 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,geset)
   (ATL_CINT M, ATL_CINT N, const SCALAR alpha, const SCALAR beta,
    TYPE *A, ATL_CINT lda)
/*
 * Sets main diagonal to beta, rest of matrix to alpha
 */
{

   ATL_INT j;
   ATL_CINT MN = Mmin(M,N);
   #ifdef TCPLX
      ATL_CINT lda2 = lda+lda;
   #else
      #define lda2 lda
   #endif

#ifdef TCPLX
   if (*alpha == *beta && alpha[1] == beta[1])
#else
   if (alpha == beta)
#endif
   {
      for (j=0; j < N; j++, A += lda2)
         Mjoin(PATL,set)(M, alpha, A, 1);
      return;
   }
   for (j=0; j < MN; j++, A += lda2)
   {
      if (j)
         Mjoin(PATL,set)(j, alpha, A, 1);
      #ifdef TCPLX
         A[j+j] = *beta;
         A[j+j+1] = beta[1];
      #else
         A[j] = beta;
      #endif
      if (M-j-1)
         Mjoin(PATL,set)(M-j-1, alpha, A+((j+1) SHIFT), 1);
   }
   for (; j < N; j++, A += lda2)
      Mjoin(PATL,set)(M, alpha, A, 1);
}
#ifndef TCPLX
   #undef lda2
#endif
