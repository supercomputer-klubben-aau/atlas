/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"

void ATL_UGEMV(ATL_CINT M, ATL_CINT N, const TYPE *A, ATL_CINT lda,
               const TYPE *X, TYPE *Y)
/*
 *  y = [0,1]*y + A*x, A is MxN,  len(X) = N, len(Y) = M
 */
{
   ATL_CINT lda2 = lda+lda, N2 = N + N;
   ATL_INT j;

   #ifdef BETA0
      Mjoin(PATL,cpsc)(M, X, A, 1, Y, 1);
      A += lda2;
      j = 2;
   #else
      j=0;
   #endif
   for (; j < N2; j += 2, A += lda2)
      Mjoin(PATL,axpy)(M, X+j, A, 1, Y, 1);
}
