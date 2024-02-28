/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"

void Mjoin(PATL,gediff)
   (const int M0, const int N, const TYPE *A, const int lda0,
    const TYPE *B, const int ldb0, TYPE *C, const int ldc0)
{
   const int M = M0 SHIFT, lda = lda0 SHIFT, ldb = ldb0 SHIFT, ldc = ldc0 SHIFT;
   int i, j;

   for (j=0; j < N; j++)
   {
      for (i=0; i != M; i++) C[i] = A[i] - B[i];
      A += lda;
      B += ldb;
      C += ldc;
   }
}
