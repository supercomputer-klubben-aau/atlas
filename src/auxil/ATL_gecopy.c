/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,gecopy)
   (const int M0, const int N, const TYPE *A, const int lda,
    TYPE *C, const int ldc)
/*
 * C <- A;  copy backwards, so 1st cols of C are retained in cache for reuse
 */
{
   int i, j;
   #ifdef TREAL
      #define M M0
      const int incA = lda+lda, incC = ldc+ldc;
   #else
      const int M = M0<<1, incA = lda<<2, incC = ldc<<2;
   #endif
   const int n = N>>1;
   const TYPE *A0, *A1;
   TYPE *C0, *C1;

   A0 = A + (lda SHIFT)*(N-2);
   A1 = A0 + (lda SHIFT);
   C0 = C + (ldc SHIFT)*(N-2);
   C1 = C0 + (ldc SHIFT);
   for (j=n; j; j--, A0 -= incA, A1 -= incA, C0 -= incC, C1 -= incC)
   {
      for (i=M-1; i >= 0; i--)
      {
         C0[i] = A0[i];
         C1[i] = A1[i];
      }
   }
   if (N-n-n)
      for (i=M-1; i >= 0; i--)
         C[i] = A[i];
}
