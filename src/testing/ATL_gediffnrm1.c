/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"
TYPE Mjoin(PATL,gediffnrm1)
   (const int M, const int N, const TYPE *A, const int lda,
    const TYPE *B, const int ldb)
/*
 * Calculates the 1-norm of (A-B)
 */
{
   const int lda2 = lda SHIFT, ldb2 = ldb SHIFT;
   const int M2 = M SHIFT;
   int i, j;
   TYPE max=0.0, t0;

   for (j=0; j < N; j++)
   {
      t0 = ATL_rzero;
      for (i=0; i != M2; i++) t0 += Mabs(A[i] - B[i]);
      if (t0 != t0)
         return(t0);
      if (t0 > max) max = t0;
      A += lda2;
      B += ldb2;
   }
   return(max);
}
