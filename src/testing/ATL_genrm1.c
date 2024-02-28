/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,genrm1)(const int M, const int N, const TYPE *A, const int lda)
/*
 * Calculates the 1-norm of a general rectangular matrix
 */
{
   const int lda2 = lda SHIFT;
   int j;
   TYPE max=0.0, t0;

   for (j=0; j < N; j++)
   {
      t0 = Mjoin(PATL,asum)(M, A, 1);
      if (t0 != t0)
         return(t0);
      if (t0 > max) max = t0;
      A += lda2;
   }
   return(max);
}
