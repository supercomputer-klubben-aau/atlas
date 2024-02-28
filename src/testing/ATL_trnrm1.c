/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,trnrm1)(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                        const int N, const TYPE *A, const int lda)
/*
 * Calculates the 1-norm of a general rectangular matrix
 */
{
   const int incA = (Uplo == AtlasUpper ? (lda SHIFT) : ((lda+1)SHIFT));
   const int ioff = (Diag == AtlasNonUnit ? 1 : 0);
   int j;
   TYPE max=0.0, t0;

   if (Uplo == AtlasUpper)
   {
      for (j=0; j < N; j++)
      {
         t0 = Mjoin(PATL,asum)(j+ioff, A, 1);
         if (t0 != t0)
            return(t0);
         if (Diag == AtlasUnit) t0 += ATL_rone;
         if (t0 > max) max = t0;
         A += incA;
      }
   }
   else
   {
      if (Diag == AtlasUnit) A += 1 SHIFT;
      for (j=N; j; j--)
      {
         t0 = Mjoin(PATL,asum)(j+ioff-1, A, 1);
         if (t0 != t0)
            return(t0);
         if (Diag == AtlasUnit) t0 += ATL_rone;
         if (t0 > max) max = t0;
         A += incA;
      }
   }
   return(max);
}
