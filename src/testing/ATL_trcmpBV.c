/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bitvec.h"
/*
 * Returns a bitvec where set bits are errors above tolerance tol,
 * For Upper, only diagonal & above are checked
 * For Lower, only diagonal & below are checked
 * Entries are ordered row-major for ease of printing
 * Combine with ATL_print2dBV for pictoral error report.
 */
void *Mjoin(PATL,trcmpBV)
   (int verb, double tol, const enum ATLAS_UPLO Uplo, int M, int N,
    const TYPE *A, int lda, const TYPE *B, int ldb)
{
   ATL_BV_t *bv;
   int i, j;
   size_t lda2=lda SHIFT, ldb2=ldb SHIFT;

   bv = ATL_NewBV(M*N);
   if (tol < 0.0)
      tol = -tol;
   for (j=0; j < N; j++, A += lda2, B += ldb2)
   {
      const int MM=(Uplo == AtlasLower) ? M:j;
      for (i=(Uplo==AtlasLower)?j:0; i < MM; i++)
      {
         #ifdef TCPLX
            const int I = i+i;
            double diff = A[I] - B[I], idiff = A[I+1] - B[I+1];
            if (diff < 0.0)
               diff = -diff;
            if (idiff < 0.0)
               idiff = 0.0;
            if (diff > tol || idiff > tol)
            {
               if (verb > 1)
                  printf("A(%d,%d)=[%e,%e];  expected=[%e,%e]\n", i, j,
                         B[I], B[I+1], A[I], A[I+1]);
               ATL_SetBitBV(bv, i*N+j);
            }
         #else
            double diff = A[i] - B[i];
            if (diff < 0.0)
               diff = -diff;

            if (diff > tol)
            {
               if (verb > 1)
                  printf("A(%d,%d)=%e;  expected=%e, diff=%e\n",
                         i, j, B[i], A[i], diff);
               ATL_SetBitBV(bv, i*N+j);
            }
         #endif

      }
   }
   return(bv);
}
