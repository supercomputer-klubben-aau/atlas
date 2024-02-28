/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#define FILLCONST -2560000000.0
void Mjoin(PATL,trgen)(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                       const int N, TYPE *A, const int lda0, const int seed)
{
   const int M = N SHIFT, lda = lda0 SHIFT;
   int i, j;

   dumb_seed(seed);
   Mjoin(PATL,gefillgap)(N, N, A, lda0);
   if (Uplo == AtlasUpper)
   {
      for (j=0; j != N; j++)
      {
         for (i=0; i != (j SHIFT); i++) A[i] = dumb_rand();
         if (Diag == AtlasNonUnit)
         {
            A[i++] = dumb_rand();
            #ifdef TCPLX
               A[i++] = dumb_rand();
            #endif
         }
         for (; i < M; i++) A[i] = FILLCONST;
         A += lda;
      }
   }
   else
   {
      for (j=0; j != N; j++)
      {
         for (i=0; i != (j SHIFT); i++) A[i] = FILLCONST;
         if (Diag == AtlasNonUnit)
         {
            A[i++] = dumb_rand();
            #ifdef TCPLX
               A[i++] = dumb_rand();
            #endif
         }
         for (; i != M; i++) A[i] = dumb_rand();
         A += lda;
      }
   }
}
