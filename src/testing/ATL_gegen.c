/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#define FILLCONST -2560000000.0

void Mjoin(PATL,gefillgap)
   (const int M, const int N, TYPE *A, const int lda0)
{
   const int lda=lda0 SHIFT, n = lda0 - M;
   int j;

   if (n)
   {
      A += M SHIFT;
      for (j=0; j < N; j++, A += lda)
         Mjoin(PATLU,set)(n SHIFT, FILLCONST, A, 1);
   }
}

int Mjoin(PATL,gechkgap)
   (const int M0, const int N, TYPE *A, const int lda0)
{
   const int M = M0 SHIFT, lda=lda0 SHIFT, n = lda0 - M0;
   int i, j, OVERWRITES=0;
   if (n)
   {
      for (j=0; j < N; j++)
      {
         for (i=M; i < lda; i++)
         {
            if (A[j*lda+i] != FILLCONST)
            {
               fprintf(stderr, "   Overwrite in lda gap, A(%d,%d) = %f!!\n",
                       i, j, A[j*lda+i]);
               OVERWRITES++;
            }
         }
      }
   }
   return(OVERWRITES);
}

void Mjoin(PATL,gegen)
   (const int M0, const int N, TYPE *A, const int lda0, const int seed)
{
   const int M = M0 SHIFT, lda = lda0 SHIFT;
   int i, j;

   dumb_seed(seed);
   Mjoin(PATL,gefillgap)(M0, N, A, lda0);
   for (j=N; j; j--)
   {
      for (i=0; i != M; i++) A[i] = dumb_rand();
      A += lda;
   }
}
