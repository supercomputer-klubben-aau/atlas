/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"

void Mjoin(PATL,geprint)
   (char *mat, const int M, const int N, const TYPE *A, const int lda0)
{
   const int lda = lda0 SHIFT;
   int i, j;

   printf("\n%s = \n",mat);
   for (i=0; i != M; i++)
   {
      #ifdef TREAL
         for (j=0; j != N; j++) printf("%f  ",A[i+j*lda]);
      #else
         for (j=0; j != N; j++)
            printf("(%f,%f)  ",A[2*i+j*lda], A[1+2*i+j*lda]);
      #endif
      printf("\n");
   }
}
