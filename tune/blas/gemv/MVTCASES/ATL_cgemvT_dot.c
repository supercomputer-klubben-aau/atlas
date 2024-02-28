/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"

void ATL_UGEMV(ATL_CINT M, ATL_CINT N, const TYPE *A, ATL_CINT lda,
               const TYPE *X, TYPE *Y)
/*
 *  y = [0,1]*y + A*x, A is MxN, storing the transpose of the matrix
 */
{
   TYPE ry, iy;
   ATL_CINT lda2 = lda+lda;
   ATL_INT j;

   for (j=0; j < N; j++, A += lda2, Y += 2)
   {
      #ifdef BETA0
         Mjoin(PATL,dotu_sub)(M, A, 1, X, 1, Y);
      #else
         #ifdef __clang__  /* workaround for clang error */
            TYPE dot[2];
            Mjoin(PATL,dotu_sub)(M, A, 1, X, 1, dot);
            *Y += *dot;
            Y[1] += dot[1];
         #else
            ry = *Y; iy = Y[1];
            Mjoin(PATL,dotu_sub)(M, A, 1, X, 1, Y);
            *Y += ry;
            Y[1] += iy;
         #endif
      #endif
   }
}
