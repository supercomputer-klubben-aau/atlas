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
   TYPE y0;
   ATL_INT j;

   for (j=0; j < N; j++, A += lda)
   {
      y0 = Mjoin(PATL,dot)(M, A, 1, X, 1);
      #ifdef BETA0
         *Y++ = y0;
      #else
         *Y++ += y0;
      #endif
   }
}
