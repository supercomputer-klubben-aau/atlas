/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010, 2009, 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lvl2.h"

void ATL_UGERK(ATL_CINT M, ATL_CINT N, const TYPE *X, const TYPE *Y,
               TYPE *A, ATL_CINT lda)
{
   ATL_CINT lda2 = lda<<1;
   const TYPE *stY = Y + N+N;
   TYPE y[2];
   do
   {
      *y = *Y;
      #ifdef Conj_
         y[1] = -Y[1];
      #else
         y[1] = Y[1];
      #endif
      Mjoin(PATL,axpy)(M, y, X, 1, A, 1);
      Y += 2;
      A += lda2;
   }
   while (Y != stY);
}
