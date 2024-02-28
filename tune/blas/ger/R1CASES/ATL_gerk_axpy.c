/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009, 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lvl2.h"
#include "atlas_prefetch.h"

void ATL_UGERK
   (ATL_CINT M, ATL_CINT N, const TYPE *X, const TYPE *Y, TYPE *A, ATL_CINT lda)
{
   const TYPE *stY = Y + N;
   if (M > 8)
   {
      do
      {
         Mjoin(PATL,axpy)(M, *Y, X, 1, A, 1);
         Y++;
         A += lda;
      }
      while (Y != stY);
   }
   else
      Mjoin(PATL,gerk_Mlt16)(M, N, ATL_rone, X, 1, Y, 1, A, lda);
}
