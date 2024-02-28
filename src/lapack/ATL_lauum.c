/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */
#include "atlas_lapack.h"
void ATL_lauum(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
               const int N, TYPE *A, const int lda)
{
   if (N > 0)
   {
      if (Order == AtlasColMajor)
      {
         if (Uplo == AtlasUpper) ATL_lauumCU(N, A, lda);
         else ATL_lauumCL(N, A, lda);
      }
      else
      {
         if (Uplo == AtlasUpper) ATL_lauumRU(N, A, lda);
         else ATL_lauumRL(N, A, lda);
      }
   }
}
