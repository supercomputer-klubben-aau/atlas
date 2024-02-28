/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_lapack.h"

int ATL_potrf(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
              const int N, TYPE *A, const int lda0)
{
   int ierr=0;
   size_t lda = lda0;
   if (N)
   {
      if (Order == AtlasColMajor)
      {
         if (Uplo == AtlasUpper) ierr = ATL_potrfU(N, A, lda);
         else ierr = ATL_potrfL(N, A, lda);
      }
      else
      {
      #ifdef TREAL
         if (Uplo == AtlasUpper) ierr = ATL_potrfL(N, A, lda);
         else ierr = ATL_potrfU(N, A, lda);
      #else
         if (Uplo == AtlasUpper) ierr = Mjoin(PATL,potrfRU)(N, A, lda);
         else ierr = Mjoin(PATL,potrfRL)(N, A, lda);
      #endif
      }
   }
   return(ierr);
}
