/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_lapack.h"

#ifdef TCPLX
   #define MyTrans CblasConjTrans
#else
   #define MyTrans CblasTrans
#endif

void ATL_potrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, const int NRHS, const TYPE *A, const int lda,
               TYPE *B, const int ldb)
{
   #ifdef TCPLX
      int j;
      const int ldb2 = ldb+ldb;
      const TYPE one[2] = {ATL_rone, ATL_rzero};
   #else
      #define one ATL_rone
   #endif

   if (!N || !NRHS) return;
   if (Order == CblasColMajor)
   {
/*
 *    Solve X = inv(U) * inv(U') * B
 */
      if (Uplo == AtlasUpper)
      {
         cblas_trsm(Order, CblasLeft, CblasUpper, MyTrans, CblasNonUnit,
                    N, NRHS, one, A, lda, B, ldb);
         cblas_trsm(Order, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
                    N, NRHS, one, A, lda, B, ldb);
      }
/*
 *    Solve X = inv(L') * inv(L) * B
 */
      else
      {
         cblas_trsm(Order, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                    N, NRHS, one, A, lda, B, ldb);
         cblas_trsm(Order, CblasLeft, CblasLower, MyTrans, CblasNonUnit,
                    N, NRHS, one, A, lda, B, ldb);
      }
   }
/*
 * For row-major, remember we have x' and b', so we must transpose usual
 * equations
 */
   else
   {
      #ifdef TCPLX
         for (j=0; j < NRHS; j++)
            Mjoin(PATLU,scal)(N, -1.0, B+j*ldb2+1, 2);
      #endif
/*
 *    solve x^T = b^T * inv(U) * inv(U^T)
 *    conj( x^H = b^H * inv(U) * inv(U^H) )  (complex)
 */
      if (Uplo == CblasUpper)
      {
         cblas_trsm(Order, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
                    NRHS, N, one, A, lda, B, ldb);
         cblas_trsm(Order, CblasRight, CblasUpper, MyTrans, CblasNonUnit,
                    NRHS, N, one, A, lda, B, ldb);
      }
/*
 *    solve x^T = b^T * inv(L^T) * inv(L)
 *    conj( x^H = b^H * inv(L^H) * inv(L) )  (complex)
 */
      else
      {
         cblas_trsm(Order, CblasRight, CblasLower, MyTrans, CblasNonUnit,
                    NRHS, N, one, A, lda, B, ldb);
         cblas_trsm(Order, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit,
                    NRHS, N, one, A, lda, B, ldb);
      }
      #ifdef TCPLX
         for (j=0; j < NRHS; j++)
            Mjoin(PATLU,scal)(N, -1.0, B+j*ldb2+1, 2);
      #endif
   }
}
