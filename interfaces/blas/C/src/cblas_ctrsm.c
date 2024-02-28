/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias3.h"
#endif
#include "atlas_level3.h"
#include "cblas.h"


void cblas_ctrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void * alpha, const void *A, const int lda,
                 void *B, const int ldb)
{
   enum CBLAS_SIDE side;
   enum CBLAS_UPLO uplo;
   int info=2000;

#ifndef NoCblasErrorChecks
   if (Order == CblasColMajor)
   {
      if (Side == CblasLeft)
      {
         if ( (lda < M) || (lda < 1) )
            info = cblas_errprn(10, info,"lda must be >= MAX(M,1): lda=%d M=%d",
                                lda, M);
      }
      else if (Side == CblasRight)
      {
         if ( (lda < N) || (lda < 1) )
            info = cblas_errprn(10, info,"lda must be >= MAX(N,1): lda=%d N=%d",
                                lda, N);
      }
      else info = cblas_errprn(2, info,
                               "SIDE must be %d or %d, but is set to %d",
                               CblasRight, CblasLeft, Side);
      if ( (ldb < M) || (ldb < 1) )
         info = cblas_errprn(12, info, "ldb must be >= MAX(M,1): ldb=%d M=%d",
                             ldb, M);
   }
   else if (Order == CblasRowMajor)
   {
      if (Side == CblasLeft)
      {
         if ( (lda < M) || (lda < 1) )
            info = cblas_errprn(10, info,"lda must be >= MAX(M,1): lda=%d M=%d",
                                lda, M);
      }
      else if (Side == CblasRight)
      {
         if ( (lda < N) || (lda < 1) )
            info = cblas_errprn(10, info,"lda must be >= MAX(N,1): lda=%d N=%d",
                                lda, N);
      }
      else info = cblas_errprn(2, info,
                               "SIDE must be %d or %d, but is set to %d",
                               CblasRight, CblasLeft, Side);
      if ( (ldb < N) || (ldb < 1) )
         info = cblas_errprn(12, info, "ldb must be >= MAX(N,1): ldb=%d N=%d",
                             ldb, N);
   }
   else
      info = cblas_errprn(1, info, "Order must be %d or %d, but is set to %d",
                          CblasRowMajor, CblasColMajor, Order);

   if (Uplo != CblasUpper && Uplo != CblasLower)
      info = cblas_errprn(3, info, "UPLO must be %d or %d, but is set to %d",
                          CblasUpper, CblasLower, Uplo);

   if (TA != AtlasNoTrans && TA != AtlasTrans && TA != AtlasConjTrans)
      info = cblas_errprn(4, info,
                          "TransA must be %d, %d or %d, but is set to %d",
                          CblasNoTrans, CblasTrans, CblasConjTrans, TA);

   if (Diag != CblasUnit && Diag != CblasNonUnit)
      info = cblas_errprn(5, info, "UPLO must be %d or %d, but is set to %d",
                          CblasUnit, CblasNonUnit, Diag);

   if (M < 0) info = cblas_errprn(6, info,
                     "M cannot be less than zero; it is set to %d.", M);
   if (N < 0) info = cblas_errprn(7, info,
                     "N cannot be less than zero; it is set to %d.", N);


   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_ctrsm", "");
      return;
   }
#endif

   if (Order == CblasColMajor)
      ATL_ctrsm(Side, Uplo, TA, Diag, M, N, alpha, A, lda, B, ldb);
   else
   {
      if (Side == CblasLeft) side = CblasRight;
      else side = CblasLeft;
      if (Uplo == CblasUpper) uplo = CblasLower;
      else uplo = CblasUpper;
      ATL_ctrsm(side, uplo, TA, Diag, N, M, alpha, A, lda, B, ldb);
   }
}
