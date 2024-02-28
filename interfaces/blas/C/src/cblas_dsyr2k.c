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


void cblas_dsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double  alpha, const double *A, const int lda,
                  const double *B, const int ldb, const double  beta,
                  double *C, const int ldc)
{
   enum CBLAS_UPLO uplo;
   enum CBLAS_TRANSPOSE trans;

#ifndef NoCblasErrorChecks
   int info = 2000;
   if (Uplo != CblasUpper && Uplo != CblasLower)
      info = cblas_errprn(2, info, "UPLO must be %d or %d, but is set to %d",
                               CblasUpper, CblasLower, Uplo);
   if (N < 0) info = cblas_errprn(4, info,
                     "N cannot be less than zero; it is set to %d.", N);
   if (K < 0) info = cblas_errprn(5, info,
                     "K cannot be less than zero; it is set to %d.", K);
   if (Order == CblasColMajor)
   {
      if (Trans == AtlasNoTrans)
      {
         if ( (lda < N) || (lda < 1) )
            info = cblas_errprn(8, info, "lda must be >= MAX(N,1): lda=%d N=%d",
                                lda, N);
         if ( (ldb < N) || (ldb < 1) )
            info = cblas_errprn(10, info,"ldb must be >= MAX(N,1): ldb=%d N=%d",
                                ldb, N);
      }
      else
      {
         if (Trans != AtlasTrans && Trans != AtlasConjTrans)
            info = cblas_errprn(3, info,
                                "Trans must be %d or %d, but is set to %d",
                                CblasNoTrans, CblasTrans, Trans);
         if ( (lda < K) || (lda < 1) )
            info = cblas_errprn(8, info, "lda must be >= MAX(K,1): lda=%d K=%d",
                                lda, K);
         if ( (ldb < K) || (ldb < 1) )
            info = cblas_errprn(10, info,"ldb must be >= MAX(K,1): ldb=%d K=%d",
                                ldb, K);
      }
   }
   else if (Order == CblasRowMajor)
   {
      if (Trans == AtlasNoTrans)
      {
         if ( (lda < K) || (lda < 1) )
            info = cblas_errprn(8, info, "lda must be >= MAX(K,1): lda=%d K=%d",
                                lda, K);
         if ( (ldb < K) || (ldb < 1) )
            info = cblas_errprn(10, info,"ldb must be >= MAX(K,1): ldb=%d K=%d",
                                ldb, K);
      }
      else
      {
         if (Trans != AtlasTrans && Trans != AtlasConjTrans)
            info = cblas_errprn(3, info,
                                "Trans must be %d or %d, but is set to %d",
                                CblasNoTrans, CblasTrans, Trans);
         if ( (lda < N) || (lda < 1) )
            info = cblas_errprn(8, info, "lda must be >= MAX(N,1): lda=%d N=%d",
                                lda, N);
         if ( (ldb < N) || (ldb < 1) )
            info = cblas_errprn(10, info,"ldb must be >= MAX(N,1): ldb=%d N=%d",
                                ldb, N);
      }
   }
   else info = cblas_errprn(1, info, "Order must be %d or %d, but is set to %d",
                            CblasRowMajor, CblasColMajor, Order);

   if ( (ldc < N) || (ldc < 1) )
      info = cblas_errprn(13, info, "ldc must be >= MAX(N,1): ldc=%d N=%d",
                          ldc, N);
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_dsyr2k", "");
      return;
   }
#endif

   if (Order == CblasColMajor)
      ATL_dsyr2k(Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   else
   {
      if (Uplo == CblasUpper) uplo = CblasLower;
      else uplo = CblasUpper;
      if (Trans == CblasNoTrans) trans = CblasTrans;
      else trans = CblasNoTrans;
      ATL_dsyr2k(uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
}
