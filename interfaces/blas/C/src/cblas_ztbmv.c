/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define DCPLX
#include "atlas_misc.h"
#include "cblas.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias2.h"
#endif
#include "atlas_level2.h"

void cblas_ztbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX)
{
   int info = 2000;
   enum CBLAS_UPLO uplo;
   enum CBLAS_TRANSPOSE ta;
   double *x = X;

#ifndef NoCblasErrorChecks
   if (Order != CblasColMajor && Order != CblasRowMajor)
      info = cblas_errprn(1, info, "Order must be %d or %d, but is set to %d",
                          CblasRowMajor, CblasColMajor, Order);
   if (Uplo != CblasUpper && Uplo != CblasLower)
      info = cblas_errprn(2, info, "UPLO must be %d or %d, but is set to %d",
                          CblasUpper, CblasLower, Uplo);
   if (TA != CblasNoTrans && TA != CblasTrans && TA != CblasConjTrans)
      info = cblas_errprn(3, info,
                          "TransA must be %d, %d or %d, but is set to %d",
                          CblasNoTrans, CblasTrans, CblasConjTrans, TA);
   if (Diag != CblasUnit && Diag != CblasNonUnit)
      info = cblas_errprn(4, info, "DIAG must be %d or %d, but is set to %d",
                          CblasUnit, CblasNonUnit, Diag);

   if (N < 0) info = cblas_errprn(5, info,
                        "N cannot be less than zero; is set to %d.", N);
   if (K < 0)
      info = cblas_errprn(6, info, "Valid K: 0 < K < N; K=%d, N=%d.", K, N);
   if (lda < K+1)
      info = cblas_errprn(8, info, "lda must be >= K+1: lda=%d K=%d", lda, K);
   if (!incX)
      info = cblas_errprn(10, info, "incX cannot be zero; is set to %d.", incX);
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_ztbmv", "");
      return;
   }
#endif
   if (incX < 0) x += (1-N)*incX<<1;
   if (Order == CblasColMajor)
      ATL_ztbmv(Uplo, TA, Diag, N, K, A, lda, x, incX);
   else
   {
      uplo = ( (Uplo == CblasUpper) ? CblasLower : CblasUpper );
      if (TA == CblasNoTrans) ta = CblasTrans;
      else if (TA == CblasConjTrans) ta = AtlasConj;
      else ta = CblasNoTrans;
      ATL_ztbmv(uplo, ta, Diag, N, K, A, lda, x, incX);
   }
}
