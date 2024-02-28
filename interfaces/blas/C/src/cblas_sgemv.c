/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define SREAL
#include "atlas_misc.h"
#include "cblas.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias2.h"
#endif
#include "atlas_level2.h"

void cblas_sgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TA,
                 const int M, const int N, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY)
{
   int info = 2000;
   #define x X
   #define y Y

#ifndef NoCblasErrorChecks
   if (TA != CblasNoTrans && TA != CblasTrans && TA != CblasConjTrans)
      info = cblas_errprn(2, info,
                          "TransA must be %d, %d or %d, but is set to %d",
                          CblasNoTrans, CblasTrans, CblasConjTrans, TA);

   if (M < 0) info = cblas_errprn(3, info,
                        "M cannot be less than zero; is set to %d.", M);
   if (N < 0) info = cblas_errprn(4, info,
                        "N cannot be less than zero; is set to %d.", N);
   if (!incX) info = cblas_errprn(9, info,
                                  "incX cannot be zero; is set to %d.", incX);
   if (!incY) info = cblas_errprn(12, info,
                                  "incY cannot be zero; is set to %d.", incY);
   if (Order == CblasColMajor)
   {
      if (lda < M || lda < 1)
         info = cblas_errprn(7, info, "lda must be >= MAX(M,1): lda=%d M=%d",
                             lda, M);
   }
   else if (Order == CblasRowMajor)
   {
      if (lda < N || lda < 1)
         info = cblas_errprn(7, info, "lda must be >= MAX(N,1): lda=%d N=%d",
                             lda, N);
   }
   else
      info = cblas_errprn(1, info, "Order must be %d or %d, but is set to %d",
                          CblasRowMajor, CblasColMajor, Order);
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_sgemv", "");
      return;
   }
#endif
   if (TA == AtlasNoTrans)
   {
      if (incX < 0) x += (1-N)*incX;
      if (incY < 0) y += (1-M)*incY;
   }
   else
   {
      if (incX < 0) x += (1-M)*incX;
      if (incY < 0) y += (1-N)*incY;
   }
   if (Order == CblasColMajor)
      ATL_sgemv(TA, M, N, alpha, A, lda, x, incX, beta, y, incY);
   else
   {
      if (TA == CblasNoTrans)
         ATL_sgemv(CblasTrans, N, M, alpha, A, lda, x, incX, beta, y, incY);
      else
         ATL_sgemv(CblasNoTrans, N, M, alpha, A, lda, x, incX, beta, y, incY);
   }
}
