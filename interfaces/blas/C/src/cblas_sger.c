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

void cblas_sger (const enum CBLAS_ORDER Order, const int M, const int N,
                 const float alpha, const float *X, const int incX,
                 const float *Y, const int incY, float *A, const int lda)
{
   int info = 2000;
   #define x X
   #define y Y

#ifndef NoCblasErrorChecks
   if (M < 0) info = cblas_errprn(2, info,
                        "M cannot be less than zero; is set to %d.", M);
   if (N < 0) info = cblas_errprn(3, info,
                        "N cannot be less than zero; is set to %d.", N);
   if (!incX) info = cblas_errprn(6, info,
                                  "incX cannot be zero; is set to %d.", incX);
   if (!incY) info = cblas_errprn(8, info,
                                  "incY cannot be zero; is set to %d.", incY);
   if (Order == CblasColMajor)
   {
      if (lda < M || lda < 1)
         info = cblas_errprn(10, info, "lda must be >= MAX(M,1): lda=%d M=%d",
                             lda, M);
   }
   else if (Order == CblasRowMajor)
   {
      if (lda < N || lda < 1)
         info = cblas_errprn(10, info, "lda must be >= MAX(N,1): lda=%d M=%d",
                             lda, N);
   }
   else
      info = cblas_errprn(1, info, "Order must be %d or %d, but is set to %d",
                          CblasRowMajor, CblasColMajor, Order);
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_sger", "");
      return;
   }
#endif

   if (incX < 0) x += (1-M)*incX;
   if (incY < 0) y += (1-N)*incY;

   if (Order == CblasColMajor)
      ATL_sger(M, N, alpha, x, incX, y, incY, A, lda);
   else
      ATL_sger(N, M, alpha, y, incY, x, incX, A, lda);
}
