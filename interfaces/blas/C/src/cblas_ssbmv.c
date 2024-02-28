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

void cblas_ssbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY)
{
   int info = 2000;
   const enum CBLAS_UPLO ruplo = (Uplo == CblasUpper) ? CblasLower : CblasUpper;
   #define x X
   #define y Y

#ifndef NoCblasErrorChecks
   if (Order != CblasColMajor && Order != CblasRowMajor)
      info = cblas_errprn(1, info, "Order must be %d or %d, but is set to %d",
                          CblasRowMajor, CblasColMajor, Order);
   if (Uplo != CblasUpper && Uplo != CblasLower)
      info = cblas_errprn(2, info,
                          "Uplo must be %d or %d, but is set to %d",
                          CblasUpper, CblasLower, Uplo);

   if (N < 0) info = cblas_errprn(3, info,
                        "N cannot be less than zero; is set to %d.", N);
   if (K < 0)
      info = cblas_errprn(4, info, "Valid K: 0 < K < N; K=%d, N=%d.", K, N);
   if (lda < K+1) info = cblas_errprn(7, info,
      "lda cannot be less than K+1;  K=%d, lda=%d\n", K, lda);
   if (!incX) info = cblas_errprn(9, info,
                                  "incX cannot be zero; is set to %d.", incX);
   if (!incY) info = cblas_errprn(12, info,
                                  "incY cannot be zero; is set to %d.", incY);
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_ssbmv", "");
      return;
   }
#endif

   if (incX < 0) x += (1-N)*incX;
   if (incY < 0) y += (1-N)*incY;
   if (Order == CblasColMajor)
      ATL_ssbmv(Uplo, N, K, alpha, A, lda, x, incX, beta, y, incY);
   else
      ATL_ssbmv(ruplo, N, K, alpha, A, lda, x, incX, beta, y, incY);
}
