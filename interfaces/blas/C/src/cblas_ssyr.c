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

void cblas_ssyr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha,
                const float *X, const int incX, float *A, const int lda)
{
   int info = 2000;
   #define x X

#ifndef NoCblasErrorChecks
   if (Order != CblasColMajor && Order != CblasRowMajor)
      info = cblas_errprn(1, info, "Order must be %d or %d, but is set to %d",
                          CblasRowMajor, CblasColMajor, Order);
   if (Uplo != CblasUpper && Uplo != CblasLower)
      info = cblas_errprn(2, info, "UPLO must be %d or %d, but is set to %d",
                          CblasUpper, CblasLower, Uplo);
   if (N < 0) info = cblas_errprn(3, info,
                        "N cannot be less than zero; is set to %d.", N);
   if (!incX) info = cblas_errprn(6, info,
                                  "incX cannot be zero; is set to %d.", incX);
   if (lda < N || lda < 1)
      info = cblas_errprn(8, info, "lda must be >= MAX(N,1): lda=%d N=%d",
                          lda, N);
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_ssyr", "");
      return;
   }
#endif

   if (incX < 0) x += (1-N)*incX;

   if (Order == CblasColMajor)
      ATL_ssyr(Uplo, N, alpha, x, incX, A, lda);
   else
      ATL_ssyr(( (Uplo == CblasUpper) ? CblasLower : CblasUpper ),
                    N, alpha, x, incX, A, lda);
}
