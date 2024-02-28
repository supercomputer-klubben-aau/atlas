/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define DREAL
#include "atlas_misc.h"
#include "cblas.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias2.h"
#endif
#include "atlas_level2.h"

void cblas_dsyr2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha,
                 const double *X, const int incX,
                 const double *Y, const int incY, double *A, const int lda)
{
   int info = 2000;
   #define x X
   #define y Y

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
   if (!incY) info = cblas_errprn(8, info,
                                  "incY cannot be zero; is set to %d.", incY);
   if (lda < N || lda < 1)
      info = cblas_errprn(10, info, "lda must be >= MAX(N,1): lda=%d N=%d",
                          lda, N);
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_dsyr2", "");
      return;
   }
#endif

   if (incX < 0) x += (1-N)*incX;
   if (incY < 0) y += (1-N)*incY;

   if (Order == CblasColMajor)
      ATL_dsyr2(Uplo, N, alpha, x, incX, y, incY, A, lda);
   else
      ATL_dsyr2(( (Uplo == CblasUpper) ? CblasLower : CblasUpper ),
                     N, alpha, y, incY, x, incX, A, lda);
}
