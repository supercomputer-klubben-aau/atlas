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

void cblas_zher(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha,
                const void *X, const int incX, void *A, const int lda)
{
   int info = 2000;
   void *vx;
   double one[2] = {ATL_rone, ATL_rzero};
   double *x0;
   const double *x=X;

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
      cblas_xerbla(info, "cblas_zher", "");
      return;
   }
#endif

   if (incX < 0) x += (1-N)*incX<<1;

   if (Order == CblasColMajor)
      ATL_zher(Uplo, N, alpha, x, incX, A, lda);
   else if (alpha != ATL_rzero)
   {
      vx = malloc(ATL_Cachelen + ATL_MulBySize(N));
      ATL_assert(vx);
      x0 = ATL_AlignPtr(vx);
      ATL_zmoveConj(N, one, x, incX, x0, 1);
      ATL_zher(( (Uplo == CblasUpper) ? CblasLower : CblasUpper ),
               N, alpha, x0, 1, A, lda);
      free(vx);
   }
   else
      ATL_zher(( (Uplo == CblasUpper) ? CblasLower : CblasUpper ),
               N, ATL_rzero, x, incX, A, lda);
}
