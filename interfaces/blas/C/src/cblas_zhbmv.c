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

void cblas_zhbmv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY)
{
   int info = 2000;
   const enum CBLAS_UPLO ruplo = (Uplo == CblasUpper) ? CblasLower : CblasUpper;
   void *vx;
   double *X0, *x = (double*) X;
   double *y = Y;
   const double *alp=alpha;
   const double *bet=beta;
   double calpha[2], cbeta[2];
   const double one[2] = {ATL_rone, ATL_rzero};
   calpha[0] = *alp;
   calpha[1] = -alp[1];
   cbeta[0] = *bet;
   cbeta[1] = -bet[1];

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
      cblas_xerbla(info, "cblas_zhbmv", "");
      return;
   }
#endif

   if (incX < 0) x += (1-N)*incX<<1;
   if (incY < 0) y += (1-N)*incY<<1;
   if (Order == CblasColMajor)
      ATL_zhbmv(Uplo, N, K, alpha, A, lda, x, incX, beta, y, incY);
   else
   {
      vx = malloc(ATL_Cachelen + 2*N*sizeof(double));
      ATL_assert(vx);
      X0 = x;
      x = ATL_AlignPtr(vx);
      ATL_zmoveConj(N, calpha, X0, incX, x, 1);
      if (*bet != ATL_rzero || bet[1] != ATL_rzero)
      {
         ATL_zscalConj(N, cbeta, y, incY);
         ATL_zhbmv(ruplo, N, K, one, A, lda, x, 1, one, y, incY);
      }
      else ATL_zhbmv(ruplo, N, K, one, A, lda, x, 1, beta, y, incY);
      free(vx);
      ATL_zscalConj(N, one, y, incY);
   }
}
