/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define SCPLX
#include "atlas_misc.h"
#include "cblas.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias2.h"
#endif
#include "atlas_level2.h"

void cblas_chpr2(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha,
                 const void *X, const int incX,
                 const void *Y, const int incY, void *A)
{
   int info = 2000;
   void *vx, *vy;
   float *x0, *y0;
   const float *x=X, *y=Y, *alp=alpha;
   const float one[2]={ATL_rone, ATL_rzero};

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
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_chpr2", "");
      return;
   }
#endif

   if (incX < 0) x += (1-N)*incX<<1;
   if (incY < 0) y += (1-N)*incY<<1;

   if (Order == CblasColMajor)
      ATL_chpr2(Uplo, N, alpha, x, incX, y, incY, A);
   else if (alp[0] != ATL_rzero || alp[1] != ATL_rzero)
   {
      vx = malloc(ATL_Cachelen + ATL_MulBySize(N));
      vy = malloc(ATL_Cachelen + ATL_MulBySize(N));
      ATL_assert(vx != NULL && vy != NULL);
      x0 = ATL_AlignPtr(vx);
      y0 = ATL_AlignPtr(vy);
      ATL_cmoveConj(N, alpha, y, incY, y0, 1);
      ATL_ccopyConj(N, x, incX, x0, 1);
      ATL_chpr2(( (Uplo == CblasUpper) ? CblasLower : CblasUpper ),
                N, one, y0, 1, x0, 1, A);
      free(vx);
      free(vy);
   }
   else ATL_chpr2(( (Uplo == CblasUpper) ? CblasLower : CblasUpper ),
                  N, alpha, y, incY, x, incX, A);
}
