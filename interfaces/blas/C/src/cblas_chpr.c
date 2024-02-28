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

void cblas_chpr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha,
                const void *X, const int incX, void *A)
{
   int info = 2000;
   void *vx;
   float one[2] = {ATL_rone, ATL_rzero};
   float *x0;
   const float *x=X;

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
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_chpr", "");
      return;
   }
#endif

   if (incX < 0) x += (1-N)*incX<<1;

   if (Order == CblasColMajor)
      ATL_chpr(Uplo, N, alpha, x, incX, A);
   else if (alpha != ATL_rzero)
   {
      vx = malloc(ATL_Cachelen + ATL_MulBySize(N));
      ATL_assert(vx);
      x0 = ATL_AlignPtr(vx);
      ATL_cmoveConj(N, one, x, incX, x0, 1);
      ATL_chpr(( (Uplo == CblasUpper) ? CblasLower : CblasUpper ),
               N, alpha, x0, 1, A);
      free(vx);
   }
   else
      ATL_chpr(( (Uplo == CblasUpper) ? CblasLower : CblasUpper ),
               N, ATL_rzero, x, incX, A);
}
