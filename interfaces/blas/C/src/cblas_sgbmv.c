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

void cblas_sgbmv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TA,
                 const int M, const int N, const int KL, const int KU,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX,
                 const float beta, float *Y, const int incY)
{
   int info = 2000;
   #define x X
   #define y Y

#ifndef NoCblasErrorChecks
   if (Order != CblasRowMajor && Order != CblasColMajor)
      info = cblas_errprn(1, info, "Order must be %d or %d, but is set to %d",
                          CblasRowMajor, CblasColMajor, Order);
   if (TA != CblasNoTrans && TA != CblasTrans && TA != CblasConjTrans)
      info = cblas_errprn(2, info,
                          "TransA must be %d, %d or %d, but is set to %d",
                          CblasNoTrans, CblasTrans, CblasConjTrans, TA);
   if (M < 0) info = cblas_errprn(3, info,
                        "M cannot be less than zero; is set to %d.", M);
   if (N < 0) info = cblas_errprn(4, info,
                        "N cannot be less than zero; is set to %d.", N);
   if (KL < 0) info = cblas_errprn(5, info,
                         "KL cannot be less than zero; is set to %d.", KL);
   if (KU < 0) info = cblas_errprn(6, info,
                         "KU cannot be less than zero; is set to %d.", KU);
   if (lda < (KL+KU+1))
      info = cblas_errprn(9, info, "lda must be >= KU+KL+1: lda=%d KU+KL+1=%d",
                          lda, KU+KL+1);
   if (!incX) info = cblas_errprn(11, info,
                                  "incX cannot be zero; is set to %d.", incX);
   if (!incY) info = cblas_errprn(14, info,
                                  "incY cannot be zero; is set to %d.", incY);
   if (info != 2000)
   {
      cblas_xerbla(info, "cblas_sgbmv", "");
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
      ATL_sgbmv(TA, M, N, KL, KU, alpha, A, lda, x, incX, beta, y, incY);
   else
   {
      if (TA == CblasNoTrans)
         ATL_sgbmv(CblasTrans, N, M, KU, KL, alpha, A, lda, x, incX,
                   beta, y, incY);
      else ATL_sgbmv(CblasNoTrans, N, M, KU, KL, alpha, A, lda, x, incX,
                     beta, y, incY);
   }
}
