/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#define DREAL
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"

int clapack_dpotrf(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, double *A, const int lda)
{
   int ierr=0;
   if (Order != CblasRowMajor && Order != CblasColMajor)
   {
      ierr = -1;
      cblas_xerbla(1, "clapack_dpotrf",
                   "Order must be %d or %d, but is set to %d\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (Uplo != CblasUpper && Uplo != CblasLower)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_dpotrf",
                   "Uplo must be %d or %d, but is set to %d\n",
                   CblasUpper, CblasLower, Uplo);
   }
   if (N < 0)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_dpotrf",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (lda < N || lda < 1)
   {
      ierr = -5;
      cblas_xerbla(5, "clapack_dpotrf",
                   "lda must be >= MAX(N,1): lda=%d N=%d\n", lda, N);
   }
   if (!ierr) ierr = ATL_dpotrf(Order, Uplo, N, A, lda);
   return(ierr);
}
