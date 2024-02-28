/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 Peter Soendergaard
 * Code contributers : Peter Soendergaard, R. Clint Whaley
 */
#define SCPLX
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"

int clapack_ctrtri(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const enum ATLAS_DIAG Diag, const int N,
                   void *A, const int lda)
{
   int ierr;
   if (Order != CblasRowMajor && Order != CblasColMajor)
   {
      ierr = -1;
      cblas_xerbla(1, "clapack_ctrtri",
                   "Order must be %d or %d, but is set to %d\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (Uplo != CblasUpper && Uplo != CblasLower)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_ctrtri",
                   "Uplo must be %d or %d, but is set to %d\n",
                   CblasUpper, CblasLower, Uplo);
   }
   if (Diag != CblasUnit && Diag != CblasNonUnit)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_ctrtri",
                   "Diag must be %d or %d, but is set to %d\n",
                   CblasNonUnit, CblasUnit, Diag);
   }
   if (N < 0)
   {
      ierr = -4;
      cblas_xerbla(4, "clapack_ctrtri",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (lda < N || lda < 1)
   {
      ierr = -6;
      cblas_xerbla(6, "clapack_ctrtri",
                   "lda must be >= MAX(N,1): lda=%d N=%d\n", lda, N);
   }
   if (ierr) ierr = ATL_ctrtri(Order, Uplo, Diag, N, A, lda);
   return(ierr);
}
