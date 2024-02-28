/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#define DCPLX
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"

int clapack_zposv(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                  const int N, const int NRHS, void *A, const int lda,
                  void *B, const int ldb)
{
   int ierr = 0;

   if (Order != CblasRowMajor && Order != CblasColMajor)
   {
      ierr = -1;
      cblas_xerbla(1, "clapack_zposv",
                   "Order must be %d or %d, but is set to %d\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (Uplo != CblasUpper && Uplo != CblasLower)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_zposv",
                   "Uplo must be %d or %d, but is set to %d\n",
                   CblasUpper, CblasLower, Uplo);
   }
   if (N < 0)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_zposv",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (NRHS < 0)
   {
      ierr = -4;
      cblas_xerbla(4, "clapack_zposv",
                   "NRHS cannot be less than zero 0,; is set to %d.\n", NRHS);
   }
   if (lda < N || lda < 1)
   {
      ierr = -6;
      cblas_xerbla(6, "clapack_zposv",
                   "lda must be >= MAX(N,1): lda=%d N=%d\n", lda, N);
   }
   if (ldb < N || ldb < 1)
   {
      ierr = -8;
      cblas_xerbla(8, "clapack_zposv",
                   "ldb must be >= MAX(N,1): ldb=%d N=%d\n", ldb, N);
   }
   if (!ierr) ierr = ATL_zpotrf(Order, Uplo, N, A, lda);
   if (!ierr) ATL_zpotrs(Order, Uplo, N, NRHS, A, lda, B, ldb);
   return(ierr);
}
