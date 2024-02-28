/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#define SREAL
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"

int clapack_sgetrs
   (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const float *A, const int lda,
    const int *ipiv, float *B, const int ldb)
{
   int ierr=0;
   if (Order != CblasRowMajor && Order != CblasColMajor)
   {
      ierr = -1;
      cblas_xerbla(1, "clapack_sgetrs",
                   "Order must be %d or %d, but is set to %d\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (Trans != CblasNoTrans && Trans != CblasTrans && Trans != CblasConjTrans)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_sgetrs",
                   "Trans must be %d, %d, or %d, but is set to %d\n",
                   CblasNoTrans, CblasTrans, CblasConjTrans);
   }
   if (N < 0)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_sgetrs",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (NRHS < 0)
   {
      ierr = -4;
      cblas_xerbla(4, "clapack_sgetrs",
                   "NRHS cannot be less than zero 0,; is set to %d.\n", NRHS);
   }
   if (lda < N || lda < 1)
   {
      ierr = -6;
      cblas_xerbla(6, "clapack_sgetrs",
                   "lda must be >= MAX(N,1): lda=%d N=%d\n", lda, N);
   }
   if (ldb < N || ldb < 1)
   {
      ierr = -9;
      cblas_xerbla(9, "clapack_sgetrs",
                   "ldb must be >= MAX(N,1): lda=%d N=%d\n", lda, N);
   }
   if (!ierr) ATL_sgetrs(Order, Trans, N, NRHS, A, lda, ipiv, B, ldb);
   return(ierr);
}
