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

int clapack_zgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                  void *A, const int lda, int *ipiv,
                  void *B, const int ldb)
/*
 * clapack_zgesv computes the solution to a system of linear equations
 *   A * X = B,
 * where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
 * The LU factorization used to factor A is dependent on the Order parameter,
 * as detailed in the leading comments of clapack_zgetrf.
 * The factored form of A is then used solve the system of equations A * X = B.
 * A is overwritten with the appropriate LU factorization, and B, which
 * contains B on input, is overwritten with the solution X on output.
 */
{
   int ierr = 0;

   if (Order != CblasRowMajor && Order != CblasColMajor)
   {
      ierr = -1;
      cblas_xerbla(1, "clapack_zgesv",
                   "Order must be %d or %d, but is set to %d.\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (N < 0)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_zgesv",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (NRHS < 0)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_zgesv",
                   "NRHS cannot be less than zero 0,; is set to %d.\n", NRHS);
   }
   if (lda < N || lda < 1)
   {
      ierr = -5;
      cblas_xerbla(5, "clapack_zgesv",
                   "lda must be >= MAX(N,1): lda=%d N=%d\n", lda, N);
   }
   if (ldb < N || ldb < 1)
   {
      ierr = -8;
      cblas_xerbla(8, "clapack_zgesv",
                   "ldb must be >= MAX(N,1): ldb=%d N=%d\n", ldb, N);

   }
   if (!ierr) ierr = ATL_zgetrf(Order, N, N, A, lda, ipiv);
   if (!ierr) ATL_zgetrs(Order, CblasNoTrans, N, NRHS, A, lda, ipiv, B, ldb);
   return(ierr);
}
