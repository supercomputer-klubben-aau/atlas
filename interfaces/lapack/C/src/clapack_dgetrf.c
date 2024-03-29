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

int clapack_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   double *A, const int lda, int *ipiv)
/*
 * Computes one of two LU factorizations based on the setting of the Order
 * parameter, as follows:
 * ----------------------------------------------------------------------------
 *                       Order == CblasColMajor
 * Column-major factorization of form
 *   A = P * L * U
 * where P is a row-permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapazoidal if M > N), and U is upper triangular
 * (upper trapazoidal if M < N).
 *
 * ----------------------------------------------------------------------------
 *                       Order == CblasRowMajor
 * Row-major factorization of form
 *   A = P * L * U
 * where P is a column-permutation matrix, L is lower triangular (lower
 * trapazoidal if M > N), and U is upper triangular with unit diagonals (upper
 * trapazoidal if M < N).
 *
 * ============================================================================
 * Let IERR be the return value of the function:
 *    If IERR == 0, successful exit.
 *    If (IERR < 0) the -IERR argument had an illegal value
 *    If (IERR > 0 && Order == CblasColMajor)
 *       U(i-1,i-1) is exactly zero.  The factorization has been completed,
 *       but the factor U is exactly singular, and division by zero will
 *       occur if it is used to solve a system of equations.
 *    If (IERR > 0 && Order == CblasRowMajor)
 *       L(i-1,i-1) is exactly zero.  The factorization has been completed,
 *       but the factor L is exactly singular, and division by zero will
 *       occur if it is used to solve a system of equations.
 */
{
   int ierr=0;
   if (Order != CblasRowMajor && Order != CblasColMajor)
   {
      ierr = -1;
      cblas_xerbla(1, "clapack_dgetrf",
                   "Order must be %d or %d, but is set to %d\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (M < 0)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_dgetrf",
                   "M cannot be less than zero 0,; is set to %d.\n", M);
   }
   if (N < 0)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_dgetrf",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (Order == CblasColMajor)
   {
      if (lda < M || lda < 1)
      {
         ierr = -6;
         cblas_xerbla(6, "clapack_dgetrf",
                      "lda must be >= MAX(M,1): lda=%d M=%d\n", lda, M);
      }
   }
   else
   {
      if (lda < N || lda < 1)
      {
         ierr = -6;
         cblas_xerbla(6, "clapack_dgetrf",
                      "lda must be >= MAX(N,1): lda=%d N=%d\n", lda, N);
      }
   }
   if (!ierr) ierr = ATL_getrf(Order, M, N, A, lda, ipiv);
   return(ierr);
}

