/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lapack.h"

void ATL_getrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
               const int N, const int NRHS, const TYPE *A, const int lda,
               const int *ipiv, TYPE *B, const int ldb)
/*
 * OK, this pivoting crap is tricky.  The trick is, when we pivot columns
 * of the matrix, this effects X but not B, and when we pivot rows, this
 * effects B, but not X.  So, must never attempt to apply a Pr
 * (row permutation matrix) to X or a Pc to B.
 */
{
   enum CBLAS_DIAG Lunit, Uunit;
   #ifdef TREAL
      #define one ATL_rone
   #else
      const TYPE one[2] = {ATL_rone, ATL_rzero};
   #endif

   if (!N || !NRHS) return;

   if (Order == CblasColMajor)
   {
/*
 *    A*X = B.  Since we have pivoted A by Pr (PA=LU), we pivot B by Pr,
 *    **and this does not effect X at all**, so we solve
 *    X = inv(U)*inv(L)*(Pr * B)
 */
      if (Trans == CblasNoTrans)
      {
         ATL_laswp(NRHS, B, ldb, 0, N, ipiv, 1);
         cblas_trsm(Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                    N, NRHS, one, A, lda, B, ldb);
         cblas_trsm(Order, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
                    N, NRHS, one, A, lda, B, ldb);
      }
/*
 *    trans(L*U = PA)  ==>  U' L' = A' P, so P is Pc, and does not effect B,
 *    U' L' Pc X = B  ==> Pc X = inv(L') * inv(U') * B, but we want
 *    X, not Pc X, so we apply inv(Pc) after doing these steps.
 */
      else
      {
         cblas_trsm(Order, CblasLeft, CblasUpper, Trans, CblasNonUnit,
                    N, NRHS, one, A, lda, B, ldb);
         cblas_trsm(Order, CblasLeft, CblasLower, Trans, CblasUnit,
                    N, NRHS, one, A, lda, B, ldb);
         ATL_laswp(NRHS, B, ldb, 0, N, ipiv, -1);
      }
   }
/*
 * For row-major arrays, we actually have X^T and B^T, so must tranpose
 * both sides of equation, so what we are solving is:  X' * A' = B'
 */
   else
   {
/*
 *    A = LU*inv(Pc), X' * (LU*inv(Pc))' = B'  ==>  X' * inv(Pc) * U' * L' = B'
 *    X' * inv(Pc) = U' * L' * B', so apply inv(Pc) after solves.
 */
      if (Trans == CblasNoTrans)
      {
         cblas_trsm(Order, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
                    NRHS, N, one, A, lda, B, ldb);
         cblas_trsm(Order, CblasRight, CblasUpper, CblasTrans, CblasUnit,
                    NRHS, N, one, A, lda, B, ldb);
         ATL_laswp(NRHS, B, ldb, 0, N, ipiv, -1);
      }
/*
 *    A' = (LU*inv(Pc))', but Pc is on rows of non-trans matrix, so:
 *    X' * (inv(Pr)*L*U) = B'
 *    X' = (Pr * B') * inv(U) * inv(L)
 *    NOTE: this case is untested
 */
      else
      {
         ATL_laswp(NRHS, B, ldb, 0, N, ipiv, 1);
         cblas_trsm(Order, CblasRight, CblasUpper, CblasNoTrans, CblasUnit,
                    NRHS, N, one, A, lda, B, ldb);
         cblas_trsm(Order, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit,
                    NRHS, N, one, A, lda, B, ldb);
      }
   }
}
