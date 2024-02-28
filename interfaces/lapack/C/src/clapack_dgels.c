/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#define DREAL
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"

#define mytrans CblasTrans
int clapack_dgels(const enum CBLAS_ORDER Order,
                  const enum CBLAS_TRANSPOSE TA,
                  ATL_CINT M, ATL_CINT N, ATL_CINT NRHS,
                  double *A, ATL_CINT lda, double *B, const int ldb)
/*
 *  GELS solves overdetermined or underdetermined linear systems
 *  involving an M-by-N matrix A, or its conjugate-transpose, using a QR
 *  or LQ factorization of A.  It is assumed that A has full rank.
 */
{
   int ierr = 0;

   if (Order != CblasRowMajor && Order != CblasColMajor)
   {
      ierr = -1;
      cblas_xerbla(1, "clapack_dgesv",
                   "Order must be %d or %d, but is set to %d.\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (TA != AtlasNoTrans && TA != mytrans)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_dgesv",
                   "Trans must be %d or %d, but is set to %d.\n",
                   CblasNoTrans, mytrans, TA);
   }
   if (M < 0)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_dgesv",
                   "M cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (N < 0)
   {
      ierr = -4;
      cblas_xerbla(4, "clapack_dgesv",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (NRHS < 0)
   {
      ierr = -5;
      cblas_xerbla(5, "clapack_dgesv",
                   "NRHS cannot be less than zero 0,; is set to %d.\n", NRHS);
   }
   if (lda < M || lda < 1)
   {
      ierr = -7;
      cblas_xerbla(7, "clapack_dgesv",
                   "lda must be >= MAX(M,1): lda=%d M=%d\n", lda, M);
   }
   if (ldb < Mmax(M,N) || ldb < 1)
   {
      ierr = -9;
      cblas_xerbla(9, "clapack_dgesv",
                   "ldb must be >= MAX(M,N,1): ldb=%d M=%d N=%d\n", ldb, M, N);

   }
   if (Order == CblasColMajor)
      ierr = ATL_dgels(TA, M, N, NRHS, A, lda, B, ldb, NULL, 0);
   else  /* row-major array */
   {
      enum CBLAS_TRANSPOSE TAr = (TA == AtlasNoTrans) ? mytrans : CblasNoTrans;
      ierr = ATL_dgels(TAr, N, M, NRHS, A, lda, B, ldb, NULL, 0);
   }
   return(ierr);
}
#undef mytrans
