/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#define DCPLX
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"

#define mytrans CblasConjTrans
int clapack_zgels(const enum CBLAS_ORDER Order,
                  const enum CBLAS_TRANSPOSE TA,
                  ATL_CINT M, ATL_CINT N, ATL_CINT NRHS,
                  void *A, ATL_CINT lda, void *B, const int ldb)
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
      cblas_xerbla(1, "clapack_zgesv",
                   "Order must be %d or %d, but is set to %d.\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (TA != AtlasNoTrans && TA != mytrans)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_zgesv",
                   "Trans must be %d or %d, but is set to %d.\n",
                   CblasNoTrans, mytrans, TA);
   }
   if (M < 0)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_zgesv",
                   "M cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (N < 0)
   {
      ierr = -4;
      cblas_xerbla(4, "clapack_zgesv",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (NRHS < 0)
   {
      ierr = -5;
      cblas_xerbla(5, "clapack_zgesv",
                   "NRHS cannot be less than zero 0,; is set to %d.\n", NRHS);
   }
   if (lda < M || lda < 1)
   {
      ierr = -7;
      cblas_xerbla(7, "clapack_zgesv",
                   "lda must be >= MAX(M,1): lda=%d M=%d\n", lda, M);
   }
   if (ldb < Mmax(M,N) || ldb < 1)
   {
      ierr = -9;
      cblas_xerbla(9, "clapack_zgesv",
                   "ldb must be >= MAX(M,N,1): ldb=%d M=%d N=%d\n", ldb, M, N);

   }
   if (Order == CblasColMajor)
      ierr = ATL_zgels(TA, M, N, NRHS, A, lda, B, ldb, NULL, 0);
   else  /* row-major array */
   {
      enum CBLAS_TRANSPOSE TAr = (TA == AtlasNoTrans) ? mytrans : CblasNoTrans;
      TYPE *a=((TYPE*)A)+1;
      ATL_CINT lda2 = lda+lda;
      ATL_INT j;
/*
 *    Complex takes only the conjugate tranpose, so conjugate entries before
 *    calling with reversed transpose setting
 */
      for (j=0; j < N; j++, a += lda2)
         Mjoin(PATLU,scal)(N, ATL_rnone, a, 2);
      ierr = ATL_zgels(TAr, N, M, NRHS, A, lda, B, ldb, NULL, 0);
   }
   return(ierr);
}
#undef mytrans
