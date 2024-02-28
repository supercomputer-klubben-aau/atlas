/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2014, 2009 R. Clint Whaley
 */
#define SREAL
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"


int clapack_sgelqf
   (const enum CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
    float *A, ATL_CINT lda, float *TAU)
{
   int ierr=0;
   if (Order != CblasRowMajor && Order != CblasColMajor)
   {
      ierr = -1;
      cblas_xerbla(1, "clapack_sgelqf",
                   "Order must be %d or %d, but is set to %d\n",
                   CblasRowMajor, CblasColMajor, Order);
   }
   if (M < 0)
   {
      ierr = -2;
      cblas_xerbla(2, "clapack_sgelqf",
                   "M cannot be less than zero 0,; is set to %d.\n", M);
   }
   if (N < 0)
   {
      ierr = -3;
      cblas_xerbla(3, "clapack_sgelqf",
                   "N cannot be less than zero 0,; is set to %d.\n", N);
   }
   if (Order == CblasColMajor)
   {
      if (lda < M || lda < 1)
      {
         ierr = -5;
         cblas_xerbla(5, "clapack_sgelqf",
                      "lda must be >= MAX(M,1): lda=%d M=%d\n", lda, M);
      }
   }
   else
   {
      if (lda < N || lda < 1)
      {
         ierr = -5;
         cblas_xerbla(5, "clapack_sgelqf",
                      "lda must be >= MAX(N,1): lda=%d N=%d\n", lda, N);
      }
   }
   if (ierr)
      return(ierr);
   if (Order == CblasColMajor)
      return(ATL_sgelqf(M, N, A, lda, TAU, NULL, 0));
   else
      return(ATL_sgeqrf(N, M, A, lda, TAU, NULL, 0));
}