/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2007 R. Clint Whaley
 */
#include "atlas_misc.h"

#ifdef SREAL
void ATL_dsgecollapse(const int M, const int N, double *A,
                      const int lda, const int ldc)
/*
 * Copies double precision array C into a float array in-place.
 */
{
   float *C = (float*)A;
   int i, j;

   ATL_assert(ldc <= 2*lda);
   for (j=0; j < N; j++, A += lda, C += ldc)
   {
      for (i=0; i < M; i++)
         C[i] = A[i];
   }
}
#elif defined(DREAL)
void ATL_qdgecollapse(const int M, const int N, ATL_QTYPE *A,
                      const int lda, const int ldc)
/*
 * Copies double precision array C into a float array in-place.
 */
{
   double *C = (double*)A;
   int i, j;

   ATL_assert(ldc <= 2*lda);
   for (j=0; j < N; j++, A += lda, C += ldc)
   {
      for (i=0; i < M; i++)
         C[i] = A[i];
   }
}
#elif defined(SCPLX)
void ATL_zcgecollapse(const int M, const int N, double *A,
                      const int lda, const int ldc)
{
   ATL_dsgecollapse(M+M, N, A, lda+lda, ldc+ldc);
}
#else
void ATL_ezgecollapse(const int M, const int N, ATL_QTYPE *A,
                      const int lda, const int ldc)
{
   ATL_qdgecollapse(M+M, N, A, lda+lda, ldc+ldc);
}
#endif
