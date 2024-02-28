/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2007 R. Clint Whaley
 */
#include "atlas_misc.h"

#ifdef SREAL
void ATL_dstrcollapse(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                      const int N, double *A, const int lda, const int ldc)
#elif defined(SCPLX)
void ATL_zctrcollapse(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                      const int N, double *A, const int lda0, const int ldc0)
#elif defined(DCPLX)
void ATL_eztrcollapse(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                      const int N, ATL_QTYPE *A, const int lda0, const int ldc0)
#else
void ATL_qdtrcollapse(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                      const int N, ATL_QTYPE *A, const int lda, const int ldc)
#endif
/*
 * Translates upper/lower triangle from double to single precision
 */
{
#if defined(SREAL) || defined(SCPLX)
   float *C = (float*)A;
#else
   double *C = (double*)A;
#endif
#ifdef TCPLX
   const int lda = lda0+lda0, ldc = ldc0+ldc0;
#endif
   int i, j, st;

   ATL_assert(ldc <= 2*lda);
   if (Uplo == AtlasUpper)
   {
      for (j=0; j < N; j++, C += ldc, A += lda)
      {
         st = (Diag == AtlasUnit) ? j-1 : j;
         for (i=0; i < st; i++)
         #ifdef TREAL
            C[i] = A[i];
         #else
         {
            C[i+i] = A[i+i];
            C[i+i+1] = A[i+i+1];
         }
         #endif
      }
   }
   else
   {
      for (j=0; j < N; j++, C += ldc, A += lda)
      {
         for (i=(Diag == AtlasUnit) ? j+1 : j; i < N; i++)
         #ifdef TREAL
            C[i] = A[i];
         #else
         {
            C[i+i] = A[i+i];
            C[i+i+1] = A[i+i+1];
         }
         #endif
      }
   }
}
