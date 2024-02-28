/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_lapack.h"
#include "atlas_lamch.h"
int Mjoin(PATL,trtrs)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
    const enum ATLAS_DIAG Diag, ATL_CINT N, ATL_CINT NRHS,
    const TYPE *A, ATL_CINT lda, TYPE *B, ATL_CINT ldb)
/*
 * Checks for singularity, and then solves system:
 *   A * X = B or A^T * X = B
 * where A is a triangular matrix (as indicated by Uplo).
 * RETURNS :
 *   0 : successful exit
 *  <0 : argument #(-return) had illegal value (start counting from 1)
 *  >0 : diag elt # (1st elt 1) was zero, so A is singular
 */
{
   #ifdef TCPLX
      TYPE one[2] = {ATL_rone, ATL_rzero};
      ATL_CINT N2=N+N;
   #else
      #define one ATL_rone
   #endif
   ATL_CINT ldap1 = (lda+1)SHIFT;
   ATL_INT i;
/*
 * Zero on diagonal means singular triangular matrix
 */
   if (Diag != AtlasUnit)
   {
      #ifdef TCPLX
         for (i=0; i < N2; i += 2, A += ldap1)
            if (SCALAR_IS_ZERO(A))
               return((i>>1)+1);
      #else
         for (i=0; i < N; i++, A += ldap1)
            if (*A == ATL_rzero)
               return(i+1);
      #endif
      A -= ldap1*N;
   }
   cblas_trsm(CblasColMajor, AtlasLeft, Uplo, TA, Diag, N, NRHS, one,
              A, lda, B, ldb);
   return(0);
}
