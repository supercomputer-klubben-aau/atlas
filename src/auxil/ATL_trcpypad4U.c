/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

/*
 * c = upper(A), and pads A with identity to make N a multiple 4.
 * NOTE: padding is at TOP of upper triangular matrix (TRSM starts at bottom)!
 * This makes c safe for use with rank-4 update that assumes N%4 == 0.
 * Sets diagonal to 1.0 if Diag==Unit, else inverts diagonal entries.
 * This routine is used for trsmKL_rk4 (ATLAS/src/blas/level3/kernel).
 */
void Mjoin(PATL,trcpypad4U)
(
   enum ATLAS_DIAG Diag,
   ATL_CINT N,                  /* size of triangular matrix A */
   const TYPE *A,               /* NxN upper triangular matrix */
   ATL_CINT lda,                /* leading dim of A */
   TYPE *c,                     /* N4xN4 cpy of A, padded to N4 with I */
   ATL_CINT ldc                 /* leading dim of A */
)
{
   ATL_UINT i, j;
   ATL_CUINT N4=(N+3)&(~3);  /* N4 = CEIL(N/4)*4 */
   const int mr = N4-N;

   for (j=0; j < mr; j++, c += ldc)
   {
      for (i=0; i < N4; i++)
         c[i] = ATL_rzero;
      c[j] = ATL_rone;
   }
   for (; j < N4; j++, c += ldc, A += lda)
   {
      for (i=0; i < mr; i++)
         c[i] = ATL_rzero;
      for (; i < j; i++)
         c[i] = A[i-mr];
      c[j] = (Diag == AtlasNonUnit) ? 1.0 / A[j-mr] : ATL_rone;
   }
}
