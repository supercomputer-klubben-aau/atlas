/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

/*
 * c = lower(A), and pads A with identity to make N a multiple 4.  This is
 * used to make c safe for use with rank-4 update that assumes N mul of 4.
 * Sets diagonal to 1.0 if Diag==Unit, else inverts diagonal entried.
 * This routine is used for trsmKL_rk4 (ATLAS/src/blas/level3/kernel).
 */
void Mjoin(PATL,trcpypad4L)
(
   enum ATLAS_DIAG Diag,
   ATL_CINT N,                  /* size of triangular matrix A */
   const TYPE *A,               /* NxN lower triangular matrix */
   ATL_CINT lda,                /* leading dim of A */
   TYPE *c,                     /* N4xN4 cpy of A, padded to N4 with I */
   ATL_CINT ldc                 /* leading dim of A */
)
{
   ATL_UINT i, j;
   ATL_CUINT N4=(N+3)&(~3);  /* N4 = CEIL(N/4)*4 */

   for (j=0; j < N; j++, c += ldc, A += lda)
   {
      if (Diag == AtlasUnit)
         c[j] = 1.0;
      else
         c[j] = 1.0 / A[j];
      for (i=j+1; i < N; i++)
         c[i] = A[i];
      for (; i < N4; i++)
         c[i] = ATL_rzero;
   }
   for (; j < N4; j++, c += ldc)
   {
      for (i=0; i < N4; i++)
         c[i] = ATL_rzero;
      c[j] = ATL_rone;
   }
}
