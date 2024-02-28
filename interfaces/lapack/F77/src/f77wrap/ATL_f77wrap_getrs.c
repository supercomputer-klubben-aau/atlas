/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_GETRS(const F77_INTEGER *ITRAN, const F77_INTEGER *N,
                   const F77_INTEGER *NRHS, const TYPE *A,
                   const F77_INTEGER *lda, const F77_INTEGER *ipiv0,
                   TYPE *B, const F77_INTEGER *ldb ATL_STRLEN_1)
{
   const int n = *N;
   int i, *ipiv;
   ipiv = malloc(n*sizeof(int));
   ATL_assert(ipiv);
   for (i=0; i != n; i++) ipiv[i] = ipiv0[i] - 1;

   ATL_getrs(AtlasColMajor, *ITRAN, *N, *NRHS, A, *lda, ipiv, B, *ldb);

   free(ipiv);
}
