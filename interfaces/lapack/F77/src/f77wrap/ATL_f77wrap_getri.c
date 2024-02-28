/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_GETRI(const F77_INTEGER *N, TYPE *A, const F77_INTEGER *lda,
                   const F77_INTEGER *ipiv0, TYPE *work, F77_INTEGER *lwork,
                   F77_INTEGER *info)
{
   const int n = *N;
   int *ipiv=NULL;
   int i, lwrk = *lwork;

   if (lwrk != -1)
   {
      ipiv = malloc(n*sizeof(int));
      ATL_assert(ipiv);
      for (i=0; i != n; i++) ipiv[i] = ipiv0[i] - 1;
   }
   *info = ATL_getri(AtlasColMajor, *N, A, *lda, ipiv, work, &lwrk);
   if (work) *work = lwrk;
   else if (*lwork == -1)
      ATL_xerbla(5, __FILE__,
                 "For workspace query, workspace cannot be NULL\n");
   if (ipiv) free(ipiv);
}
