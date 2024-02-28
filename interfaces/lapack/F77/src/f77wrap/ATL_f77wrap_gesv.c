/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_GESV(const F77_INTEGER *N, const F77_INTEGER *NRHS,
                  TYPE *A, const F77_INTEGER *lda, F77_INTEGER *ipiv0,
                  TYPE *B, const F77_INTEGER *ldb, F77_INTEGER *info)
{
   const int n = *N;
   int i;
   #ifdef ATL_FunkyInts
      int *ipiv;
      ipiv = malloc(n*sizeof(int));
      ATL_assert(ipiv);
   #else
      #define ipiv ipiv0
   #endif

   *info = ATL_getrf(AtlasColMajor, *N, *N, A, *lda, ipiv);
   if (*info == 0)
      ATL_getrs(AtlasColMajor, AtlasNoTrans, *N, *NRHS, A, *lda, ipiv, B, *ldb);
   #ifdef ATL_FunkyInts
      for (i=0; i != n; i++) ipiv0[i] = ipiv[i] + 1;
      free(ipiv);
   #else
      for (i=0; i != n; i++) ipiv[i]++;
   #endif
}
