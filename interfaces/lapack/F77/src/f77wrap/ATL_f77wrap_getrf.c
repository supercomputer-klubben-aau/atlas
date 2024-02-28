/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_GETRF(const F77_INTEGER *M, const F77_INTEGER *N,
                   TYPE *A, const F77_INTEGER *lda, F77_INTEGER *ipiv0,
                   F77_INTEGER *info)
{
   const int MN = Mmin(*M,*N);
   int i;
   #ifdef ATL_FunkyInts
      int *ipiv;
      ipiv = malloc(MN*sizeof(int));
      ATL_assert(ipiv);
   #else
      #define ipiv ipiv0
   #endif
   *info = ATL_getrf(AtlasColMajor, *M, *N, A, *lda, ipiv);
   #ifdef ATL_FunkyInts
      for (i=0; i != MN; i++) ipiv0[i] = ipiv[i] + 1;
      free(ipiv);
   #else
      for (i=0; i != MN; i++) ipiv[i]++;
   #endif
}
