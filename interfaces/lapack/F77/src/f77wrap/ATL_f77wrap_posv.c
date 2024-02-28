/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"

void F77WRAP_POSV(const F77_INTEGER *iuplo, const F77_INTEGER *N,
                  const F77_INTEGER *NRHS, TYPE *A, const F77_INTEGER *lda,
                  TYPE *B, const F77_INTEGER *ldb, F77_INTEGER *info
                  ATL_STRLEN_1)
{
   *info = ATL_potrf(AtlasColMajor, *iuplo, *N, A, *lda);
   if (*info == 0)
      ATL_potrs(AtlasColMajor, *iuplo, *N, *NRHS, A, *lda, B, *ldb);
}
