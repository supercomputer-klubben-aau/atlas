/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"

void F77WRAP_GELS(const F77_INTEGER *itrans, const F77_INTEGER *M,
                  const F77_INTEGER *N, const F77_INTEGER *NRHS,
                  TYPE *A, const F77_INTEGER *lda, TYPE *B,
                  const F77_INTEGER *ldb, TYPE *work, const F77_INTEGER *lwork,
                  F77_INTEGER *info)
{
   *info = ATL_gels(*itrans, *M, *N, *NRHS, A, *lda, B, *ldb, work, *lwork);
}
