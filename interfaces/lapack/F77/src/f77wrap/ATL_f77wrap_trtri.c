/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_TRTRI(const F77_INTEGER *iuplo, const F77_INTEGER *idiag,
                   const F77_INTEGER *N, TYPE *A,
                   const F77_INTEGER *lda, F77_INTEGER *info
                   ATL_STRLEN_1)
{
   *info = ATL_trtri(AtlasColMajor, *iuplo, *idiag, *N, A, *lda);
}
