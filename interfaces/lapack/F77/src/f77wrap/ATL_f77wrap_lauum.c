/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_LAUUM(const F77_INTEGER *iuplo, const F77_INTEGER *N,
                   TYPE *A, const F77_INTEGER *lda, F77_INTEGER *info
                   ATL_STRLEN_1)
{
   *info = 0;
   ATL_lauum(AtlasColMajor, *iuplo, *N, A, *lda);
}
