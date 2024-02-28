/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"

void F77WRAP_GEQRF(const F77_INTEGER *M, const F77_INTEGER *N,
                   TYPE *A, const F77_INTEGER *lda, TYPE *TAU, TYPE *work,
                   const F77_INTEGER *ldw, F77_INTEGER *info)
{
   *info = ATL_geqrf(*M, *N, A, *lda, TAU, work, *ldw);
}
