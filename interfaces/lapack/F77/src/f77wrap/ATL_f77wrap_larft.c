/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_LARFT
   (const F77_INTEGER *idirect, const F77_INTEGER *istore,
    const F77_INTEGER *N, const F77_INTEGER *K,
    TYPE *V, const F77_INTEGER *ldv, const TYPE *TAU,
    TYPE *T, const F77_INTEGER *ldt)
{
   ATL_larft(*idirect, *istore, *N, *K, V, *ldv, TAU, T, *ldt);
}
