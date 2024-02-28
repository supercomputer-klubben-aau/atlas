/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_LARFB
   (const F77_INTEGER *iside, const F77_INTEGER *itrans,
    const F77_INTEGER *idirect, const F77_INTEGER *istore,
    const F77_INTEGER *M, const F77_INTEGER *N, const F77_INTEGER *K,
    TYPE *V, const F77_INTEGER *ldv, TYPE *T, const F77_INTEGER *ldt,
    TYPE *C, const F77_INTEGER *ldc, TYPE *work, const F77_INTEGER *ldwork)
{
   ATL_larfb(*iside, *itrans, *idirect, *istore, *M, *N, *K, V, *ldv, T, *ldt,
             C, *ldc, work, *ldwork);
}
