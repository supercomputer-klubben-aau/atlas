/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_f77blas.h"
#include "atlas_f77wrap.h"
void F77set(const F77_INTEGER *N, const TYPE *alpha, TYPE *X,
            const F77_INTEGER *incX)
{
   int incx=(*incX);
   if (incx < 0) incx = -incx;
   #ifdef TREAL
      Mjoin(PATL,set)(*N, *alpha, X, incx);
   #else
      Mjoin(PATL,set)(*N, alpha, X, incx);
   #endif
}
