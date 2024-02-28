/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define SCPLX
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias1.h"
#endif
#include "atlas_level1.h"
#include "cblas.h"

void cblas_cscal(const int N, const void *alpha, void *X, const int incX)
{
   if (N > 0 && incX > 0)
      ATL_cscal(N, alpha, X, incX);
}
