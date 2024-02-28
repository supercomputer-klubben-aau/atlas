/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define DCPLX
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias1.h"
#endif
#include "atlas_level1.h"
#include "cblas.h"

double cblas_dzasum(const int N, const void *X, const int incX)
{
   if (N > 0 && incX > 0)
      return(ATL_dzasum(N, X, incX));
   return(0.0);
}
