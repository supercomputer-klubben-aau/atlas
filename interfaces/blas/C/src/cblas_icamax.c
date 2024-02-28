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

int cblas_icamax(const int N, const void *X, const int incX)
{
   if (N > 0 && incX > 0)
     return(ATL_icamax(N, X, incX));
   return(0);
}
