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

void cblas_csscal(const int N, const float alpha, void *X, const int incX)
{
   float al[2];
   if (N > 0 && incX > 0)
   {
      al[0] = alpha;
      al[1] = ATL_rzero;
      if (incX < 0) ATL_cscal(N, al, X, -incX);
      else ATL_cscal(N, al, X, incX);
   }
}
