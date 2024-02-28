/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define DREAL
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias1.h"
#endif
#include "atlas_level1.h"
#include "cblas.h"

void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY)
{
   if (N > 0)
   {
      if (incX < 0)
      {
         if (incY < 0) ATL_daxpy(N, alpha, X, -incX, Y, -incY);
         else ATL_daxpy(N, alpha, X+(1-N)*incX, incX, Y, incY);
      }
      else if (incY < 0) ATL_daxpy(N, alpha, X+(N-1)*incX, -incX, Y, -incY);
      else ATL_daxpy(N, alpha, X, incX, Y, incY);
   }
}
