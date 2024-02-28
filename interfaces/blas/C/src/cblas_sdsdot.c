/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#define SREAL
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias1.h"
#endif
#include "atlas_level1.h"
#include "cblas.h"

float cblas_sdsdot(const int N, const float alpha, const float *X,
                   const int incX, const float *Y, const int incY)
{
   if (N > 0)
   {
      if (incX < 0)
      {
         if (incY < 0) return(ATL_sdsdot(N, alpha, X, -incX, Y, -incY));
         else return(ATL_sdsdot(N, alpha, X+(1-N)*incX, incX, Y, incY));
      }
      else if (incY < 0)
         return(ATL_sdsdot(N, alpha, X+(N-1)*incX, -incX, Y, -incY));
      else return(ATL_sdsdot(N, alpha, X, incX, Y, incY));
   }
   else return(0.0f);
}
