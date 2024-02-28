/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#include "atlas_misc.h"
#include "atlas_level1.h"

float ATL_sdsdot(const int N, const float alpha, const float *X,
                 const int incX, const float *Y, const int incY)
{
   int i;
   double dot = alpha;
   for (i=N; i; i--, X += incX, Y += incY) dot += (double)(*X) * (double) (*Y);
   return(dot);
}
