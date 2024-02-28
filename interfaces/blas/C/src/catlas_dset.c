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

void catlas_dset(const int N, const double alpha, double *X, const int incX)
{
   ATL_dset( N, alpha, X, (incX >= 0 ? incX : -incX) );
}
