/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77rotg )
(
   TYPE                      * A,
#ifdef TREAL
   TYPE                      * B,
#else
   const SCALAR              B,
#endif
   TYPE                      * C,
   TYPE                      * S
)
{
#ifdef TREAL
   TYPE                      * b = B;
#else
   TYPE                      b[2];

   *b = *B; b[1] = B[1];
#endif

   F77rotg( A, b, C, S );
}
