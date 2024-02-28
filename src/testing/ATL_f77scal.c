/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77scal )
(
   const int                 N,
   const SCALAR              ALPHA,
   TYPE                      * X,
   const int                 INCX
)
{
   const F77_INTEGER         F77N = N, F77incx = Mabs(INCX);
#ifdef TCPLX
   TYPE                      alpha[2];

   *alpha   = *ALPHA;
   alpha[1] = ALPHA[1];
#else
   TYPE                      alpha = ALPHA;
#endif

   if( INCX < 0 ) X -= ( ( 1 - N ) * INCX ) SHIFT;

   F77scal( &F77N, SADD alpha, X, &F77incx );
}
