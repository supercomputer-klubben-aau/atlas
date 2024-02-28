/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, Mjoin( UPR, f77scal ) )
(
   const int                 N,
   const TYPE                ALPHA,
   TYPE                      * X,
   const int                 INCX
)
{
   const F77_INTEGER         F77N = N, F77incx = Mabs(INCX);
   TYPE                      alpha = ALPHA;

   if( INCX < 0 ) X -= ( ( 1 - N ) * INCX ) SHIFT;

   F77rscal( &F77N, &alpha, X, &F77incx );
}
