/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

int Mjoin( PATL, f77amax )
(
   const int                 N,
   const TYPE                * X,
   const int                 INCX
)
{
   const F77_INTEGER         F77N = N, F77incx = Mabs(INCX);
   int                       imax = 0;

   if( N > 0 )
   {
      if( INCX < 0 ) X -= ( ( 1 - N ) * INCX ) SHIFT;
      imax = F77amax( &F77N, X, &F77incx ) - 1;
   }
   return( imax );
}
