/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

#ifdef TREAL
TYPE Mjoin( PATL, f77nrm2 )
#else
TYPE Mjoin( PATLU, Mjoin( PRE, f77nrm2 ) )
#endif
(
   const int                 N,
   const TYPE                * X,
   const int                 INCX
)
{
   TYPE                      nrm2;
   const F77_INTEGER         F77N=N, F77incx = Mabs(INCX);
   if( INCX < 0 ) X -= ( ( 1 - N ) * INCX ) SHIFT;

   F77nrm2( &F77N, X, &F77incx, &nrm2 );

   return( nrm2 );
}
