/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

#ifdef TREAL
TYPE Mjoin( PATL, f77asum )
#else
TYPE Mjoin( PATLU, Mjoin( PRE, f77asum ) )
#endif
(
   const int                 N,
   const TYPE                * X,
   const int                 INCX
)
{
   TYPE                      asum;
   const F77_INTEGER         F77N = N, F77incx = Mabs(INCX);
   if( INCX < 0 ) X -= ( ( 1 - N ) * INCX ) SHIFT;

   F77asum( &F77N, X, &F77incx, &asum );

   return( asum );
}
