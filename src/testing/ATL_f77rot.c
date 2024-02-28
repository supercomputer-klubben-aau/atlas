/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

#ifdef TREAL
void Mjoin( PATL, f77rot )
#else
void Mjoin( PATL, Mjoin( UPR, f77rot ) )
#endif
(
   const int                 N,
   TYPE                      * X,
   const int                 INCX,
   TYPE                      * Y,
   const int                 INCY,
#ifdef TREAL
   const SCALAR              C,
   const SCALAR              S
#else
   const TYPE                C,
   const TYPE                S
#endif
)
{
#ifdef ATL_FunkyInts
   const F77_INTEGER         F77N = N, F77incx = INCX, F77incy = INCY;
#else
   #define F77N              N
   #define F77incx           INCX
   #define F77incy           INCY
#endif
   TYPE                      c = C, s = S;

   if( INCX < 0 ) X -= ( ( 1 - N ) * INCX ) SHIFT;
   if( INCY < 0 ) Y -= ( ( 1 - N ) * INCY ) SHIFT;

#ifdef TREAL
   F77rot( &F77N, X, &F77incx, Y, &F77incy, SADD c, SADD s );
#else
   F77rot( &F77N, X, &F77incx, Y, &F77incy, &c,     &s     );
#endif
}
