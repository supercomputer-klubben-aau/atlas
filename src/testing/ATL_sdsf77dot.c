/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

float ATL_sdsf77dot
(
   const int                  N,
   const float                B,
   const float                * X,
   const int                  INCX,
   const float                * Y,
   const int                  INCY
)
{
   float                      b = B, dot;
#ifdef ATL_FunkyInts
   const F77_INTEGER          F77N = N, F77incx = INCX, F77incy = INCY;
#else
   #define F77N               N
   #define F77incx            INCX
   #define F77incy            INCY
#endif

   if( INCX < 0 ) X -= ( 1 - N ) * INCX;
   if( INCY < 0 ) Y -= ( 1 - N ) * INCY;

   F77sdsdot( &F77N, &b, X, &F77incx, Y, &F77incy, &dot );

   return( dot );
}
