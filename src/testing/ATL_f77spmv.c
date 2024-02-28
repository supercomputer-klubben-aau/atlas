/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77spmv )
(
   const enum ATLAS_UPLO     UPLO,
   const int                 N,
   const SCALAR              ALPHA,
   const TYPE                * A,
   const TYPE                * X,
   const int                 INCX,
   const SCALAR              BETA,
   TYPE                      * Y,
   const int                 INCY
)
{
#if defined( StringSunStyle )
   #if defined( ATL_FunkyInts )
   F77_INTEGER               ONE = 1;
   #else
   int                       ONE = 1;
   #endif
#elif defined( StringStructVal ) || defined( StringStructPtr )
   F77_CHAR                  fuplo;
#elif defined( StringCrayStyle )
   F77_CHAR                  fuplo;
#endif

   char                      cuplo;

#ifdef ATL_FunkyInts
   const F77_INTEGER         F77N    = N,
                             F77incx = INCX, F77incy = INCY;
#else
   #define F77N              N
   #define F77incx           INCX
   #define F77incy           INCY
#endif

#ifdef TCPLX
   TYPE                      alpha[2], beta[2];

   *alpha = *ALPHA; alpha[1] = ALPHA[1];
   *beta  = *BETA;  beta [1] = BETA [1];
#else
   TYPE                      alpha = ALPHA, beta = BETA;
#endif

   if( UPLO == AtlasUpper ) cuplo = 'U';
   else                     cuplo = 'L';

   if( INCX < 0 ) X -= ( ( 1 - N ) * INCX ) SHIFT;
   if( INCY < 0 ) Y -= ( ( 1 - N ) * INCY ) SHIFT;

#if   defined( StringSunStyle  )
   F77spmv( &cuplo, &F77N,        SADD alpha, A,          X, &F77incx,
            SADD beta, Y, &F77incy, ONE );
#elif defined( StringCrayStyle )
   fuplo = ATL_C2F_TransChar( cuplo );
   F77spmv( fuplo,  &F77N,        SADD alpha, A,          X, &F77incx,
            SADD beta, Y, &F77incy );
#elif defined( StringStructVal )
   fuplo.len = 1; fuplo.cp = &cuplo;
   F77spmv( fuplo,  &F77N,        SADD alpha, A,          X, &F77incx,
            SADD beta, Y, &F77incy );
#elif defined( StringStructPtr )
   fuplo.len = 1; fuplo.cp = &cuplo;
   F77spmv( &fuplo, &F77N,        SADD alpha, A,          X, &F77incx,
            SADD beta, Y, &F77incy );
#else
   (void) fprintf( stderr, "\n\nF77/C interface not defined!!\n\n" );
   exit( -1 );
#endif
}
