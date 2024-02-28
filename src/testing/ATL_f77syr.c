/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77syr )
(
   const enum ATLAS_UPLO     UPLO,
   const int                 N,
   const SCALAR              ALPHA,
   const TYPE                * X,
   const int                 INCX,
   TYPE                      * A,
   const int                 LDA
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
   const F77_INTEGER         F77N = N, F77lda  = LDA, F77incx = INCX;
#else
   #define F77N              N
   #define F77lda            LDA
   #define F77incx           INCX
#endif

#ifdef TCPLX
   TYPE                      alpha[2];

   *alpha = *ALPHA; alpha[1] = ALPHA[1];
#else
   TYPE                      alpha = ALPHA;
#endif

   if( UPLO == AtlasUpper ) cuplo = 'U';
   else                     cuplo = 'L';

   if( INCX < 0 ) X -= ( ( 1 - N ) * INCX ) SHIFT;

#if   defined( StringSunStyle  )
   F77syr( &cuplo, &F77N, SADD alpha, X, &F77incx, A, &F77lda, ONE );
#elif defined( StringCrayStyle )
   fuplo = ATL_C2F_TransChar( cuplo );
   F77syr( fuplo,  &F77N, SADD alpha, X, &F77incx, A, &F77lda );
#elif defined( StringStructVal )
   fuplo.len = 1; fuplo.cp = &cuplo;
   F77syr( fuplo,  &F77N, SADD alpha, X, &F77incx, A, &F77lda );
#elif defined( StringStructPtr )
   fuplo.len = 1; fuplo.cp = &cuplo;
   F77syr( &fuplo, &F77N, SADD alpha, X, &F77incx, A, &F77lda );
#else
   (void) fprintf( stderr, "\n\nF77/C interface not defined!!\n\n" );
   exit( -1 );
#endif
}
