/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 * Code contributers : Antoine P. Petitet, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77gbmv )
(
   const enum ATLAS_TRANS    Trans,
   const int                 M,
   const int                 N,
   const int                 KL,
   const int                 KU,
   const SCALAR              ALPHA,
   const TYPE                * A,
   const int                 LDA,
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
   F77_CHAR                  ftran;
#elif defined( StringCrayStyle )
   F77_CHAR                  ftran;
#endif

   char                      ctran;

#ifdef ATL_FunkyInts
   const F77_INTEGER         F77M    = M,   F77N    = N,
                             F77KL   = KL,  F77KU   = KU,
                             F77lda  = LDA, F77incx = INCX, F77incy = INCY;
#else
   #define F77M M
   #define F77N N
   #define F77KL             KL
   #define F77KU             KU
   #define F77lda            LDA
   #define F77incx           INCX
   #define F77incy           INCY
#endif
   int                       nX, nY;

#ifdef TCPLX
   TYPE                      alpha[2], beta[2];

   *alpha = *ALPHA; alpha[1] = ALPHA[1];
   *beta  = *BETA;  beta [1] = BETA [1];
#else
   TYPE                      alpha = ALPHA, beta = BETA;
#endif

   if( Trans == AtlasNoTrans )
   {
      nX = N;
      nY = M;
      ctran = 'N';
   }
   else
   {
      nX = M;
      nY = N;
      if( Trans == AtlasTrans ) ctran = 'T';
      else                      ctran = 'C';
   }

   if( INCX < 0 ) X -= ( ( 1 - nX ) * INCX ) SHIFT;
   if( INCY < 0 ) Y -= ( ( 1 - nY ) * INCY ) SHIFT;

#if   defined( StringSunStyle  )
   F77gbmv( &ctran, &F77M, &F77N, &F77KL, &F77KU, SADD alpha, A, &F77lda,
            X, &F77incx, SADD beta, Y, &F77incy, ONE );
#elif defined( StringCrayStyle )
   ftran = ATL_C2F_TransChar( ctran );
   F77gbmv( ftran,  &F77M, &F77N, &F77KL, &F77KU, SADD alpha, A, &F77lda,
            X, &F77incx, SADD beta, Y, &F77incy );
#elif defined( StringStructVal )
   ftran.len = 1; ftran.cp = &ctran;
   F77gbmv( ftran,  &F77M, &F77N, &F77KL, &F77KU, SADD alpha, A, &F77lda,
            X, &F77incx, SADD beta, Y, &F77incy );
#elif defined( StringStructPtr )
   ftran.len = 1; ftran.cp = &ctran;
   F77gbmv( &ftran, &F77M, &F77N, &F77KL, &F77KU, SADD alpha, A, &F77lda,
            X, &F77incx, SADD beta, Y, &F77incy );
#else
   (void) fprintf( stderr, "\n\nF77/C interface not defined!!\n\n" );
   exit( -1 );
#endif
}
