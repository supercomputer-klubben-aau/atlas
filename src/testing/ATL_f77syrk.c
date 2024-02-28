/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 * Code contributers : Antoine P. Petitet, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77syrk )
(
   const enum ATLAS_UPLO     UPLO,
   const enum ATLAS_TRANS    TRANS,
   const int                 N,
   const int                 K,
   const SCALAR              ALPHA,
   const TYPE                * A,
   const int                 LDA,
   const SCALAR              BETA,
   TYPE                      * C,
   const int                 LDC
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
   F77_CHAR                  ftran;
#elif defined( StringCrayStyle )
   F77_CHAR                  fuplo;
   F77_CHAR                  ftran;
#endif

   char                      cuplo, ctran;

#ifdef ATL_FunkyInts
   const F77_INTEGER         F77N   = N,   F77K   = K,
                             F77lda = LDA, F77ldc = LDC;
#else
   #define F77N              N
   #define F77K              K
   #define F77lda            LDA
   #define F77ldc            LDC
#endif

#ifdef TCPLX
   TYPE                      alpha[2], beta [2];

   *alpha = *ALPHA; alpha[1] = ALPHA[1];
   *beta  = *BETA;  beta [1] = BETA [1];
#else
   TYPE                      alpha = ALPHA, beta = BETA;
#endif

   if(      UPLO  == AtlasLower   ) cuplo = 'L';
   else                             cuplo = 'U';

#ifdef TREAL
   if(      TRANS == AtlasNoTrans ) ctran = 'N';
   else if( TRANS == AtlasTrans   ) ctran = 'T';
   else                             ctran = 'C';
#else
   if(      TRANS == AtlasNoTrans ) ctran = 'N';
   else                             ctran = 'T';
#endif

#if   defined( StringSunStyle  )
   F77syrk( &cuplo, &ctran, &F77N, &F77K, SADD alpha, A, &F77lda,
            SADD beta, C, &F77ldc, ONE, ONE );
#elif defined( StringCrayStyle )
   fuplo = ATL_C2F_TransChar( cuplo );
   ftran = ATL_C2F_TransChar( ctran );
   F77syrk( fuplo,  ftran,  &F77N, &F77K, SADD alpha, A, &F77lda,
            SADD beta, C, &F77ldc );
#elif defined( StringStructVal )
   fuplo.len = 1; fuplo.cp = &cuplo;
   ftran.len = 1; ftran.cp = &ctran;
   F77syrk( fuplo,  ftran,  &F77N, &F77K, SADD alpha, A, &F77lda,
            SADD beta, C, &F77ldc );
#elif defined( StringStructPtr )
   fuplo.len = 1; fuplo.cp = &cuplo;
   ftran.len = 1; ftran.cp = &ctran;
   F77syrk( &fuplo, &ftran, &F77N, &F77K, SADD alpha, A, &F77lda,
            SADD beta, C, &F77ldc );
#else
   (void) fprintf( stderr, "\n\nF77/C interface not defined!!\n\n" );
   exit( -1 );
#endif
}
