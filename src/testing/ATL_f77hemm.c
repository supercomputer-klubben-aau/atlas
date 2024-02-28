/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 * Code contributers : Antoine P. Petitet, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77hemm )
(
   const enum ATLAS_SIDE     SIDE,
   const enum ATLAS_UPLO     UPLO,
   const int                 M,
   const int                 N,
   const SCALAR              ALPHA,
   const TYPE                * A,
   const int                 LDA,
   const TYPE                * B,
   const int                 LDB,
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
   F77_CHAR                  fside;
   F77_CHAR                  fuplo;
#elif defined( StringCrayStyle )
   F77_CHAR                  fside;
   F77_CHAR                  fuplo;
#endif

   char                      cside, cuplo;

#ifdef ATL_FunkyInts
   const F77_INTEGER         F77M   = M,   F77N   = N,
                             F77lda = LDA, F77ldb = LDB, F77ldc = LDC;
#else
   #define F77M              M
   #define F77N              N
   #define F77lda            LDA
   #define F77ldb            LDB
   #define F77ldc            LDC
#endif

#ifdef TCPLX
   TYPE                      alpha[2], beta[2];

   *alpha = *ALPHA; alpha[1] = ALPHA[1];
   *beta  = *BETA;  beta [1] = BETA [1];
#else
   TYPE                      alpha = ALPHA, beta = BETA;
#endif

   if( SIDE == AtlasRight ) cside = 'R';
   else                     cside = 'L';

   if( UPLO == AtlasLower ) cuplo = 'L';
   else                     cuplo = 'U';

#if   defined( StringSunStyle  )
   F77hemm( &cside, &cuplo, &F77M, &F77N, SADD alpha, A, &F77lda,
            B, &F77ldb, SADD beta, C, &F77ldc, ONE, ONE );
#elif defined( StringCrayStyle )
   fside = ATL_C2F_TransChar( cside );
   fuplo = ATL_C2F_TransChar( cuplo );
   F77hemm( fside,  fuplo,  &F77M, &F77N, SADD alpha, A, &F77lda,
            B, &F77ldb, SADD beta, C, &F77ldc );
#elif defined( StringStructVal )
   fside.len = 1; fside.cp = &cside;
   fuplo.len = 1; fuplo.cp = &cuplo;
   F77hemm( fside,  fuplo,  &F77M, &F77N, SADD alpha, A, &F77lda,
            B, &F77ldb, SADD beta, C, &F77ldc );
#elif defined( StringStructPtr )
   fside.len = 1; fside.cp = &cside;
   fuplo.len = 1; fuplo.cp = &cuplo;
   F77hemm( &fside, &fuplo, &F77M, &F77N, SADD alpha, A, &F77lda,
            B, &F77ldb, SADD beta, C, &F77ldc );
#else
   (void) fprintf( stderr, "\n\nF77/C interface not defined!!\n\n" );
   exit( -1 );
#endif
}
