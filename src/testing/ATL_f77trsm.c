/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 * Code contributers : Antoine P. Petitet, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77trsm )
(
   const enum ATLAS_SIDE     SIDE,
   const enum ATLAS_UPLO     UPLO,
   const enum ATLAS_TRANS    TRANS,
   const enum ATLAS_DIAG     DIAG,
   const int                 M,
   const int                 N,
   const SCALAR              ALPHA,
   const TYPE                * A,
   const int                 LDA,
   TYPE                      * B,
   const int                 LDB
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
   F77_CHAR                  ftran;
   F77_CHAR                  fdiag;
#elif defined( StringCrayStyle )
   F77_CHAR                  fside;
   F77_CHAR                  fuplo;
   F77_CHAR                  ftran;
   F77_CHAR                  fdiag;
#endif

   char                      cside, cuplo, ctran, cdiag;

#ifdef ATL_FunkyInts
   const F77_INTEGER         F77M   = M,   F77N   = N,
                             F77lda = LDA, F77ldb = LDB;
#else
   #define F77M              M
   #define F77N              N
   #define F77lda            LDA
   #define F77ldb            LDB
#endif

#ifdef TCPLX
   TYPE                      alpha[2];

   *alpha = *ALPHA; alpha[1] = ALPHA[1];
#else
   TYPE                      alpha = ALPHA;
#endif

   if(      TRANS == AtlasNoTrans ) ctran = 'N';
   else if( TRANS == AtlasTrans   ) ctran = 'T';
   else                             ctran = 'C';

   if(       SIDE  == AtlasRight  ) cside = 'R';
   else                             cside = 'L';

   if(       UPLO  == AtlasLower  ) cuplo = 'L';
   else                             cuplo = 'U';

   if(       DIAG  == AtlasUnit   ) cdiag = 'U';
   else                             cdiag = 'N';

#if   defined( StringSunStyle  )
   F77trsm( &cside, &cuplo, &ctran, &cdiag, &F77M, &F77N, SADD alpha,
            A, &F77lda, B, &F77ldb, ONE, ONE, ONE, ONE );
#elif defined( StringCrayStyle )
   fside = ATL_C2F_TransChar( cside );
   fuplo = ATL_C2F_TransChar( cuplo );
   ftran = ATL_C2F_TransChar( ctran );
   fdiag = ATL_C2F_TransChar( cdiag );
   F77trsm( fside,  fuplo,  ftran,  fdiag,  &F77M, &F77N, SADD alpha,
            A, &F77lda, B, &F77ldb );
#elif defined( StringStructVal )
   fside.len = 1; fside.cp = &cside;
   fuplo.len = 1; fuplo.cp = &cuplo;
   ftran.len = 1; ftran.cp = &ctran;
   fdiag.len = 1; fdiag.cp = &cdiag;
   F77trsm( fside,  fuplo,  ftran,  fdiag,  &F77M, &F77N, SADD alpha,
            A, &F77lda, B, &F77ldb );
#elif defined( StringStructPtr )
   fside.len = 1; fside.cp = &cside;
   fuplo.len = 1; fuplo.cp = &cuplo;
   ftran.len = 1; ftran.cp = &ctran;
   fdiag.len = 1; fdiag.cp = &cdiag;
   F77trsm( &fside, &fuplo, &ftran, &fdiag, &F77M, &F77N, SADD alpha,
            A, &F77lda, B, &F77ldb );
#else
   (void) fprintf( stderr, "\n\nF77/C interface not defined!!\n\n" );
   exit( -1 );
#endif
}
