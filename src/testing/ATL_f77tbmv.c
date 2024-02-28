/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 * Code contributers : Antoine P. Petitet, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77tbmv )
(
   const enum ATLAS_UPLO     UPLO,
   const enum ATLAS_TRANS    TRANS,
   const enum ATLAS_DIAG     DIAG,
   const int                 N,
   const int                 K,
   const TYPE                * A,
   const int                 LDA,
   TYPE                      * X,
   const int                 INCX
)
{
#if defined( StringSunStyle )
   #if defined( ATL_FunkyInts )
   F77_INTEGER               ONE = 1;
   #else
   int                       ONE = 1;
   #endif
#elif defined( StringStructVal ) || defined( StringStructPtr )
   F77_CHAR                  fuplo, ftran, fdiag;
#elif defined( StringCrayStyle )
   F77_CHAR                  fuplo, ftran, fdiag;
#endif

   char                      cuplo, ctran, cdiag;

#ifdef ATL_FunkyInts
   const F77_INTEGER         F77N = N, F77K = K, F77lda = LDA, F77incx = INCX;
#else
   #define F77N              N
   #define F77K              K
   #define F77lda            LDA
   #define F77incx           INCX
#endif

   if(      UPLO  == AtlasUpper   ) cuplo = 'U';
   else                             cuplo = 'L';

   if(      DIAG  == AtlasNonUnit ) cdiag = 'N';
   else                             cdiag = 'U';

   if(      TRANS == AtlasNoTrans ) ctran = 'N';
   else if( TRANS == AtlasTrans   ) ctran = 'T';
   else                             ctran = 'C';

   if( INCX < 0 ) X -= ( ( 1 - N ) ) * INCX SHIFT;

#if defined(StringSunStyle)
   F77tbmv( &cuplo, &ctran, &cdiag, &F77N, &F77K, A, &F77lda, X, &F77incx,
            ONE, ONE, ONE );
#elif defined(StringCrayStyle)
   ftran = ATL_C2F_TransChar( ctran );
   fdiag = ATL_C2F_TransChar( cdiag );
   fuplo = ATL_C2F_TransChar( cuplo );
   F77tbmv( fuplo,  ftran,  fdiag,  &F77N, &F77K, A, &F77lda, X, &F77incx );
#elif defined(StringStructVal)
   fuplo.len = 1; fuplo.cp = &cuplo;
   ftran.len = 1; ftran.cp = &ctran;
   fdiag.len = 1; fdiag.cp = &cdiag;
   F77tbmv( fuplo,  ftran,  fdiag,  &F77N, &F77K, A, &F77lda, X, &F77incx );
#elif defined(StringStructPtr)
   fuplo.len = 1; fuplo.cp = &cuplo;
   F77tbmv( &fuplo, &ftran, &fdiag, &F77N, &F77K, A, &F77lda, X, &F77incx );
   ftran.len = 1; ftran.cp = &ctran;
   fdiag.len = 1; fdiag.cp = &cdiag;
#else
   (void) fprintf( stderr, "\n\nF77/C interface not defined!!\n\n" );
   exit(-1);
#endif
}
