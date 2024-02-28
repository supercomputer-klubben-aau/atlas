/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,tbnrm1)(const enum ATLAS_UPLO UPLO, const enum ATLAS_DIAG DIAG,
                        const int N, const int K, const TYPE *A, const int LDA)
/*
 * Calculates the 1-norm of a triangular band matrix
 */
{
   int i, i0, i1, iaij, j, jaj, l, lda2 = ( LDA SHIFT );
   TYPE max=0.0, t0;

   if( UPLO == AtlasUpper )
   {
      for( j = 0, jaj  = 0; j < N; j++, jaj += lda2 )
      {
         l       = K - j;
         i0      = ( j - K > 0 ? j - K : 0 );

         t0 = ATL_rzero;
         for( i = i0, iaij = ((l+i0) SHIFT)+jaj; i < j; i++, iaij += (1 SHIFT) )
         {
#ifdef TREAL
            t0 += Mabs( A[iaij] );
#else
            t0 += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
#endif
            if (t0 != t0)
               return(t0);
         }
         if( DIAG == AtlasNonUnit ) t0 += ATL_rone;

         if (t0 > max) max = t0;
      }
   }
   else
   {
      for( j = N-1, jaj = (N-1)*lda2; j >= 0; j--, jaj -= lda2 )
      {
         t0 = ATL_rzero;
         if( DIAG == AtlasNonUnit ) t0 += ATL_rone;

         i1   = ( N - 1 > j + K ? j + K : N - 1 );
         for( i  = j+1, iaij = (1 SHIFT)+jaj; i <= i1; i++, iaij += (1 SHIFT) )
         {
#ifdef TREAL
            t0 += Mabs( A[iaij] );
#else
            t0 += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
#endif
            if (t0 != t0)
               return(t0);
         }
         if (t0 > max) max = t0;
      }
   }
   return(max);
}
