/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,tpnrm1)(const enum ATLAS_UPLO UPLO, const enum ATLAS_DIAG DIAG,
                        const int N, const TYPE *A)
/*
 * Calculates the 1-norm of a triangular packed matrix
 */
{
   int i, iaij, j;
   TYPE max=0.0, t0;

   if( UPLO == AtlasUpper )
   {
      for( j = 0, iaij= 0; j < N; j++ )
      {
         t0 = ATL_rzero;
         for( i = 0; i < j; i++, iaij += (1 SHIFT) )
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
         iaij += (1 SHIFT);
      }
   }
   else
   {
      for( j = N-1, iaij = ((((N-1)*(N+2)) >> 1) SHIFT); j >= 0; j-- )
      {
         t0 = ATL_rzero;
         if( DIAG == AtlasNonUnit ) t0 += ATL_rone;
         iaij += (1 SHIFT);
         for( i = j+1; i < N; i++, iaij += (1 SHIFT) )
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

         iaij -= ( ( N - j ) << (1 SHIFT) ) + (1 SHIFT);
      }
   }
   return( max );
}
