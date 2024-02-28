/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,hbnrm)
(const enum ATLAS_UPLO UPLO, const int N, const int K,
 const TYPE *A, const int LDA)
{
   int                        i, i0, i1, iaij, iy, j, jaj, ky = 0, l,
                              lda2 = (LDA SHIFT);
   TYPE max=ATL_rzero, t0, * work= NULL;

   if( N <= 0 ) return( ATL_rzero );

   work = (TYPE *)malloc( N * sizeof( TYPE ) );
   if( work == NULL )
   {fprintf( stderr, "mem alloc failed in [sb,hb]nrm, bye ...\n" ); exit( 1 );}
   else { for( i = 0; i < N; i++ ) work[i] = ATL_rzero; }

   if( UPLO == AtlasUpper )
   {
      for( j = 0, jaj = 0; j < N; j++, jaj += lda2 )
      {
         t0      = ATL_rzero;

         l     = K - j;
         i0    = ( j - K > 0 ? j - K : 0 );

         for( i = i0, iaij  = ((l+i0) SHIFT)+jaj, iy = ky;
              i < j;  i++, iaij += (1 SHIFT), iy += 1 )
         {
            work[iy] += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
            t0       += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
         }
         work[j] += Mabs( A[iaij] ) + t0;

         if( j >= K ) { ky += 1; }
      }
   }
   else
   {
      for( j = 0, jaj = 0; j < N; j++, jaj += lda2 )
      {
         t0      = ATL_rzero;
         work[j] = Mabs( A[jaj] );
         i1     = ( N - 1 > j + K ? j + K : N - 1 );
         for( i = j+1, iaij = (1 SHIFT)+jaj; i <= i1; i++,
              iaij += (1 SHIFT) )
         {
            work[i] += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
            t0      += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
         }
         work[j] += t0;
      }
   }

   max = work[0];
   for( j = 1; j < N; j++ )
   {
      const TYPE t0=work[j];
      if (t0 != t0 || max < t0)
         max = t0;
   }

   if( work ) free( work );

   return( max );
}
