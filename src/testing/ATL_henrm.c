/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,henrm)
(const enum ATLAS_UPLO UPLO, const int N, const TYPE *A, const int LDA)
{
   int                        i, iaij, j, jaj, lda2 = ( LDA SHIFT ),
                              ldap12 = ( ( LDA + 1 ) SHIFT );
   TYPE max=ATL_rzero, t0, * work= NULL;

   if( N <= 0 ) return( ATL_rzero );

   work = (TYPE *)malloc( N * sizeof( TYPE ) );
   if( work == NULL )
   {fprintf( stderr, "mem alloc failed in [sy,he]nrm, bye ...\n" ); exit( 1 );}
   else { for( i = 0; i < N; i++ ) work[i] = ATL_rzero; }

   if( UPLO == AtlasUpper )
   {
      for( j = 0, jaj = 0; j < N; j++, jaj += lda2 )
      {
         work[j] = t0 = ATL_rzero;

         for( i = 0, iaij = jaj; i < j; i++, iaij += (1 SHIFT) )
         {
            work[i] += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
            t0      += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
         }
         work[j] += Mabs( A[iaij] ) + t0;
      }
   }
   else
   {
      for( j = 0, jaj = 0; j < N; j++, jaj += ldap12 )
      {
         t0      = ATL_rzero;
         work[j] = Mabs( A[jaj] );
         for( i = j+1, iaij = jaj+(1 SHIFT); i < N; i++, iaij += (1 SHIFT) )
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
