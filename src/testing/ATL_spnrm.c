/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,spnrm)
(const enum ATLAS_UPLO UPLO, const int N, const TYPE *A)
{
   int                        i, iaij, j;
   TYPE max=ATL_rzero, t0, * work= NULL;

   if( N <= 0 ) return( ATL_rzero );

   work = (TYPE *)malloc( N * sizeof( TYPE ) );
   if( work == NULL )
   {fprintf( stderr, "mem alloc failed in [sp,hp]nrm, bye ...\n" ); exit( 1 );}
   else { for( i = 0; i < N; i++ ) work[i] = ATL_rzero; }

   if( UPLO == AtlasUpper )
   {
      for( j = 0, iaij = 0; j < N; j++ )
      {
         work[j] = t0 = ATL_rzero;

         for( i = 0; i < j; i++, iaij += (1 SHIFT) )
         {
#ifdef TREAL
            work[i] += Mabs( A[iaij] );
            t0      += Mabs( A[iaij] );
#else
            work[i] += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
            t0      += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
#endif
         }
#ifdef TREAL
         work[j] += Mabs( A[iaij] ) + t0;
#else
         work[j] += Mabs( A[iaij] ) + Mabs( A[iaij+1] ) + t0;
#endif
         iaij    += (1 SHIFT);
      }
   }
   else
   {
      for( j = 0, iaij = 0; j < N; j++ )
      {
         t0      = ATL_rzero;
#ifdef TREAL
         work[j] = Mabs( A[iaij] );
#else
         work[j] = Mabs( A[iaij] ) + Mabs( A[iaij+1] );
#endif

         iaij    += (1 SHIFT);
         for( i = j+1; i < N; i++, iaij += (1 SHIFT) )
         {
#ifdef TREAL
            work[i] += Mabs( A[iaij] );
            t0      += Mabs( A[iaij] );
#else
            work[i] += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
            t0      += Mabs( A[iaij] ) + Mabs( A[iaij+1] );
#endif
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
