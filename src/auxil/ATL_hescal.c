/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"

void Mjoin( PATL, hescal )
(
   const enum ATLAS_UPLO      UPLO,
   const int                  M,
   const int                  N,
   const TYPE                 ALPHA,
   TYPE                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * ATL_hescal  scales a (trapezoidal)  Hermitian  m-by-n matrix A by the
 * real scalar alpha.  The imaginary parts of the diagonal elements of A
 * need not be set on input,  they are assumed to be zero,  and  on exit
 * they are set to zero.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   int                        i, incA, j, mn;
   TYPE                       * a_j;
/* ..
 * .. Executable Statements ..
 *
 */
   if( UPLO == AtlasLower )
   {
      incA = ( ( LDA - M + 1 ) SHIFT ); mn = Mmin( M, N );

      if(      ALPHA == ATL_rzero )
      {
         for( j = 0; j < mn; j++ )
         {
            for( i = j; i < M; i++, A += 2 ) { *A = A[1] = ATL_rzero; }
            A += incA;
            incA += 2;
         }
      }
      else if( ALPHA != ATL_rone )
      {
         for( j = 0; j < mn; j++ )
         {
            *A *= ALPHA; A[1] = ATL_rzero; A += 2;

            for( i = j+1; i < M; i++, A += 2 ) { *A *= ALPHA; A[1] *= ALPHA; }
            A += incA;
            incA += 2;
         }
      }
   }
   else
   {
      incA = ( LDA SHIFT );

      if(      ALPHA == ATL_rzero )
      {
         for( j = 0, mn = M - N; j < N; j++, mn++, A += incA )
         {
            a_j = A;
            for( i = 0; i <= mn; i++, a_j += 2 )
            {
               *a_j = a_j[1] = ATL_rzero;
            }
         }
      }
      else if( ALPHA != ATL_rone )
      {
         for( j = 0, mn = M - N; j < N; j++, mn++, A += incA )
         {
            a_j = A;
            for( i = 0; i < mn; i++, a_j += 2 )
            {
               *a_j   *= ALPHA;
               a_j[1] *= ALPHA;
            }
            *a_j *= ALPHA; a_j[1] = ATL_rzero;
         }
      }
   }
}
