/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"

TYPE Mjoin(PATL,gbnrm1)(const int M, const int N, const int KL, const int KU,
const TYPE *A, const int LDA)
/*
 * Calculates the 1-norm of a general band rectangular matrix
 */
{
   int i, i0, i1, iaij, j, jaj, k, lda2 = ( LDA SHIFT );
   TYPE max=ATL_rzero, t0;

   for( j = 0, jaj = 0; j < N; j++, jaj += lda2 )
   {
      k  = KU - j;
      i0 = ( j - KU > 0 ? j - KU : 0 );
      i1 = ( M - 1 > j + KL ? j + KL : M - 1 );

      t0 =  ATL_rzero;
      for( i = i0, iaij = ((k+i0) SHIFT)+jaj; i <= i1; i++, iaij += (1 SHIFT) )
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
   return(max);
}
