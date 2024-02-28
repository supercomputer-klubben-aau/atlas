/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 Peter Soendergaard
 */
#include "atlas_lapack.h"
#include "atlas_lvl3.h"

int ATL_trtri(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
	      const enum ATLAS_DIAG Diag, const int N, TYPE *A, const int lda)
{
   const int ldap1 = (lda+1)SHIFT;
   int i;

   if (N > 0)
   {
/*
 *    Check for singularity if nonunit
 */
      if (Diag == AtlasNonUnit)
      {
         for (i=0; i != N; i++, A += ldap1)
         {
            #ifdef TREAL
               if (*A == ATL_rzero) return(i+1);
            #else
               if (*A == ATL_rzero && A[1] == ATL_rzero) return(i+1);
            #endif
         }
         A -= N*ldap1;
      }
      if (Uplo == AtlasUpper)
      {
         if (Order == AtlasColMajor) return(ATL_trtriCU(Diag, N, A, lda));
         else return(ATL_trtriRU(Diag, N, A, lda));
      }
      else
      {
         if (Order == AtlasColMajor) return(ATL_trtriCL(Diag, N, A, lda));
         else return(ATL_trtriRL(Diag, N, A, lda));
      }
   }
   return(0);
}
