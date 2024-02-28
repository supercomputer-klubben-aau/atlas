/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#define SCPLX
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"

int clapack_cgetri(const enum CBLAS_ORDER Order, const int N, void *A,
                   const int lda, const int *ipiv)
{
   int ierr=0, lwrk;
   void *vp;

   lwrk = Mjoin(PATL,laGetB)(N, 0, N, 0);
   if (lwrk <= N) lwrk *= N;
   else lwrk = N*N;
   vp = malloc(ATL_Cachelen + ATL_MulBySize(lwrk));
   if (vp)
   {
      ierr = ATL_getri(Order, N, A, lda, ipiv, ATL_AlignPtr(vp), &lwrk);
      free(vp);
   }
   else
   {
      cblas_xerbla(7, "clapack_cgetri",
                   "Cannot allocate workspace of %d\n", lwrk);
      return(-7);
   }
   return(ierr);
}
