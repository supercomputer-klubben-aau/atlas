/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */
#include "atlas_lapack.h"
#include "atlas_lvl3.h"
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))
int ATL_getri(const enum CBLAS_ORDER Order, const int N, TYPE *A, const int lda,
              const int *ipiv, TYPE *wrk, int *lwrk)
{
   int ierr=0;
   if (*lwrk != -1)
   {
      if (Order == AtlasRowMajor)
         ierr = ATL_getriR(N, A, lda, ipiv, wrk, *lwrk);
      else ierr = ATL_getriC(N, A, lda, ipiv, wrk, *lwrk);
   }
   else *lwrk = N*ATL_VWipgen_98KB;
   return(ierr);
}
