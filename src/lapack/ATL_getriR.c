/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */
#include "atlas_lapack.h"
#include "atlas_lvl3.h"
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))

static void trcpzeroU(const int M, const int N, TYPE *U, const int ldu,
                      TYPE *C, const int ldc)
/*
 * Copies an upper row-major array from U, zeroing U, U is unit, so
 * diagonal is not copied
 */
{
   const int ldu2 = ldu SHIFT, ldc2 = ldc SHIFT, N2 = N SHIFT;
   int i, j;

   for (i=0; i != M; i++)
   {
      for (j=(i+1)SHIFT; j < N2; j++)
      {
         C[j] = U[j];
         U[j] = ATL_rzero;
      }
      C += ldc2;
      U += ldu2;
   }
}
int ATL_getriR(const int N, TYPE *A, const int lda, const int *ipiv,
               TYPE *wrk, const int lwrk)
{
   int jb, nb, I, ndown, iret;
   const int lda2 = lda SHIFT;
   #ifdef TREAL
      const TYPE one=ATL_rone, none=ATL_rnone;
   #else
      const TYPE one[2]={ATL_rone,ATL_rzero}, none[2]={ATL_rnone, ATL_rzero};
   #endif

   iret = ATL_trtri(CblasRowMajor, CblasLower, CblasNonUnit, N, A, lda);
   if (!iret && N > 1)
   {
/*
 *    Find largest NB we can use with our provided workspace
 */
    jb = lwrk / N;
    nb = Mjoin(PATL,laGetB)(N, 0, N, 0);
    if (jb < nb)
    {
       if (jb >= ATL_VWipgen_100LCMMN)
          nb = (jb/ATL_VWipgen_100LCMMN)*ATL_VWipgen_100LCMMN;
       else if (jb >= 4)
          nb = 4;
       else
          nb = jb;
    }
    if (!nb) return(-6);  /* need at least 1 row of workspace */
/*
 *    Only first iteration will have partial block, unroll it
 */
      jb = N - (N/nb)*nb;
      if (!jb) jb = nb;
      I = N - jb;
      A += lda2*I;
      trcpzeroU(jb, jb, A+(I SHIFT), lda, wrk, jb);
      cblas_trsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
                 jb, N, one, wrk, jb, A, lda);
      if (I)
      {
         do
         {
            I -= nb;
            A -= nb*lda2;
            ndown = N-I;
            trcpzeroU(nb, ndown, A+(I SHIFT), lda, wrk, ndown);
            cblas_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nb, N,
                       ndown-nb, none, wrk+(nb SHIFT), ndown, A+nb*lda2, lda,
                       one, A, lda);
            cblas_trsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans,
                       CblasUnit, nb, N, one, wrk, ndown, A, lda);
         }
         while(I);
      }
/*
 *    Apply row interchanges
 */
      for (I=N-2; I >= 0; I--)
      {
         jb = ipiv[I];
         if (jb != I) cblas_swap(N, A+I*lda2, 1, A+jb*lda2, 1);
      }
   }
   return(iret);
}
