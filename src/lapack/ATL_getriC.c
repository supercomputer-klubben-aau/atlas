/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */
#include "atlas_lapack.h"
#include "atlas_lvl3.h"
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))

static void trcpzeroL(const int M, const int N, TYPE *L, const int ldl,
                      TYPE *C, const int ldc)
/*
 * Copies lower triangle from L, replacing with zeros
 */
{
   const int M2 = M SHIFT, ldl2 = ldl SHIFT, ldc2 = ldc SHIFT;
   int i, j;

   for (j=0; j != N; j++)
   {
      for (i=((j+1)SHIFT); i < M2; i++)
      {
         C[i] = L[i];
         L[i] = ATL_rzero;
      }
      C += ldc2;
      L += ldl2;
   }
}

int ATL_getriC(const int N, TYPE *A, const int lda, const int *ipiv,
               TYPE *wrk, const int lwrk)
{
   const int lda2 = lda SHIFT;
   int J, jb, nb, nright, iret;
   TYPE *A0 = A;
   #ifdef TREAL
      const TYPE one=ATL_rone, none=ATL_rnone;
   #else
      const TYPE one[2]={ATL_rone,ATL_rzero}, none[2]={ATL_rnone, ATL_rzero};
   #endif

   iret = ATL_trtri(CblasColMajor, CblasUpper, CblasNonUnit, N, A, lda);
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
      if (!nb) return(-6);  /* need at least 1 col of workspace */
/*
 *    Only first iteration will have partial block, unroll it
 */
      jb = N - (N/nb)*nb;
      if (!jb) jb = nb;
      J = N - jb;
      A += lda2*J;
      trcpzeroL(jb, jb, A+(J SHIFT), lda, wrk, jb);
      cblas_trsm(CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit,
                 N, jb, one, wrk, jb, A, lda);
      if (J)
      {
         do
         {
            J -= nb;
            A -= nb*lda2;
            nright = N-J;
            trcpzeroL(nright, nb, A+(J SHIFT), lda, wrk, nright);
            cblas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, nb,
                       nright-nb, none, A+nb*lda2, lda, wrk+(nb SHIFT), nright,
                       one, A, lda);
            cblas_trsm(CblasColMajor, CblasRight, CblasLower, CblasNoTrans,
                       CblasUnit, N, nb, one, wrk, nright, A, lda);
         }
         while(J);
      }
/*
 *    Apply column interchanges
 */
      for (J=N-2; J >= 0; J--)
      {
         jb = ipiv[J];
         if (jb != J) cblas_swap(N, A+J*lda2, 1, A+jb*lda2, 1);
      }
   }
   return(iret);
}
