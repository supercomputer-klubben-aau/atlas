/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2003 R. Clint Whaley
 */
/*
 * HERK actually uses real gemm, so use real blocking factors
 */
#ifdef SCPLX
   #include "atlas_samm_sum.h"
#elif defined(DCPLX)
   #include "atlas_damm_sum.h"
#endif
#include "atlas_misc.h"
#define ATL_NOAMM 1
#include "atlas_lvl3.h"
#undef ATL_NOAMM
#include "atlas_level3.h"
#include "atlas_level1.h"
#include "atlas_lapack.h"
#include <math.h>

#define ATL_potrfRL Mjoin(PATL,potrfRL)
int ATL_potrfRL(const int N, TYPE *A, const int lda)
{
   TYPE *An, *Ar;
   int Nleft, Nright, ierr;
   static const TYPE ONE[2] = {ATL_rone, ATL_rzero};
   const size_t lda2=lda+lda;

   if (N > 1)
   {
      Nleft = N >> 1;
      #ifdef ATL_VWipgen_100LCMMN
         if (Nleft > ATL_VWipgen_100LCMMN<<1)
            Nleft = (Nleft/ATL_VWipgen_100LCMMN)*ATL_VWipgen_100LCMMN;
      #endif
      Nright = N - Nleft;
      ierr = ATL_potrfRL(Nleft, A, lda);
      if (!ierr)
      {
         Ar = A + Nleft * lda2;
         An = Ar + Nleft+Nleft;
         cblas_trsm(CblasRowMajor, CblasRight, CblasLower, CblasConjTrans,
                    CblasNonUnit, Nright, Nleft, ONE, A, lda, Ar, lda);
         cblas_herk(CblasRowMajor, CblasLower, CblasNoTrans, Nright, Nleft,
                    ATL_rnone, Ar, lda, ATL_rone, An, lda);
         ierr = ATL_potrfRL(Nright, An, lda);
         if (ierr) return(ierr+Nleft);
      }
      else return(ierr);
   }
   else if (N == 1)
   {
      if (*A > ATL_rzero)
      {
         *A = sqrt(*A);
         A[1] = ATL_rzero;
      }
      else return(1);
   }
   return(0);
}
