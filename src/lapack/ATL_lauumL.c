/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */
#include "atlas_lapack.h"
#ifdef TREAL
   #define my_syrk cblas_syrk
   #define my_trans CblasTrans
#else
   #define my_syrk cblas_herk
   #define my_trans CblasConjTrans
#endif

#ifdef RowMajor_
   #define MyOrder CblasRowMajor
   #define ATL_lauumU Mjoin(PATL,lauumRU)
   #define ATL_lauumL Mjoin(PATL,lauumRL)
#else
   #define MyOrder CblasColMajor
   #define ATL_lauumU Mjoin(PATL,lauumCU)
   #define ATL_lauumL Mjoin(PATL,lauumCL)
#endif

void ATL_lauumL(const int N, TYPE *A, const int lda)
{
   int Nleft, Nright;
   #ifdef TREAL
      const TYPE one=ATL_rone;
   #else
      const TYPE one[2]={ATL_rone, ATL_rzero};
   #endif
   TYPE *G, *U0=A, *U1;

   if (N > 1)
   {
      Nleft = N >> 1;
      #ifdef ATL_VWipgen_100LCMMN
         if (Nleft > ATL_VWipgen_100LCMMN)
            Nleft = (Nleft/ATL_VWipgen_100LCMMN)*ATL_VWipgen_100LCMMN;
      #endif
      Nright = N - Nleft;
      #ifdef RowMajor_
         G = A + Nleft*(lda SHIFT);
         U1 = G + (Nleft SHIFT);
      #else
         G = A + (Nleft SHIFT);
         U1 = G + Nleft*(lda SHIFT);
      #endif
      ATL_lauumL(Nleft, U0, lda);
      my_syrk(MyOrder, CblasLower, my_trans, Nleft, Nright, ATL_rone,
              G, lda, ATL_rone, U0, lda);
      cblas_trmm(MyOrder, CblasLeft, CblasLower, my_trans, CblasNonUnit,
                 Nright, Nleft, one, U1, lda, G, lda);
      ATL_lauumL(Nright, U1, lda);
   }
   else *A = *A * *A;
}
