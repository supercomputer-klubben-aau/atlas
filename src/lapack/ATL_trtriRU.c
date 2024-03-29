/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 Peter Soendergaard
 */
#include "atlas_lapack.h"
#include "atlas_lvl3.h"

#define REAL_RECURSE_LIMIT 4

#ifdef TREAL

#define SAFE_INVERT(A_) *(A_) = ATL_rone / *(A_);
#define TRI_MUL(A_, B_) \
{ \
     *(B_) = (*(A_)) * (*(B_)); \
}
#define TRI_NEG(A_) \
{ \
     *(A_) = - *(A_); \
}
#else

#define SAFE_INVERT(A_) Mjoin(PATL,cplxinvert)(1, A_, 1, A_, 1);
#define TRI_MUL(A_, B_) \
{ \
    (B_)[0] = ((A_)[0])*((B_)[0])-((A_)[1])*((B_)[1]); \
    (B_)[1] = ((A_)[1])*((B_)[0])+((A_)[0])*((B_)[1]); \
}
#define TRI_NEG(A_) \
{ \
     (A_)[0] = - (A_)[0]; \
     (A_)[1] = - (A_)[1]; \
}
#endif

static int ATL_trtriRU_4(const enum ATLAS_DIAG Diag, TYPE *A, const int lda)
{
    TYPE *pA0=A, *pA1=A+lda, *pA2=A+2*lda, *pA3=A+3*lda;
    TYPE A01=pA0[1], A02=pA0[2], A03=pA0[3];
    TYPE             A12=pA1[2], A13=pA1[3];
    TYPE                         A23=pA2[3];
    TYPE tmp;

    if (Diag == AtlasNonUnit)
    {
       SAFE_INVERT(pA0);
       SAFE_INVERT(pA1+1);
       SAFE_INVERT(pA2+2);
       SAFE_INVERT(pA3+3);

       pA0[1] = -A01*pA1[1]*pA0[0];
       pA1[2] = -A12*pA2[2]*pA1[1];
       pA2[3] = -A23*pA3[3]*pA2[2];

       pA0[2] = -(A01*pA1[2]+A02*pA2[2])*pA0[0];
       pA1[3] = -(A12*pA2[3]+A13*pA3[3])*pA1[1];

       pA0[3] = -(A01*pA1[3]+A02*pA2[3]+A03*pA3[3])*pA0[0];
    }
    else
    {
       pA0[1] = -A01;
       pA1[2] = -A12;
       pA2[3] = -A23;

       pA0[2] = -(A01*pA1[2]+A02);
       pA1[3] = -(A12*pA2[3]+A13);

       pA0[3] = -(A01*pA1[3]+A02*pA2[3]+A03);
    }
    return(0);
}


static int ATL_trtriRU_3(const enum ATLAS_DIAG Diag, TYPE *A, const int lda)
{

    TYPE *pA0=A, *pA1=A+lda, *pA2=A+2*lda;
    TYPE A01=pA0[1], A02=pA0[2];
    TYPE             A12=pA1[2];

    TYPE *B01 = pA0+1;
    TYPE *B02 = pA0+2;
    TYPE *B12 = pA1+2;

    TYPE tmp;

    if (Diag == AtlasNonUnit)
    {
       SAFE_INVERT(pA0);
       SAFE_INVERT(pA1+1);
       SAFE_INVERT(pA2+2);
       *B01 = -A01*pA1[1]*pA0[0];
       *B12 = -A12*pA2[2]*pA1[1];
       *B02 = -(A01*(*B12)+A02*pA2[2])*pA0[0];
    }
    else
    {
       *B01 = -A01;
       *B12 = -A12;
       *B02 = -(A01*(*B12)+A02);
    }
    return(0);
}

int ATL_trtriRU(const enum ATLAS_DIAG Diag, const int N, TYPE *A, const int lda)
{
  int ierr = 0;

   TYPE *Age, *Atr;
   TYPE tmp;
   int Nleft, Nright;
   #ifdef TREAL
      #define one ATL_rone
      #define mone -ATL_rone
      #define none ATL_rnone
   #else
      static const TYPE one[2] = {ATL_rone, ATL_rzero};
      static const TYPE mone[2] = {-ATL_rone, ATL_rzero};
      static const TYPE none[2] = {ATL_rnone, ATL_rzero};
   #endif

#ifdef TREAL
   if (N > REAL_RECURSE_LIMIT)
#else
   if (N > 1)
#endif
   {
      Nleft = N >> 1;
      #ifdef ATL_VWipgen_100LCMMN
         if (Nleft > ATL_VWipgen_100LCMMN)
            Nleft = (Nleft/ATL_VWipgen_100LCMMN)*ATL_VWipgen_100LCMMN;
      #endif
      Nright = N - Nleft;

      Age = A + (Nleft SHIFT);
      Atr = A + (Nleft * (lda+1) SHIFT);

      cblas_trsm(AtlasRowMajor, AtlasRight, AtlasUpper, AtlasNoTrans, Diag,
                  Nleft, Nright, one, Atr, lda, Age, lda);

      cblas_trsm(AtlasRowMajor, AtlasLeft, AtlasUpper, AtlasNoTrans, Diag,
                  Nleft, Nright, mone, A, lda, Age, lda);

      ierr = ATL_trtriRU(Diag, Nleft, A, lda);
      if (ierr!=0) return(ierr);
      ierr = ATL_trtriRU(Diag, Nright, Atr, lda);
      if (ierr!=0) return(ierr+Nleft);
   }
   else
   {
#ifdef TREAL
     if (N==4) return(ATL_trtriRU_4(Diag,A,lda));
     else if (N==3) return(ATL_trtriRU_3(Diag,A,lda));
     else
       if ( N == 2)
       {
         if (Diag == AtlasNonUnit)
         {
            SAFE_INVERT(A);
            SAFE_INVERT(A+((lda+1) SHIFT));
            TRI_MUL(A,A+(1 SHIFT));
            TRI_MUL(A+((lda+1) SHIFT),A+(1 SHIFT));
         }
         TRI_NEG(A+(1 SHIFT));

       }
       else
#endif
     {
       if (Diag == AtlasNonUnit) SAFE_INVERT(A);
     }
   }

   return(ierr);

}

