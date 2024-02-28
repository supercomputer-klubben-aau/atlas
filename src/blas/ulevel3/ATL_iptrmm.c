/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Majedul Sujon
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_lvl3.h"

#ifdef Right_
   #include Mstr(Mjoin(ATLAS_PRE,trmmRU_view.h))
   #define INCR 1
   #define ATL_iptrmm Mjoin(PATL,iptrmmR)
   #define ATL_trmm_ST Mjoin(PATL,trmmR_ST)
   #define ATL_trmm_alloc Mjoin(PATL,trmmR_alloc)
   #define ATL_tminfoLN Mjoin(PATL,tminfoR_LN)
   #define ATL_tminfoLT Mjoin(PATL,tminfoR_LT)
   #define ATL_tminfoUN Mjoin(PATL,tminfoR_UN)
   #define ATL_tminfoUT Mjoin(PATL,tminfoR_UT)
   #ifdef TCPLX
      #define ATL_tminfoLC Mjoin(PATL,tminfoR_LC)
      #define ATL_tminfoUC Mjoin(PATL,tminfoR_UC)
   #endif
   #define ATL_trmm_LNUT Mjoin(PATL,trmmR_LNUT)
   #define ATL_trmm_LTUN Mjoin(PATL,trmmR_LTUN)
#else
   #include Mstr(Mjoin(ATLAS_PRE,trmmLU_view.h))
   #define INCR ldb
   #define ATL_iptrmm Mjoin(PATL,iptrmmL)
   #define ATL_trmm_ST Mjoin(PATL,trmmL_ST)
   #define ATL_trmm_alloc Mjoin(PATL,trmmL_alloc)
   #define ATL_tminfoLN Mjoin(PATL,tminfoL_LN)
   #define ATL_tminfoLT Mjoin(PATL,tminfoL_LT)
   #define ATL_tminfoUN Mjoin(PATL,tminfoL_UN)
   #define ATL_tminfoUT Mjoin(PATL,tminfoL_UT)
   #ifdef TCPLX
      #define ATL_tminfoLC Mjoin(PATL,tminfoL_LC)
      #define ATL_tminfoUC Mjoin(PATL,tminfoL_UC)
   #endif
   #define ATL_trmm_LNUT Mjoin(PATL,trmmL_LNUT)
   #define ATL_trmm_LTUN Mjoin(PATL,trmmL_LTUN)
#endif
int ATL_iptrmm
(
   ATL_UINT bv,  /* 0:Right, 1:Upper, 2:TransA, 3: Conj, 4:NonUnit */
   ATL_CSZT               NT,  /* number of diagonals in triangle */
   ATL_CSZT               NR,  /* number of right hand sides */
   const SCALAR           alpha,
   const TYPE             *A,
   ATL_CSZT               lda,
   TYPE                   *B,
   ATL_CSZT               ldb
)
{
/*
 * Handle degerated cases
 */
   if (NT == 1)
   {
      if (bv&16)   /* if non-unit, work to do, else nothing to do at M=1! */
      #ifdef TCPLX
      {
         TYPE scl[2], tmp;
         *scl = *A;
         tmp = A[1];
         scl[1] = (bv&8) ? -tmp : tmp;
         if (!SCALAR_IS_ONE(alpha))
         {
            TYPE ra, ia;
            ra = scl[0] * alpha[0] - scl[1] * alpha[1];
            ia = scl[1] * alpha[0] + scl[0] * alpha[1];
            scl[0] = ra; scl[1] = ia;
         }
         Mjoin(PATL,scal)(NR, scl, B, INCR);
      }
      #else
         Mjoin(PATL,scal)(NR, alpha * *A, B, INCR);
      #endif
      else if (!SCALAR_IS_ONE(alpha))  /* unit may still need to scale! */
         Mjoin(PATL,scal)(NR, alpha, B, INCR);
      return(0);
   }
/*
 * if NT (KB) is not greater than MAX_KB, apply trmm_ST code
 */
#ifdef Right_
   else if (NT <= ATL_VIEW_MAX_KB && NR <= ATL_VIEW_MAX_MB)
#else
   else if (NT <= ATL_VIEW_MAX_KB && NR <= ATL_VIEW_MAX_NB)
#endif
      return(ATL_trmm_ST(bv, NT, NR, alpha, A, lda, B, ldb));
/*
 * NT (KB) is large, so apply outer loop for trmm
 */
   else
   {
      int res;
      unsigned int ctu; /* Conj Trans Upper*/
      ATL_SZT Tsz, Rsz, Csz; /* workspace size for TRMM */
      TYPE *BW, *L, *R, *w;
      void *vp=NULL;
      ipinfo_t gip;
      tminfo_t si;
   #ifdef TCPLX
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
      TYPE ZERO[2] = {ATL_rzero, ATL_rzero};
   #else
      #define ONE ATL_rone
      #define ZERO ATL_rzero
   #endif
      enum ATLAS_TRANS TA;
      #ifdef TCPLX
         if ((bv&12)==12)
            TA = AtlasConjTrans;
         else
      #endif
      TA = (bv&4) ? AtlasTrans:AtlasNoTrans;
/*
 *    Get ipinfo for the GEMM
 */
      if (bv & 1) /* Right*/
         Mjoin(PATL,ipnekInfo)(&gip, AtlasNoTrans, TA, NR, NT, NT, ldb, lda, ldb,
                               alpha, ZERO);
      else  /* Left */
         Mjoin(PATL,ipmekInfo)(&gip, TA, AtlasNoTrans, NT,NR,NT, lda, ldb, ldb,
                               alpha, ZERO);
/*
 *    select tminfo and one of four variation of outer loops
 */
      void (*tminfo) (tminfo_t *ip, ipinfo_t *gip, const enum ATLAS_DIAG Diag,
            ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda,
            ATL_CSZT ldb);
      int (*ammtrmm)(ipinfo_t *ip, tminfo_t *si,
            ATL_CSZT M, ATL_CSZT N, const SCALAR alpha, const TYPE *A,
            ATL_CSZT lda, TYPE *X, ATL_CSZT ldx, ATL_CSZT Tsz, ATL_CSZT Rsz,
            ATL_CSZT Csz, TYPE *tbw, TYPE *L, TYPE *R, TYPE *w );

      ctu = (bv>>1)&7; /* masked three bits responsible for Conj Trans Upper*/
      if (!ctu)
         tminfo = ATL_tminfoLN;
      else if (ctu == 2)
         tminfo = ATL_tminfoLT;
      #ifdef TCPLX
         else if (ctu == 6)
            tminfo = ATL_tminfoLC;
      #endif
      else if (ctu == 1)
         tminfo = ATL_tminfoUN;
      else if (ctu == 3)
         tminfo = ATL_tminfoUT;
      #ifdef TCPLX
      else if (ctu == 7)
         tminfo = ATL_tminfoUC;
      #endif
      else /* unknown trmm case! */
           ATL_assert(0);
      if (!ctu || ctu==3 || ctu == 7)
         ammtrmm = ATL_trmm_LNUT;
      else
         ammtrmm = ATL_trmm_LTUN;
/*
 *    Get TRMM info
 */
      if (bv&1)
         tminfo(&si, &gip, (bv&16)?AtlasNonUnit:AtlasUnit, NR, NT, alpha, lda,
               ldb);
      else
         tminfo(&si, &gip, (bv&16)?AtlasNonUnit:AtlasUnit, NT, NR, alpha, lda,
               ldb);
/*
 *    allocate workspace for trmm outer loop
 */
      vp = ATL_trmm_alloc(&gip, &si, NT, &Tsz, &Rsz, &Csz, &BW, &L, &R, &w);
      if (!vp)
         return(1); /* not enough memory, reucrse again  */
/*
 *    apply trmm outer-loop
 */
      if (bv & 1)
         res = ammtrmm(&gip, &si, NR, NT, alpha, A, lda, B, ldb, Tsz, Rsz, Csz,
               BW, L, R, w);
      else
         res = ammtrmm(&gip, &si, NT, NR, alpha, A, lda, B, ldb, Tsz, Rsz, Csz,
               BW, L, R, w);

      free(vp);
      return(res);
   }
   return(0);
}
