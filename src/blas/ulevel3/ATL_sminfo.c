/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(ATLAS_UPR,amm_kern.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_a1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_a1.h))
#ifdef TCPLX
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_a1.h))
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_a1.h))
#endif
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_aN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_aN.h))
#ifdef TCPLX
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_aN.h))
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_aN.h))
#endif

#define ATL_DECL_ 1
#ifdef Right_
   #ifdef Upper_
      #error "Upper not supported"
   #else
      #ifdef Trans_
         #include Mstr(Mjoin(ATLAS_PRE,utrsmR_LT.h))
         #include Mstr(Mjoin(ATLAS_PRE,trsmRT_view.h))
         #define ATL_sminfo Mjoin(PATL,sminfoR_LT)
         #define INCA lda
         #define BV (1|4)
      #else
         #include Mstr(Mjoin(ATLAS_PRE,utrsmR_LN.h))
         #include Mstr(Mjoin(ATLAS_PRE,trsmRN_view.h))
         #define ATL_sminfo Mjoin(PATL,sminfoR_LN)
         #define INCA 1
         #define BV 1
      #endif
   #endif
#else  /* Left */
   #ifdef Upper_
      #error "Upper not supported"
   #else
      #ifdef Trans_
         #include Mstr(Mjoin(ATLAS_PRE,utrsmL_LT.h))
         #include Mstr(Mjoin(ATLAS_PRE,trsmLT_view.h))
         #define ATL_sminfo Mjoin(PATL,sminfoL_LT)
         #define INCA 1
         #define BV 4
      #else
         #include Mstr(Mjoin(ATLAS_PRE,utrsmL_LN.h))
         #include Mstr(Mjoin(ATLAS_PRE,trsmLN_view.h))
         #define ATL_sminfo Mjoin(PATL,sminfoL_LN)
         #define INCA lda
         #define BV 0
      #endif
   #endif
#endif
#undef ATL_DECL_
   /* bv: 0:Right, 1:Upper, 2:TransA, 3: Conj, 4:NonUnit */
int ATL_sminfo
   (sminfo_t *ip, ATL_CUINT bv, ATL_CSZT N, ATL_CSZT R,
    const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb)
{
   int iv=ATL_VIEW_BEST_IDX, imm, icpA, icpB, icpC, k;
   int ALLT=0;
   float spf;
   #ifndef TCPLX
      ATL_iptr_t iptmp;
   #endif
   ATL_cparr_t *cp = (ATL_cparr_t *)&(ip->b2blk);
/*
 * Really need to see if its better to do gemm+trsm in this case, but for
 * now just assume not
 */
   if (N > ATL_VIEW_BEST_KB) /* inner product case */
   {
      int nb;
      for (; iv; iv--)
      {
         nb = ATL_GetViewKB(iv);
         if (N >= nb)
            break;
      }
      if (N > nb && iv < ATL_VIEW_BEST_IDX)
         iv++;
   }
   ALLT = ATL_trsm_allT(iv);
   ip->bv = (ALLT<<5) | bv;
   ATL_GetViewInfo(iv, spf, k, ip->rb, ip->kb, imm, icpA, icpB, k,k,k,k);
   icpA += icpA;
   icpB += icpB;
   #ifdef TCPLX
      ATL_AMM_info(imm, &ip->amm_b0, &ip->amm_b1, &ip->amm_bn, ip->mmflg,
                   ip->mu, ip->nu, ip->ku, k, k, k, k);
   #else
      ATL_AMM_info(imm, &ip->amm_b0, &ip->amm_b1, &iptmp, ip->mmflg,
                   ip->mu, ip->nu, ip->ku, k, k, k, k);
   #endif
   ip->utrsm = ATL_findutrsm(ALLT, ip->mu, ip->nu);
   if (ALLT)
   {
      k = ip->mu;
      ip->mu = ip->nu;
      ip->nu = k;
   }
   #ifdef Right_  /* LT_UN or LN_UT case */
      *cp = (ALLT) ? ATL_CpyFromATg_a1[icpB]:ATL_CpyFromATg_a1[icpA];
   #else
      *cp = (ALLT) ? ATL_CpyFromANg_a1[icpA]:ATL_CpyFromANg_a1[icpB];
   #endif
   cp = (ATL_cparr_t *)&(ip->a2blk);
   #ifdef Trans_  /* LT or UN case */
      if (bv&2) /* UN case */
      {
         #ifdef Right_
            ip->incA = 1 SHIFT;
         #else
            ip->incA = lda SHIFT;
         #endif
         #ifdef TCPLX
            if ((bv&8))  /* UC case! */
               #ifdef Right_
                  *cp  = (ALLT)?ATL_CpyFromAHg_aN[icpA]:ATL_CpyFromAHg_aN[icpB];
               #else
                  *cp  = (ALLT)?ATL_CpyFromACg_aN[icpB]:ATL_CpyFromACg_aN[icpA];
               #endif
            else
         #endif
         #ifdef Right_
            *cp = (ALLT) ? ATL_CpyFromANg_aN[icpA]:ATL_CpyFromANg_aN[icpB];
         #else
            *cp = (ALLT) ? ATL_CpyFromATg_aN[icpB]:ATL_CpyFromATg_aN[icpA];
         #endif
      }
      else      /* LT case */
      {
         #ifdef Right_
            ip->incA = lda SHIFT;
         #else
            ip->incA = 1 SHIFT;
         #endif
         #ifdef TCPLX
            if (bv&8) /* LH case! */
               #ifdef Right_
                  *cp  = (ALLT)?ATL_CpyFromAHg_aN[icpA]:ATL_CpyFromAHg_aN[icpB];
               #else
                  *cp  = (ALLT)?ATL_CpyFromACg_aN[icpB]:ATL_CpyFromACg_aN[icpA];
               #endif
            else
         #endif
         #ifdef Right_
            *cp = (ALLT) ? ATL_CpyFromATg_aN[icpA]:ATL_CpyFromATg_aN[icpB];
         #else
            *cp = (ALLT) ? ATL_CpyFromANg_aN[icpB]:ATL_CpyFromANg_aN[icpA];
         #endif
      }
   #else /* LN or UT case, these cases use incA backwards (negative) */
      if (bv&2)  /* UT case */
      {
         #ifdef Right_
            ip->incA = lda SHIFT;
         #else
            ip->incA = 1 SHIFT;
         #endif
         #ifdef TCPLX
            if (bv&8) /* UH case! */
               #ifdef Right_
                  *cp  = (ALLT)?ATL_CpyFromAHg_aN[icpA]:ATL_CpyFromAHg_aN[icpB];
               #else
                  *cp  = (ALLT)?ATL_CpyFromACg_aN[icpB]:ATL_CpyFromACg_aN[icpA];
               #endif
            else
         #endif
         #ifdef Right_
            *cp = (ALLT) ? ATL_CpyFromATg_aN[icpA]:ATL_CpyFromATg_aN[icpB];
         #else
            *cp = (ALLT) ? ATL_CpyFromANg_aN[icpB]:ATL_CpyFromANg_aN[icpA];
         #endif
      }
      else /* LN case */
      {
         #ifdef Right_
            ip->incA = 1 SHIFT;
         #else
            ip->incA = lda SHIFT;
         #endif
         #ifdef TCPLX
            if (bv&8) /* LC case! */
               #ifdef Right_
                  *cp  = (ALLT)?ATL_CpyFromAHg_aN[icpA]:ATL_CpyFromAHg_aN[icpB];
               #else
                  *cp  = (ALLT)?ATL_CpyFromACg_aN[icpB]:ATL_CpyFromACg_aN[icpA];
               #endif
            else
         #endif
         #ifdef Right_
            *cp = (ALLT) ? ATL_CpyFromANg_aN[icpA]:ATL_CpyFromANg_aN[icpB];
         #else
            *cp = (ALLT) ? ATL_CpyFromATg_aN[icpB]:ATL_CpyFromATg_aN[icpA];
         #endif
      }
   #endif
   return(iv);
}
