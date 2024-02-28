/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Majedul Sujon
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_lvl3.h"

#include Mstr(Mjoin(ATLAS_UPR,amm_kern.h))

/* cases: LN, UT, LT, UN ... LC/UC included */
#ifdef Right_
   #ifdef Upper_
      #ifdef Trans_
#include Mstr(Mjoin(ATLAS_PRE,utrmmR_LN.h)) /* LN & UT share same header file*/
         #ifdef Conj_
            #define ATL_tminfo Mjoin(PATL,tminfoR_UC)
         #else
            #define ATL_tminfo Mjoin(PATL,tminfoR_UT)
         #endif
      #else
#include Mstr(Mjoin(ATLAS_PRE,utrmmR_LT.h)) /* LT & UN share same header file*/
         #define ATL_tminfo Mjoin(PATL,tminfoR_UN)
      #endif
   #else /* lower */
      #ifdef Trans_
#include Mstr(Mjoin(ATLAS_PRE,utrmmR_LT.h))
         #ifdef Conj_
            #define ATL_tminfo Mjoin(PATL,tminfoR_LC)
         #else
            #define ATL_tminfo Mjoin(PATL,tminfoR_LT)
         #endif
      #else
#include Mstr(Mjoin(ATLAS_PRE,utrmmR_LN.h))
         #define ATL_tminfo Mjoin(PATL,tminfoR_LN)
      #endif
   #endif
#else /* left side */
   #ifdef Upper_
      #ifdef Trans_
#include Mstr(Mjoin(ATLAS_PRE,utrmmL_LN.h)) /* LN & UT share same header file*/
         #ifdef Conj_
            #define ATL_tminfo Mjoin(PATL,tminfoL_UC)
         #else
            #define ATL_tminfo Mjoin(PATL,tminfoL_UT)
         #endif
      #else
#include Mstr(Mjoin(ATLAS_PRE,utrmmL_LT.h)) /* LT & UN share same header file*/
         #define ATL_tminfo Mjoin(PATL,tminfoL_UN)
      #endif
   #else /* lower */
      #ifdef Trans_
#include Mstr(Mjoin(ATLAS_PRE,utrmmL_LT.h))
         #ifdef Conj_
            #define ATL_tminfo Mjoin(PATL,tminfoL_LC)
         #else
            #define ATL_tminfo Mjoin(PATL,tminfoL_LT)
         #endif
      #else
#include Mstr(Mjoin(ATLAS_PRE,utrmmL_LN.h))
         #define ATL_tminfo Mjoin(PATL,tminfoL_LN)
      #endif
   #endif
#endif

void ATL_tminfo
(
   tminfo_t *ip,
   ipinfo_t *gip,
   const enum ATLAS_DIAG Diag,
   ATL_CSZT M,
   ATL_CSZT N,
   const SCALAR alpha,
   ATL_CSZT lda,
   ATL_CSZT ldb
)
{
#ifdef TCPLX
   TYPE *alpL, *alpR, *alpC;
#else
   TYPE alpL, alpR, alpC;
#endif
   tcm2am_t *t2a = (tcm2am_t *)&(ip->t2blk);
   cm2am_t *r2a = (cm2am_t *)&(ip->r2blk);
   ablk2cmat_t *blk2c = (ablk2cmat_t *)&(ip->blk2c);

   ip->mu = ATL_TRMMK_MU;
   ip->nu = ATL_TRMMK_NU;
   ip->ku = ATL_TRMMK_KU;
   ip->vlen = ATL_TRMMK_VLEN;
   ip->kvec = ATL_TRMMK_KVEC;
   ip->kb = ATL_TRMM_KB;
   ip->flg = 0;
   #ifdef Right_
      #ifdef Trans_
         ip->incA = lda SHIFT;
      #else
         ip->incA = 1 SHIFT;
      #endif
   #else
      #ifdef Trans_
         ip->incA = 1 SHIFT;
      #else
         ip->incA = lda SHIFT;
      #endif
   #endif
/*
 * TRMM Kernel
 */
   ip->amm_b0 = Mjoin(PATL,trmmK_b0);
   ip->amm_b1 =  Mjoin(PATL,trmmK_b1);
   #ifdef TCPLX
      ip->amm_bn = Mjoin(PATL,trmmK_bn);
   #endif
/*
 * check whether the kernel matches with gemm
 */
   if (gip->mu == ip->mu && gip->nu == ip->nu
         && gip->ku == ip->ku && gip->vlen == ip->vlen)
   {
/*
 *    check mdim of both kernels
 */
      if (ip->kvec == ATL_AMM_KMAJOR(ATL_AMM_GetFLAG(gip->idxK)))
         ip->flg |= 1;
   }
/*
 * TRMM copy routines....
 */
   if (ip->flg & 1) /* trmm kernel similar to gemm */
   {
      alpC = gip->alpC;
      #ifdef Right_
         alpL = gip->alpB;
         alpR = gip->alpA;
      #else
         alpL = gip->alpA;
         alpR = gip->alpB;
      #endif
/*
 *    l2a
 */
      #ifdef Upper_
         if (SCALAR_IS_ONE(alpL))
         {
            if (Diag == AtlasUnit)
               *t2a = Mjoin(PATL,trmm_a2blk_Up_diagU_a1);
            else
               *t2a = Mjoin(PATL,trmm_a2blk_Up_a1);
         }
         else
         {
            if (Diag == AtlasUnit)
               *t2a = SCALAR_IS_NONE(alpL)? Mjoin(PATL,trmm_a2blk_Up_diagU_aN):
                  Mjoin(PATL,trmm_a2blk_Up_diagU_aX);
            else
               *t2a = SCALAR_IS_NONE(alpL)? Mjoin(PATL,trmm_a2blk_Up_aN):
                  Mjoin(PATL,trmm_a2blk_Up_aX);
         }
      #else /* lower */
         if (SCALAR_IS_ONE(alpL))
         {
            if (Diag == AtlasUnit)
               *t2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_a1);
            else
               *t2a = Mjoin(PATL,trmm_a2blk_Lo_a1);
         }
         else
         {
            if (Diag == AtlasUnit)
               *t2a = SCALAR_IS_NONE(alpL)? Mjoin(PATL,trmm_a2blk_Lo_diagU_aN):
                  Mjoin(PATL,trmm_a2blk_Lo_diagU_aX);
            else
               *t2a = SCALAR_IS_NONE(alpL)? Mjoin(PATL,trmm_a2blk_Lo_aN):
                  Mjoin(PATL,trmm_a2blk_Lo_aX);
         }
      #endif
/*
 *    r2a: using gemm's alpha distribution
 */
      if (SCALAR_IS_ONE(alpR))
         *r2a = Mjoin(PATL,trmm_b2blk_a1);
      else
         *r2a = SCALAR_IS_NONE(alpR) ? Mjoin(PATL,trmm_b2blk_aN) :
            Mjoin(PATL,trmm_b2blk_aX);
/*
 *    blk2c: using gemm's alpha distribution
 *    NOTE: beta is always zero for trmm
 */
      if (SCALAR_IS_ONE(alpC))
         *blk2c = Mjoin(PATL,trmm_blk2c_a1b0);
      else
         *blk2c = SCALAR_IS_NONE(alpR) ? Mjoin(PATL,trmm_blk2c_aNb0) :
            Mjoin(PATL,trmm_blk2c_aXb0);
   }
/*
 * TRMM uses different kernel, no copies can be shared
 * so optimize copy routine considering TRMM
 */
   else
   {
      *blk2c = Mjoin(PATL,trmm_blk2c_a1b0);
      if (SCALAR_IS_ONE(alpha))
      {
         #ifdef Upper_
            if (Diag == AtlasUnit)
               *t2a = Mjoin(PATL,trmm_a2blk_Up_diagU_a1);
            else
               *t2a = Mjoin(PATL,trmm_a2blk_Up_a1);
         #else
            if (Diag == AtlasUnit)
               *t2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_a1);
            else
               *t2a = Mjoin(PATL,trmm_a2blk_Lo_a1);
         #endif
         *r2a = Mjoin(PATL,trmm_b2blk_a1);
      }
      else  /* alpha is not one */
      {
   #ifdef Right_
      if (N > (M>>1)) /* multiply alpha in tri-copy*/
   #else
      if (M > (N>>1)) /* multiply alpha in tri-copy*/
   #endif
      {
         *r2a = Mjoin(PATL,trmm_b2blk_a1);
         if (SCALAR_IS_NONE(alpha))
         {
            #ifdef Upper_
               if (Diag == AtlasUnit)
                  *t2a = Mjoin(PATL,trmm_a2blk_Up_diagU_aN);
               else
                  *t2a = Mjoin(PATL,trmm_a2blk_Up_aN);
            #else
               if (Diag == AtlasUnit)
                  *t2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_aN);
               else
                  *t2a = Mjoin(PATL,trmm_a2blk_Lo_aN);
            #endif
         }
         else /* alphaX, since alpha can't be 0 here */
         {
            #ifdef Upper_
               if (Diag == AtlasUnit)
                  *t2a = Mjoin(PATL,trmm_a2blk_Up_diagU_aX);
               else
                  *t2a = Mjoin(PATL,trmm_a2blk_Up_aX);
            #else
               if (Diag == AtlasUnit)
                  *t2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_aX);
               else
                  *t2a = Mjoin(PATL,trmm_a2blk_Lo_aX);
            #endif
         }
      }
      else /* multiply alpha in gemm copy */
      {
            #ifdef Upper_
               if (Diag == AtlasUnit)
                  *t2a = Mjoin(PATL,trmm_a2blk_Up_diagU_a1);
               else
                  *t2a = Mjoin(PATL,trmm_a2blk_Up_a1);
            #else
               if (Diag == AtlasUnit)
                  *t2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_a1);
               else
                  *t2a = Mjoin(PATL,trmm_a2blk_Lo_a1);
            #endif
            if (SCALAR_IS_NONE(alpha))
               *r2a = Mjoin(PATL,trmm_b2blk_aN);
            else
               *r2a = Mjoin(PATL,trmm_b2blk_aX);
         }
      }
   }
}
