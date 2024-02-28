/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Majedul Sujon
 * Code contributers : Majedul Sujon, Rakib Hasan
 */
#include "atlas_misc.h"
/* cases: LN, UT, LT, UN ... LC/UC included */
#ifdef Upper_
   #ifdef Trans_
#include Mstr(Mjoin(ATLAS_PRE,utrmmL_LN.h)) /* LN & UT share same header file*/
      #ifdef Conj_
         #define ATL_utrmmL Mjoin(PATL,utrmmL_UC)
         #define LUC 1
      #else
         #define ATL_utrmmL Mjoin(PATL,utrmmL_UT)
         #define LUT 1
      #endif
   #else
#include Mstr(Mjoin(ATLAS_PRE,utrmmL_LT.h)) /* LT & UN share same header file*/
      #define ATL_utrmmL Mjoin(PATL,utrmmL_UN)
      #define LUN 1
   #endif
#else /* lower */
   #ifdef Trans_
#include Mstr(Mjoin(ATLAS_PRE,utrmmL_LT.h))
      #ifdef Conj_
         #define ATL_utrmmL Mjoin(PATL,utrmmL_LC)
         #define LLC 1
      #else
         #define ATL_utrmmL Mjoin(PATL,utrmmL_LT)
         #define LLT 1
      #endif
   #else
#include Mstr(Mjoin(ATLAS_PRE,utrmmL_LN.h))
      #define ATL_utrmmL Mjoin(PATL,utrmmL_LN)
      #define LLN 1
   #endif
#endif

int ATL_utrmmL
(
   const enum ATLAS_DIAG Diag,
   ATL_CSZT M,         /* it's actually M */
   ATL_CSZT N,         /* it's actually N*/
   const SCALAR alpha,
   const TYPE *A,      /* A is always tri-matrix: here left one */
   ATL_CSZT lda,
   TYPE *X,            /* it's B: the output and full matrix*/
   ATL_CSZT ldx
)
{
   void *vp=NULL;
   TYPE *pt, *pr, *pc;
   ATL_SZT sz, szC, szT, szR;
   ATL_SZT szFull, szPan, szCorner;
   int tmu, tku, tuu, ntfu, tKr;
   const int MU = ATL_TRMMK_MU,
             NU = ATL_TRMMK_NU,
             KU = ATL_TRMMK_KU,
             VLEN = ATL_TRMMK_VLEN;
   ATL_CSZT K = ((M+KU-1)/KU)*KU;
   ATL_SZT nmu = (M+MU-1)/MU, nnu=(N+NU-1)/NU;
   ATL_SZT mb=nmu*MU, nb=nnu*NU;
   ammkern_t trmmK_b0 = Mjoin(PATL,trmmK_b0);
   ammkern_t trmmK_b1 =  Mjoin(PATL,trmmK_b1);
   #ifdef TCPLX
      ammkern_t trmmK_bn = Mjoin(PATL,trmmK_bn);
      TYPE ONE[2] = {ATL_rone, ATL_rzero}, ZERO[2] = {ATL_rzero, ATL_rzero};
   #else
      #define ONE ATL_rone
      #define ZERO ATL_rzero
      #define NONE ATL_rnone
   #endif
   tcm2am_t l2a;
   cm2am_t r2a;
   ablk2cmat_t blk2c = Mjoin(PATL,trmm_blk2c_a1b0);
/*
 * Selecting copy routines...
 * multiplying with alpha in copy routines:
 *    MxM/2, MXN
 */
   if (SCALAR_IS_ONE(alpha))
   {
      #ifdef Upper_
         if (Diag == AtlasUnit)
            l2a = Mjoin(PATL,trmm_a2blk_Up_diagU_a1);
         else
            l2a = Mjoin(PATL,trmm_a2blk_Up_a1);
      #else
         if (Diag == AtlasUnit)
            l2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_a1);
         else
            l2a = Mjoin(PATL,trmm_a2blk_Lo_a1);
      #endif
      r2a = Mjoin(PATL,trmm_b2blk_a1);
   }
   else
   {
      if (N > M/2) /* multiply alpha in tri-copy*/
      {
         r2a = Mjoin(PATL,trmm_b2blk_a1);
         if (SCALAR_IS_NONE(alpha))
         {
         #ifdef Upper_
            if (Diag == AtlasUnit)
               l2a = Mjoin(PATL,trmm_a2blk_Up_diagU_aN);
            else
               l2a = Mjoin(PATL,trmm_a2blk_Up_aN);
         #else
            if (Diag == AtlasUnit)
               l2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_aN);
            else
               l2a = Mjoin(PATL,trmm_a2blk_Lo_aN);
         #endif
         }
         else /* alphaX, since alpha can't be 0 here */
         {
         #ifdef Upper_
            if (Diag == AtlasUnit)
               l2a = Mjoin(PATL,trmm_a2blk_Up_diagU_aX);
            else
               l2a = Mjoin(PATL,trmm_a2blk_Up_aX);
         #else
            if (Diag == AtlasUnit)
               l2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_aX);
            else
               l2a = Mjoin(PATL,trmm_a2blk_Lo_aX);
         #endif
         }
      }
      else /* multiply alpha in gemm copy */
      {
         #ifdef Upper_
            if (Diag == AtlasUnit)
               l2a = Mjoin(PATL,trmm_a2blk_Up_diagU_a1);
            else
               l2a = Mjoin(PATL,trmm_a2blk_Up_a1);
         #else
            if (Diag == AtlasUnit)
               l2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_a1);
            else
               l2a = Mjoin(PATL,trmm_a2blk_Lo_a1);
         #endif
         if (SCALAR_IS_NONE(alpha))
            r2a = Mjoin(PATL,trmm_b2blk_aN);
         else
            r2a = Mjoin(PATL,trmm_b2blk_aX);
      }
   }
/*
 * calculate and allocate memory for workspace
 * NOTE: This calculation is used from Rakib's TRMM tester: trmmKtst.c
 */
#if EXTRA_SPACE
   szT = mb * K;
#else
   tmu = MU, tku = KU;
   tuu = Mmax(tmu, tku);
   /* first do full blocks */
   ntfu = M/tuu;
   tKr = M - ntfu*tuu;
   szFull = (ntfu*(ntfu+1)/2) * tuu * tuu;
   #if ( (!defined(Upper_) && !defined(Trans_)) || \
         (defined(Upper_) && defined(Trans_)) )
      szPan = ((tKr+tmu-1)/tmu)*ntfu * tuu * tmu;
   #else
      szPan = ((tKr+tku-1)/tku)*ntfu * tuu * tku;
   #endif
   /* finally the partial block */
   szCorner = ((tKr + tmu - 1)/tmu) * ((tKr + tku - 1)/tku) * tmu * tku;
   szT = szFull + szPan + szCorner;
#endif
   szR = nb*K;
   szC = (((MU*NU+VLEN-1)/VLEN)*VLEN)*nmu*nnu;
   sz = ATL_MulBySize(szT + MU*KU + szR + NU*KU + szC + (MU+MU)*NU +
                      3*ATL_Cachelen);
   if (sz < ATL_MaxMalloc)
      vp = malloc(sz);
   if (!vp) return(2);
   pt = ATL_AlignPtr(vp);
   pr = pt + (szT SHIFT);
   pr = ATL_AlignPtr(pr);
   pc = pr + (szR SHIFT);
   pc = ATL_AlignPtr(pc);
   #ifdef TCPLX
      l2a(M, alpha, A, lda, pt, pt+szT);
      r2a(M, N, alpha, X, ldx, pr, pr+szR);
      trmmK_b0(nmu, nnu, K, pt+szT, pr+szR, pc,     pt,     pr+szR, pc+szC);
      trmmK_b0(nmu, nnu, K, pt,     pr+szR, pc+szC, pt,     pr,     pc);
      trmmK_bn(nmu, nnu, K, pt,     pr,     pc,     pt+szT, pr,     pc+szC);
      trmmK_b1(nmu, nnu, K, pt+szT, pr,     pc+szC, pt+szT, pr+szR, pc);
      blk2c(M, N, ONE, pc, pc+szC, ZERO, X, ldx);
   #else
      l2a(M, alpha, A, lda, pt);
      r2a(M, N, alpha, X, ldx, pr);
      trmmK_b0(nmu, nnu, K, pt, pr, pc, pt, pr, pc);
      blk2c(M, N, ONE, pc, ZERO, X, ldx);
   #endif
   free(vp);
   return(0);
}
