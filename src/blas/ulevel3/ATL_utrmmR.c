/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Majedul Sujon
 * Code contributers : Majedul Sujon, Rakib Hasan
 */
#include "atlas_misc.h"

/* cases: LN, UT, LT, UN ... LC/UC included */
#ifdef Upper_
   #ifdef Trans_
#include Mstr(Mjoin(ATLAS_PRE,utrmmR_LN.h)) /* LN & UT share same header */
      #ifdef Conj_
         #define ATL_utrmmR Mjoin(PATL,utrmmR_UC)
         #define RUC 1
      #else
         #define ATL_utrmmR Mjoin(PATL,utrmmR_UT)
         #define RUT 1
      #endif
   #else
#include Mstr(Mjoin(ATLAS_PRE,utrmmR_LT.h)) /* LT & UN share same header */
      #define ATL_utrmmR Mjoin(PATL,utrmmR_UN)
      #define RUN 1
   #endif
#else /* lower */
   #ifdef Trans_
#include Mstr(Mjoin(ATLAS_PRE,utrmmR_LT.h)) /* LT & UN share same header */
      #ifdef Conj_
         #define ATL_utrmmR Mjoin(PATL,utrmmR_LC)
         #define RLC 1
      #else
         #define ATL_utrmmR Mjoin(PATL,utrmmR_LT)
         #define RLT 1
      #endif
   #else
#include Mstr(Mjoin(ATLAS_PRE,utrmmR_LN.h)) /* LN & UT share same header */
      #define ATL_utrmmR Mjoin(PATL,utrmmR_LN)
      #define RLN 1
   #endif
#endif

int ATL_utrmmR
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
   int tnu, tku, tuu, ntfu, tKr;
   const int MU = ATL_TRMMK_MU,
             NU = ATL_TRMMK_NU,
             KU = ATL_TRMMK_KU,
             VLEN = ATL_TRMMK_VLEN;
   ATL_CSZT K = ((N+KU-1)/KU)*KU;
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
   cm2am_t l2a;
   tcm2am_t r2a;
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
            r2a = Mjoin(PATL,trmm_a2blk_Up_diagU_a1);
         else
            r2a = Mjoin(PATL,trmm_a2blk_Up_a1);
      #else
         if (Diag == AtlasUnit)
            r2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_a1);
         else
            r2a = Mjoin(PATL,trmm_a2blk_Lo_a1);
      #endif
      l2a = Mjoin(PATL,trmm_b2blk_a1);
   }
   else
   {
      if (M > N/2) /* multiply alpha in tri-copy*/
      {
         l2a = Mjoin(PATL,trmm_b2blk_a1);
         if (SCALAR_IS_NONE(alpha))
         {
         #ifdef Upper_
            if (Diag == AtlasUnit)
               r2a = Mjoin(PATL,trmm_a2blk_Up_diagU_aN);
            else
               r2a = Mjoin(PATL,trmm_a2blk_Up_aN);
         #else
            if (Diag == AtlasUnit)
               r2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_aN);
            else
               r2a = Mjoin(PATL,trmm_a2blk_Lo_aN);
         #endif
         }
         else /* alphaX, since alpha can't be 0 here */
         {
         #ifdef Upper_
            if (Diag == AtlasUnit)
               r2a = Mjoin(PATL,trmm_a2blk_Up_diagU_aX);
            else
               r2a = Mjoin(PATL,trmm_a2blk_Up_aX);
         #else
            if (Diag == AtlasUnit)
               r2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_aX);
            else
               r2a = Mjoin(PATL,trmm_a2blk_Lo_aX);
         #endif
         }
      }
      else /* multiply alpha in gemm copy */
      {
         #ifdef Upper_
            if (Diag == AtlasUnit)
               r2a = Mjoin(PATL,trmm_a2blk_Up_diagU_a1);
            else
               r2a = Mjoin(PATL,trmm_a2blk_Up_a1);
         #else
            if (Diag == AtlasUnit)
               r2a = Mjoin(PATL,trmm_a2blk_Lo_diagU_a1);
            else
               r2a = Mjoin(PATL,trmm_a2blk_Lo_a1);
         #endif
         if (SCALAR_IS_NONE(alpha))
            l2a = Mjoin(PATL,trmm_b2blk_aN);
         else
            l2a = Mjoin(PATL,trmm_b2blk_aX);

      }
   }
/*
 * calculate and allocate memory for workspace
 * NOTE: This calculation is used from Rakib's TRMM tester: trmmKtst.c
 */
   tnu = NU, tku = KU;
   tuu = Mmax(tnu, tku);
   /* first do full blocks */
   ntfu = N/tuu;
   tKr = N - ntfu*tuu;
   szFull = (ntfu*(ntfu+1)/2) * tuu * tuu;
   /* now the partial panel */
   #if ( (!defined(Upper_) && !defined(Trans_)) || \
            (defined(Upper_) && defined(Trans_)) )
      szPan = ((tKr+tku-1)/tku)*ntfu * tuu * tku;
   #else
      szPan = ((tKr+tnu-1)/tnu)*ntfu * tuu * tnu;
   #endif
   /* finally the partial block */
   szCorner = ((tKr + tnu - 1)/tnu) * ((tKr + tku - 1)/tku) * tnu * tku;
   szT = szFull + szPan + szCorner;
   szR = mb*K;
   szC = (((MU*NU+VLEN-1)/VLEN)*VLEN)*nmu*nnu;
   sz = ATL_MulBySize(szT + NU*KU + szR + MU*KU + szC + (MU+MU)*NU +
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
      l2a(N, M, alpha, X, ldx, pr, pr+szR);
      r2a(N, alpha, A, lda, pt, pt+szT);
      trmmK_b0(nmu, nnu, K, pr+szR, pt+szT, pc,     pr,     pt+szT, pc+szC);
      trmmK_b0(nmu, nnu, K, pr,     pt+szT, pc+szC, pr,     pt,     pc);
      trmmK_bn(nmu, nnu, K, pr,     pt,     pc,     pr+szR, pt,     pc+szC);
      trmmK_b1(nmu, nnu, K, pr+szR, pt,     pc+szC, pr+szR, pt+szT, pc);
      blk2c(M, N, ONE, pc, pc+szC, ZERO, X, ldx);
   #else
      l2a(N, M, alpha, X, ldx, pr);
      r2a(N, alpha, A, lda, pt);
      trmmK_b0(nmu, nnu, K, pr, pt, pc, pr, pt, pc);
      blk2c(M, N, ATL_rone, pc, ATL_rzero, X, ldx);
   #endif
   free(vp);
   return(0);
}
