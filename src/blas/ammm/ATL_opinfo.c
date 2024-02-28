#include "atlas_cache.h"
#include "atlas_amm.h"
#ifdef TREAL
   #define  ATL_DECL_ 1
#endif
#include Mstr(Mjoin(ATLAS_UPR,amm_kern.h))
#ifndef TREAL
   #define  ATL_DECL_ 1
#endif
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_aX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_aX.h))
#ifdef TCPLX
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_aX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_aX.h))
#endif
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_aXbX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_aXbN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_aXb1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_aXb0.h))

#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_aN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_aN.h))
#ifdef TCPLX
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_aN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_aN.h))
#endif
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_aNbX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_aNbN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_aNb1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_aNb0.h))

#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_a1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_a1.h))
#ifdef TCPLX
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_a1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_a1.h))
#endif
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1bX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1bN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1b1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1b0.h))


void Mjoin(PATL,opinfo) /* pop opinfo_t using atlas_opgen_view K-3 idx */
   (opinfo_t *out, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    ATL_CSZT M, ATL_CSZT N, ATL_CSZT K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta)
{
   ammkern_t amm1, amm0, ammN;
   ATL_cparr_t *cp = (ATL_cparr_t *)&out->a2blk;
   double spmfMM, spmfCP;
   size_t nfmblks, npmblks, nfnblks, npnblks;
   unsigned int icpA, icpB, icpC, imm, k, mu, nu, ku, vlen, flgK;
   unsigned int sz, nmu, nnu, pnmu, pnnu, mb, pmb, nb, pnb, MB, NB, kb;
   const unsigned int iv=K-3;
   #ifdef TCPLX
      static const TYPE CONE[2] = {ATL_rone, ATL_rzero};
   #else
      ATL_iptr_t iptmp;
      #define CONE ATL_rone
   #endif

   ATL_assert(iv < ATL_VWopgen_NCASES);
   ATL_GetVWopgenInfo(iv, spmfMM, MB, NB, k, imm, icpA, icpB, icpC, k, k, k);
   ATL_assert(imm < ATL_AMM_NCASES);
   icpA += icpA;
   icpB += icpB;
   icpC += icpC;
   #ifdef TCPLX
      ATL_AMM_info(imm, &out->amm_b0, &out->amm_b1, &out->amm_bn, flgK,
                   mu, nu, ku, vlen, k, k, k);
   #else
      ATL_AMM_info(imm, &out->amm_b0, &iptmp, &iptmp, flgK,
                   mu, nu, ku, vlen, k, k, k);
   #endif

   if (M > MB)
   {
      mb = MB;
      nfmblks = M / MB;
      nmu = MB / mu;
      out->mF = pmb = M - nfmblks*MB;
      pnmu = (pmb+mu-1)/mu;
      pmb = pnmu * mu;
      npmblks = (pmb) ? 1 : 0;
   }
   else
   {
      npmblks = 1;
      pnmu = nmu = (M+mu-1)/mu;
      pmb = mb = nmu * mu;
      nfmblks = 0;
      out->mF = M;
   }
   if (N > NB)
   {
      nb = NB;
      nfnblks = N / NB;
      nnu = NB / nu;
      out->nF = pnb = N - nfnblks*NB;
      pnnu = (pnb+nu-1)/nu;
      pnb = pnnu * nu;
      npnblks = (pnb) ? 1 : 0;
   }
   else
   {
      npnblks = 1;
      pnnu = nnu = (N+nu-1)/nu;
      pnb = nb = nnu * nu;
      nfnblks = 0;
      out->nF = N;
   }

   kb = ((K+ku-1)/ku)*ku;
   out->nfmblks = nfmblks;
   out->npmblks = npmblks;
   out->mb = mb;
   out->pmb = pmb;
   out->nfnblks = nfnblks;
   out->npnblks = npnblks;
   out->nb = nb;
   out->pnb = pnb;
   out->mu = mu;
   out->nu = nu;
   out->ku = ku;
   out->vlen = vlen;
   out->lda = lda;
   out->ldb = ldb;
   out->ldc = ldc;

   out->nmu = nmu;
   out->nmuF = out->pnmu = pnmu;
   out->nnu = nnu;
   out->nnuF = out->pnnu = pnnu;

   sz = ((mu*nu+vlen-1)/vlen)*vlen;
   sz *= (out->nmu) ? out->nmu : out->pnmu;
   sz *= (out->nnu) ? out->nnu : out->pnnu;
   out->szC = sz;
   out->exsz = (mu*nu)<<1;
   out->idx = iv;
   out->kb = K;
   out->KB = kb;
   out->lda = lda;
   out->ldb = ldb;
   out->ldc = ldc;
   out->pszA = pmb*kb;
   out->szA = kb * mb;
   out->pszB = kb*pnb;
   out->szB = kb * nb;
   if (vlen > 1)
   {
      out->szB = ((out->szB+vlen-1)/vlen)*vlen;
      out->szA = ((out->szA+vlen-1)/vlen)*vlen;
      out->pszB = ((out->pszB+vlen-1)/vlen)*vlen;
      out->pszA = ((out->pszA+vlen-1)/vlen)*vlen;
   }
   out->alpA = out->alpB = CONE;
   out->beta = beta;
   #ifdef TCPLX
      out->ONE = CONE;
   if (TA == AtlasNoTrans || TA == AtlasConj)
   #else
   if (TA == AtlasNoTrans)
   #endif
   {
      out->incAm  = mb SHIFT;
      out->pincAm = pmb SHIFT;
   }
   else
   {
      out->incAm = lda*(mb SHIFT);
      out->pincAm = lda*(pmb SHIFT);
   }
   #ifdef TCPLX
   if (TB == AtlasNoTrans || TB == AtlasConj)
   #else
   if (TB == AtlasNoTrans)
   #endif
   {
      out->incBn = ldb*(nb SHIFT);
      out->pincBn = ldb*(pnb SHIFT);
   }
   else
   {
      out->incBn  = nb SHIFT;
      out->pincBn = pnb SHIFT;
   }
/*
 * Once we have ATL_ammmN written, put alpha on A if there's only 1 blk of
 * of B.  For now, always put on B and use ATL_ammmM
 */
   #if 0
   if (nnblks > 1 && nmblks < 2) /* If we only have row panel of C */
   {                             /* apply alpha to A */
      #ifdef TCPLX
         if (TB == AtlasNoTrans)
            out->b2blk = ATL_AMM_BN2BLK_a1[idx];
         else if (TB == AtlasTrans)
            out->b2blk = ATL_AMM_BT2BLK_a1[idx];
         else
            out->b2blk = (TB == AtlasConjTrans) ?
               ATL_AMM_BH2BLK_a1[idx] : ATL_AMM_BC2BLK_a1[idx];

         if (SCALAR_IS_ONE(alpha))
         {
            if (TA == AtlasNoTrans)
               out->a2blk = ATL_AMM_AT2BLK_a1[idx];
            else if (TA == AtlasTrans)
               out->a2blk = ATL_AMM_AN2BLK_a1[idx];
            else
               out->a2blk = (TA == AtlasConjTrans) ?
                  ATL_AMM_AC2BLK_a1[idx] : ATL_AMM_AH2BLK_a1[idx];
         }
         else if (SCALAR_IS_NONE(alpha))
         {
            if (TA == AtlasNoTrans)
               out->a2blk = ATL_AMM_AT2BLK_an[idx];
            else if (TA == AtlasTrans)
               out->a2blk = ATL_AMM_AN2BLK_an[idx];
            else
               out->a2blk = (TA == AtlasConjTrans) ?
                  ATL_AMM_AC2BLK_an[idx] : ATL_AMM_AH2BLK_an[idx];
         }
         else
         {
            if (TA == AtlasNoTrans)
               out->a2blk = ATL_AMM_AT2BLK_aX[idx];
            else if (TA == AtlasTrans)
               out->a2blk = ATL_AMM_AN2BLK_aX[idx];
            else
               out->a2blk = (TA == AtlasConjTrans) ?
                       ATL_AMM_AC2BLK_aX[idx] : ATL_AMM_AH2BLK_aX[idx];
         }
      #else
         out->b2blk = (TB == AtlasNoTrans) ?
                      ATL_AMM_BN2BLK_a1[idx] : ATL_AMM_BT2BLK_a1[idx];
         if (SCALAR_IS_ONE(alpha))
            out->a2blk = (TA == AtlasNoTrans) ?
                      ATL_AMM_AT2BLK_a1[idx] : ATL_AMM_AN2BLK_a1[idx];
         else if (SCALAR_IS_NONE(alpha))
            out->a2blk = (TA == AtlasNoTrans) ?
                      ATL_AMM_AT2BLK_an[idx] : ATL_AMM_AN2BLK_an[idx];
         else
            out->a2blk = (TA == AtlasNoTrans) ?
                      ATL_AMM_AT2BLK_aX[idx] : ATL_AMM_AN2BLK_aX[idx];
      #endif
   }
   else  /* apply alpha to B */
   #endif
   {
      out->alpB = alpha;
      #ifdef TCPLX
         if (TA == AtlasNoTrans)
            *cp = ATL_CpyFromATg_a1[icpA];
         else if (TA == AtlasTrans)
            *cp = ATL_CpyFromANg_a1[icpA];
         else
            *cp = (TA == AtlasConjTrans) ?
                  ATL_CpyFromACg_a1[icpA]: ATL_CpyFromAHg_a1[icpA];

         cp = (ATL_cparr_t *)&out->b2blk;
         if (SCALAR_IS_ONE(alpha))
         {
            if (TB == AtlasNoTrans)
               *cp = ATL_CpyFromANg_a1[icpB];
            else if (TB == AtlasTrans)
               *cp = ATL_CpyFromATg_a1[icpB];
            else
               *cp = (TB == AtlasConjTrans) ?
                     ATL_CpyFromAHg_a1[icpB]: ATL_CpyFromACg_a1[icpB];
         }
         else if (SCALAR_IS_NONE(alpha))
         {
            if (TB == AtlasNoTrans)
               *cp = ATL_CpyFromANg_aN[icpB];
            else if (TB == AtlasTrans)
               *cp = ATL_CpyFromATg_aN[icpB];
            else
               *cp = (TB == AtlasConjTrans) ?
                     ATL_CpyFromAHg_aN[icpB]: ATL_CpyFromACg_aN[icpB];
         }
         else
         {
            if (TB == AtlasNoTrans)
               *cp = ATL_CpyFromANg_aX[icpB];
            else if (TB == AtlasTrans)
               *cp = ATL_CpyFromATg_aX[icpB];
            else
               *cp = (TB == AtlasConjTrans) ?
                     ATL_CpyFromAHg_aX[icpB]: ATL_CpyFromACg_aN[icpB];
         }
      #else
         *cp = (TA == AtlasNoTrans) ?
               ATL_CpyFromATg_a1[icpA]:ATL_CpyFromANg_a1[icpA];
         cp = (ATL_cparr_t *) &out->b2blk;
         if (SCALAR_IS_ONE(alpha))
            *cp = (TB == AtlasNoTrans) ?
                  ATL_CpyFromANg_a1[icpB]:ATL_CpyFromATg_a1[icpB];
         else if (SCALAR_IS_NONE(alpha))
            *cp = (TB == AtlasNoTrans) ?
                  ATL_CpyFromANg_aN[icpB]:ATL_CpyFromATg_aN[icpB];
         else
            *cp = (TB == AtlasNoTrans) ?
                  ATL_CpyFromANg_aX[icpB]:ATL_CpyFromATg_aX[icpB];
      #endif
   }
   cp = (ATL_cparr_t *) &out->blk2C;
   if (SCALAR_IS_NONE(beta))
      *cp = ATL_CpyIntoCNg_a1bN[icpC];
   else if (SCALAR_IS_ONE(beta))
      *cp = ATL_CpyIntoCNg_a1b1[icpC];
   else
      *cp = (SCALAR_IS_ZERO(beta)) ?
            ATL_CpyIntoCNg_a1b0[icpC] : ATL_CpyIntoCNg_a1bX[icpC];
}

