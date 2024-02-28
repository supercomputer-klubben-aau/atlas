#include "atlas_cache.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(ATLAS_UPR,amm_kern.h))
#define  ATL_DECL_ 1
#include Mstr(Mjoin(ATLAS_PRE,ipmek_view.h))
#undef ATL_DECL_
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

void Mjoin(PATL,ipmekInfoPop)  /* populate ip wt given params */
(
   ipinfo_t *ip,        /* output */
   int idx,             /* what amm kernel index to use */
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   size_t M,
   size_t N,
   size_t K,
   size_t lda,
   size_t ldb,
   size_t ldc,
   const SCALAR alpha,
   const SCALAR beta,
   size_t nfmblks,
   size_t npmblks,
   ATL_UINT mb,
   ATL_UINT pmb,
   size_t nfnblks,
   size_t npnblks,
   ATL_UINT nb,
   ATL_UINT pnb
)
{
   ATL_cparr_t *cp = (ATL_cparr_t *)&(ip->a2blk);
   #ifdef TCPLX
      static const TYPE CONE[2] = {ATL_rone, ATL_rzero};
   #else
      ATL_iptr_t iptmp;
      #define CONE ATL_rone
   #endif
   ATL_UINT vlen, mu, nu, ku, kb, nkblks, KB0, kb0, szC;
   ATL_UINT imm, ik1, icpA, icpB, icpC, flgK, kbmin, k;
   float spmf;

   ip->idx = idx;
   ATL_GetVWipmekInfo(idx, spmf, k, k, kb, imm, icpA, icpB, icpC, k,
                      kbmin, k);
   ip->idxK = imm;
   icpA += icpA;
   icpB += icpB;
   icpC += icpC;
   #ifdef TCPLX
      ATL_AMM_info(imm, &ip->amm_b0, &ip->amm_b1, &ip->amm_bn, flgK,
                   mu, nu, ku, vlen, k, k, ik1);
   #else
      ATL_AMM_info(imm, &ip->amm_b0, &ip->amm_b1, &iptmp, flgK,
                   mu, nu, ku, vlen, k, k, ik1);
   #endif
   ip->mu = mu;
   ip->nu = nu;
   ip->ku = ku;
   ip->kb = kb;
   ip->vlen = vlen;
   ip->lda = lda;
   ip->ldb = ldb;
   ip->ldc = ldc;

   ip->nfmblks = nfmblks;
   ip->npmblks = npmblks;
   ip->nfnblks = nfnblks;
   ip->npnblks = npnblks;
   ip->mb = mb;
   ip->pmb = pmb;
   ip->nmu = mb / mu;
   ip->pnmu = pmb / mu;
   if (npmblks)
      ip->mF = M - nfmblks*mb - (npmblks-1)*pmb;
   else
      ip->mF = M - (nfmblks-1)*mb;
   ip->nmuF = (ip->mF+mu-1) / mu;

   ip->nb = nb;
   ip->pnb = pnb;
   ip->nnu = nb / nu;
   ip->pnnu = pnb / nu;
   if (npnblks)
      ip->nF = N - nfnblks*nb - (npnblks-1)*pnb;
   else
      ip->nF = N - (nfnblks-1)*nb;
   ip->nnuF = (ip->nF+nu-1) / nu;

/*
 * Compute K remainder block, and how it affects kernel to use
 */
   nkblks = K / kb;
   KB0 = kb0 = K - nkblks*kb;
   ip->ammK1_b0 = ip->amm_b0;  /* hope: normal kerns for K-clean */
   #ifdef TCPLX
      ip->ONE = CONE;
   #endif
   if (kb0)  /* K not a multiple of kb, need initial block */
   {
      if (ATL_AMM_KMAJOR(flgK))
      {
         KB0 = ((kb0+vlen-1)/vlen)*vlen; /* k-maj req K mul of vlen */
         if (KB0 != kb)
         {
            if (ATL_AMM_KRUNTIME(flgK))
            {
               if ((kbmin && KB0 < kbmin) || (KB0/ku)*ku != KB0)
                  ip->ammK1_b0 = NULL;
            }
            else
               ip->ammK1_b0 = NULL;
         }
      }
      else if (ATL_AMM_KRUNTIME(flgK))
      {
         KB0 = ((kb0+ku-1)/ku)*ku;
         if ((kbmin && kb0 < kbmin) || KB0-kb0 > ku)
            ip->ammK1_b0 = NULL;
      }
      else /* compile-time K needs K1 clean whenever KB0 != kb0 */
         ip->ammK1_b0 = NULL;
      if (!ip->ammK1_b0) /* must use ik1 kerns instead! */
      {
         ATL_assert(ik1);
         ik1--;
         #ifdef TCPLX
            ATL_AMM_kinfo(ik1, &ip->ammK1_b0, &ip->ammK1_b1, &ip->ammK1_bn);
         #else
            ATL_AMM_kinfo(ik1, &ip->ammK1_b0, &ip->ammK1_b1, &iptmp);
         #endif
         KB0 = ATL_AMM_GetKU(ik1);
         KB0 = ((kb0+KB0-1)/KB0)*KB0;
      }
      else
      {
         ip->ammK1_b1 = ip->amm_b1;
         #ifdef TCPLX
            ip->ammK1_bn = ip->amm_bn;
         #endif
      }
   }
   else
   {
      ip->ammK1_b1 = ip->amm_b1;
      #ifdef TCPLX
         ip->ammK1_bn = ip->amm_bn;
      #endif
      kb0 = KB0 = kb;
      nkblks--;
      ATL_assert(nkblks >= 0);
   }
   ip->KB0 = KB0;
   ip->kb0 = kb0;
   ip->exsz = (mu+mu)*nu;
   szC = ((mu*nu+vlen-1)/vlen)*vlen;
   szC *= ip->nnu * ip->nmu;
   ip->szC = szC;


   ip->alpA = ip->alpB = ip->alpC = CONE;
   if (M < N)  /* alpha goes on A or C */
   {
      if (M < K)
         ip->alpC = alpha;
      else
         ip->alpA = alpha;
   }
   else if (N < K)
      ip->alpC = alpha;
   else
      ip->alpB = alpha;
   if (SCALAR_IS_NONE(ip->alpA))
   {
      #ifdef TCPLX
      if (TA == AtlasConjTrans)
         *cp = ATL_CpyFromACg_aN[icpA];
      else if (TA == AtlasConj)
         *cp = ATL_CpyFromAHg_aN[icpA];
      else
      #endif
         *cp = (TA == AtlasNoTrans) ?
               ATL_CpyFromATg_aN[icpA] : ATL_CpyFromANg_aN[icpA];
   }
   else if (SCALAR_IS_ONE(ip->alpA))
   {
      #ifdef TCPLX
      if (TA == AtlasConjTrans)
         *cp = ATL_CpyFromACg_a1[icpA];
      else if (TA == AtlasConj)
         *cp = ATL_CpyFromAHg_a1[icpA];
      else
      #endif
         *cp = (TA == AtlasNoTrans) ?
               ATL_CpyFromATg_a1[icpA] : ATL_CpyFromANg_a1[icpA];
   }
   else  /* alphaA = X */
   {
      #ifdef TCPLX
      if (TA == AtlasConjTrans)
         *cp = ATL_CpyFromACg_aX[icpA];
      else if (TA == AtlasConj)
         *cp = ATL_CpyFromAHg_aX[icpA];
      else
      #endif
         *cp = (TA == AtlasNoTrans) ?
               ATL_CpyFromATg_aX[icpA] : ATL_CpyFromANg_aX[icpA];
   }

   cp = (ATL_cparr_t *) &(ip->b2blk);
   if (SCALAR_IS_NONE(ip->alpB))
   {
      #ifdef TCPLX
      if (TB == AtlasConjTrans)
         *cp = ATL_CpyFromAHg_aN[icpB];
      else if (TB == AtlasConj)
         *cp = ATL_CpyFromACg_aN[icpB];
      else
      #endif
         *cp = (TB == AtlasNoTrans) ?
            ATL_CpyFromANg_aN[icpB] : ATL_CpyFromATg_aN[icpB];
   }
   else if (SCALAR_IS_ONE(ip->alpB))
   {
      #ifdef TCPLX
      if (TB == AtlasConjTrans)
         *cp = ATL_CpyFromAHg_a1[icpB];
      else if (TB == AtlasConj)
         *cp = ATL_CpyFromACg_a1[icpB];
      else
      #endif
         *cp = (TB == AtlasNoTrans) ?
            ATL_CpyFromANg_a1[icpB] : ATL_CpyFromATg_a1[icpB];
   }
   else  /* alphaB = X */
   {
      #ifdef TCPLX
      if (TB == AtlasConjTrans)
         *cp = ATL_CpyFromAHg_aX[icpB];
      else if (TB == AtlasConj)
         *cp = ATL_CpyFromACg_aX[icpB];
      else
      #endif
         *cp = (TB == AtlasNoTrans) ?
            ATL_CpyFromANg_aX[icpB] : ATL_CpyFromATg_aX[icpB];
   }

   cp = (ATL_cparr_t *) &(ip->blk2c_b1);
   if (SCALAR_IS_ONE(ip->alpC))
   {
      *cp = ATL_CpyIntoCNg_a1b1[icpC];
      cp = (ATL_cparr_t *) &(ip->blk2c);
      if (SCALAR_IS_NONE(beta))
         *cp = ATL_CpyIntoCNg_a1bN[icpC];
      else if (SCALAR_IS_ONE(beta))
         ip->blk2c = ip->blk2c_b1;
      else
         *cp = (SCALAR_IS_ZERO(beta)) ?
            ATL_CpyIntoCNg_a1b0[icpC] : ATL_CpyIntoCNg_a1bX[icpC];
   }
   else if (SCALAR_IS_NONE(ip->alpC))
   {
      *cp = ATL_CpyIntoCNg_aNb1[icpC];
      cp = (ATL_cparr_t *) &(ip->blk2c);
      if (SCALAR_IS_NONE(beta))
         *cp = ATL_CpyIntoCNg_aNbN[icpC];
      else if (SCALAR_IS_ONE(beta))
         ip->blk2c = ip->blk2c_b1;
      else
         *cp = (SCALAR_IS_ZERO(beta)) ?
            ATL_CpyIntoCNg_aNb0[icpC] : ATL_CpyIntoCNg_aNbX[icpC];
   }
   else /* alphaC = X */
   {
      *cp = ATL_CpyIntoCNg_aXb1[icpC];
      cp = (ATL_cparr_t *) &(ip->blk2c);
      if (SCALAR_IS_NONE(beta))
         *cp = ATL_CpyIntoCNg_aXbN[icpC];
      else if (SCALAR_IS_ONE(beta))
         ip->blk2c = ip->blk2c_b1;
      else
         *cp = (SCALAR_IS_ZERO(beta)) ?
            ATL_CpyIntoCNg_aXb0[icpC] : ATL_CpyIntoCNg_aXbX[icpC];
   }

   ip->nfkblks = nkblks;
   if (IS_COLMAJ(TA))
   {
      ip->incAk  = (kb SHIFT)*lda;
      ip->incAm  = (mb SHIFT);
      ip->pincAm = (pmb SHIFT);
   }
   else
   {
      ip->incAk  = (kb SHIFT);
      ip->incAm  = (mb SHIFT)*lda;
      ip->pincAm = (pmb SHIFT)*lda;
   }
   if (IS_COLMAJ(TB))
   {
      ip->incBk  = (kb SHIFT);
      ip->incBn  = (nb SHIFT)*ldb;
      ip->pincBn = (pnb SHIFT)*ldb;
   }
   else
   {
      ip->incBk  = (kb SHIFT)*ldb;
      ip->incBn  = (nb SHIFT);
      ip->pincBn = (pnb SHIFT);
   }
   ip->szA = mb*kb;
   ip->szB = kb*nb;
   ip->pszA = pmb*kb;
   ip->pszB = kb*pnb;
   if (vlen > 1)
   {
      ip->szB = ((ip->szB+vlen-1)/vlen)*vlen;
      ip->szA = ((ip->szA+vlen-1)/vlen)*vlen;
      ip->pszB = ((ip->pszB+vlen-1)/vlen)*vlen;
      ip->pszA = ((ip->pszA+vlen-1)/vlen)*vlen;
   }
}
#ifndef TCPLX
   #undef CONE
#endif
void Mjoin(PATL,ipmekInfo)
(
   ipinfo_t *ip,        /* output */
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   size_t M,
   size_t N,
   size_t K,
   size_t lda,
   size_t ldb,
   size_t ldc,
   const SCALAR alpha,
   const SCALAR beta
)
{
   size_t nfmblks, npmblks, nfnblks, npnblks;
   const ATL_iptr_t *ipp;
   ATL_UINT mb, pmb, nb, pnb, mu, nu, nmu, pnmu, nnu, pnnu;
   ATL_UINT MB, NB, KB, k, imm, idx=ATL_VIEW_BEST_IDX;
   float spf;

   #if ATL_VIEW_NCASES > 1
   if (K < ATL_VIEW_BEST_KB)
   {
      for (idx=0; idx < ATL_VIEW_LAST_IDX; idx++)
         if (ATL_GetViewKB(idx) >= K)
            break;
   }
   #endif
   ATL_GetViewInfo(idx, spf, MB, NB, KB, imm, k, k, k, k, k, k);
   ipp = ATL_AMM_ML + ATL_AMM_Idx2Entry(imm) + 3;
   ATL_AMM_iinfo(ipp, k, mu, nu, k, k, k, k, k);
/*
 * Consider expanding M/N block to cover whole C if best block factor is
 * already close in size, to avoid having one large block, and one tiny
 */
   if (NB+NB > N && MB+MB > M)
   {
      mb = ((M+mu-1)/mu)*mu;
      nb = ((N+nu-1)/nu)*nu;
      if (ATL_MulBySize(mb*nb + 2*(mb+nb)*KB) <= LLPC_SZ)
      {
         MB = mb;
         NB = nb;
         KB = MB; /* KB must be equal to MB for mek */
      }
   }
   if (MB >= M)
   {
      npmblks = 1;
      nfmblks = 0;
      pnmu = nmu = (M+mu-1)/mu;
      mb = pmb = nmu * mu;
   }
   else
   {
      mb = MB;
      nfmblks = M / MB;
      nmu = MB / mu;
      pmb = M - nfmblks*MB;
      if (pmb)
      {
         npmblks = 1;
         pnmu = (pmb+mu-1) / mu;
         pmb = pnmu * mu;
      }
      else
         npmblks = pnmu = pmb = 0;
   }
   if (NB >= N)
   {
      npnblks = 1;
      nfnblks = 0;
      pnnu = nnu = (N+nu-1)/nu;
      nb = pnb = nnu * nu;
   }
   else
   {
      nb = NB;
      nfnblks = N / NB;
      nnu = NB / nu;
      pnb = N - nfnblks*NB;
      if (pnb)
      {
         npnblks = 1;
         pnnu = (pnb+nu-1) / nu;
         pnb = pnnu * nu;
      }
      else
         npnblks = pnnu = pnb = 0;
   }
   Mjoin(PATL,ipmekInfoPop)(ip, idx, TA, TB, M, N, K, lda, ldb, ldc, alpha,
                             beta, nfmblks, npmblks, mb, pmb, nfnblks,
                             npnblks, nb, pnb);
}
