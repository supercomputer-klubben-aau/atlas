#include "atlas_cache.h"
#include "atlas_amm.h"
#define ATL_WANT_ILCM
#include "atlas_iopt.h"
#include Mstr(Mjoin(ATLAS_UPR,amm_kern.h))
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1bX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1bN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1b1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1b0.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_aX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_aX.h))
#ifdef TCPLX
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_aX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_aX.h))
#endif

#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_aN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_aN.h))
#ifdef TCPLX
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_aN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_aN.h))
#endif

#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_a1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_a1.h))
#ifdef TCPLX
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_a1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_a1.h))
#endif


int Mjoin(PATL,opsyrkInfo) /* pop opinfo_t using atlas_opgen_view K-3 idx */
   (opinfo_t *out, int flag, enum ATLAS_TRANS TA,
    ATL_CSZT N, ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta)
/*
 * RETURNS: 0 if opsyrk can do this, problem, non-zero if it can't.
 *          If non-zero is returned, out may not be fully populated.
 */
{
   ammkern_t amm1, amm0, ammN;
   ATL_cparr_t *cp = (ATL_cparr_t *)&out->a2blk;
   double spmfMM, spmfCP;
   const int CONJ=flag&1;
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

   if (iv >= ATL_VWopgen_NCASES)
      return(1);
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

   if (MB != NB)   /* need to make square */
   {
      int U;
      U = ATL_iLCM(mu, nu);  /* square NB must be multiple of U */
      if (N < U)
        return(2);
      NB = (MB+NB)>>1;       /* average B keeps cache usage similar */
      NB = Mmin(NB, N);
      if (NB <= U)
         NB = U;
      else if (NB >= 3*U)
         NB = (NB/U)*U;
      else
         NB = ((NB+U-1)/U)*U;
   }
   if (N > NB)  /* if we have more than 1 block, set everything up */
   {
      mb = nb = NB;
      nfnblks = N / NB;
      nmu = NB / mu;
      nnu = NB / nu;
      out->mF = out->nF = pnb = N - nfnblks*NB;
      pnmu = (pnb+mu-1)/mu;
      pnnu = (pnb+nu-1)/nu;
      pmb = pnmu * mu;
      pnb = pnnu * nu;
      npnblks = (pnb) ? 1 : 0;
   }
   else   /* else only one, perhaps partial, block */
   {
      npnblks = 1;
      pnnu = nnu = (N+nu-1)/nu;
      pnmu = nmu = (N+mu-1)/mu;
      pnb = nb = nnu * nu;
      pmb = mb = pnmu * mu;
      nfnblks = 0;
      out->mF = out->nF = N;
   }

   kb = ((K+ku-1)/ku)*ku;

   out->nfmblks = out->nfnblks = nfnblks;
   out->npmblks = out->npnblks = npnblks;
   out->mb = mb;
   out->nb = nb;
   out->pmb = pmb;
   out->pnb = pnb;
   out->mu = mu;
   out->nu = nu;
   out->ku = ku;
   out->vlen = vlen;
   out->ldb = out->lda = lda;
   out->ldc = ldc;

   out->nmu = nmu;
   out->nnu = nnu;
   out->nmuF = out->pnmu = pnmu;
   out->nnuF = out->pnnu = pnnu;

   sz = ((mu*nu+vlen-1)/vlen)*vlen;
   sz *= (nnu) ? nmu*nnu : pnmu*pnnu;
   out->szC = sz;
   out->exsz = (mu*nu)<<1;
   out->idx = iv;
   out->kb = K;
   out->KB = kb;
   out->pszA = pmb*kb;
   out->pszB = pnb*kb;
   out->szB = out->szA = kb * nb;
   if (vlen > 1)
   {
      out->szB = ((out->szB+vlen-1)/vlen)*vlen;
      out->szA = ((out->szA+vlen-1)/vlen)*vlen;
      out->pszB = ((out->pszB+vlen-1)/vlen)*vlen;
      out->pszA = ((out->pszA+vlen-1)/vlen)*vlen;
   }
   out->alpA = out->alpB = CONE;
   out->beta = beta;
   out->alpB = alpha;
   #ifdef TCPLX
      out->ONE = CONE;
   if (TA == AtlasNoTrans || TA == AtlasConj)
   #else
   if (TA == AtlasNoTrans)
   #endif
   {  /* A noTrans, B Trans */
      *cp = ATL_CpyFromATg_a1[icpA];
      cp = (ATL_cparr_t *)&out->b2blk;
      #ifdef TCPLX
         if (SCALAR_IS_NONE(alpha))
            *cp = (CONJ) ? ATL_CpyFromAHg_aN[icpB] : ATL_CpyFromATg_aN[icpB];
         else if (SCALAR_IS_ONE(alpha))
            *cp = (CONJ) ? ATL_CpyFromAHg_a1[icpB] : ATL_CpyFromATg_a1[icpB];
         else
            *cp = (CONJ) ? ATL_CpyFromAHg_aX[icpB] : ATL_CpyFromATg_aX[icpB];
      #else
         if (SCALAR_IS_NONE(alpha))
            *cp = ATL_CpyFromATg_aN[icpB];
         else
            *cp = (SCALAR_IS_ONE(alpha))? ATL_CpyFromATg_a1[icpB]
                                         :ATL_CpyFromATg_aX[icpB];
      #endif
      out->incAm  = mb SHIFT;
      out->pincAm = pmb SHIFT;
      out->incBn  = nb SHIFT;
      out->pincBn = pnb SHIFT;
   }
   else  /* A (Conj)transpose, B no trans */
   {
      #ifdef TCPLX
         *cp = CONJ ? ATL_CpyFromACg_a1[icpA]:ATL_CpyFromANg_a1[icpA];
      #else
         *cp = ATL_CpyFromANg_a1[icpA];
      #endif
      cp = (ATL_cparr_t *)&out->b2blk;
      if (SCALAR_IS_NONE(alpha))
         *cp = ATL_CpyFromANg_aN[icpB];
      else
         *cp = (SCALAR_IS_ONE(alpha)) ?  ATL_CpyFromANg_a1[icpB]
                                       : ATL_CpyFromANg_aX[icpB];
      out->incAm = lda*(mb SHIFT);
      out->pincAm = lda*(pmb SHIFT);
      out->incBn = lda*(nb SHIFT);
      out->pincBn = lda*(pnb SHIFT);
   }
   cp = (ATL_cparr_t *) &out->blk2C;
   if (SCALAR_IS_NONE(beta))
      *cp = ATL_CpyIntoCNg_a1bN[icpC];
   else if (SCALAR_IS_ONE(beta))
      *cp = ATL_CpyIntoCNg_a1b1[icpC];
   else
      *cp = (SCALAR_IS_ZERO(beta)) ?
            ATL_CpyIntoCNg_a1b0[icpC] : ATL_CpyIntoCNg_a1bX[icpC];
   return(0);
}
