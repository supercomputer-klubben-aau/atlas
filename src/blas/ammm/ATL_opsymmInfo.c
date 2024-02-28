#include "atlas_cache.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(ATLAS_UPR,amm_kern.h))
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1bX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1bN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1b1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,IntoCNg_a1b0.h))
#ifdef TCPLX
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_aX.h))
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_aX.h))
#endif
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_a1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_a1.h))

int Mjoin(PATL,opsymmInfo) /* pop opinfo_t wt atlas_opgen_view K-3 */
   (opinfo_t *op, ATL_UINT bv,
    ATL_CSZT M, ATL_CSZT N, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta)
/*
 * bv: 0:Left, 1:Upper, 2:HEMM
 * RETURNS: 0 if A fits in cache for outer-product, else non-zero.
 */
{
   ATL_cparr_t *cp = (ATL_cparr_t *)&op->blk2C;
   ATL_iptr_t szSet;
   double spmfMM;
   #ifdef TCPLX
      static const TYPE CONE[2] = {ATL_rone, ATL_rzero};
   #else
      #define CONE ATL_rone
   #endif
   ATL_CUINT K = (bv&1) ? M:N, iv=K-3;
   ATL_UINT MB, NB, KB, mb, nb, k, imm, icpA, icpB, icpC, mu, nu, ku, nmu, nnu;

   if (iv >= ATL_VWopgen_NCASES)
         return(1);
   ATL_GetVWopgenInfo(iv, spmfMM, MB, NB, k, imm, icpA, icpB, icpC, k, k, k);
   ATL_assert(imm < ATL_AMM_NCASES);
   icpA += icpA;
   icpB += icpB;
   icpC += icpC;
   #ifdef TCPLX
      ATL_AMM_info(imm, &op->amm_b0, &op->amm_b1, &op->amm_bn, k,
                   mu, nu, ku, op->vlen, k, k, k);
   #else
      ATL_AMM_info(imm, &op->amm_b0, &szSet, &szSet, k,
                   mu, nu, ku, op->vlen, k, k, k);
   #endif
   op->idx = iv;
   op->mu = mu;
   op->nu = nu;
   op->ku = ku;
   op->exsz = (mu*nu)<<1;
   op->ldc = ldc;
   op->alpA = op->alpB = CONE;
   cp = (ATL_cparr_t *) &op->blk2C;
   if (SCALAR_IS_NONE(beta))
      *cp = ATL_CpyIntoCNg_a1bN[icpC];
   else if (SCALAR_IS_ONE(beta))
      *cp = ATL_CpyIntoCNg_a1b1[icpC];
   else
      *cp = (SCALAR_IS_ZERO(beta)) ?
         ATL_CpyIntoCNg_a1b0[icpC] : ATL_CpyIntoCNg_a1bX[icpC];
/*
 * If M==K ask at least A fit in LLPC
 */
   KB = ((K+ku-1)/ku)*ku;
   if (bv&1)  /* Left, M=K, symmetric matrix A is gemm's A */
   {
      nmu = (M+mu-1)/mu;
      mb = nmu*mu;
      if (mb > MB) /* must shrink NB to keep working set same */
      {
         ATL_UINT t;
         szSet = (MB+NB)*KB + MB*NB;
         t = mb*KB;
         if (szSet < t)
            return(2);
         nb = (szSet - t) / (KB+mb);
         if (nb < nu)
            return(2);
      }
      else
         nb = NB;
      nnu = nb / nu;
      nb = nnu*nu;
      op->lda = lda;
      op->ldb = ldb;
      op->kb = M;
/*
 *    For complex, hermitian makes it so we need worry about transpose on the
 *    herimitian matrix, and we also cannot apply a complex alpha with an
 *    imaginary component to A, since it won't appear in known-zero imag diag
 */
      cp = (ATL_cparr_t *)&op->a2blk;
      #ifdef TCPLX
         if (bv&4)  /* hermitian, not symmetric! */
         {
            *cp = ATL_CpyFromATg_a1[icpA];
            cp = (ATL_cparr_t *)&op->b2blk;
            if (alpha[1] == ATL_rzero)
               *cp = ATL_CpyFromANg_a1[icpB];
            else
            {
               *cp = ATL_CpyFromANg_aX[icpB];
               op->alpB = alpha;
            }
         }
         else
      #endif
      {
         *cp = ATL_CpyFromANg_a1[icpA];
         cp = (ATL_cparr_t *)&op->b2blk;
         *cp = ATL_CpyFromANg_a1[icpB];
      }
   }
   else      /* Right, N=K, symm's symmetric A is gemm's B */
   {
      nnu = (N+nu-1)/nu;
      nb = nnu*nu;
      if (nb > NB) /* must shrink MB to keep working set same */
      {
         ATL_CUINT t = nb*KB;
         szSet = (MB+NB)*KB + MB*NB;
         if (szSet < t)
            return(2);
         mb = (szSet - t) / (KB+nb);
         if (mb < mu)
            return(2);
      }
      else
         mb = MB;
      nmu = mb / mu;
      mb = nmu*mu;
      op->pincBn = op->incBn = 0;  /* not used for 1 blk opsymm */
      op->lda = ldb;
      op->ldb = lda;
      op->kb = N;
      cp = (ATL_cparr_t *)&op->b2blk;
      #ifdef TCPLX
         if (bv&4)  /* hermitian, not symmetric! */
         {
            *cp = ATL_CpyFromANg_a1[icpB];
            cp = (ATL_cparr_t *)&op->a2blk;
            if (alpha[1] == ATL_rzero)
               *cp = ATL_CpyFromATg_a1[icpA];
            else
            {
               *cp = ATL_CpyFromATg_aX[icpA];
               op->alpA = alpha;
            }
         }
         else
      #endif
      {
         *cp = ATL_CpyFromANg_a1[icpB];
         cp = (ATL_cparr_t *)&op->a2blk;
         *cp = ATL_CpyFromATg_a1[icpA];
      }
   }
   op->KB = KB;
   op->beta = beta;
   #ifdef TCPLX
      op->ONE = CONE;
   #endif
   if (M > mb)
   {
      op->mb = mb;
      op->nmu = nmu;
      op->nfmblks = M / mb;
      op->mF = M - mb*op->nfmblks;
      if (op->mF)
      {
         op->nmuF = (op->mF+mu-1)/mu;
         op->pmb = op->nmuF * mu;
         op->npmblks = 1;
      }
      else
         op->nmuF = op->pmb = op->npmblks = 0;
      op->npmblks = (op->pmb) ? 1 : 0;
      op->pnmu = op->nmuF;
   }
   else
   {
      nmu = (M+mu-1)/mu;
      op->npmblks = 1;
      op->nfmblks = 0;
      op->nmuF = op->pnmu = op->nmu = nmu;
      op->mb = op->pmb = nmu*mu;
      op->mF = M;
   }
   if (N > nb)
   {
      op->nb = nb;
      op->nnu = nnu;
      op->nfnblks = N / nb;
      op->nF = N - nb*op->nfnblks;
      if (op->nF)
      {
         op->nnuF = (op->nF+nu-1)/nu;
         op->pnb = op->nnuF * nu;
         op->npnblks = 1;
      }
      else
         op->nnuF = op->pnb = op->npnblks = 0;
      op->npnblks = (op->pnb) ? 1 : 0;
      op->pnnu = op->nnuF;
   }
   else
   {
      nnu = (N+nu-1)/nu;
      op->npnblks = 1;
      op->nfnblks = 0;
      op->nnuF = op->pnnu = op->nnu = nnu;
      op->nb = op->pnb = nnu*nu;
      op->nF = N;
   }
   op->szA = op->mb * op->KB;
   op->szB = op->nb * op->KB;
   op->pszA = op->pmb * op->KB;
   op->pszB = op->pnb * op->KB;
   op->szC = op->mb * op->nb;
   k = op->vlen;
   k = ((mu*nu+k-1)/k)*k;
   k *= nmu * nnu;
   op->szC = k;
   if (bv&1)  /* Left,  M=K, symm matrix A is gemm's A */
   {
      op->incBn = op->nb * (ldb SHIFT);
      op->pincBn = op->pnb * (ldb SHIFT);
      op->incAm = op->pincAm = 0;  /* not used for opsymm */
   }
   else       /* Right, N=K, symm matrix A is gemm's B */
   {
      op->incAm = op->mb SHIFT;
      op->pincAm = op->pmb SHIFT;
      op->incBn = op->pincBn = 0;  /* not used for opsymm */
   }
   if (op->vlen > 1)
   {
      ATL_CUINT vlen=op->vlen;
      op->szB = ((op->szB+vlen-1)/vlen)*vlen;
      op->szA = ((op->szA+vlen-1)/vlen)*vlen;
      op->pszB = ((op->pszB+vlen-1)/vlen)*vlen;
      op->pszA = ((op->pszA+vlen-1)/vlen)*vlen;
   }
   #if 0
   printf("symmB=(%u,%u,%u),(%u,%u,%u), OB=(%u,%u)\n", op->mb, op->nb, op->KB,
          op->pmb, op->pnb, op->kb, MB, NB);
   #endif
   return(0);
}
