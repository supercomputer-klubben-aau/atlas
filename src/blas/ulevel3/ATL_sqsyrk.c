/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_level2.h"
#include "atlas_level3.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_sqsyrk.h))
/*
 * Service routine, particularly for parallel.  Takes its blocking from ip
 * (assuming that is what is being used below diagonal blocks
 * flag bits, meaning if set (opposite if unset):
 * 0/1: C is upper
 * 1/2: TA == AtlasNoTrans
 * 2/4: use beta=0 syrk kernel (else use beta=1)
 *
 * If (Uplo==Upper && sy2blk)
 *    (1) flag&1 == 1; (2) wU non-NULL; (3) sy2blk is beta=0.
 * ==> wU can be aliased wt wS if you don't need correct wS on output
 */
#ifdef Conj_
   #define syrkBlk Mjoin(PATL,sqherkBlk)
   #define syrk_K  Mjoin(PATL,sqherk_KLoop)
   #define syrk  Mjoin(PATL,sqherk)
#else
   #define syrkBlk Mjoin(PATL,sqsyrkBlk)
   #define syrk_K  Mjoin(PATL,sqsyrk_KLoop)
   #define syrk  Mjoin(PATL,sqsyrk)
#endif
/*#define MEM_DEBUG 1*/

void syrkBlk
(
   ipinfo_t *ip,
   int flag,      /* bitvec: 0: set means C is upper, 1: set TA==AtlasNoTrans */
   size_t d,      /* which global diagonal blk of C is being computed */
   size_t k,      /* which K block of A is being computed */
   const TYPE *A, /* if non-NULL, base A ptr to copy */
   cm2am_t sy2blk,/* copy A to syrk storage */
   ablk2cmat_t blk2c, /* if non-NULL sy storage to C storage copy func */
   const SCALAR beta, /* only needed if blk2c non-NULL */
   TYPE *C,       /* if blk2c && non-NULL, which C to write to */
   TYPE *rS,      /* real portion of wS (unused for real routines) */
   TYPE *wS,      /* space to store syrk A; */
   TYPE *wA,      /* if non-NULL, ip-based A workspace */
   TYPE *wAn,     /* next A wrkspc to be prefetched */
   TYPE *wB,      /* if non-NULL ip-based At workspace */
   TYPE *wBn,     /* next B wrkspc to be prefetched */
   TYPE *rC,      /* real portion of wC (unused for real routines) */
   TYPE *wC,      /* if non-NULL: ptr to syrk-storage C wrkspc */
   TYPE *wU       /* reflection space for Upper if non-NULL */
)
{
   ATL_CUINT kb = (k != ip->nfkblks) ? ip->kb : ip->kb0;
   ATL_UINT nb, kbS, nnu;
   size_t nfblks = ip->nfnblks;
   #ifdef TCPLX
      TYPE *rA, *rB;
   #endif
   if (d == nfblks + ip->npnblks - 1)
   {
      nb = ip->nF;
      #ifdef TCPLX
         if (d < nfblks)
         {
            rA = wA + ip->szA;
            rB = wB + ip->szB;
         }
         else
         {
            rA = wA + ip->pszA;
            rB = wB + ip->pszB;
         }
      #endif
   }
   else if (d < nfblks)
   {
      nb = ip->nb;
      #ifdef TCPLX
         rA = wA + ip->szA;
         rB = wB + ip->szB;
      #endif
   }
   else
   {
      nb = ip->pnb;
      #ifdef TCPLX
         rA = wA + ip->pszA;
         rB = wB + ip->pszB;
      #endif
   }
   nnu = (nb+ATL_SYRKK_NU-1)/ATL_SYRKK_NU;
   kbS = ((kb+ATL_SYRKK_KU-1)/ATL_SYRKK_KU)*ATL_SYRKK_KU;
   if (A)  /* want to copy input array */
   {
      A = IdxA_ip(ip, A, d, k);
      if (wA)  /* want to copy A to gemm storage too! */
      {
         #ifdef TCPLX
            ip->a2blk(kb, nb, ip->alpA, A, ip->lda, rA, wA);
         #else
            ip->a2blk(kb, nb, ip->alpA, A, ip->lda, wA);
         #endif
         if (ip->a2blk == sy2blk)
         {
            wS = wA;
            #ifdef TCPLX
               rS = rA;
            #endif
         }
      }
      if (wB)  /* want to copy At to gemm storage too! */
      {
         #ifdef TCPLX
            ip->b2blk(kb, nb, ip->alpB, A, ip->ldb, rB, wB);
         #else
            ip->b2blk(kb, nb, ip->alpB, A, ip->ldb, wB);
         #endif
         if (ip->b2blk == sy2blk)
         {
            wS = wB;
            #ifdef TCPLX
               rS = rB;
            #endif
         }
      }
      if (wS != wA && wS != wB)
      {
         #ifdef TCPLX
            sy2blk(kb, nb, ip->ONE, A, ip->lda, rS, wS);
         #else
            sy2blk(kb, nb, ATL_rone, A, ip->lda, wS);
         #endif
      }
   }
   if (wC && wS)  /* want to compute SYRK on this block into wC */
   {
      #ifdef TCPLX
         #ifdef Conj_
            TYPE *crA=(flag&2)?rS:wS, *ciA=(flag&2)?wS:rS;
            if (flag&4)
            {
               Mjoin(PATL,amsyrkK_b0)(nnu, nnu, kbS, wS, wS, rC, crA, ciA, wC);
               Mjoin(PATL,amsyrkK_b0)(nnu, nnu, kbS, crA, ciA, wC, rS, rS, rC);
            }
            else
            {
               Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kbS, wS, wS, rC, crA, ciA, wC);
               Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kbS, crA, ciA, wC, rS, rS, rC);
            }
            Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kbS, rS, rS, rC, ciA, crA, wC);
            Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kbS, ciA, crA, wC, wS, wS, rC);
         #else
            if (flag&4)
            {
               Mjoin(PATL,amsyrkK_b0)(nnu, nnu, kbS, wS, wS, rC, rS, wS, wC);
               Mjoin(PATL,amsyrkK_b0)(nnu, nnu, kbS, rS, wS, wC, rS, rS, rC);
            }
            else
            {
               Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kbS, wS, wS, rC, rS, wS, wC);
               Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kbS, rS, wS, wC, rS, rS, rC);
            }
            Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kbS, rS, rS, rC, wS, rS, wC);
            Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kbS, wS, rS, wC, wS, wS, rC);
         #endif
      #else
         if (flag&4)
            Mjoin(PATL,amsyrkK_b0)(nnu, nnu, kbS, wS, wS, wC, wS, wS, wC);
         else
            Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kbS, wS, wS, wC, wS, wS, wC);
      #endif
   }
   if (blk2c)
   {
      const size_t ldc=ip->ldc;
      #ifdef TCPLX
         const TYPE *alp = ip->alpC;
      #else
         TYPE alp = ip->alpC;
      #endif
      if (ip->alpA != ip->alpB)
         alp = (ip->alpA != ip->alpC) ? ip->alpA : ip->alpB;
      if (C)
      {
         C = IdxC_ip(ip, C, d, d);
         #ifdef TCPLX
            blk2c(nb, nb, alp, rC, wC, beta, C, ldc);
            #ifdef Conj_  /* must zero complex part of diagonal! */
               Mjoin(PATLU,zero)(nb, C+1, (ldc+1)SHIFT);
            #endif
         #else
            blk2c(nb, nb, alp, wC, beta, C, ldc);
         #endif
      }
   }
   #if 0
   if ((flag&1) && C && wU)  /* need to reflect from Lower to Upper */
   {
      size_t ldc2=(ip->ldc)SHIFT;
      #ifdef Conj_
         for (k=0; k < nb; k++, C += ldc2, wU += 2)
         {
            Mjoin(PATL,axpbyConj)(k+1, ip->ONE, wU, nb, beta, C, 1);
            C[((k+1)<<1)-1] = ATL_rzero;
         }
      #else
         if (SCALAR_IS_ZERO(beta))
            for (k=0; k < nb; k++, wU += (1 SHIFT), C += ldc2)
               Mjoin(PATL,copy)(k+1, wU, nb, C, 1);
         else
            for (k=0; k < nb; k++, wU += (1 SHIFT), C += ldc2)
            #ifdef TCPLX
               Mjoin(PATL,axpby)(k+1, ip->ONE, wU, nb, beta, C, 1);
            #else
               Mjoin(PATL,axpby)(k+1, ATL_rone, wU, nb, beta, C, 1);
            #endif
      #endif
   }
   #endif
}

/*
 * This helper function computes one block of C by looping over the K dim.
 * flag bits, meaning if set (opposite if unset):
 * 0/1: C is upper
 * 1/2: TA == AtlasNoTrans
 */
void syrk_K  /* inner-product based syrk/herk loop over K loop */
(
   ipinfo_t *ip,
   int flag,      /* bitvec: 0: set means C is upper, 1: set TA==AtlasNoTrans */
   size_t d,      /* which global diagonal blk of C is being computed */
   const TYPE *A, /* if non-NULL, base A ptr to copy */
   cm2am_t sy2blk,/* copy A to syrk storage */
   ablk2cmat_t blk2c, /* if non-NULL sy storage to C storage copy func */
   const SCALAR beta, /* only needed if blk2c non-NULL */
   TYPE *C,       /* if blk2c non-NULL, which C to write to, else ignored */
   TYPE *rS,      /* real portion of wS (unused in real routines) */
   TYPE *wS,      /* space to store syrk A; */
   TYPE *wA,      /* if non-NULL, ip-based A workspace */
   TYPE *wB,      /* if non-NULL ip-based At workspace */
   TYPE *rC,      /* real portion of wC (unused in real routines) */
   TYPE *wC,      /* if non-NULL: ptr to syrk-storage C wrkspc */
   TYPE *wU       /* NBxNB wrkspc needed for Upper C storage & blk2c != NULL */
)
{
   const size_t nfkblks = ip->nfkblks;
   size_t k;
   ATL_UINT szA, szB;
   if (d < ip->nfnblks)
      szA = szB = ip->szA;
   else
   {
      szA = ip->pszA;
      szB = ip->pszB;
   }
   #ifdef TCPLX
      szA += szA;
      szB += szB;
   #endif
/*
 * For first K block, use beta=0 kernel to init wC
 */
   if (nfkblks)
   {
      TYPE *wAn=(wA)? wA+szA:NULL, *wBn=(wB) ? wB+szB : NULL;
      syrkBlk(ip, flag|4, d, 0, A, sy2blk, NULL, beta, NULL, rS, wS,
              wA, wAn, wB, wBn, rC, wC, wU);
      wA = wAn;
      wB = wBn;
   }
   else /* this is first & last block! */
   {
      syrkBlk(ip, flag|4, d, 0, A, sy2blk, blk2c, beta, C, rS, wS,
              wA, wA, wB, wB, rC, wC, wU);
      return;
   }
/*
 * Handle all blocks except first (handled above) & last (handled below)
 */
   for (k=1; k < nfkblks; k++)
   {
      TYPE *wAn=(wA)? wA+szA:NULL, *wBn=(wB) ? wB+szB : NULL;
      syrkBlk(ip, flag, d, k, A, sy2blk, NULL, beta, NULL, rS, wS,
              wA, wAn, wB, wBn, rC, wC, wU);
      wA = wAn;
      wB = wBn;
   }
/*
 * Last block actually writes to C
 */
   syrkBlk(ip, flag, d, k, A, sy2blk, blk2c, beta, C, rS, wS,
           wA, wA, wB, wB, rC, wC, wU);
}

int syrk
(
   const enum ATLAS_UPLO  Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT  N,
   ATL_CSZT K,
   #ifdef Conj_
   const TYPE ralpha,
   #else
   const SCALAR alpha,
   #endif
   const TYPE *A,
   ATL_CSZT lda,
   #ifdef Conj_
   const TYPE rbeta,
   #else
   const SCALAR beta,
   #endif
   TYPE *C,
   ATL_CSZT ldc
)
{
   size_t sz, szA, szB, szC, szS, nnblks, extra;
   void *vp=NULL;
   TYPE *wA, *wB, *wC, *wS, *wCs, *rC, *rCs, *rS;
   double timG;
   int nb, nbS, flg, idx;
   int ierr;
   cm2am_t sy2blk;
   ablk2cmat_t blk2sy, blk2c;
   ipinfo_t ip;
   #ifdef Conj_
      TYPE alpha[2]={ralpha, ATL_rzero}, beta[2]={rbeta, ATL_rzero};
      const enum ATLAS_TRANS TB=(TA==AtlasNoTrans)?AtlasConjTrans:AtlasNoTrans;
   #else
      const enum ATLAS_TRANS TB = (TA==AtlasNoTrans) ? AtlasTrans:AtlasNoTrans;
   #endif
   Mjoin(PATL,ipmenInfo)(&ip, TA, TB, N,N,K, lda, lda, ldc, alpha, beta);
   #if 0
   printf("D=(%u,%u), B=(%u,%u,%u), b=(%u,%u,%u), nB=(%u,%u,%u)\n",
          (unsigned int)N,(unsigned int)K, ip.mb, ip.nb, ip.kb,
          ip.pmb, ip.pnb, ip.kb0, ip.nfmblks, ip.nfnblks, ip.nfkblks);
   #endif
/*
 * Will eventually need syrk timed for all square blocks to select best case.
 * For now, just pretend syrk time doesn't matter
 */
/*
 * Need space for only one column-panel of At
 */
   szB = ip.nfnblks ? ip.szA : ip.pszA;
   szB *= (ip.nfkblks+1);
/*
 * A needs entire matrix minus one row/col panel
 */
   szA = szB * (ip.nfnblks + ip.npnblks - 1);
   nb = (ip.nfnblks) ? ip.nb : ip.pnb;
   nbS = (nb+ATL_SYRKK_NU-1)/ATL_SYRKK_NU;
   szC = ((nbS+1)*nbS)>>1;  /* only need lower tri blks, not full nnu*nnu */
   nbS *= ATL_SYRKK_NU;
   szC *= ((ATL_SYRKK_NU*ATL_SYRKK_NU+ATL_SYRKK_VLEN-1)/ATL_SYRKK_VLEN)
          * ATL_SYRKK_VLEN;
   extra = ip.exsz;
   extra = Mmax(extra, (ATL_SYRKK_NU+ATL_SYRKK_NU)*ATL_SYRKK_NU);
   szS = ((ip.kb + ATL_SYRKK_KU-1)/ATL_SYRKK_KU)*ATL_SYRKK_KU;
   szS *= nbS;
   sz = Mmax(ip.szC,szC);
   sz = ATL_MulBySize(sz + szS + szA+szB + extra) + 5*ATL_Cachelen;
   if (sz < ATL_MaxMalloc)
      vp = malloc(sz);
   if (!vp)
      return(1);  /* keep recurring, can't malloc space! */
   wA = ATL_AlignPtr(vp);
   wB = wA + (szA SHIFT);
   wB = ATL_AlignPtr(wB);
   wS = wB + (szB SHIFT);
   wS = ATL_AlignPtr(wS);
   wC = wS + (szS SHIFT);
   wC = ATL_AlignPtr(wC);
   if (!szA)
      wA = NULL;
   wCs = wC;
   #ifdef TCPLX
      rC = wC + ip.szC;
      rCs = wC + szC;
      rS = wS + szS;
   #else
      rCs = rC = wC;
      rS = wS;
   #endif
/*
 * ============================================================================
 * First, we compute diagonals of C, and in the process we will copy A/A^T
 * for use in computing non-diaginal blocks using inner-product amm.
 * ============================================================================
 */
   sy2blk = IS_COLMAJ(TA) ? Mjoin(PATL,a2blk_syrkT) : Mjoin(PATL,a2blk_syrkN);

   if (Uplo == AtlasLower)
   {
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNb1);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1b1):Mjoin(PATL,SyrkIntoC_aXb1);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNbN);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1bN):Mjoin(PATL,SyrkIntoC_aXbN);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNb0);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1b0):Mjoin(PATL,SyrkIntoC_aXb0);
      }
      else
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNbX);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1bX):Mjoin(PATL,SyrkIntoC_aXbX);
      }
   }
   else
   {
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNb1_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1b1_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXb1_L2UT);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNbN_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1bN_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXbN_L2UT);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNb0_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1b0_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXb0_L2UT);
      }
      else
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNbX_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1bX_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXbX_L2UT);
      }
   }
   nnblks = ip.nfnblks + ip.npnblks;
   flg = (TA == AtlasNoTrans) ? 2 : 0;
/*
 * Upper doesn't need to copy first row panel of A or last col panel of At
 * to GEMM storage.  If C only block, should call syrk_IP instead!
 */
   if (Uplo == AtlasUpper)
   {
      flg |= 1;
      syrk_K(&ip, flg, 0, A, sy2blk, blk2sy, beta, C, rS, wS, wB, NULL,
             rCs, wC, wCs);
      if (N > nb)
      {
         const size_t incC = nb*((ldc+1)SHIFT);
         size_t i;
/*
 *       Compute all C blks within this rowpan of C, copy rest of At
 */
         blk2c = ip.blk2c;
         Mjoin(PATL,iploopsNK)(&ip, 0, 1, NULL, A, C, 11, wB, wA, rC, wC,
                               beta, blk2c);
/*
 *       Loop over all col-pans of C, excepting first & last;
 *       A already fully copied, B will be copied by syrk_Kloop call.
 */
         for (i=1; i < nnblks-1; i++)
         {
            syrk_K(&ip, flg, i, A, sy2blk, blk2sy, beta, C, rS, wS, wB, NULL,
                   rCs, wC, wCs);
            wA += szB SHIFT;
            Mjoin(PATL,iploopsNK)(&ip, i, i+1, NULL, NULL, C+i*(nb SHIFT), 11,
                                  wB, wA, rC, wC, beta, blk2c);
         }
/*
 *       Last colpan is only 1 diag blk, so don't copy A or B for gemm
 */
         syrk_K(&ip, flg, i, A, sy2blk, blk2sy, beta, C, rS, wS, NULL, NULL,
                rCs, wC, wCs);
      }
   }
/*
 * Lower doesn't need to copy first row panel of A or last col panel of At
 * to GEMM storage.  If C only block, should call syrk_IP instead!
 */
   else
   {
      const size_t incC = (ldc SHIFT) * nb;
      size_t j;
      TYPE *c;
/*
 *    Compute first diag block, copying At to wB for use by all of this colpan
 *    of C's gemm computation
 */
      syrk_K(&ip, flg, 0, A, sy2blk, blk2sy, beta, C, rS, wS, NULL, wB,
             rCs, wC, wCs);
      if (N > nb)
      {
/*
 *       Compute all C blks within this colpan of C, copying rest of A
 */
         blk2c = ip.blk2c;
         Mjoin(PATL,iploopsMK)(&ip, 1, 0, A, NULL, C, 7, wA, wB, rC, wC,
                               beta, blk2c);
/*
 *       Loop over all col-pans of C, excepting first & last;
 *       A already fully copied, B will be copied by syrk_Kloop call.
 */
         c = C + incC;
         for (j=1; j < nnblks-1; j++)
         {
            syrk_K(&ip, flg, j, A, sy2blk, blk2sy, beta, C, rS, wS,
                   NULL, wB, rCs, wC, wCs);
            wA += szB SHIFT;
            Mjoin(PATL,iploopsMK)(&ip, j+1, j, NULL, NULL, c, 7, wA, wB, rC, wC,
                                  beta, blk2c);
            c += incC;
         }
/*
 *       Last colpan is only 1 diag blk, so don't copy A or B for gemm
 */
         syrk_K(&ip, flg, j, A, sy2blk, blk2sy, beta, C, rS, wS,
                NULL, NULL, rCs, wC, wCs);
      }
   }
   free(vp);
   return(0);
}
