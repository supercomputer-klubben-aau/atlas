/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_level2.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_sqsyrk.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_umsyrk.h))
#ifdef Conj_
   #define sqsyrkBlk Mjoin(PATL,sqherkBlk_OP)
   #define umsyrkBlk Mjoin(PATL,umherkBlk_OP)
   #define opsqsyrk Mjoin(PATL,opsqherk)
   #define opumsyrk Mjoin(PATL,opumherk)
   #define opsyrk Mjoin(PATL,opherk)
#else
   #define sqsyrkBlk Mjoin(PATL,sqsyrkBlk_OP)
   #define umsyrkBlk Mjoin(PATL,umsyrkBlk_OP)
   #define opsqsyrk Mjoin(PATL,opsqsyrk)
   #define opumsyrk Mjoin(PATL,opumsyrk)
   #define opsyrk Mjoin(PATL,opsyrk)
#endif
/*
 * Indexes both A and C from base ptrs according to d (diagonal blk)
 */
void sqsyrkBlk
(
   opinfo_t *op,
   int flag,      /* bitvec: 0: set means C is upper, 1: set TA==AtlasNoTrans */
   size_t d,      /* which global diagonal blk of C is being computed */
   const TYPE *A, /* if non-NULL, base A ptr to copy */
   cm2am_t sy2blk,/* copy A to syrk storage */
   ablk2cmat_t blk2c, /* if non-NULL sy storage to C storage copy func */
   TYPE *C,       /* if blk2c non-NULL, which C to write to, else ignored */
   TYPE *rS,      /* real ptr (unused for real types) */
   TYPE *wS,      /* space to store syrk A; */
   TYPE *wA,      /* if non-NULL, op-based A workspace */
   TYPE *wAn,     /* next A wrkspc to be prefetched */
   TYPE *wB,      /* if non-NULL op-based At workspace */
   TYPE *wBn,     /* next B wrkspc to be prefetched */
   TYPE *rC,      /* real ptr (unused for real routs) */
   TYPE *wC,      /* if non-NULL: ptr to syrk-storage C wrkspc */
   TYPE *wU       /* NBxNB wrkspc needed for Upper C storage & blk2c != NULL */
)
{
   ATL_CSZT lda = op->lda, nfblks = op->nfnblks;
   ATL_CUINT KB = op->KB, kb = op->kb;
   ATL_UINT nb, kbS, nnu;
   #ifdef TCPLX
      ATL_UINT szC, szA;
      TYPE *rA, *rB;
   #endif
   const SCALAR beta=op->beta;
   if (d == nfblks + op->npnblks - 1)  /* last block is SYRK only */
   {
      nb = op->nF;
      nb = (nb) ? nb : op->nb;
      #ifdef TCPLX
         if (d < nfblks)
         {
            rA = wA + op->szA;
            rB = wB + op->szB;
         }
         else
         {
            rA = wA + op->pszA;
            rB = wB + op->pszB;
         }
      #endif
   }
   else if (d < nfblks)
   {
      nb = op->nb;
      #ifdef TCPLX
         rA = wA + op->szA;
         rB = wB + op->szB;
      #endif
   }
   else
   {
      nb = op->pnb;
      #ifdef TCPLX
         rA = wA + op->pszA;
         rB = wB + op->pszB;
      #endif
   }
   nnu = (nb+ATL_SQSYRKK_NU-1)/ATL_SQSYRKK_NU;
   kbS = ((kb+ATL_SQSYRKK_KU-1)/ATL_SQSYRKK_KU)*ATL_SQSYRKK_KU;
   if (A)  /* want to copy input array */
   {
/*
 *    Move A ptr to d'th block
 */
      if (d)
      {
         size_t n = Mmin(d, nfblks);
         A += n*op->incAm;
         n = d - n;  /* # of partial blocks remaining in d */
         A += n*op->pincAm;
      }
      if (wA)  /* want to copy A to gemm storage too! */
      {
         #ifdef TCPLX
            op->a2blk(kb, nb, op->alpA, A, lda, rA, wA);
         #else
            op->a2blk(kb, nb, op->alpA, A, lda, wA);
         #endif
         if (op->a2blk == sy2blk)
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
            op->b2blk(kb, nb, op->alpB, A, lda, rB, wB);
         #else
            op->b2blk(kb, nb, op->alpB, A, lda, wB);
         #endif
         if (op->b2blk == sy2blk)
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
            sy2blk(kb, nb, op->ONE, A, lda, rS, wS);
         #else
            sy2blk(kb, nb, ATL_rone, A, lda, wS);
         #endif
      }
   }
   if (wC)  /* want to compute SYRK on this block into wC */
   {
      #ifdef TCPLX
         #ifdef Conj_
            TYPE *crA=(flag&2)?rS:wS, *ciA=(flag&2)?wS:rS;
            Mjoin(PATL,sqsyrkK_b0)(nnu, nnu, kbS, wS, wS, rC, crA, ciA, wC);
            Mjoin(PATL,sqsyrkK_b0)(nnu, nnu, kbS, crA, ciA, wC, rS, rS, rC);
            Mjoin(PATL,sqsyrkK_b1)(nnu, nnu, kbS, rS, rS, rC, ciA, crA, wC);
            Mjoin(PATL,sqsyrkK_bn)(nnu, nnu, kbS, ciA, crA, wC, wS, wS, rC);
         #else
            Mjoin(PATL,sqsyrkK_b0)(nnu, nnu, kbS, wS, wS, rC, rS, wS, wC);
            Mjoin(PATL,sqsyrkK_b0)(nnu, nnu, kbS, rS, wS, wC, rS, rS, rC);
            Mjoin(PATL,sqsyrkK_bn)(nnu, nnu, kbS, rS, rS, rC, wS, rS, wC);
            Mjoin(PATL,sqsyrkK_b1)(nnu, nnu, kbS, wS, rS, wC, wS, wS, rC);
         #endif
      #else
         Mjoin(PATL,sqsyrkK_b0)(nnu, nnu, kbS, wS, wS, wC, wS, wS, wC);
      #endif
   }
   if (blk2c)
   {
      const size_t ldc=op->ldc;
      int k;
      #ifdef TCPLX
         const TYPE *alp = (op->alpA == op->ONE) ? op->alpB : op->alpA;
      #else
         TYPE alp = (op->alpA == ATL_rone) ? op->alpB : op->alpA;
      #endif
      C += d*op->nb*((ldc+1)SHIFT);
      if (flag&1)  /* Upper matrix */
      {
         #ifdef TCPLX
            blk2c(nb, nb, alp, rC, wC, beta, C, ldc);
            #ifdef Conj_  /* must zero complex part of diagonal! */
               Mjoin(PATLU,zero)(nb, C+1, (ldc+1)SHIFT);
            #endif
         #else
            blk2c(nb, nb, alp, wC, beta, C, ldc);
         #endif
      }
      else
      {
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
}
int opsqsyrk
(
   opinfo_t *op,
   const enum ATLAS_UPLO  Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT  N,
   ATL_CSZT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   size_t sz, szA, szC, szU, szS, nnblks, szAblk;
   void *vp;
   TYPE *wA, *wB, *wC, *wS, *wU, *rC, *rCs, *rS;
   int nb, nbS, flg, idx, extra;
   cm2am_t sy2blk;
   ablk2cmat_t blk2sy, blk2c;

   nnblks = op->nfnblks + op->npnblks;
   sy2blk = IS_COLMAJ(TA) ? Mjoin(PATL,a2blk_sqsyrkT) : Mjoin(PATL,a2blk_sqsyrkN);

   if (Uplo == AtlasLower)
   {
      szU = 0;
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,sqSyrkIntoC_aNb1);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,sqSyrkIntoC_a1b1):Mjoin(PATL,sqSyrkIntoC_aXb1);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,sqSyrkIntoC_aNbN);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,sqSyrkIntoC_a1bN):Mjoin(PATL,sqSyrkIntoC_aXbN);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,sqSyrkIntoC_aNb0);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,sqSyrkIntoC_a1b0):Mjoin(PATL,sqSyrkIntoC_aXb0);
      }
      else
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,sqSyrkIntoC_aNbX);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,sqSyrkIntoC_a1bX):Mjoin(PATL,sqSyrkIntoC_aXbX);
      }
   }
   else  /* C is Upper */
   {
      szU = 0;
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,sqSyrkIntoC_aNb1_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                        Mjoin(PATL,sqSyrkIntoC_a1b1_L2UT)
                        :Mjoin(PATL,sqSyrkIntoC_aXb1_L2UT);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,sqSyrkIntoC_aNbN_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                        Mjoin(PATL,sqSyrkIntoC_a1bN_L2UT)
                        :Mjoin(PATL,sqSyrkIntoC_aXbN_L2UT);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,sqSyrkIntoC_aNb0_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                        Mjoin(PATL,sqSyrkIntoC_a1b0_L2UT)
                        :Mjoin(PATL,sqSyrkIntoC_aXb0_L2UT);
      }
      else
      {
         if (SCALAR_IS_NONE(alpha))
            blk2sy = Mjoin(PATL,sqSyrkIntoC_aNbX_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(alpha) ?
                        Mjoin(PATL,sqSyrkIntoC_a1bX_L2UT)
                        :Mjoin(PATL,sqSyrkIntoC_aXbX_L2UT);
      }
   }
   extra = (ATL_SQSYRKK_NU+ATL_SQSYRKK_NU)*ATL_SQSYRKK_NU;
   flg = (TA == AtlasNoTrans) ? 2 : 0;
   if (TA == AtlasUpper)
   {
      flg |= 1;
      extra -= Mmin(extra, szU);
   }
   if (nnblks == 1)  /* we've got a 1 block of SYRK only! */
   {
      nbS = (N+ATL_SQSYRKK_NU-1)/ATL_SQSYRKK_NU;
      szC = ((nbS+1)*nbS)>>1;  /* only need lower tri blks, not full nnu*nnu */
      nbS *= ATL_SQSYRKK_NU;
      szC *= ((ATL_SQSYRKK_NU*ATL_SQSYRKK_NU+ATL_SQSYRKK_VLEN-1)/ATL_SQSYRKK_VLEN)
             * ATL_SQSYRKK_VLEN;
      #if ATL_SQSYRKK_KVEC > 1
         szS = ((K+ATL_SQSYRKK_KVEC-1)/ATL_SQSYRKK_KVEC)*ATL_SQSYRKK_KVEC;
         szS *= nbS;
      #else
         szS = nbS * K;
      #endif
      szU = Mmax(op->szC, szU);  /* Majedul: don't need anymore */
      sz = ATL_MulBySize(szU + szC + szS + extra) + 3*ATL_Cachelen;
      vp = malloc(sz);
      if (!vp)
         return(1);
      wS = ATL_AlignPtr(vp);
      wC = wS + (szS SHIFT);
      wC = ATL_AlignPtr(wC);
      wU = wC + (szC SHIFT);
      wU = ATL_AlignPtr(wU);
      #ifdef TCPLX
         rC = wC + szC;
         rS = wS + szS;
      #else
         rC = wC;
         rS = wS;
      #endif
      flg |= (Uplo == AtlasUpper) ? 1 : 0;
      sqsyrkBlk(op, flg, 0, A, sy2blk, blk2sy, C, rS, wS,
              NULL, NULL, NULL, NULL, rC, wC, wU);
      free(vp);
      return(0);
   }
/*
 * If we reach here, we have at rank-K SYRK update requiring both SYRK & GEMM
 * Since nnblks > 1, nfnblks > 1 as well.
 */
   nbS = (op->nb+ATL_SQSYRKK_NU-1)/ATL_SQSYRKK_NU;
   szC = ((nbS+1)*nbS)>>1;  /* only need lower tri blks, not full nnu*nnu */
   nbS *= ATL_SQSYRKK_NU;
   szC *= ((ATL_SQSYRKK_NU*ATL_SQSYRKK_NU+ATL_SQSYRKK_VLEN-1)/ATL_SQSYRKK_VLEN)
          * ATL_SQSYRKK_VLEN;
   #if ATL_SQSYRKK_KVEC > 1
      szS = ((K+ATL_SQSYRKK_KVEC-1)/ATL_SQSYRKK_KVEC)*ATL_SQSYRKK_KVEC;
      szS *= nbS;
   #else
      szS = nbS * K;
   #endif
   szU = Mmax(op->szC, szU);  /* Majedul: no need anymore */
   szAblk = op->szA;
   szA = szAblk * (nnblks-1);
   /*
    * FIXED: wC is used both for syrk and gemm.
    * We need to allocate Mmax(szC, op->szC) for wC
    */
   sz = Mmax(szC, op->szC);
   sz = ATL_MulBySize(szA+szAblk + szS + sz + szU + extra) + 5*ATL_Cachelen;
   vp = malloc(sz);
   if (!vp)
      return(1);
   wA = ATL_AlignPtr(vp);
   wB = wA + (szA SHIFT);
   wB = ATL_AlignPtr(wB);
   wS = wB + (szAblk SHIFT);
   wS = ATL_AlignPtr(wS);
   wC = wS + (szS SHIFT);
   wC = ATL_AlignPtr(wC);
   wU = wC + (szC SHIFT);
   wU = ATL_AlignPtr(wU);
   #ifdef TCPLX
      rC = wC + op->szC;
      rCs = wC + szC;
      rS = wS + szS;
   #else
      rC = wC;
      rCs = wC;
      rS = wS + szS;
   #endif
   if (Uplo == AtlasLower)
   {
      size_t j;
/*
 *    Do first diagonal block, don't copy A blk to GEMM storage since it is
 *    used only for this diagonal (SYRK)
 */
      sqsyrkBlk(op, flg, 0, A, sy2blk, blk2sy, C, rS, wS, NULL, NULL, wB, wB,
              rCs, wC, wU);
/*
 *    Do rank-K update on ge blks beneath diag, cop->ing entire A at same time
 */
      Mjoin(PATL,oploopsM)(op, 1, 0, A, NULL, C, 1, wA, wB, rC, wC);
/*
 *    For remaining column panels of C, sqsyrkBlk copies B, reusues wA
 */
      for (j=1; j < nnblks-1; j++)
      {
         sqsyrkBlk(op, flg, j, A, sy2blk, blk2sy, C, rS, wS, NULL, NULL, wB, wB,
                 rCs, wC, wU);
         wA += (szAblk SHIFT);
         Mjoin(PATL,oploopsM)(op, j+1, j, NULL, NULL, C, 1, wA, wB, rC, wC);
      }
/*
 *    Last col panel is only one diagonal
 */
      sqsyrkBlk(op, flg, j, A, sy2blk, blk2sy, C, rS, wS, NULL, NULL, NULL, NULL,
              rCs, wC, wU);
   }
   else
   {
      size_t i;
/*
 *    Do first diagonal block, don't copy A^T blk to GEMM storage since it is
 *    used only for this diagonal (SYRK)
 */
      flg |= 1;
      sqsyrkBlk(op, flg, 0, A, sy2blk, blk2sy, C, rS, wS, wB, wB, NULL, NULL,
              rCs, wC, wU);
/*
 *    Do rank-K update on ge blks beneath diag, copying entire A^T at same time
 */
      Mjoin(PATL,oploopsN)(op, 0, 1, NULL, A, C, 2, wB, wA, rC, wC);
/*
 *    For remaining column panels of C, sqsyrkBlk copies A, reuses A^T
 */
      for (i=1; i < nnblks-1; i++)
      {
         sqsyrkBlk(op, flg, i, A, sy2blk, blk2sy, C, rS, wS, wB, wB, NULL, NULL,
                 rCs, wC, wU);
         wA += (szAblk SHIFT);
         Mjoin(PATL,oploopsN)(op, i, i+1, NULL, NULL, C, 2, wB, wA, rC, wC);
      }
/*
 *    Last col panel is only one diagonal
 */
      sqsyrkBlk(op, flg, i, A, sy2blk, blk2sy, C, rS, wS, NULL, NULL, NULL, NULL,
              rCs, wC, wU);
   }
   free(vp);
   return(0);
}
void umsyrkBlk
(
   opinfo_t *op,
   int flag,      /* bitvec: 0: set means C is upper, 1: set TA==AtlasNoTrans */
   size_t d,      /* which global diagonal blk of C is being computed */
   const TYPE *A, /* if non-NULL, base A ptr to copy in wA */
   const TYPE *B, /* if non-NULL, base A ptr to copy in wB */
   ablk2cmat_t blk2c, /* if non-NULL sy storage to C storage copy func */
   TYPE *C,       /* if blk2c non-NULL, which C to write to, else ignored */
   TYPE *wA,      /* if non-NULL, op-based A workspace */
   TYPE *wAn,     /* next A wrkspc to be prefetched */
   TYPE *wB,      /* if non-NULL op-based At workspace */
   TYPE *wBn,     /* next B wrkspc to be prefetched */
   TYPE *rC,      /* real ptr (unused for real routs) */
   TYPE *wC      /* if non-NULL: ptr to syrk-storage C wrkspc */
)
{
   ATL_CSZT lda = op->lda, nfblks = op->nfnblks;
   ATL_CSZT ldb = op->ldb;
   ATL_CUINT KB = op->KB, kb = op->kb;
   ATL_UINT nb, kbS, nnu, nmu;
#ifdef TCPLX
   ATL_UINT szC, szA;
   TYPE *rA, *rB;
   TYPE *iA=wA, *iB=wB; /* may not need ... wA, wB can be used !!!*/
#endif
   const SCALAR beta=op->beta;
   TYPE *wS;
   /*assert(op->mb==op->nb);*/
   if (d == nfblks + op->npnblks - 1) /* last syrk only block*/
   {
      nb = op->nF;
      nb = (nb) ? nb : op->nb;
      #ifdef TCPLX
         if (d < nfblks)
         {
            rA = wA + op->szA;
            rB = wB + op->szB;
         }
         else
         {
            rA = wA + op->pszA;
            rB = wB + op->pszB;
         }
      #endif
   }
   else if (d < nfblks)
   {
      nb = op->nb;
      #ifdef TCPLX
         rA = wA + op->szA;
         rB = wB + op->szB;
      #endif
   }
   else
   {
      nb = op->pnb;
      #ifdef TCPLX
         rA = wA + op->pszA;
         rB = wB + op->pszB;
      #endif
   }
/*
 * calc nmu and nnu. one must be multiple of other
 */
   if (op->mu > op->nu)
   {
      nmu = (nb+op->mu-1)/op->mu;
      nnu = nmu * (op->mu/op->nu); /* NU must be multiple of MU */
   }
   else
   {
      nnu = (nb+op->nu-1)/op->nu;
      nmu = nnu * (op->nu/op->mu); /* NU must be multiple of MU */
   }
/*
 * we can use KB since we are using same kernel as gemm
 */
   kbS = op->KB;

   if (A)  /* want to copy input array in wA */
   {
/*
 *    Move A ptr to d'th block
 */
      if (d)
      {
         size_t n = Mmin(d, nfblks);
         A += n*op->incAm;
         n = d - n;  /* # of partial blocks remaining in d */
         A += n*op->pincAm;
      }
      #ifdef TCPLX
         op->a2blk(kb, nb, op->alpA, A, lda, rA, wA);
      #else
         op->a2blk(kb, nb, op->alpA, A, lda, wA);
      #endif
   }
   if (B) /* want to copy input array in wB */
   {
      if (d)
      {
         size_t n = Mmin(d, nfblks);
         B += n*op->incBn;
         n = d - n;  /* # of partial blocks remaining in d */
         B += n*op->pincBn;
      }

      #ifdef TCPLX
         op->b2blk(kb, nb, op->alpB, B, ldb, rB, wB);
      #else
         op->b2blk(kb, nb, op->alpB, B, ldb, wB);
      #endif
   }
   if (wC)  /* want to compute SYRK on this block into wC */
   {
/*
 *    NOTE: no need to specially handle herk/conj, since we are using gemm's
 *    copy routines and similar kernels
 */
      #ifdef TCPLX
         Mjoin(PATL,umsyrkK_b0)(nmu, nnu, kbS, wA, wB, rC, rA, wB, wC);
         Mjoin(PATL,umsyrkK_b0)(nmu, nnu, kbS, rA, wB, wC, rA, rB, rC);
         Mjoin(PATL,umsyrkK_bn)(nmu, nnu, kbS, rA, rB, rC, wA, rB, wC);
         Mjoin(PATL,umsyrkK_b1)(nmu, nnu, kbS, wA, rB, wC, wA, wB, wC);
      #else
         Mjoin(PATL,umsyrkK_b0)(nmu, nnu, kbS, wA, wB, wC, wA, wB, wC);
      #endif
   }
   if (blk2c)
   {
      const size_t ldc=op->ldc;
      int k;
      #ifdef TCPLX
         const TYPE *alp = (op->alpA == op->ONE) ? op->alpB : op->alpA;
      #else
         TYPE alp = (op->alpA == ATL_rone) ? op->alpB : op->alpA;
      #endif
      C += d*op->nb*((ldc+1)SHIFT);
      if (flag&1)  /* Upper matrix */
      {
         #ifdef TCPLX
            blk2c(nb, nb, alp, rC, wC, beta, C, ldc);
            #ifdef Conj_  /* must zero complex part of diagonal! */
               Mjoin(PATLU,zero)(nb, C+1, (ldc+1)SHIFT);
            #endif
         #else
            blk2c(nb, nb, alp, wC, beta, C, ldc);
         #endif
      }
      else
      {
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
}

int opumsyrk
(
   opinfo_t *op,
   const enum ATLAS_UPLO  Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT  N,
   ATL_CSZT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   size_t sz, szA, szC, szU, szS, nnblks, szAblk;
   void *vp;
   TYPE *wA, *wB, *wC, *wS, *wU, *rC, *rCs, *rS;
   int nb, nbS, flg, idx, extra;
   cm2am_t sy2blk;
   ablk2cmat_t blk2sy, blk2c;
   int smnu;

   nnblks = op->nfnblks + op->npnblks;
   smnu = Mmax(op->mu, op->nu);           /* max unrolling */
/*
 * select c-copy
 */
   if (Uplo == AtlasLower)
   {
      if (SCALAR_IS_ONE(beta))
         blk2sy = Mjoin(PATL,umSyrkIntoC_a1b1);
      else if (SCALAR_IS_NONE(beta))
         blk2sy = Mjoin(PATL,umSyrkIntoC_a1bN);
      else if (SCALAR_IS_ZERO(beta))
         blk2sy = Mjoin(PATL,umSyrkIntoC_a1b0);
      else
         blk2sy = Mjoin(PATL,umSyrkIntoC_a1bX);
   }
   else  /* C is Upper */
   {
      if (SCALAR_IS_ONE(beta))
         blk2sy = Mjoin(PATL,umSyrkIntoC_a1b1_L2UT);
      else if (SCALAR_IS_NONE(beta))
         blk2sy = Mjoin(PATL,umSyrkIntoC_a1bN_L2UT);
      else if (SCALAR_IS_ZERO(beta))
         blk2sy = Mjoin(PATL,umSyrkIntoC_a1b0_L2UT);
      else
         blk2sy = Mjoin(PATL,umSyrkIntoC_a1bX_L2UT);
   }
   flg = (TA == AtlasNoTrans) ? 2 : 0;
   if (TA == AtlasUpper)
      flg |= 1;
   if (nnblks == 1)  /* we've got a 1 block of SYRK only! */
   {
      extra = (smnu+smnu)*smnu;
      nbS = (N+smnu-1)/smnu;
      szC = ((nbS+1)*nbS)>>1;
      nbS *= smnu;
      szC *= ((smnu*smnu+op->ku-1)/op->ku)*op->ku;
      szS = nbS * op->KB;         /* KB is already rounded with ku */
      szA = szAblk = szS;
      sz = ATL_MulBySize(szA+szAblk+szC+extra) + 3*ATL_Cachelen;
      vp = malloc(sz);
      if (!vp)
         return(1);
      wA = ATL_AlignPtr(vp);
      wB = wA + (szA SHIFT);
      wB = ATL_AlignPtr(wB);
      wC = wB + (szAblk SHIFT);
      wC = ATL_AlignPtr(wC);
      #ifdef TCPLX
         rC = wC + szC;
      #else
         rC = wC;
      #endif
      flg |= (Uplo == AtlasUpper) ? 1 : 0;
      umsyrkBlk(op, flg, 0, A, A, blk2sy, C, wA, wA, wB, wB, rC, wC);
      free(vp);
      return(0);
   }
/*
 * If we reach here, we have at rank-K SYRK update requiring both SYRK & GEMM
 * Since nnblks > 1, nfnblks > 1 as well.
 */
/*
 * syrk needs less space than gemm, so we can reuse same wC space
 */
   szC = op->szC;
   szAblk = op->szA;
   szA = szAblk * (nnblks-1);
   extra = op->exsz;           /* why needed?? */
   sz = ATL_MulBySize(szA + szAblk + szC + extra) + 3*ATL_Cachelen;

   vp = malloc(sz);
   if (!vp)
      return(1);
   wA = ATL_AlignPtr(vp);
   wB = wA + (szA SHIFT);
   wB = ATL_AlignPtr(wB);
   wC = wB + (szAblk SHIFT);
   wC = ATL_AlignPtr(wC);
   #ifdef TCPLX
      rC = wC + op->szC;
   #else
      rC = wC;
   #endif
   if (Uplo == AtlasLower)
   {
      size_t j;
/*
 *    SYRK uses same copy routines as GEMM. copy of B can be reused in GEMM
 *    of the column panels
 */
      umsyrkBlk(op, flg, 0, A, A, blk2sy, C, wA, wA, wB, wB, rC, wC);
/*
 *    Do rank-K update on ge blks beneath diag, copying entire A at same time
 */
      Mjoin(PATL,oploopsM)(op, 1, 0, A, NULL, C, 1, wA, wB, rC, wC);
/*
 *    For remaining column panels of C, syrkBlk copies B, re-uses wA
 */
      for (j=1; j < nnblks-1; j++)
      {
/*
 *       ****** big idea:
 *       1. Resue A from previous wA calculated in oploopsM last time.
 *          After the first call of oploopsM the whole wA is populated
 *       2. compute wB which can be reused in subsequent oploopsM
 */
         umsyrkBlk(op, flg, j, NULL, A, blk2sy, C, wA,wA, wB,wB, rC, wC);
         wA += (szAblk SHIFT);
         Mjoin(PATL,oploopsM)(op, j+1, j, NULL, NULL, C, 1, wA, wB, rC, wC);
      }
/*
 *    Last col panel is only one diagonal
 */
      umsyrkBlk(op, flg, j, NULL, A, blk2sy, C, wA, wA, wB, wB, rC, wC);
   }
   else
   {
      size_t i;
/*
 *    GEMM will reuse A*T, Syrk can reuse B later... wB is wA
 */
      flg |= 1;
      umsyrkBlk(op, flg, 0, A, A, blk2sy, C, wB, wB, wA, wA, rC, wC);
/*
 *    Do rank-K update on ge blks beneath diag, copying entire A^T at same time
 */
      Mjoin(PATL,oploopsN)(op, 0, 1, NULL, A, C, 2, wB, wA, rC, wC);
/*
 *    For remaining column panels of C, syrkBlk copies A, reuses A^T
 */
      for (i=1; i < nnblks-1; i++)
      {
         umsyrkBlk(op, flg, i, A, NULL, blk2sy, C, wB,wB, wA,wA, rC, wC);
         wA += (szAblk SHIFT);
         Mjoin(PATL,oploopsN)(op, i, i+1, NULL, NULL, C, 2, wB, wA, rC, wC);
      }
/*
 *    Last col panel is only one diagonal
 */
      umsyrkBlk(op, flg, i, A, NULL, blk2sy, C, wB, wB, wA, wA, rC, wC);
   }
   free(vp);
   return(0);
}

int opsyrk
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
   int i, ismul;
   int ierr;
   #ifdef Conj_
      TYPE alpha[2]={ralpha, ATL_rzero}, beta[2]={rbeta, ATL_rzero};
      const enum ATLAS_TRANS TB=(TA==AtlasNoTrans)?AtlasConjTrans:AtlasNoTrans;
   #else
      const enum ATLAS_TRANS TB = (TA==AtlasNoTrans) ? AtlasTrans:AtlasNoTrans;
   #endif
   opinfo_t op;
/*
 * find opsyrkinfo
 */
   #ifdef Conj_
      if (Mjoin(PATL,opsyrkInfo)(&op, 1, TA, N, K, lda, ldc, alpha, beta))
   #else
      if (Mjoin(PATL,opsyrkInfo)(&op, 0, TA, N, K, lda, ldc, alpha, beta))
   #endif
         return(2);
   #if 0
   fprintf(stderr, "D=(%u,%u), B=(%u,%u,%u), pB=(%u,%u,%u)\n",
           (unsigned int)N, (unsigned int)K, op.mb, op.nb, op.KB,
           op.pmb, op.pnb, op.kb);
   #endif
/*
 * mu,nu,ku, must match with amm_umsyrk.h to apply umsyrk until we manage
 * to have opgenview for UM
 */
   ismul = 1;
   if (op.mu != ATL_UMSYRKK_MU || op.nu != ATL_UMSYRKK_NU
         || op.ku !=ATL_UMSYRKK_KU )
      ismul = 0;
/*
 * try to use umsyrk first if possible
 */
   if (ismul)
      ierr = opumsyrk(&op, Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
   else
      ierr = opsqsyrk(&op, Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
   return(ierr);
}
