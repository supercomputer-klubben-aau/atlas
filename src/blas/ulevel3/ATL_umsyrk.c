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
#include Mstr(Mjoin(ATLAS_PRE,amm_umsyrk.h))
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
   #define syrkBlk Mjoin(PATL,umherkBlk)
   #define syrk_K  Mjoin(PATL,umherk_KLoop)
   #define syrk  Mjoin(PATL,umherk)
#else
   #define syrkBlk Mjoin(PATL,umsyrkBlk)
   #define syrk_K  Mjoin(PATL,umsyrk_KLoop)
   #define syrk  Mjoin(PATL,umsyrk)
#endif
/*#define MEM_DEBUG 1*/

void syrkBlk
(
   ipinfo_t *ip,
   int flag,      /* bitvec: 0: set means C is upper, 1: set TA==AtlasNoTrans*/
   size_t d,      /* which global diagonal blk of C is being computed */
   size_t k,      /* which K block of A is being computed */
   const TYPE *A, /* if non-NULL, base A ptr to copy */
   const TYPE *B, /* if non-NULL, base A ptr to copy */
   ablk2cmat_t blk2c, /* if non-NULL sy storage to C storage copy func */
   const SCALAR beta, /* only needed if blk2c non-NULL */
   TYPE *C,       /* if blk2c && non-NULL, which C to write to */
   TYPE *wA,      /* if non-NULL, ip-based A workspace */
   TYPE *wAn,     /* next A wrkspc to be prefetched */
   TYPE *wB,      /* if non-NULL ip-based At workspace */
   TYPE *wBn,     /* next B wrkspc to be prefetched */
   TYPE *rC,      /* real portion of wC (unused for real routines) */
   TYPE *wC      /* if non-NULL: ptr to syrk-storage C wrkspc */
)
{
   ATL_CUINT kb = (k != ip->nfkblks) ? ip->kb : ip->kb0;
   ATL_UINT mb, nb, kbS, nnu, nmu;
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
/*
 * calc nmu and nnu. one must be multiple of other
 */
   if (ip->mu > ip->nu)
   {
      nmu = (nb+ip->mu-1)/ip->mu;
      nnu = nmu * (ip->mu/ip->nu); /* NU must be multiple of MU */
   }
   else
   {
      nnu = (nb+ip->nu-1)/ip->nu;
      nmu = nnu * (ip->nu/ip->mu); /* NU must be multiple of MU */
   }
   kbS = ((kb+ip->ku-1)/ip->ku)*ip->ku; /* same as KB0 */
/*
 * copy A
 */
   if (A)  /* want to copy input array */
   {
      A = IdxA_ip(ip, A, d, k);
      #ifdef TCPLX
         ip->a2blk(kb, nb, ip->alpA, A, ip->lda, rA, wA);
      #else
         ip->a2blk(kb, nb, ip->alpA, A, ip->lda, wA);
      #endif
   }
   if (B)
   {
      B = IdxB_ip(ip, B, k, d);
      #ifdef TCPLX
         ip->b2blk(kb, nb, ip->alpB, B, ip->ldb, rB, wB);
      #else
         ip->b2blk(kb, nb, ip->alpB, B, ip->ldb, wB);
      #endif
   }
   if (wC)  /* want to compute SYRK on this block into wC */
   {
      #ifdef TCPLX
         if (flag&4)
         {
            Mjoin(PATL,amsyrkK_b0)(nmu, nnu, kbS, wA, wB, rC, rA, wB, wC);
            Mjoin(PATL,amsyrkK_b0)(nmu, nnu, kbS, rA, wB, wC, rA, rB, rC);
         }
         else
         {
            Mjoin(PATL,amsyrkK_bn)(nmu, nnu, kbS, wA, wB, rC, rA, wB, wC);
            Mjoin(PATL,amsyrkK_b1)(nmu, nnu, kbS, rA, wB, wC, rA, rB, rC);
         }
         Mjoin(PATL,amsyrkK_bn)(nmu, nnu, kbS, rA, rB, rC, wA, rB, wC);
         Mjoin(PATL,amsyrkK_b1)(nmu, nnu, kbS, wA, rB, wC, wA, wB, wC);
      #else
         if (flag&4)
            Mjoin(PATL,amsyrkK_b0)(nmu, nnu, kbS, wA, wB, wC, wA, wB, wC);
         else
            Mjoin(PATL,amsyrkK_b1)(nmu, nnu, kbS, wA, wB, wC, wA, wB, wC);
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
}

void syrk_K  /* inner-product based syrk/herk loop over K loop */
(
   ipinfo_t *ip,
   int flag,      /* bitvec: 0: set means C is upper, 1: set TA==AtlasNoTrans*/
   size_t d,      /* which global diagonal blk of C is being computed */
   const TYPE *A, /* if non-NULL, base A ptr to copy */
   const TYPE *B, /* if non-NULL, base B ptr to copy */
   ablk2cmat_t blk2c, /* if non-NULL sy storage to C storage copy func */
   const SCALAR beta, /* only needed if blk2c non-NULL */
   TYPE *C,       /* if blk2c non-NULL, which C to write to, else ignored */
   TYPE *wA,      /* if non-NULL, ip-based A workspace */
   TYPE *wB,      /* if non-NULL ip-based At workspace */
   TYPE *rC,      /* real portion of wC (unused in real routines) */
   TYPE *wC       /* if non-NULL: ptr to syrk-storage C wrkspc */
)
{
   const size_t nfkblks = ip->nfkblks;
   size_t k;
   ATL_UINT szA, szB;
   if (d < ip->nfnblks)
   {
      /*szA = szB = ip->szA;*/
      szA = ip->szA;
      szB = ip->szB;
   }
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
      syrkBlk(ip, flag|4, d, 0, A, B, NULL, beta, NULL, wA,wAn, wB,wBn,
                 rC, wC);
      wA = wAn;
      wB = wBn;
   }
   else /* this is first & last block! */
   {
      syrkBlk(ip, flag|4, d, 0, A, B, blk2c, beta, C, wA, wA, wB, wB,
                rC, wC);
      return;
   }
/*
 * Handle all blocks except first (handled above) & last (handled below)
 */
   for (k=1; k < nfkblks; k++)
   {
      TYPE *wAn=(wA)? wA+szA:NULL, *wBn=(wB) ? wB+szB : NULL;
      syrkBlk(ip, flag, d, k, A, B, NULL, beta, NULL, wA, wAn, wB, wBn,
                rC, wC);
      wA = wAn;
      wB = wBn;
   }
/*
 * Last block actually writes to C
 */
   syrkBlk(ip, flag, d, k, A, B, blk2c, beta, C, wA,wA, wB,wB, rC,wC);
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
   int k, ierr;
   ablk2cmat_t blk2sy, blk2c;
   ipinfo_t ip;
   #ifdef Conj_
      TYPE alpha[2]={ralpha, ATL_rzero}, beta[2]={rbeta, ATL_rzero};
      const enum ATLAS_TRANS TB=(TA==AtlasNoTrans)?AtlasConjTrans:AtlasNoTrans;
   #else
      const enum ATLAS_TRANS TB = (TA==AtlasNoTrans) ? AtlasTrans:AtlasNoTrans;
   #endif
#ifdef MEM_DEBUG
   void *vA=NULL, *vB=NULL, *vC=NULL;
#endif
   Mjoin(PATL,ipmenUMInfo)(&ip, TA, TB, N,N,K, lda, lda, ldc, alpha, beta);
   nb = (ip.nfnblks) ? ip.nb : ip.pnb;
   nnblks = ip.nfnblks + ip.npnblks;
/*
 * Main idea:
 * 1. with syrk+gemm of full blks, syrk can always reuse the workspace of gemm
 *    no need to allocate extra space for syrk
 * 2. syrk only: need to allocate only considering the syrk
 */
/*
 * Need space for only one column-panel of At
 */
   szB = ip.nfnblks ? ip.szA : ip.pszA;
   szB *= (ip.nfkblks+1);
/*
 * space calc for single syrk block
 */
   if (nnblks == 1)    /* single blk in n-dimension: syrk only*/
   {
      int smnu=(ip.mu > ip.nu)?ip.mu:ip.nu; /* max of mu, nu*/
      nbS = (nb+smnu-1)/smnu;
      szC = ((nbS+1)*nbS)>>1;
      nbS *= smnu;
/*
 *    considering mu and nu as one is multiple of other, the size of subblock
 *    is smnu*smnu (rounded with ku).
 */
      szC *= ((smnu*smnu+ip.ku-1)/ip.ku)*ip.ku;
      szS =  ((ip.kb+ip.ku-1)/ip.ku)*ip.ku;
      szA = nbS * szS * (ip.nfkblks + 1);
      extra = (smnu+smnu)*smnu;  /* may not need */
   }
   else
   {
/*
 *    A needs entire matrix minus one row/col panel
 */
      szA = szB * (ip.nfnblks + ip.npnblks - 1);
      szC = ip.szC;
      extra = ip.exsz;
   }
#ifdef MEM_DEBUG
/*
 * DEBUG: to debug, I allocate memory individually.... so that I can use
 * valgrind to find out of memory access error
 */
   sz = ATL_MulBySize(szA) + 2*ATL_Cachelen;
   if (sz < ATL_MaxMalloc)
      vA = malloc(sz);
   if (!vA)
      return(1);  /* keep recurring, can't malloc space! */
   wA = ATL_AlignPtr(vA);

   sz = ATL_MulBySize(szB) + 2*ATL_Cachelen;
   if (sz < ATL_MaxMalloc)
      vB = malloc(sz);
   if (!vB)
      return(1);  /* keep recurring, can't malloc space! */
   wB = ATL_AlignPtr(vB);

   sz = ATL_MulBySize(szC + extra) + 2*ATL_Cachelen;
   if (sz < ATL_MaxMalloc)
      vC = malloc(sz);
   if (!vC)
      return(1);  /* keep recurring, can't malloc space! */
   wC = ATL_AlignPtr(vC);
#else
/*
 * use it for timing to reduce mulitple malloc.
 */
   sz = ATL_MulBySize(szC + szA+szB + extra) + 5*ATL_Cachelen;
   if (sz < ATL_MaxMalloc)
      vp = malloc(sz);
   if (!vp)
      return(1);  /* keep recurring, can't malloc space! */
   wA = ATL_AlignPtr(vp);
   wB = wA + (szA SHIFT);
   wB = ATL_AlignPtr(wB);
   wC = wB + (szB SHIFT);
   wC = ATL_AlignPtr(wC);
#endif
   #ifdef TCPLX
      rC = wC + szC;
   #else
      rC = wC;
   #endif
/*
 * ============================================================================
 * First, we compute diagonals of C, and in the process we will copy A/A^T
 * for use in computing non-diaginal blocks using inner-product amm.
 * ============================================================================
 */
   flg = (TA == AtlasNoTrans) ? 2 : 0;
/*
 * NOTE: alpha may be already applied on A-copy. So, use ip.alpC instead of
 * alpha to select appropriate copy routines
 */
   if (Uplo == AtlasLower)
   {
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(ip.alpC))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNb1);
         else
            blk2sy = SCALAR_IS_ONE(ip.alpC) ?
                     Mjoin(PATL,SyrkIntoC_a1b1):Mjoin(PATL,SyrkIntoC_aXb1);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(ip.alpC))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNbN);
         else
            blk2sy = SCALAR_IS_ONE(ip.alpC) ?
                     Mjoin(PATL,SyrkIntoC_a1bN):Mjoin(PATL,SyrkIntoC_aXbN);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(ip.alpC))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNb0);
         else
            blk2sy = SCALAR_IS_ONE(ip.alpC) ?
                     Mjoin(PATL,SyrkIntoC_a1b0):Mjoin(PATL,SyrkIntoC_aXb0);
      }
      else
      {
         if (SCALAR_IS_NONE(ip.alpC))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNbX);
         else
            blk2sy = SCALAR_IS_ONE(ip.alpC) ?
                     Mjoin(PATL,SyrkIntoC_a1bX):Mjoin(PATL,SyrkIntoC_aXbX);
      }
   }
   else
   {
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(ip.alpC))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNb1_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(ip.alpC) ?
                     Mjoin(PATL,SyrkIntoC_a1b1_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXb1_L2UT);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(ip.alpC))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNbN_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(ip.alpC) ?
                     Mjoin(PATL,SyrkIntoC_a1bN_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXbN_L2UT);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(ip.alpC))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNb0_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(ip.alpC) ?
                     Mjoin(PATL,SyrkIntoC_a1b0_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXb0_L2UT);
      }
      else
      {
         if (SCALAR_IS_NONE(ip.alpC))
            blk2sy = Mjoin(PATL,SyrkIntoC_aNbX_L2UT);
         else
            blk2sy = SCALAR_IS_ONE(ip.alpC) ?
                     Mjoin(PATL,SyrkIntoC_a1bX_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXbX_L2UT);
      }
   }
/*
 * syrk-upper
 */
   if (Uplo == AtlasUpper)
   {
      flg |= 1;
/*
 *    1st syrk: copies both formats, no reuse
 */
      syrk_K(&ip, flg, 0, A, A, blk2sy, beta, C, wB, wA, rC, wC);
      if (N > nb)
      {
         const size_t incC = nb*((ldc+1)SHIFT);
         size_t i;
/*
 *       Compute all C blks within this rowpan of C, copies B (At), reuse A
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
/*
 *          syrk copies A and reuses B copy
 */
            syrk_K(&ip, flg, i, A, NULL, blk2sy, beta, C, wB, wA,
                          rC, wC);
            wA += szB SHIFT;
/*
 *          gemm can reuse both the copies
 */
            Mjoin(PATL,iploopsNK)(&ip, i, i+1, NULL, NULL, C+i*(nb SHIFT), 11,
                                  wB, wA, rC, wC, beta, blk2c);
         }
/*
 *       Last colpan is only 1 diag blk
 */
         syrk_K(&ip, flg, i, A, NULL, blk2sy, beta, C, wB, wA,
                       rC, wC);
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
 *    Compute first diag block, copies both
 */
      syrk_K(&ip, flg, 0, A, A, blk2sy, beta, C, wA, wB, rC, wC);
      if (N > nb)
      {
/*
 *       Compute all C blks within this colpan of C, copying rest of A, reuse B
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
/*
 *          syrk reuse A and copies B
 */
            syrk_K(&ip, flg, j, NULL, A, blk2sy, beta, C, wA, wB,
                          rC, wC);
            wA += szB SHIFT;
/*
 *          gemm reuses both the copies
 */
            Mjoin(PATL,iploopsMK)(&ip, j+1, j, NULL, NULL, c, 7, wA, wB, rC, wC,
                                  beta, blk2c);
            c += incC;
         }
/*
 *       Last colpan is only 1 diag blk, still reuses A
 */
         syrk_K(&ip, flg, j, NULL, A, blk2sy, beta, C, wA, wB,
                       rC, wC);
      }
   }
#ifdef MEM_DEBUG
   free(vA);
   free(vB);
   free(vC);
#else
   free(vp);
#endif
   return(0);
}
