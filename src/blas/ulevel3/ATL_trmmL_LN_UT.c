/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Rakib Hasan
 * Code contributers : Rakib Hasan, Majedul Sujon
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_reflvl3.h"

#define trmmBlk Mjoin(PATL,trmmLLNBlk)
#define trmmL_LNUT Mjoin(PATL,trmmL_LNUT)

void trmmBlk
(
   tminfo_t *si,       /* trmm info */
   int flag,            /* 1: c-workspace shared with gemm */
   ATL_CSZT M,
   ATL_CSZT N,
   const SCALAR alpha,
   const TYPE *A,       /* if non-NULL, base A ptr to copy */
   ATL_CSZT lda,
   TYPE *B,             /* if non-NULL, base B ptr to copy */
   ATL_CSZT ldb,
   TYPE *C,             /* if non-NULL, base C ptr to copy the result */
   ATL_CSZT ldc,
   TYPE *wA,            /* workspace for A */
   TYPE *wB,            /* workspace for B */
   TYPE *rC,            /* real portion of wC (unused for real routines) */
   TYPE *wC,            /* workspace for C */
   ATL_UINT szA,        /* workspace size of iA */
   ATL_UINT szB         /* workspace size of iB, may vary based on sharing */
)
{
   const int mu = si->mu;
   const int nu = si->nu;
   const int ku = si->ku;
   ATL_CSZT nmu = (M + mu - 1) / mu;
   ATL_CSZT nnu = (N + nu - 1) / nu;
   ATL_CSZT K = ((M + ku - 1) / ku) * ku;
   tcm2am_t t2a = si->t2blk;
   cm2am_t r2a = si->r2blk;
   ablk2cmat_t blk2c = si->blk2c;
   ammkern_t trmmK_b0 = si->amm_b0;
   ammkern_t trmmK_b1 = si->amm_b1;
   #ifdef TCPLX
      ammkern_t trmmK_bn = si->amm_bn;
      TYPE *rA = wA + szA;
      TYPE *rB = wB + szB;
      TYPE ZERO[2] = {ATL_rzero, ATL_rzero};
   #else
      #define ZERO ATL_rzero
   #endif
/*
 * Appropriate copy routines with alpha has been selected from tminfo
 * so, it's safe to pass alpha through all routines
 */
   #ifdef TCPLX
      if (A) t2a(M, alpha, A, lda, rA, wA);
      if (B) r2a(M, N, alpha, B, ldb, rB, wB);
      if (flag&1) /* shared with gemm, so accumulate the result */
      {
         trmmK_bn(nmu, nnu, K, wA, wB, rC, rA, wB, wC);
         trmmK_b1(nmu, nnu, K, rA, wB, wC, rA, rB, rC);
      }
      else
      {
         trmmK_b0(nmu, nnu, K, wA, wB, rC, rA, wB, wC);
         trmmK_b0(nmu, nnu, K, rA, wB, wC, rA, rB, rC);
      }
      trmmK_bn(nmu, nnu, K, rA, rB, rC, wA, rB, wC);
      trmmK_b1(nmu, nnu, K, wA, rB, wC, wA, wB, rC);
      if (C) blk2c(M, N, alpha, rC, wC, ZERO, C, ldc);
   #else
      if (A) t2a(M, alpha, A, lda, wA);
      if (B) r2a(M, N, alpha, B, ldb, wB);
      if (flag&1) /* shared with gemm, so accumulate the result */
         trmmK_b1(nmu, nnu, K, wA, wB, rC, wA, wB, rC);
      else
         trmmK_b0(nmu, nnu, K, wA, wB, rC, wA, wB, rC);
      if (C) blk2c(M, N, alpha, rC, ZERO, C, ldc);
   #endif
}

int trmmL_LNUT
(
   ipinfo_t *ip,   /* ipinfo for gemm */
   tminfo_t *si,  /* ipinfo for trmm */
   ATL_CSZT M,
   ATL_CSZT N,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   TYPE *X,
   ATL_CSZT ldx,
   ATL_CSZT Tsz,   /* workspace size of single A Block for TRMM */
   ATL_CSZT Rsz,   /* workspace size of single B Block for TRMM */
   ATL_CSZT Csz,   /* workspace size of result for TRMM, 0 if shared with gemm*/
   TYPE *tbw,      /* workspace ptr for TRMM */
   TYPE *L,        /* workspace ptr for whole A matrix */
   TYPE *RW,       /* workspace ptr for B panel */
   TYPE *w         /* workspace ptr for C in gemm (may be shared with trmm) */
)
{
   #ifdef TCPLX
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
      TYPE ZERO[2] = {ATL_rzero, ATL_rzero};
   #else
      #define ONE ATL_rone
      #define ZERO ATL_rzero
   #endif
   ATL_SZT i, j;
   const int ainc = si->incA;
   int tincb;
   unsigned int ib, jb, nfkblks, Mr;
   const int MV = 3; /* move A & B */
   TYPE *x, *l, *wL, *tcw;
   TYPE *rC, *wC, *wT, *rT;
   unsigned int NB, szA, szC, nmblks, incb;
   /*const int MB = ip->kb;*/
   const int MB = ip->mb;
   int flag = si->flg;

   #ifdef DEBUG
   if (ip->mb != ip->kb)
   {
/*
 *    NOTE: MB may not be equal to KB when we have only one block
 */
      fprintf(stderr, "MB is not equal to KB!!!\n");
   }
   #endif

   szA = ip->szA;
   szC = ip->szC;
   #ifdef TCPLX
      wC = w;
      rC = wC + szC;
   #else
      wC = rC = w;
   #endif
   if (flag&1) /* trmm shares c-space with gemm */
   {
      wT = wC;
      rT = rC;
   }
   else /* not shared with gemm */
   {
      tcw = tbw + (Rsz SHIFT);
      #ifdef TCPLX
         wT = tcw;
         rT = tcw + Csz;
      #else
         wT = rT = tcw;
      #endif
   }
   NB = ip->nb;
   nmblks = ip->nfmblks + ip->npmblks;
   nfkblks = (M-MB-1) / MB;
   Mr = M - (M/MB)*MB;
   if (!Mr) Mr = MB;
/*
 * NOTE: Gemm will always be called on full K-blks,
 * only TRMM applied on partial k-blks. So, make KB0 and kb0 equal to kb
 */
   ip->KB0 = ip->kb0 = ip->kb;
/*
 * for each col panel Bcpan in B
 */
   for (j=0, jb=0, x=X; j < N; j += NB, jb++, x += (NB SHIFT)*ldx)
   {
      int mb, nb;
      const int DoCopyA = !j;
      TYPE *a, *b, *Ac = ((TYPE*)A), *xc = x, *xb, *rb;

      nb = Mmin((N-j), NB);
      mb = Mmin(Mr, M);

      incb = nb < NB ? ip->pszB : ip->szB; /* size of B workspace for panel */
      ip->nfkblks = nfkblks; /* reset nfkblks for the new panel */
      Ac += ((M - mb) SHIFT) * (lda+1);  /* A --> point to last diag blk */
      xc += ((M - mb) SHIFT);            /* B --> last blk of the panel */
/*
 *    for each block Bblk of Bcpan from bottom to top
 */
      for (i=M-mb, ib=1, l=L; i > 0; i -= MB, ib++)
      {
         a = DoCopyA ? (Ac - i*ainc) : NULL;
         b = (ib == 1) ? x : NULL;
/*
 *       gemm call: store tmp result to workspace, not copied out.
 *       rC,wC has the result, not copied out to X
 */
         Mjoin(PATL,iploopsK)(ip, nmblks-ib, jb, a, b, xc, MV, l, RW, rC, wC,
               ZERO, NULL);
         l += ((ip->nfkblks+1) SHIFT) * szA; /* point to A workspace for TRMM */
         wL = l;
         a = DoCopyA ? Ac : NULL; /* when need to copy, point to diagonal blk */
/*
 *       share B copies with gemm if both kernels use same mu,nu,ku,mdim and
 *       gemm already did the copy
 *       share workspace for C as well to minimize extra blk2c copy
 */
         if ((flag&1) && (ib!=1))
         {
            xb = NULL;
            rb = RW + ( ((i/MB) SHIFT)*incb);
            tincb = incb;
         }
         else /* can't share copies with gemm */
         {
            xb = xc;
            rb = tbw;
            tincb = Rsz;
         }
         trmmBlk(si, si->flg, mb, nb, alpha, a, lda, xb, ldx, xc, ldx, wL, rb,
                 rT, wT, Tsz, tincb);
/*
 *       if no sharing, copy out w with blk2c_b1
 *       NOTE: TRMM always copies the result to X
 */
         if (!(flag&1)) /* no sharing means need extra blk2c copy */
         {
            #ifdef TCPLX
               ip->blk2c_b1(mb, nb, ip->alpC, rC, wC, ONE, xc, ldx);
            #else
               ip->blk2c_b1(mb, nb, ip->alpC, rC, ONE, xc, ldx);
            #endif
         }
         l += (Tsz SHIFT); /* skip A-copy space for TRMM */
         mb = MB; /* for 1st iter, mb = mr, next on, it's MB */
         Ac -= (mb SHIFT) * (lda+1); /* point to next up diagonal */
         xc -= (mb SHIFT);
         ip->nfkblks--; /* reduce number of k blocks for next gemm */
      }
      wL = l;
      a = DoCopyA ? Ac : NULL;
/*
 *    share B copies with gemm if possible
 */
      if ((flag&1) && M > mb) /* inner loop execuites */
      {
         xb = NULL;
         /*rb = RW + ( ((i/MB) SHIFT)*incb) ;*/
         rb = RW; /* always first block at this point */
         tincb = incb;
      }
      else
      {
         xb = xc;
         rb = tbw;
         tincb = Rsz;
      }
/*
 *    TRMM only block, don't accumulate the result
 */
      trmmBlk(si, 0, mb, nb, alpha, a, lda, xb, ldx, xc, ldx, wL, rb, rT,
            wT, Tsz, tincb);
   }
   return(0);
}
