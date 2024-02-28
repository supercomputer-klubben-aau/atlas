/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Rakib Hasan
 * Code contributers : Rakib Hasan, Majedul Sujon
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_reflvl3.h"

#define trmmBlk Mjoin(PATL,trmmRLNBlk)
#define trmmR_LNUT Mjoin(PATL,trmmR_LNUT)

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
   ATL_CSZT szA,        /* workspace size of iA */
   ATL_CSZT szB         /* workspace size of iB, may vary based on sharing */
)
{
   const int mu = si->mu;
   const int nu = si->nu;
   const int ku = si->ku;
   ATL_CSZT nmu = (M + mu - 1) / mu;
   ATL_CSZT nnu = (N + nu - 1) / nu;
   ATL_CSZT K = ((N + ku - 1) / ku) * ku;
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
 * Appropriate copies with alpha are selected from tminfo.. safe to pass alpha
 */
   #ifdef TCPLX
      if (A) t2a(N, alpha, A, lda, rA, wA);
      if (B) r2a(N, M, alpha, B, ldc, rB, wB);
      if (flag&1)
      {
         trmmK_bn(nmu, nnu, K, wB, wA, rC, rB, wA, wC);
         trmmK_b1(nmu, nnu, K, rB, wA, wC, rB, rA, rC);
      }
      else
      {
         trmmK_b0(nmu, nnu, K, wB, wA, rC, rB, wA, wC);
         trmmK_b0(nmu, nnu, K, rB, wA, wC, rB, rA, rC);
      }
      trmmK_bn(nmu, nnu, K, rB, rA, rC, wB, rA, wC);
      trmmK_b1(nmu, nnu, K, wB, rA, wC, wB, wA, rC);
      if (C) blk2c(M, N, alpha, rC, wC, ZERO, C, ldc);
   #else
      if (A) t2a(N, alpha, A, lda, wA);
      if (B) r2a(N, M, alpha, B, ldb, wB);
      if (flag&1)
         trmmK_b1(nmu, nnu, K, wB, wA, rC, wB , wA, rC);
      else
         trmmK_b0(nmu, nnu, K, wB, wA, rC, wB , wA, rC);
      if (C) blk2c(M, N, alpha, rC, ZERO, C, ldc);
   #endif
}

int trmmR_LNUT
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
   unsigned int ib, jb, nfkblks, Nr;
   const int MV = 3; /* move A & B */
   TYPE *x, *l, *wL, *tcw;
   TYPE *rC, *wC, *rT, *wT;
   unsigned int MB, szB, szC, nnblks, inca;
   const int NB = ip->nb; /* nb & kb should be same */
   int tincb;
   int flag = si->flg;

   #ifdef DEBUG
   if(ip->nb != ip->kb)
   {
/*
 *       NOTE: NB may not be equal to KB when we have only one block
 */
         fprintf(stderr, "NB is not equal to KB!!!\n");
   }
   #endif
   szB = ip->szB;
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
   MB = ip->mb;
   nnblks = ip->nfnblks + ip->npnblks;
   nfkblks = (N-NB-1) / NB;
   Nr = N - (N/NB)*NB;
   if (!Nr) Nr = NB;
/*
 * NOTE: Gemm will always be called on full K-blks,
 * only TRMM applied on partial k-blks. So, make KB0 and kb0 equal to kb
 */
   ip->KB0 = ip->kb0 = ip->kb;

   for (i=0, ib=0, x=X; i < M; i += MB, ib++, x += (MB SHIFT))
   {
      int nb, mb;
      const int DoCopyA = !i;
      TYPE *a, *b, *d, *Ac = ((TYPE*)A), *xc = x, *rb, *xb;
      mb = Mmin((M-i), MB);
      nb = Mmin(Nr, N);
      inca = mb < MB ? ip->pszA : ip->szA; /* szA actually is for B*/
      ip->nfkblks = nfkblks; /* reset nfkblks for the new panel */

      for (j=nb, jb=1, l=L; j < N; j += NB, jb++)
      {
         a = DoCopyA ? (Ac+ nb*ainc) : NULL;
         b = (jb == 1) ? (xc + (nb SHIFT)*ldx) : NULL;
/*
 *       gemm call: store tmp result to workspace, not copied out.
 *       rC,wC has the result, not copied out to X
 *       NOTE: since remainder is kept at the beginnig of N dimension,
 *       it represents last blk in rowpan with respect to gemm (iploops)
 */
         Mjoin(PATL,iploopsK)(ip, ib, nnblks-jb, b, a, xc, MV,
               RW+(jb-1)*(inca SHIFT), l, rC, wC, ZERO, NULL);
         l += ((ip->nfkblks+1) SHIFT) * szB;
         wL = l;
         a = DoCopyA ? Ac : NULL;
/*
 *       share B copies with gemm if both kernels use same mu,nu,ku,mdim and
 *       gemm already did the copy
 *       share workspace for C as well to minimize extra blk2c copy
 */
         if ((flag&1) && (jb!=1))
         {
            xb = NULL;
            rb = RW + (jb-2)* (inca SHIFT) ;
            tincb = inca;
         }
         else
         {
            xb = xc;
            rb = tbw;
            tincb = Rsz;
         }
         trmmBlk(si, flag, mb, nb, alpha, a, lda, xb, ldx, xc, ldx, wL, rb, rT,
               wT, Tsz, tincb);
/*
 *       if no sharing, copy out w with blk2c_b1
 *       NOTE: TRMM always copies the result to X
 */
         if (!(flag&1))
         {
            #ifdef TCPLX
               ip->blk2c_b1(mb, nb, ip->alpC, rC, wC, ONE, xc, ldx);
            #else
               ip->blk2c_b1(mb, nb, ip->alpC, rC, ONE, xc, ldx);
            #endif
         }
         l += (Tsz SHIFT);
         Ac += (nb SHIFT) * (lda+1);
         xc += (nb SHIFT) * ldx;
         nb = NB;
         ip->nfkblks--; /* reduce number of k blocks for next gemm */
      }
      wL = l;
      a = DoCopyA ? Ac : NULL;
/*
 *    share B copies with gemm if possible
 */
      if ((flag&1) && (N>nb))
      {
            xb = NULL;
            rb = RW + (jb-2)* (inca SHIFT) ;
            tincb = inca; /* inca actually for B*/
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
      trmmBlk(si, 0, mb, nb, alpha, a, lda, xb, ldx, xc, ldx, wL, rb, rT, wT,
            Tsz, tincb);
   }
   return(0);
}
