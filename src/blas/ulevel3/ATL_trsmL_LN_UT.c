/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Rakib Hasan
 * Code contributers : Rakib Hasan, Majedul Sujon
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_reflvl3.h"
#define ATL_trsmL_LNUT Mjoin(PATL,trsmL_LNUT)

int ATL_trsmL_LNUT
(
   ipinfo_t *ip,     /* ipinfo for gemm */
   sminfo_t *si,     /* ipinfo for trsm */
   const enum ATLAS_DIAG Diag,
   ATL_CSZT N,       /* M */
   ATL_CSZT R,       /* N */
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   TYPE *X,
   ATL_CSZT ldx,
   ATL_CSZT Tsz,     /* workspace size of single A Block for TRSM */
   ATL_CSZT Rsz,     /* size of R in TRSM */
   TYPE *diag,       /* workspace ptr for diagonals */
   TYPE *L,          /* workspace ptr for A matrix used both in gemm and trsm */
   TYPE *RW,         /* workspace ptr for B panel */
   TYPE *w           /* workspace ptr for C in gemm, used in trsm as well */
)
{
   cm2am_t b2blk = NULL;
   #ifdef TCPLX
      TYPE ONE[3] = {ATL_rone, ATL_rzero, ATL_rzero}, *ZERO=ONE+1;
      TYPE NONE[2] = {ATL_rnone, ATL_rzero};
   #else
      #define ONE ATL_rone
      #define ZERO ATL_rzero
      #define NONE ATL_rnone
   #endif
   const int ainc = si->incA;
   ATL_SZT r;
   int i, ib, jb;
   const int MV = 3; /* move A & B */
   TYPE *x, *l, *rw, *wL, *wR, *wC;
   unsigned int NB, MB, szA, szC, incb;
   void (*utrsm)(sminfo_t*,const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx,
         TYPE *diag, TYPE *L, TYPE *RW, TYPE *w);

   utrsm = si->utrsm;

   szA = ip->szA;
   szC = ip->szC;
   b2blk = ip->b2blk;
   NB = ip->nb;
   MB = ip->mb;

   wR = w;
   wC = w + (Rsz SHIFT);
/*
 * NOTE: Gemm will always be called on full K-blks,
 * only TRSM applied on partial k-blks. So, make KB0 and kb0 equal to kb
 */
   ip->KB0 = ip->kb0 = ip->kb;

   for (r=0, jb=0, x=X; r < R; r += NB, jb++, x += (NB SHIFT)*ldx)
   {
      int mb, nb = R - r;
      const int DoCopy = !r;
      TYPE *a, *d, *Ac = ((TYPE*)A), *xc = x;
      nb = Mmin(nb, NB);
      mb = Mmin(MB, N);
      incb = nb < NB ? ip->pszB : ip->szB;
/*
 *    do the first triangle
 */
      wL = L;
      a = DoCopy ? Ac : NULL;
      utrsm(si, Diag, mb, nb, alpha, a, lda, xc, ldx, diag, wL, wR, wC);

      for (i=mb, ib=1, l=L+(Tsz SHIFT);
            i < N; i += MB, ib++, l+=(Tsz SHIFT))
      {
         mb = Mmin(N-i, MB);
         Ac += (MB SHIFT) * (lda+1);
         xc += (MB SHIFT);
         #ifdef TCPLX
         {
            TYPE *iR = RW+(((ib-1)*incb) SHIFT), *rR = iR + incb;
            b2blk(MB, nb, ONE, xc-(MB SHIFT), ldx, rR, iR);
         }
         #else
            b2blk(MB, nb, ONE, xc-MB, ldx, RW+(ib-1)*incb);
         #endif
         a = DoCopy ? (Ac-i*ainc) : NULL;
         ip->nfkblks = ib - 1;
         Mjoin(PATL,iploopsK)(ip, ib, jb, a, NULL, xc, MV, l, RW, w, w+szC,
               alpha, ip->blk2c);
/*
 *       now do the diagonal trsm
 */
         l += (ib SHIFT) * szA;
         d = diag + (i SHIFT);
         wL = l;
         a = DoCopy ? Ac : NULL;
         utrsm(si, Diag, mb, nb, ONE, a, lda, xc, ldx, d, wL, wR, wC);
      }
   }
   return(0);
}
