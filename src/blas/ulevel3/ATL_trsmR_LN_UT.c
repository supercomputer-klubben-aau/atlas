/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Rakib Hasan
 * Code contributers : Rakib Hasan, Majedul Sujon
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_reflvl3.h"

#define ATL_trsmR_LNUT Mjoin(PATL,trsmR_LNUT)

int ATL_trsmR_LNUT
(
   ipinfo_t *ip,     /* ipinfo for gemm */
   sminfo_t *si,     /* ipinfo for trsm */
   const enum ATLAS_DIAG Diag,
   ATL_CSZT R,
   ATL_CSZT N,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   TYPE *X,
   ATL_CSZT ldx,
   ATL_CSZT Tsz,
   ATL_CSZT Rsz,
   TYPE *diag,      /* workspace ptr for diagonals */
   TYPE *L,        /* workspace ptr for whole A matrix */
   TYPE *RW,       /* workspace ptr for B panel */
   TYPE *w         /* workspace ptr for C in gemm (may be shared with trsm) */
)
{
   cm2am_t a2blk = NULL;
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
   TYPE *x, *l, *wL, *wR, *wC;
   unsigned int szA, szB, szC;
   unsigned int MB, NB, nnblks, inca;
   void (*utrsm)(sminfo_t*,const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx,
         TYPE *diag, TYPE *L, TYPE *RW, TYPE *w);
   utrsm = si->utrsm;

   nnblks = ip->nfnblks + ip->npnblks;
   szA = ip->szA;
   szB = ip->szB;
   szC = ip->szC;
   a2blk = ip->a2blk;
   MB = ip->mb;
   NB = ip->nb;

   wR = w;
   wC = w + (Rsz SHIFT);
/*
 * NOTE: Gemm will always be called on full K-blks,
 * only TRSM applied on partial k-blks. So, make KB0 and kb0 equal to kb
 */
   ip->KB0 = ip->kb0 = ip->kb;

   for (r=0, ib=0, x=X; r < R; r += MB, ib++, x += (MB SHIFT))
   {
      int nb, mb = R - r;
      const int DoCopy = !r;
      TYPE *a, *d, *Ac = ((TYPE*)A), *xc = x;
      mb = Mmin(mb, MB);
      nb = Mmin(NB, N);
      inca = mb < MB ? ip->pszA : ip->szA;
      Ac += ((N - nb) SHIFT) * (lda+1);
      xc += ((N - nb) SHIFT) * ldx;
/*
 *    do the first triangle
 */
      wL = L;
      a = DoCopy ? Ac : NULL;
      utrsm(si, Diag, mb, nb, alpha, a, lda, xc, ldx, diag, wL, wR, wC);

      for (i=N-nb, jb=1, l=L+(Tsz SHIFT);
            i > 0; i -= NB, jb++, l+=(Tsz SHIFT))
      {
         nb = Mmin(i, NB);
         Ac -= (nb SHIFT) * (lda+1);
         xc -= (nb SHIFT) * ldx;
         #ifdef TCPLX
         {
            TYPE *iR = RW+(nnblks-jb)*(inca SHIFT), *rR = iR + inca;
            a2blk(NB, mb, ONE, xc+(nb SHIFT)*ldx, ldx, rR, iR);
         }
         #else
            a2blk(NB, mb, ONE, xc+nb*ldx, ldx, RW+(nnblks-jb)*inca);
         #endif
         a = DoCopy ? (Ac+nb*ainc) : NULL;
         ip->nfkblks = jb - 1;
         Mjoin(PATL,iploopsK)(ip, ib, jb, NULL, a, xc, MV,
               RW+(nnblks-jb)*(inca SHIFT), l, w, w+szC, alpha, ip->blk2c);
/*
 *       now do the diagonal trsm
 */
         l += (jb SHIFT) * szB;
         d = diag + ((N-i) SHIFT);
         wL = l;
         a = DoCopy ? Ac : NULL;
         utrsm(si, Diag, mb, nb, ONE, a, lda, xc, ldx, d, wL, wR, wC);
      }
   }
   return(0);
}
