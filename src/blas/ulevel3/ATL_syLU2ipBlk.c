/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#define ATL_GLOBIDX 1
#include "atlas_amm.h"
#undef ATL_GLOBIDX
/*
 * D is ldd*kb workspace, ldd >= mb
 * This guy will copy a block from col-major symmetric storage (A) to
 * access-major (in W).
 * The coordinates (IBLK,KBLK) are logical, and therefore indicate which
 * W block is being filled in.  Logical entry (I,J) is start of block,
 * with I=mb*IBLK, J=kb*KBLK.
 * If (I,J) is in the accessed triangle, we need to copy an mbxkb block from it,
 * otherwise we copy a transposed block from the accessed triangle to create
 * the reflected block.  The actual block is mbxkb, so the transposed block
 * will be kbxmb.
 *
 * This version has to double copy blocks containing the diagonal so that it
 * can work if mb and kb are not multiples.  It is possible to avoid this
 * double copy if they are multiples, but since large SYMM rarely used, and
 * this is low order cost, we don't currently bother
 */

void Mjoin(PATL,syLU2ipBlk) /* left, upper hermitian copy to amm */
(
   ipinfo_t *ip,  /* gemm info to copy to */
   ATL_CUINT bv,  /* 0: wA full mat, 1: wA k-panel */
   cm2am_t cpN,   /* copy to use in Lower triangle */
   cm2am_t cpT,   /* copy to use in upper triangle (transpose of lower) */
   ATL_iptr_t IBLK, ATL_iptr_t KBLK,  /* ip block coordinates to copy */
   const TYPE *A, /* Left, upper Hermitian matrix to copy */
   TYPE *W,       /* output array for amm A */
   TYPE *D,       /* ldd*kb workspace for handling diagonal blocks */
   ATL_iptr_t ldd /* ldd >= mb */
)
{
   unsigned int mb, kb;
   const ATL_iptr_t nfmblks = ip->nfmblks, lda=ip->lda;
   ATL_iptr_t I, K=KBLK*ip->kb, szA;
   #ifdef TCPLX
      ATL_iptr_t ldd2=ldd+ldd;
      TYPE *rw, *iw=W;

      if (bv&1)
         iw = IdxAw_ip(ip, W, IBLK, KBLK);
      else if (bv&2)
         iw = IdxAw_ip(ip, W, 0, KBLK);
   #else
      #define ldd2 ldd
      TYPE *w=W;

      if (bv&1)
         w = IdxAw_ip(ip, W, IBLK, KBLK);
      else if (bv&1)
         w = IdxAw_ip(ip, W, 0, KBLK);
   #endif
   kb = (KBLK < ip->nfkblks) ? ip->kb : ip->kb0;
   if (IBLK < nfmblks)
   {
      mb = ip->mb;
      szA = ip->szA;
      I = IBLK*ip->mb;
   }
   else
   {
      mb = (IBLK < nfmblks+ip->npmblks-1) ? ip->pmb : ip->mF;
      szA = ip->pszA;
      I = nfmblks*ip->mb + (IBLK-nfmblks)*ip->pmb;
   }
   #ifdef TCPLX
      rw = iw + szA;
   #endif
   if (K >= I) /* block starts in upper part, mb may extend to lower */
   {
      if (K >= I+mb) /* block is all upper */
      #ifdef TCPLX
         cpN(kb, mb, ip->alpA, A+((I+lda*K)SHIFT), lda, rw, iw);
      #else
         cpN(kb, mb, ip->alpA, A+I+lda*K, lda, w);
      #endif
      else           /* block has diagonal in it! */
      {
         ATL_iptr_t i, k;
         TYPE *d = D;
         /* sycp A->D */
/*
 *       Copy block to D, getting upper elts through reflection
 */
         for (k=0; k < kb; k++, d += ldd2)
         {
            const ATL_iptr_t KK=K+k;
            for (i=0; i < mb; i++)
            {
               const ATL_iptr_t II=I+i;
               const ATL_iptr_t IA = ((KK >= II) ? II+lda*KK:KK+lda*II)SHIFT;
               #ifdef TCPLX
                  const ATL_iptr_t ii=i+i;
                  d[ii] = A[IA];
                  d[ii+1] = A[IA+1];
               #else
                  d[i] = A[IA];
               #endif
            }
         }
         #if 0
            Mjoin(PATL,geprint)("Al", mb, kb, A+I+K*lda, lda);
            Mjoin(PATL,geprint)("Dn", mb, kb, D, ldd);
         #endif
         #ifdef TCPLX
            cpN(kb, mb, ip->alpA, D, ldd, rw, iw);
         #else
            cpN(kb, mb, ip->alpA, D, ldd, w);
         #endif
      }
   }
/*
 * This block actually (K,I), and block to be copied is kbxmb, so danger is
 * that this block crosses the diagonal (K+kb > I).
 */
   else if (K+kb <= I) /* logical block is wholly contained in lower portion */
   #ifdef TCPLX
      cpT(kb, mb, ip->alpA, A+((K+lda*I)<<1), lda, rw, iw);
   #else
      cpT(kb, mb, ip->alpA, A+K+lda*I, lda, w);
   #endif
   else        /* block is reflected across diagonal! */
   {
      ATL_iptr_t i, k;
      TYPE *d = D;
      #ifdef TCPLX
         ATL_CUINT kb2=kb+kb;
      #else
         #define kb2 kb
      #endif
/*
 *    Copy block to D, getting lower elts through reflection
 */
      for (k=0; k < mb; k++, d += kb2)
      {
         for (i=0; i < kb; i++)
         {
            const ATL_iptr_t KK=K+i, II=I+k;
            const ATL_iptr_t IA = ((II <= KK) ? II+lda*KK:KK+lda*II)SHIFT;
            #ifdef TCPLX
               ATL_CUINT ii=i+i;
               d[ii] = A[IA];
               d[ii+1] = A[IA+1];
            #else
               d[i] = A[IA];
            #endif
         }
      }
      #if 0
         Mjoin(PATL,geprint)("AU", kb, mb, A+K+I*lda, lda);
         Mjoin(PATL,geprint)("Dt", kb, mb, D, kb);
      #endif
      #ifdef TCPLX
         cpT(kb, mb, ip->alpA, D, kb, rw, iw);
      #else
         cpT(kb, mb, ip->alpA, D, kb, w);
      #endif
   }
   #if 0
   {
      ATL_iptr_t MM=IBLK*ip->mb+mb, KK=KBLK*ip->kb+kb, K4=Mmin(KK,4);
      char nm[64];
      sprintf(nm, "S%lu_%lu", (unsigned long)MM, (unsigned long)KK);
      Mjoin(PATL,ipprint)(stdout, ip, 1, nm, MM, K4, W);
   }
   #endif
}
