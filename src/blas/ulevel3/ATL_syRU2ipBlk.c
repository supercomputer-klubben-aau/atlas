/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#define ATL_GLOBIDX 1
#include "atlas_amm.h"
#undef ATL_GLOBIDX
/*
 * D is kb*mb workspace
 * This guy will copy a block from col-major symmetric storage (A) to
 * access-major (in W).
 * The coordinates (KBLK,JBLK) are logical, and therefore indicate which
 * W block is being filled in.  Logical entry (K,J) is start of block,
 * with K=kb*KBLK, J=nb*JBLK.
 * If (K,J) is in the upper triangle, we need to copy an kbxmb block from it,
 * otherwise we copy a transposed block from the accessed triangle to create
 * the reflected block.  The actual block is kbxnb, so the transposed block
 * will be nbxkb.
 *
 * This version has to double copy blocks containing the diagonal so that it
 * can work if nb and kb are not multiples.  It is possible to avoid this
 * double copy if they are multiples, but since large SYMM rarely used, and
 * this is low order cost, we don't currently bother
 */
void Mjoin(PATL,syRU2ipBlk) /* left, upper hermitian copy to amm */
(
   ipinfo_t *ip,  /* gemm info to copy to */
   ATL_CUINT bv,  /* bitvec: 0:W is full mat, 1:W is Kpanel, else block only */
   cm2am_t cpN,   /* copy to use in Lower triangle */
   cm2am_t cpT,   /* copy to use in upper triangle (transpose of lower) */
   ATL_iptr_t KBLK, ATL_iptr_t JBLK,  /* ip block coordinates to copy */
   const TYPE *B, /* Right, upper Hermitian/symmetric matrix to copy */
   TYPE *W,       /* output array for amm B */
   TYPE *D        /* kb*nb workspace for handling diagonal blocks */
)
{
   unsigned int nb, kb;
   const ATL_iptr_t nfnblks = ip->nfnblks, ldb=ip->ldb;
   ATL_iptr_t J, K=KBLK*ip->kb, szB;
   #ifdef TCPLX
      TYPE *rw, *iw=W;

      if (bv&1)
         iw = IdxBw_ip(ip, W, KBLK, JBLK);
      else if (bv&2)
         iw = IdxBw_ip(ip, W, KBLK, 0);
   #else
      TYPE *w=W;

      if (bv&1)
         w = IdxBw_ip(ip, W, KBLK, JBLK);
      else if (bv&2)
         w = IdxBw_ip(ip, W, KBLK, 0);
   #endif
   kb = (KBLK < ip->nfkblks) ? ip->kb : ip->kb0;
   if (JBLK < nfnblks)
   {
      nb = ip->nb;
      szB = ip->szB;
      J = JBLK*nb;
   }
   else
   {
      nb = (JBLK < nfnblks+ip->npnblks-1) ? ip->pnb : ip->nF;
      szB = ip->pszB;
      J = nfnblks*ip->nb + (JBLK-nfnblks)*ip->pnb;
   }
   #ifdef TCPLX
      rw = iw + szB;
   #endif
   if (J >= K) /* block is in upper part, kb may extend to upper */
   {
      if (J >= K+kb) /* block is all upper */
      #ifdef TCPLX
         cpN(kb, nb, ip->alpB, B+((K+ldb*J)SHIFT), ldb, rw, iw);
      #else
         cpN(kb, nb, ip->alpB, B+K+ldb*J, ldb, w);
      #endif
      else           /* block has diagonal in it! */
      {
         ATL_iptr_t j, k;
         TYPE *d = D;
         #ifdef TCPLX
            ATL_CUINT kb2=kb+kb;
         #else
            #define kb2 kb
         #endif
/*
 *       Copy block to D, getting upper elts through reflection
 */
         for (j=0; j < nb; j++, d += kb2)
         {
            const ATL_iptr_t JJ=J+j;
            for (k=0; k < kb; k++)
            {
               const ATL_iptr_t KK=K+k;
               const ATL_iptr_t IA = ((JJ >= KK) ? KK+ldb*JJ:JJ+ldb*KK)SHIFT;
               #ifdef TCPLX
                  const ATL_iptr_t ii=k+k;
                  d[ii] = B[IA];
                  d[ii+1] = B[IA+1];
               #else
                  d[k] = B[IA];
               #endif
            }
         }
         #ifdef TCPLX
            cpN(kb, nb, ip->alpB, D, kb, rw, iw);
         #else
            cpN(kb, nb, ip->alpB, D, kb, w);
         #endif
      }
   }
/*
 * The block is actually at (J,K), and block to be copied is nbxkb, so danger
 * is that this block crosses the diagonal (J+kb > K).
 * Therefore, D is nbxkb in this case.
 */
   else if (J+nb <= K) /* logical block is wholly contained in upper portion */
   #ifdef TCPLX
      cpT(kb, nb, ip->alpB, B+((J+ldb*K)<<1), ldb, rw, iw);
   #else
      cpT(kb, nb, ip->alpB, B+J+ldb*K, ldb, w);
   #endif
   else        /* block is reflected across diagonal! */
   {
      ATL_iptr_t j, k;
      TYPE *d = D;
      #ifdef TCPLX
         ATL_CUINT nb2=nb+nb;
      #else
         #define nb2 nb
      #endif
/*
 *    Copy block to D, getting upper elts through reflection
 */
      for (j=0; j < kb; j++, d += nb2)
      {
         for (k=0; k < nb; k++)
         {
            const ATL_iptr_t KK=K+j, JJ=J+k;
            const ATL_iptr_t IA = ((JJ >= KK) ? KK+ldb*JJ:JJ+ldb*KK)SHIFT;
            #ifdef TCPLX
               ATL_CUINT ii=k+k;
               d[ii] = B[IA];
               d[ii+1] = B[IA+1];
            #else
               d[k] = B[IA];
            #endif
         }
      }
      #ifdef TCPLX
         cpT(kb, nb, ip->alpB, D, nb, rw, iw);
      #else
         cpT(kb, nb, ip->alpB, D, nb, w);
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
