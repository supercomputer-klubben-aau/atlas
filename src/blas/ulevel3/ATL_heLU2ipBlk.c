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

void Mjoin(PATL,heLU2ipBlk) /* left, upper hermitian copy to amm */
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
   ATL_iptr_t ldd2=ldd+ldd;
   TYPE *rw, *iw=W;

/* printf("B(%lu,%lu)\n", IBLK, KBLK); */
   if (bv&1)
      iw = IdxAw_ip(ip, W, IBLK, KBLK);
   else if (bv&2)
      iw = IdxAw_ip(ip, W, 0, KBLK);
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
   rw = iw + szA;
   if (K >= I) /* block starts in upper part, mb may extend to lower */
   {
      if (K > I+mb) /* block is all lower */
         cpN(kb, mb, ip->alpA, A+((I+lda*K)SHIFT), lda, rw, iw);
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
               if (KK > II)  /* Still upper, so no need to negate imag */
               {
                  const ATL_iptr_t IA=(II+lda*KK)SHIFT;
                  ATL_CUINT ii = i+i;
                  d[ii] = A[IA];
                  d[ii+1] = A[IA+1];
               }
               else if (II > KK) /* Lower, must negate imag */
               {
                  const ATL_iptr_t IA=(KK+lda*II)SHIFT;
                  ATL_CUINT ii = i+i;
                  d[ii] = A[IA];
                  d[ii+1] = -A[IA+1];
               }
               else /* II==JJ, diagonal, must zero imag */
               {
                  const ATL_iptr_t IA=II*((1+lda)SHIFT);
                  ATL_CUINT ii = i+i;
                  d[ii] = A[IA];
                  d[ii+1] = ATL_rzero;
               }
            }
         }
         cpN(kb, mb, ip->alpA, D, ldd, rw, iw);
      }
   }
/*
 * This block actually (K,I), and block to be copied is kbxmb, so danger is
 * that this block touches the diagonal (K+kb >= I).
 */
   else if (K+kb < I) /* logical block is wholly contained in lower portion */
      cpT(kb, mb, ip->alpA, A+((K+lda*I)<<1), lda, rw, iw);
   else        /* block is reflected across diagonal! */
   {
      ATL_iptr_t i, k;
      TYPE *d = D;
      ATL_CUINT kb2=kb+kb;
/*
 *    Copy block to D, getting lower elts through reflection
 */
      for (k=0; k < mb; k++, d += kb2)
      {
         for (i=0; i < kb; i++)
         {
            const ATL_iptr_t KK=K+i, II=I+k;
            if (KK > II)  /* logically in upper, so just copy */
            {
               const ATL_iptr_t IA = (II+lda*KK)SHIFT;
               ATL_CUINT ii=i+i;
               d[ii] = A[IA];
               d[ii+1] = -A[IA+1]; /* because cpT will negate! */
            }
            else if (KK < II)   /* logically in lower, so conjugate */
            {
               const ATL_iptr_t IA = (KK+lda*II)SHIFT;
               ATL_CUINT ii=i+i;
               d[ii] = A[IA];
               d[ii+1] = A[IA+1]; /* cpT will negate as needed for U */
            }
            else /* diagonal, zero imag */
            {
               const ATL_iptr_t IA = ((lda+1)SHIFT)*II;
               ATL_CUINT ii=i+i;
               d[ii] = A[IA];
               d[ii+1] = ATL_rzero;
            }
         }
      }
      cpT(kb, mb, ip->alpA, D, kb, rw, iw);
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
