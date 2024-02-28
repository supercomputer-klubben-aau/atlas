#include "atlas_misc.h"
#define ATL_GLOBIDX 1
#include "atlas_amm.h"
#undef ATL_GLOBIDX
#include Mstr(Mjoin(Mjoin(ATLAS_PRE,amm),_sum.h))
/*
 * This function does one step of matmul while optionally copying data
 * from A or B.  It assumes (rC,iC) already initialized, unless k==0,
 * in which case it is zeroed by using the _b0 kernel.
 * A,B,wA,wB should all be base pointers, and any indexing will be done here.
 *
 * bv is a bitvec with following meanings for each bit position
 * 0 : wA is 1 kpanel in size
 * 1 : wA has space for full matrix
 * 2 : wB is 1 kpanel in size
 * 3 : wB has space for full matrix
 * 4 : use beta=0, else use beta=1
 */
void Mjoin(PATL,doAmmBlk)
(
   ipinfo_t *ip,
   ATL_CUINT bv,          /* bitvec, see explanation above */
   ATL_iptr_t i,          /* which global row block of C is being computed */
   ATL_iptr_t j,          /* which global col block of C is being computed */
   ATL_iptr_t k,          /* which global K block is being multiplied */
   const TYPE *A,         /* if non-NULL, A to cpy from */
   const TYPE *B,         /* if non-NULL, B to cpy from */
   TYPE *C,               /* if non-NULL, C to write to after doing multiply */
   TYPE *wA,              /* wrkspace for block-copied A */
   TYPE *wB,              /* wrkspace for block-copied B */
   TYPE *rC,              /* ptr to real portion of C access-major workspace */
   TYPE *iC,              /* ptr to complex C workspace */
   const SCALAR beta,     /* scaling factor for C */
   const ablk2cmat_t blk2c
)
/*
 * A & B ptrs should point to start of K iteration in this rout
 */
{
   const ATL_iptr_t nfmblks=ip->nfmblks,nfnblks=ip->nfnblks,nfkblks=ip->nfkblks;
   const ATL_iptr_t nmblks=ip->npmblks+nfmblks, nnblks=ip->npnblks+nfnblks;
   ATL_iptr_t inca, incb;
   #ifdef TCPLX
      TYPE *rA, *rB;
   #endif
   ATL_UINT nmu, nnu, mb, nb;
   ATL_CUINT kb0=ip->kb0, kb=(k < nfkblks) ? ip->kb : kb0;

   if (i == nmblks-1)
   {
      nmu = ip->nmuF;
      mb = ip->mF;
      inca = (i < nfmblks) ? ip->szA : ip->pszA;
   }
   else if (i < nfmblks)
   {
      nmu = ip->nmu;
      mb = ip->mb;
      inca = ip->szA;
   }
   else
   {
      nmu = ip->pnmu;
      mb = ip->pmb;
      inca = ip->pszA;
   }

   if (j == nnblks-1)
   {
      nnu = ip->nnuF;
      nb = ip->nF;
      incb = (j < nfnblks) ? ip->szB : ip->pszB;
   }
   else if (j < nfnblks)
   {
      nnu = ip->nnu;
      nb = ip->nb;
      incb = ip->szB;
   }
   else
   {
      nnu = ip->pnnu;
      nb = ip->pnb;
      incb = ip->pszB;
   }
   if (bv&8)
      wB = IdxBw_ip(ip, wB, k, j);
   else if (bv&4)
      wB = IdxBw_ip(ip, wB, k, 0);
   #ifdef TCPLX
      rB = wB + incb;
   #endif
   if (B) /* need to copy from original B */
   {
      B = IdxB_ip(ip, B, k, j);
      #ifdef TCPLX
         ip->b2blk(kb, nb, ip->alpB, B, ip->ldb, rB, wB);
      #else
         ip->b2blk(kb, nb, ip->alpB, B, ip->ldb, wB);
      #endif
   }
   if (bv&2)
      wA = IdxAw_ip(ip, wA, i, k);
   else if (bv&1)
      wA = IdxAw_ip(ip, wA, 0, k);
   #ifdef TCPLX
      rA = wA + inca;
   #endif
   if (A) /* need to copy from original A */
   {
      A = IdxA_ip(ip, A, i, k);
      #ifdef TCPLX
         ip->a2blk(kb, mb, ip->alpA, A, ip->lda, rB, wB);
      #else
         ip->a2blk(kb, mb, ip->alpA, A, ip->lda, wB);
      #endif
   }
   if (!rC)   /* If rC NULL, user wants only A/B copy done */
     return;
   if ((k < nfkblks) | (kb0 == ip->kb))  /* usual case of full kb */
   {
      if (bv&16)    /* rC,iC not initialized, use beta=0 */
      {
         #ifdef TCPLX
            ammkern_t amm = ip->amm_b0;
            amm(nmu, nnu, kb, wA, wB, rC, rA, wB, iC);
            amm(nmu, nnu, kb, rA, wB, iC, rA, rB, rC);
            ip->amm_bn(nmu, nnu, kb, rA, rB, rC, wA, rB, iC);
            ip->amm_b1(nmu, nnu, kb, wA, rB, iC, rA+inca, rB+incb, rC);
         #else
            ip->amm_b0(nmu, nnu, kb, wA, wB, rC, wA+inca, wA+inca, rC);
         #endif
      }
      else  /* rC,iC initialized, so use beta=1 */
      {
         #ifdef TCPLX
            ammkern_t amm_b1 = ip->amm_b1, amm_bn = ip->amm_bn;
            amm_bn(nmu, nnu, kb, wA, wB, rC, rA, wB, iC);
            amm_b1(nmu, nnu, kb, rA, wB, iC, rA, rB, rC);
            amm_bn(nmu, nnu, kb, rA, rB, rC, wA, rB, iC);
            amm_b1(nmu, nnu, kb, wA, rB, iC, rA+inca, rB+incb, rC);
         #else
            ip->amm_b1(nmu, nnu, kb, wA, wB, rC, wA+inca, wB+incb, rC);
         #endif
      }
   }
   else /* final kb0 comes from end of K panel */
   {
      ATL_CUINT KB=ip->KB0;
      if (bv&16)    /* rC,iC not initialized, use beta=0 */
      {
         #ifdef TCPLX
            ammkern_t amm = ip->ammK1_b0;
            amm(nmu, nnu, KB, wA, wB, rC, rA, wB, iC);
            amm(nmu, nnu, KB, rA, wB, iC, rA, rB, rC);
            ip->ammK1_bn(nmu, nnu, KB, rA, rB, rC, wA, rB, iC);
            ip->ammK1_b1(nmu, nnu, KB, wA, rB, iC, wA, rB, rC);
         #else
            ip->ammK1_b0(nmu, nnu, KB, wA, wB, rC, wA, wB, rC);
         #endif
      }
      else  /* rC,iC initialized, so use beta=1 */
      {
         #ifdef TCPLX
            ammkern_t amm_b1 = ip->ammK1_b1, amm_bn = ip->ammK1_bn;
            amm_bn(nmu, nnu, KB, wA, wB, rC, rA, wB, iC);
            amm_b1(nmu, nnu, KB, rA, wB, iC, rA, rB, rC);
            amm_bn(nmu, nnu, KB, rA, rB, rC, wA, rB, iC);
            amm_b1(nmu, nnu, KB, wA, rB, iC, wA, rB, rC);
         #else
            ip->ammK1_b1(nmu, nnu, KB, wA, wB, rC, wA, wB, rC);
         #endif
      }
   }
   if (blk2c)
   #ifdef TCPLX
      blk2c(mb, nb, ip->alpC, rC, iC, beta, C, ip->ldc);
   #else
      blk2c(mb, nb, ip->alpC, rC, beta, C, ip->ldc);
   #endif
}

void Mjoin(PATL,iploopsK)
(
   ipinfo_t *ip,
   size_t i,              /* which global row block of C is being computed */
   size_t j,              /* which global col block of C is being computed */
   const TYPE *A,         /* if non-NULL, Portion A to cpy */
   const TYPE *B,         /* if non-NULL, Portion B to cpy, else ignored */
   TYPE *C,               /* original C to write ans to at end of Kloop */
   int MV,                /* 1: move a, 2: move b, 3: move A & B */
   TYPE *a,               /* wrkspace for block-copied A */
   TYPE *b,               /* wrkspace for block-copied B */
   TYPE *rC,              /* ptr to real portion of C access-major workspace */
   TYPE *iC,              /* ptr to complex C workspace */
   const SCALAR beta,     /* scaling factor for C */
   const ablk2cmat_t blk2c   /* C = alpC*(rC,iC) + beta*C */
)
{
   cm2am_t a2blk=ip->a2blk, b2blk=ip->b2blk;
   ammkern_t amm=ip->ammK1_b0;
   ATL_CUINT kb=ip->kb, kb0=ip->kb0, KB0=ip->KB0;
   ATL_CUINT mu=ip->mu, nu=ip->nu;
   ATL_UINT nmu, nnu, mb, nb;
   const size_t nfmblks=ip->nfmblks, nfnblks=ip->nfnblks, nfkblks=ip->nfkblks;
   const size_t nmblks=ip->npmblks+nfmblks, nnblks=ip->npnblks+nfnblks;
   const size_t lda=ip->lda, ldb=ip->ldb;
   const size_t incAk=ip->incAk, incBk=ip->incBk;
   size_t inca, incb, k;
   TYPE *an, *bn;
   #ifdef TCPLX
      ammkern_t amm_bn=ip->amm_bn, amm_b1=ip->amm_b1;
      size_t cinca, cincb;
      const TYPE *alpA=ip->alpA, *alpB=ip->alpB;
   #else
      const TYPE alpA=ip->alpA, alpB=ip->alpB;
   #endif

   if (i == nmblks-1)
   {
      nmu = ip->nmuF;
      mb = ip->mF;
      inca = (i < nfmblks) ? ip->szA : ip->pszA;
   }
   else if (i < nfmblks)
   {
      nmu = ip->nmu;
      mb = ip->mb;
      inca = ip->szA;
   }
   else
   {
      nmu = ip->pnmu;
      inca = ip->pszA;
      mb = ip->pmb;
   }

   if (j == nnblks-1)
   {
      nnu = ip->nnuF;
      nb = ip->nF;
      incb = (j < nfnblks) ? ip->szB : ip->pszB;
   }
   else if (j < nfnblks)
   {
      nnu = ip->nnu;
      incb = ip->szB;
      nb = ip->nb;
   }
   else
   {
      nnu = ip->pnnu;
      incb = ip->pszB;
      nb = ip->pnb;
   }

#ifdef TCPLX
   cinca = inca;
   cincb = incb;
   inca = (MV&1) ? inca+inca : 0;
   incb = (MV&2) ? incb+incb : 0;
/*
 * =========================================================
 * Peel first iteration to handle K remainder and use beta=0
 * =========================================================
 */
/*
 * If last K block of different size, take it from end of vector to avoid
 * screwing up alignment
 */
   if (kb0 != kb)
   {
      TYPE *iA, *rA, *iB, *rB;
      an = a; bn = b;
      iA = a + nfkblks*inca;
      rA = iA + cinca;
      if (A)
         a2blk(kb0, mb, alpA, A+incAk*nfkblks, lda, rA, iA);
      iB = b + nfkblks*incb;
      rB=iB+cincb;
      if (B)
         b2blk(kb0, nb, alpB, B+incBk*nfkblks, ldb, rB, iB);
      if (iC)
      {
         amm(nmu, nnu, KB0, iA, iB, rC, rA, iB, iC);
         amm(nmu, nnu, KB0, rA, iB, iC, rA, rB, rC);
         ip->ammK1_bn(nmu, nnu, KB0, rA, rB, rC, iA, rB, iC);
         ip->ammK1_b1(nmu, nnu, KB0, iA, rB, iC, an, bn, rC);
      }
   }
/*
 * Otherwise just take first block, so any hardware prefetch isn't confused
 */
   else
   {
      TYPE *iA=a, *rA=iA+cinca, *iB=b, *rB=iB+cincb;
      an=a+inca; bn=b+incb;
      if (A)
      {
         a2blk(kb0, mb, alpA, A, lda, rA, iA);
         A += incAk;
      }
      if (B)
      {
         b2blk(kb0, nb, alpB, B, ldb, rB, iB);
         B += incBk;
      }
      if (iC)
      {
         amm(nmu, nnu, kb, iA, iB, rC, rA, iB, iC);
         amm(nmu, nnu, kb, rA, iB, iC, rA, rB, rC);
         ip->ammK1_bn(nmu, nnu, kb, rA, rB, rC, iA, rB, iC);
         ip->ammK1_b1(nmu, nnu, kb, iA, rB, iC, an, bn, rC);
      }
   }
   a = an;
   b = bn;

   for (k=0; k < nfkblks; k++)
   {
      TYPE *iA=a, *rA=iA+cinca, *iB=b, *rB=iB+cincb;
      an=a+inca; bn=b+incb;
      if (A)
      {
         a2blk(kb, mb, alpA, A, lda, rA, iA);
         A += incAk;
      }
      if (B)
      {
         b2blk(kb, nb, alpB, B, ldb, rB, iB);
         B += incBk;
      }
      an = a + inca;
      bn = b + incb;
      if (iC)
      {
         amm_bn(nmu, nnu, kb, iA, iB, rC, rA, iB, iC);
         amm_b1(nmu, nnu, kb, rA, iB, iC, rA, rB, rC);
         amm_bn(nmu, nnu, kb, rA, rB, rC, iA, rB, iC);
         amm_b1(nmu, nnu, kb, iA, rB, iC, an, bn, rC);
      }

      a = an;
      b = bn;
   }
   if (blk2c)
      blk2c(mb, nb, ip->alpC, rC, iC, beta, C, ip->ldc);
#else
   if (!(MV&1))
      inca = 0;
   if (!(MV&2))
      incb = 0;
/*
 * =========================================================
 * Peel first iteration to handle K remainder and use beta=0
 * =========================================================
 */
/*
 * If last K block of different size, take it from end of vector to avoid
 * screwing up alignment
 */
   if (kb0 != kb)
   {
      an=a; bn=b;
      a += nfkblks*inca;
      if (A)
         a2blk(kb0, mb, alpA, A+incAk*nfkblks, lda, a);
      b += nfkblks*incb;
      if (B)
         b2blk(kb0, nb, alpB, B+incBk*nfkblks, ldb, b);
      if (rC)
         amm(nmu, nnu, KB0, a, b, rC, an, bn, rC);
   }
/*
 * Otherwise just take first block, so any hardware prefetch isn't confused
 */
   else
   {
      an=a+inca; bn=b+incb;
      if (A)
      {
         a2blk(kb0, mb, alpA, A, lda, a);
         A += incAk;
      }
      if (B)
      {
         b2blk(kb0, nb, alpB, B, ldb, b);
         B += incBk;
      }
      if (rC)
         amm(nmu, nnu, KB0, a, b, rC, an, bn, rC);
   }
   amm = ip->amm_b1;
   a = an;
   b = bn;

   for (k=0; k < nfkblks; k++)
   {
      if (A)
      {
         a2blk(kb, mb, alpA, A, lda, a);
         A += incAk;
      }
      if (B)
      {
         b2blk(kb, nb, alpB, B, ldb, b);
         B += incBk;
      }
      an = a + inca;
      bn = b + incb;
      if (rC)
         amm(nmu, nnu, kb, a, b, rC, an, bn, rC);
      a = an;
      b = bn;
   }
   if (blk2c)
      blk2c(mb, nb, ip->alpC, rC, beta, C, ip->ldc);
#endif
}

/*
 * MVS: bit pattern with following meaning by bit position:
 *   0: set: move A workspace;  unset: A work only mb*nb
 *   1: set: move B workspace;  unset: B work only kb*nb
 *   2: set: A wrkspc hold full mat;  unset: A wrkspc holds panel only
 *      if bit 0 unset, this bit has no affect
 *   3: set: B wrkspc hold full mat;  unset: B wrkspc holds panel only
 *      if bit 1 unset, this bit has no affect
 *
 * NOTE: we assume global (i,j), but this function assumes any j indexing
 *       has been rolled into the B/C ptrs.  It does only i indexing for A&C.
 *       When i|j > 0, we do *NOT* move around in the workspace as we would
 *       if we copied them ourselves.  This allows callers that don't need
 *       full workspace (eg., triangular loops) to avoid allocating it.
 */
void Mjoin(PATL,iploopsMK)
(
   ipinfo_t *ip,
   size_t i,              /* which glbl rowblk to start at (ignore blks <i) */
   size_t j,              /* which global col block of C is being computed */
   const TYPE *A,         /* if non-NULL, A to cpy, else ignored */
   const TYPE *B,         /* if non-NULL, Portion B to cpy, else ignored */
   TYPE *C,               /* original C to write ans to at end of Kloop */
   const int MVS,         /* 1:move A workspace, 2:move B wrkspc, 3: mv both */
   TYPE *a,               /* wrkspace for block-copied A */
   TYPE *b,               /* wrkspace for block-copied B */
   TYPE *rC,              /* ptr to real portion of C access-major workspace */
   TYPE *iC,              /* ptr to complex C workspace */
   const SCALAR beta,     /* scaling factor for C */
   const ablk2cmat_t blk2c   /* C = alpC*(rC,iC) + beta*C */
)
{
   const size_t nfmblks=ip->nfmblks, nmblks=ip->npmblks+nfmblks;
   const size_t BA=((MVS&1)&&(MVS&4)) ? (~0):0;
   const TYPE *pB=B;
   if (i < nfmblks)
   {
      const size_t incAm = (A) ? ip->incAm:0, incCm = (ip->mb)SHIFT;
      const size_t pansz = BA & ((ip->nfkblks + 1)*((ip->szA)SHIFT));
      if (i)
      {
         A += i * incAm;
         C += i * incCm;
      }
      for (; i < nfmblks; i++)
      {
         Mjoin(PATL,iploopsK)(ip, i, j, A, pB, C, MVS, a, b,
                              rC, iC, beta, blk2c);
         A += incAm;
         C += incCm;
         a += pansz;
         pB = NULL;  /* copy B col-panel only on first iteration */
      }
   }
   else if (i && nfmblks)
   {
      if (A)
         A += nfmblks * ip->incAm;
      C += nfmblks * ((ip->mb)SHIFT);
   }
   if (i < nmblks)
   {
      const size_t incAm = (A) ? ip->pincAm:0, incCm = (ip->pmb)SHIFT;
      const size_t pansz = BA & ((ip->nfkblks + 1)*((ip->pszA)SHIFT));
      if (i > nfmblks)
      {
         size_t k = i - nfmblks;
         A += k * incAm;
         C += k * incCm;
      }
      for (; i < nmblks; i++)
      {
         Mjoin(PATL,iploopsK)(ip, i, j, A, pB, C, MVS, a, b,
                              rC, iC, beta, blk2c);
         A += incAm;
         C += incCm;
         a += pansz;
         pB = NULL;  /* copy B col-panel only on first iteration */
      }
   }
}
/*
 * MVS: bit pattern with following meaning by bit position:
 *   0: set: move A workspace;  unset: A work only mb*nb
 *   1: set: move B workspace;  unset: B work only kb*nb
 *   2: set: A wrkspc hold full mat;  unset: A wrkspc holds panel only
 *      if bit 0 unset, this bit has no affect
 *   3: set: B wrkspc hold full mat;  unset: B wrkspc holds panel only
 *      if bit 1 unset, this bit has no affect
 *
 * NOTE: we assume global (i,j), but this function assumes any i indexing
 *       has been rolled into the A/C ptrs.  It does only j indexing for B&C.
 *       When i|j > 0, we do *NOT* move around in the workspace as we would
 *       if we copied them ourselves.  This allows callers that don't need
 *       full workspace (eg., triangular loops) to avoid allocating it.
 */
void Mjoin(PATL,iploopsNK)
(
   ipinfo_t *ip,
   size_t i,              /* which global row block of C is being computed */
   size_t j,              /* which glbl colblk to start at (ignore blks <j) */
   const TYPE *A,         /* if non-NULL, A to cpy, else ignored */
   const TYPE *B,         /* if non-NULL, Portion B to cpy, else ignored */
   TYPE *C,               /* original C to write ans to at end of Kloop */
   const int MVS,         /* 1:move A workspace, 2:move B wrkspc, 3: mv both */
   TYPE *a,               /* wrkspace for block-copied A */
   TYPE *b,               /* wrkspace for block-copied B */
   TYPE *rC,              /* ptr to real portion of C access-major workspace */
   TYPE *iC,              /* ptr to complex C workspace */
   const SCALAR beta,     /* scaling factor for C */
   const ablk2cmat_t blk2c   /* C = alpC*(rC,iC) + beta*C */
)
{
   const size_t nfnblks=ip->nfnblks, nnblks=ip->npnblks+nfnblks;
   const size_t BB=((MVS&2)&&(MVS&8)) ? (~0):0;
   const TYPE *pA=A;
   if (j < nfnblks)
   {
      const size_t incBn = (B) ? ip->incBn:0, incCn = ip->ldc*((ip->nb)SHIFT);
      const size_t pansz = BB & ((ip->nfkblks + 1)*((ip->szB)SHIFT));
      if (j)
      {
         B += j * incBn;
         C += j * incCn;
      }
      for (; j < nfnblks; j++)
      {
         Mjoin(PATL,iploopsK)(ip, i, j, pA, B, C, MVS, a, b,
                              rC, iC, beta, blk2c);
         B += incBn;
         C += incCn;
         b += pansz;
         pA = NULL;  /* copy A row-panel only on first iteration */
      }
   }
   else if (j && nfnblks)
   {
      if (B)
         B += nfnblks * ip->incBn;
      C += nfnblks * ip->ldc * ((ip->mb)SHIFT);
   }
   if (j < nnblks)
   {
      const size_t incBn = (B) ? ip->pincBn:0, incCn = ip->ldc*((ip->pnb)SHIFT);
      const size_t pansz = BB & ((ip->nfkblks + 1)*((ip->pszB)SHIFT));
      if (j > nfnblks)
      {
         size_t k = j - nfnblks;
         B += k * incBn;
         C += k * incCn;
      }
      for (; j < nnblks; j++)
      {
         Mjoin(PATL,iploopsK)(ip, i, j, pA, B, C, MVS, a, b,
                              rC, iC, beta, blk2c);
         B += incBn;
         C += incCn;
         b += pansz;
         pA = NULL;  /* copy A row-panel only on first iteration */
      }
   }
}

/*
 * A & B ptrs should always be base ptrs to this func.
 * MVS: bit pattern with following meaning by bit position:
 *   0: set: move A workspace;  unset: A work only mb*nb
 *   1: set: move B workspace;  unset: B work only kb*nb
 *   2: set: A wrkspc hold full mat;  unset: A wrkspc holds panel only
 *      if bit 0 unset, this bit has no affect
 *   3: set: B wrkspc hold full mat;  unset: B wrkspc holds panel only
 *      if bit 1 unset, this bit has no affect
 *
 */
void Mjoin(PATL,iploopsNMK)
(
   ipinfo_t *ip,
   size_t i,              /* which glbl rowblk to start at (ignore blks <i) */
   size_t j,              /* which global col block of C is being computed */
   const TYPE *A,         /* if non-NULL, A to cpy, else ignored */
   const TYPE *B,         /* if non-NULL, Portion B to cpy, else ignored */
   TYPE *C,               /* original C to write ans to at end of Kloop */
   const int MVS,         /* 1:move A workspace, 2:move B wrkspc, 3: mv both */
   TYPE *a,               /* wrkspace for block-copied A */
   TYPE *b,               /* wrkspace for block-copied B */
   TYPE *rC,              /* ptr to real portion of C access-major workspace */
   TYPE *iC,              /* ptr to complex C workspace */
   const SCALAR beta,     /* scaling factor for C */
   const ablk2cmat_t blk2c   /* C = alpC*(rC,iC) + beta*C */
)
{
   const size_t nfnblks=ip->nfnblks, nnblks=ip->npnblks+nfnblks, ldc=ip->ldc;
   const size_t BMSKA=((MVS&1)&&(MVS&4)) ? 0:(~0);
   const size_t BMSKB=((MVS&2)&&(MVS&8)) ? (~0):0;
   if (j < nfnblks)
   {
      const size_t pansz = BMSKB & ((ip->nfkblks+1)*((ip->szB)SHIFT));
      const size_t incBn = ip->incBn, incCn = ldc*((ip->nb)SHIFT);
      if (j)
      {
         B += j*incBn;
         C += j*incCn;
      }
      for (; j < nfnblks; j++)
      {
         Mjoin(PATL,iploopsMK)(ip, i, j, A, B, C, MVS, a, b,
                               rC, iC, beta, blk2c);
         A = (TYPE*)(((size_t)A) & BMSKA);
         b += pansz;
         B += incBn;
         C += incCn;
      }
   }
   else if (j && nfnblks)
   {
      B += j * ip->incBn;
      C += j * ldc * ((ip->nb)SHIFT);
   }
   if (j < nnblks)
   {
      const size_t incBn = ip->pincBn, incCn = ldc*((ip->pnb)SHIFT);
      const size_t pansz = BMSKB & ((ip->nfkblks+1)*((ip->pszB)SHIFT));
      if (j > nfnblks)
      {
         size_t k = j - nfnblks;
         B += k * incBn;
         C += k * incCn;
      }
      for(; j < nnblks; j++)
      {
         Mjoin(PATL,iploopsMK)(ip, i, j, A, B, C, MVS, a, b,
                               rC, iC, beta, blk2c);
         A = (TYPE*)(((size_t)A) & BMSKA);
         b += pansz;
         B += incBn;
         C += incCn;
      }
   }
}
/*
 * A & B ptrs should always be base ptrs to this func.
 * MVS: bit pattern with following meaning by bit position:
 *   0: set: move A workspace;  unset: A work only mb*nb
 *   1: set: move B workspace;  unset: B work only kb*nb
 *   2: set: A wrkspc hold full mat;  unset: A wrkspc holds panel only
 *      if bit 0 unset, this bit has no affect
 *   3: set: B wrkspc hold full mat;  unset: B wrkspc holds panel only
 *      if bit 1 unset, this bit has no affect
 *
 */
void Mjoin(PATL,iploopsMNK)
(
   ipinfo_t *ip,
   size_t i,              /* which glbl rowblk to start at (ignore blks <i) */
   size_t j,              /* which global col block of C is being computed */
   const TYPE *A,         /* if non-NULL, A to cpy, else ignored */
   const TYPE *B,         /* if non-NULL, Portion B to cpy, else ignored */
   TYPE *C,               /* original C to write ans to at end of Kloop */
   const int MVS,         /* 1:move A workspace, 2:move B wrkspc, 3: mv both */
   TYPE *a,               /* wrkspace for block-copied A */
   TYPE *b,               /* wrkspace for block-copied B */
   TYPE *rC,              /* ptr to real portion of C access-major workspace */
   TYPE *iC,              /* ptr to complex C workspace */
   const SCALAR beta,     /* scaling factor for C */
   const ablk2cmat_t blk2c   /* C = alpC*(rC,iC) + beta*C */
)
{
   const size_t nfmblks=ip->nfmblks, nmblks=ip->npmblks+nfmblks, ldc=ip->ldc;
   const size_t BMSKA=((MVS&1)&&(MVS&4)) ? (~0):0;
   const size_t BMSKB=((MVS&2)&&(MVS&8)) ? 0:(~0);
   if (i < nfmblks)
   {
      const size_t pansz = BMSKA & ((ip->nfkblks+1)*((ip->szA)SHIFT));
      const size_t incAm = ip->incAm, incCm = (ip->mb)SHIFT;
      if (i)
      {
         A += i*incAm;
         C += i*incCm;
      }
      for (; i < nfmblks; i++)
      {
         Mjoin(PATL,iploopsNK)(ip, i, j, A, B, C, MVS, a, b,
                               rC, iC, beta, blk2c);
         B = (TYPE*)(((size_t)B) & BMSKB);
         a += pansz;
         A += incAm;
         C += incCm;
      }
   }
   else if (i && nfmblks)
   {
      A += i * ip->incAm;
      C += i * ((ip->mb)SHIFT);
   }
   if (i < nmblks)
   {
      const size_t incAm = ip->pincAm, incCm = (ip->pmb)SHIFT;
      const size_t pansz = BMSKA & ((ip->nfkblks+1)*((ip->pszA)SHIFT));
      if (i > nfmblks)
      {
         size_t k = i - nfmblks;
         A += k * incAm;
         C += k * incCm;
      }
      for(; i < nmblks; i++)
      {
         Mjoin(PATL,iploopsNK)(ip, i, j, A, B, C, MVS, a, b,
                               rC, iC, beta, blk2c);
         B = (TYPE*)(((size_t)B) & BMSKB);
         a += pansz;
         A += incAm;
         C += incCm;
      }
   }
}
