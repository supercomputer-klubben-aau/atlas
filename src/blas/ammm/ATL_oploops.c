/*
 * These loops are used for outer product GEMM formulations,
 * where K <= ATL_rkAMM_LASTKB.
 *
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(Mjoin(ATLAS_PRE,amm),_sum.h))

void Mjoin(PATL,ammPadZero)(ATL_CUINT npan, ATL_CUINT npad, ATL_CUINT szPan,
#ifdef TCPLX
                            TYPE *rB, TYPE *iB)
#else
                            TYPE *B)
#endif
/*
 * For M- or N-vectorized kerns, zeros last npad entries in each K-panel
 */
{
   ATL_UINT k, i;
   #ifdef TCPLX
      rB += szPan - npad;
      iB += szPan - npad;
      for (k=0; k < npan; k++, rB += szPan, iB += szPan)
         for (i=0; i < npad; i++)
            rB[i] = iB[i] = ATL_rzero;
   #else
      B += szPan - npad;
      for (k=0; k < npan; k++, B += szPan)
         for (i=0; i < npad; i++)
            B[i] = ATL_rzero;
   #endif
}
/*
 * This function used primarily by threaded rank-K update cases.  It has
 * following options:
 *   If A/B i/j block can be copied by passing base A/B ptr
 *   If rC is NULL, then do only the indicated A/B copy, w/o ammm
 *   if C is NULL, then don't copy out, else do
 */
void Mjoin(PATL,opblk)
(
   opinfo_t *op,
   size_t i,      /* which global row blk of C is being computed */
   size_t j,      /* which global column blk of C is being computed */
   const TYPE *A, /* if non-NULL, base A ptr to copy */
   const TYPE *B, /* if non-NULL, base B ptr to copy */
   TYPE *C,       /* C <- alpha * A * B + beta * C */
   TYPE *pA,      /* wrkspace for block-copied A */
   TYPE *pAn,     /* next A wrkspc to be prefetched */
   TYPE *pB,      /* wrkspace for block-copied B */
   TYPE *pBn,     /* next B wrkspc to be prefetched */
   TYPE *rC,      /* if non-NULL: ptr to real part of access-major C wrkspc */
   TYPE *iC       /* ptr to cplx portion of access-major C workspace */
)
{
   const size_t nfnblks=op->nfnblks, nnblks=nfnblks+op->npnblks;
   const size_t nfmblks=op->nfmblks, nmblks=nfmblks+op->npmblks;
   size_t m, n;
   const unsigned int szA = ((i < op->nfmblks) ? op->szA:op->pszA);
   const unsigned int szB = ((j < op->nfnblks) ? op->szB:op->pszB);
   #ifdef TCPLX
      TYPE *iA=pA, *iB=pB;
      TYPE *rA = pA + szA;
      TYPE *rB = pB + szB;
   #endif
   ATL_CUINT kb = op->KB, K=op->kb;
   ATL_UINT nb, nnu, mb, nmu;

   if (j < nfnblks)
   {
      nb = op->nb;
      nnu = op->nnu;
      n = j;
   }
   else if (j != nnblks-1)
   {
      nb = op->pnb;
      nnu = op->pnnu;
      n = nfnblks;
   }
   else
   {
      nb = op->nF;
      nnu = op->nnuF;
      n = nfnblks;
   }
   j -= n;
   if (i < nfmblks)
   {
      mb = op->mb;
      nmu = op->nmu;
      m = i;
   }
   else if (i != nmblks-1)
   {
      mb = op->pmb;
      nmu = op->pnmu;
      m = nfmblks;
   }
   else
   {
      mb = op->mF;
      nmu = op->nmuF;
      m = nfmblks;
   }
   i -= m;
   if (B)
   {
      B += op->incBn * n + op->pincBn * j;
      #ifdef TCPLX
         op->b2blk(K, nb, op->alpB, B, op->ldb, rB, iB);
      #else
         op->b2blk(K, nb, op->alpB, B, op->ldb, pB);
      #endif
   }
   if (A)
   {
      A += op->incAm * m + op->pincAm * i;
      #ifdef TCPLX
         op->a2blk(K, mb, op->alpA, A, op->lda, rA, iA);
      #else
         op->a2blk(K, mb, op->alpA, A, op->lda, pA);
      #endif
   }

   if (rC)   /* want to do matrix multiplication */
   {
      #ifdef TCPLX
         ammkern_t amm_b0=op->amm_b0, amm_b1=op->amm_b1, amm_bn=op->amm_bn;
         amm_b0(nmu, nnu, kb, iA, iB, rC, rA, iB, iC);
         amm_b0(nmu, nnu, kb, rA, iB, iC, rA, rB, rC);
         amm_bn(nmu, nnu, kb, rA, rB, rC, iA, rB, iC);
         amm_b1(nmu, nnu, kb, iA, rB, iC, pAn, pBn, rC);
         if (C)
         {
            C += (op->nb * n + op->pnb * j)*((op->ldc)SHIFT) +
                 ((op->mb * m + op->pmb * i)SHIFT);
            op->blk2C(mb, nb, op->ONE, rC, iC, op->beta, C, op->ldc);
         }
      #else
         op->amm_b0(nmu, nnu, kb, pA, pB, rC, pAn, pBn, rC);
         if (C)
         {
            C += (op->nb * n + op->pnb * j)*(op->ldc) +
                 (op->mb * m + op->pmb * i);
            op->blk2C(mb, nb, ATL_rone, rC, op->beta, C, op->ldc);
         }
      #endif
   }
}

void Mjoin(PATL,oploopsM)
(
   opinfo_t *op,
   size_t i,      /* which global row blk of C is being computed */
   size_t j,      /* which global column blk of C is being computed */
   const TYPE *A, /* if non-NULL, portion of A to copy */
   const TYPE *B, /* if non-NULL, portion of B to copy */
   TYPE *C,       /* C <- alpha * A * B + beta * C */
   int MV,        /* 1: move a */
   TYPE *a,       /* wrkspace for block-copied A */
   TYPE *b,       /* wrkspace for block-copied B */
   TYPE *rC,      /* ptr to real portion of access-major C workspace */
   TYPE *iC       /* ptr to cplx portion of access-major C workspace */
)
{
   #ifdef TCPLX
      TYPE *iA=a, *rA, *iB=b, *rB=b+op->szB;
      const TYPE *alpA=op->alpA, *alpB=op->alpB, *beta=op->beta;
      const ammkern_t amm0=op->amm_b0, amm1=op->amm_b1, ammN=op->amm_bn;
   #else
      const ammkern_t amm=op->amm_b0;
   #endif
   const cm2am_t a2blk=op->a2blk, b2blk=op->b2blk;
   const ablk2cmat_t blk2C=op->blk2C;
   const size_t nfmblks=op->nfmblks, nmblks=op->npmblks+nfmblks, lda=op->lda;
   const size_t nfnblks=op->nfnblks, nnblks=op->npnblks+nfnblks;
   #ifndef TCPLX
      const TYPE alpA=op->alpA, alpB=op->alpB, beta=op->beta;
   #endif
   ATL_CUINT kb=op->kb, KB=op->KB;
   ATL_CUINT mu=op->mu, nu=op->nu;
   ATL_UINT nmu, mb, nnu=op->nnu, nb=op->nb;

/*
 * If needed, copy B block before starting operations
 */
   if (B)
   {
      ATL_UINT szB=op->szB;
      if (j)  /* need to increment B & C before copy */
      {
         const size_t incCn=(op->ldc)SHIFT;
         size_t nblks=Mmin(j, nfnblks);
         B += nblks * op->incBn;
         C += nblks*nb*incCn;
         if (j >= nfnblks)  /* need to switch to smaller nb */
         {
            szB = op->pszB;
            #ifdef TCPLX
               rB = iB + szB;
            #endif
            nb = op->pnb;
            nnu = op->pnnu;
            nblks = j - nfnblks;
            B += nblks*op->pincBn;
            C += nblks*nb*incCn;
         }
         else
            szB = op->szB;
      }
      if (j == nnblks-1 && j >= nfnblks)
      {
         szB = op->pszB;
         nb = op->nF;
         nnu = op->nnuF;
      }
      #ifdef TCPLX
         b2blk(kb, nb, alpB, B, op->ldb, rB, iB);
      #else
         b2blk(kb, nb, alpB, B, op->ldb, b);
      #endif
   }
   else   /* still need to inc C ptr & set nnu/nb */
   {
      size_t incBn=op->incBn, incCn=(op->ldc)SHIFT;
      size_t nblks=Mmin(j, nfnblks);
      C += nblks*nb*incCn;
      if (j >= nfnblks)  /* need to switch to smaller nb */
      {
         #ifdef TCPLX
            rB = iB + op->pszB;
         #endif
         nb = op->pnb;
         nnu = op->pnnu;
         nblks = j - nfnblks;
         C += nblks*nb*incCn;
         if (j == nnblks-1)
         {
            nb = op->nF;
            nnu = op->nnuF;
         }
      }
   }
/*
 * If present block is one of first large blocks, loop over large blks
 */
   if (i < nfmblks)
   {
      const size_t incAm = (A) ? op->incAm:0, incCm = (op->mb)SHIFT, m;
      ATL_UINT szA = op->szA;

      if (i)               /* if we need to skip some blocks */
      {                    /* Increment A/C ptrs as needed for i */
         if (A)
            A += i * incAm;
         C += i * incCm;
      }
      #ifdef TCPLX
         rA = iA + szA;
      #endif
      nmu = op->nmu;
      mb = op->mb;
      for (; i < nfmblks; i++)
      {
         #ifdef TCPLX
            TYPE *an=iA;
         #else
            TYPE *an=a;
         #endif
         if (A)
         {
            #ifdef TCPLX
               a2blk(kb, mb, alpA, A, lda, rA, iA);
            #else
               a2blk(kb, mb, alpA, A, lda, a);
            #endif
            A = (nmblks != 1) ? (A+incAm) : NULL;
         }
         if (MV&1)
         #ifdef TCPLX
            an = rA + szA;
         #else
            an += szA;
         #endif
         #ifdef TCPLX
            amm0(nmu, nnu, KB, iA, iB, rC, rA, iB, iC);
            amm0(nmu, nnu, KB, rA, iB, iC, rA, rB, rC);
            ammN(nmu, nnu, KB, rA, rB, rC, iA, rB, iC);
            amm1(nmu, nnu, KB, iA, rB, iC, an, iB, iC);
            blk2C(mb, nb, op->ONE, rC, iC, beta, C, op->ldc);
            iA = an;
            rA = an + szA;
         #else
            amm(nmu, nnu, KB, a, b, rC, an, b, rC);
            blk2C(mb, nb, ATL_rone, rC, beta, C, op->ldc);
            a = an;
         #endif
         C += incCm;
      }
   }
/*
 * If present block past end of large blocks, increment A/C appropriately
 */
   else if (i && nfmblks)
   {
      if (A)
         A += nfmblks * op->incAm;
      C += nfmblks * ((op->mb)SHIFT);
   }
/*
 * If there are still small blocks to do after doing all large blocks
 */
   if (i < nmblks)
   {
      const size_t incAm = (A) ? op->pincAm:0, incCm = (op->pmb)SHIFT;
      ATL_UINT szA = op->pszA;

      if (i > nfmblks)  /* I need to move past some small blocks too */
      {
         size_t nblks = i - nfmblks;
         A += nblks * incAm;
         C += nblks * incCm;
      }
      mb = op->pmb;
      nmu = op->pnmu;
      for (; i < nmblks; i++)
      {
         #ifdef TCPLX
            TYPE *an=iA;
         #else
            TYPE *an=a;
         #endif
         if (i == nmblks-1)
         {
            mb = op->mF;
            nmu = op->nmuF;
         }
         #ifdef TCPLX
            rA = iA + szA;
         #endif
         if (A)
         {
            #ifdef TCPLX
               a2blk(kb, mb, alpA, A, lda, rA, iA);
            #else
               a2blk(kb, mb, alpA, A, lda, a);
            #endif
            A = (nmblks != 1) ? (A+incAm) : NULL;
         }
         #ifdef TCPLX
            if (MV&1)
               an = rA + szA;
            amm0(nmu, nnu, KB, iA, iB, rC, rA, iB, iC);
            amm0(nmu, nnu, KB, rA, iB, iC, rA, rB, rC);
            ammN(nmu, nnu, KB, rA, rB, rC, iA, rB, iC);
            amm1(nmu, nnu, KB, iA, rB, iC, an, iB, iC);
            blk2C(mb, nb, op->ONE, rC, iC, beta, C, op->ldc);
            iA = an;
            rA = an + szA;
         #else
            if (MV&1)
               an += szA;
            amm(nmu, nnu, KB, a, b, rC, an, b, rC);
            blk2C(mb, nb, ATL_rone, rC, beta, C, op->ldc);
            a = an;
         #endif
         C += incCm;
      }
   }
}
void Mjoin(PATL,oploopsN)
(
   opinfo_t *op,
   size_t i,      /* which global row blk of C is being computed */
   size_t j,      /* which global column blk of C is being computed */
   const TYPE *A, /* if non-NULL, portion of A to copy */
   const TYPE *B, /* if non-NULL, portion of B to copy */
   TYPE *C,       /* C <- alpha * A * B + beta * C */
   int MV,        /* 2: move b */
   TYPE *a,       /* wrkspace for block-copied A */
   TYPE *b,       /* wrkspace for block-copied B */
   TYPE *rC,      /* ptr to real portion of access-major C workspace */
   TYPE *iC       /* ptr to cplx portion of access-major C workspace */
)
{
   #ifdef TCPLX
      TYPE *iA=a, *rA=a+op->szA, *iB=b, *rB;
      const TYPE *alpA=op->alpA, *alpB=op->alpB, *beta=op->beta;
      const ammkern_t amm0=op->amm_b0, amm1=op->amm_b1, ammN=op->amm_bn;
   #else
      const ammkern_t amm=op->amm_b0;
   #endif
   const cm2am_t a2blk=op->a2blk, b2blk=op->b2blk;
   const ablk2cmat_t blk2C=op->blk2C;
   const size_t nfmblks=op->nfmblks, nmblks=op->npmblks+nfmblks, ldb=op->ldb;
   const size_t nfnblks=op->nfnblks, nnblks=op->npnblks+nfnblks, ldc=op->ldc;
   #ifndef TCPLX
      const TYPE alpA=op->alpA, alpB=op->alpB, beta=op->beta;
   #endif
   ATL_CUINT kb=op->kb, KB=op->KB;
   ATL_CUINT mu=op->mu, nu=op->nu;
   ATL_UINT nmu=op->nmu, mb=op->mb, nnu, nb, szA = op->szA;

/*
 * If needed, copy A block before starting operations
 */
   if (A)
   {
      if (i)  /* need to increment A & C before copy */
      {
         size_t nblks=Mmin(i, nfmblks);
         A += nblks * op->incAm;
         C += nblks*(mb SHIFT);
         if (i >= nfmblks)  /* need to switch to smaller mb */
         {
            szA = op->pszA;
            #ifdef TCPLX
               rA = iA + szA;
            #endif
            mb = op->pmb;
            nmu = op->pnmu;
            nblks = i - nfmblks;
            A += nblks*op->pincAm;
            C += nblks*(mb SHIFT);
         }
      }
      if (i == nmblks-1 && i >= nfmblks)
      {
         szA = op->pszA;
         mb = op->mF;
         nmu = op->nmuF;
      }
      #ifdef TCPLX
         a2blk(kb, mb, alpA, A, op->lda, rA, iA);
      #else
         a2blk(kb, mb, alpA, A, op->lda, a);
      #endif
   }
   else   /* still need to inc C ptr & set nnu/nb */
   {
      size_t nblks=Mmin(i, nfmblks);
      C += nblks*(mb SHIFT);
      if (i >= nfmblks)  /* need to switch to smaller mb */
      {
         szA = op->pszA;
         #ifdef TCPLX
            rA = iA + szA;
         #endif
         mb = op->pmb;
         nmu = op->pnmu;
         nblks = i - nfmblks;
         C += nblks*(mb SHIFT);
         if (i == nmblks-1)
         {
            mb = op->mF;
            nmu = op->nmuF;
         }
      }
   }
/*
 * If present block is one of first large blocks, loop over large blks
 */
   if (j < nfnblks)
   {
      const size_t incBn = (B) ? op->incBn:0, incCn = ldc*(op->nb)SHIFT;
      ATL_UINT szB = op->szB;
      ATL_CUINT incBw = (MV&2) ? ((szB)SHIFT):0;

      if (j)               /* if we need to skip some blocks */
      {                    /* Increment B/C ptrs as needed for j */
         if (B)
            B += j * incBn;
         C += j * incCn;
      }
      #ifdef TCPLX
         rB = iB + szB;
      #endif
      nnu = op->nnu;
      nb = op->nb;
      for (; j < nfnblks; j++)
      {
         #ifdef TCPLX
            TYPE *bn=iB+incBw;
         #else
            TYPE *bn=b+incBw;
         #endif
         if (B)
         {
            #ifdef TCPLX
               b2blk(kb, nb, alpB, B, ldb, rB, iB);
            #else
               b2blk(kb, nb, alpB, B, ldb, b);
            #endif
            B = (nnblks != 1) ? (B+incBn) : NULL;
         }
         #ifdef TCPLX
            amm0(nmu, nnu, KB, iA, iB, rC, rA, iB, iC);
            amm0(nmu, nnu, KB, rA, iB, iC, rA, rB, rC);
            ammN(nmu, nnu, KB, rA, rB, rC, iA, rB, iC);
            amm1(nmu, nnu, KB, iA, rB, iC, iA, bn, iC);
            blk2C(mb, nb, op->ONE, rC, iC, beta, C, ldc);
            iB = bn;
            rB = bn + szB;
         #else
            amm(nmu, nnu, KB, a, b, rC, a, bn, rC);
            blk2C(mb, nb, ATL_rone, rC, beta, C, ldc);
            b = bn;
         #endif
         C += incCn;
      }
   }
/*
 * If present block past end of large blocks, increment B/C appropriately
 */
   else if (j && nfnblks)
   {
      if (B)
         B += nfnblks * op->incBn;
      C += nfnblks * (op->nb * ldc)SHIFT;
   }
/*
 * If there are still small blocks to do after doing all large blocks
 */
   if (j < nnblks)
   {
      const size_t incBn = (B) ? op->pincBn:0, incCn = ldc*(op->pnb)SHIFT;
      ATL_UINT szB = op->pszB;
      ATL_CUINT incBw = (MV&2) ? ((szB)SHIFT):0;

      if (j > nfnblks)  /* I need to move past some small blocks too */
      {
         size_t nblks = j - nfnblks;
         B += nblks * incBn;
         C += nblks * incCn;
      }
      nb = op->pnb;
      nnu = op->pnnu;
      for (; j < nnblks; j++)
      {
         #ifdef TCPLX
            TYPE *bn = iB + incBw;
         #else
            TYPE *bn = b + incBw;
         #endif
         if (j == nnblks-1)
         {
            nb = op->nF;
            nnu = op->nnuF;
         }
         #ifdef TCPLX
            rB = iB + szB;
         #endif
         if (B)
         {
            #ifdef TCPLX
               b2blk(kb, nb, alpB, B, ldb, rB, iB);
            #else
               b2blk(kb, nb, alpB, B, ldb, b);
            #endif
            B = (nnblks != 1) ? (B+incBn) : NULL;
         }
         #ifdef TCPLX
            amm0(nmu, nnu, KB, iA, iB, rC, rA, iB, iC);
            amm0(nmu, nnu, KB, rA, iB, iC, rA, rB, rC);
            ammN(nmu, nnu, KB, rA, rB, rC, iA, rB, iC);
            amm1(nmu, nnu, KB, iA, rB, iC, iA, bn, iC);
            blk2C(mb, nb, op->ONE, rC, iC, beta, C, op->ldc);
            iB = bn;
            rB = bn + szB;
         #else
            amm(nmu, nnu, KB, a, b, rC, a, bn, rC);
            blk2C(mb, nb, ATL_rone, rC, beta, C, op->ldc);
            b = bn;
         #endif
         C += incCn;
      }
   }
}
