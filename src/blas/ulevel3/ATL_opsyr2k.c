/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_level1.h"
#include "atlas_level2.h"
#include "atlas_kernel3.h"
#include "atlas_reflevel3.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#ifdef Conj_
   #define syr2k_OP Mjoin(PATL,opher2k)
#else
   #define syr2k_OP Mjoin(PATL,opsyr2k)
#endif

int syr2k_OP
(
   const enum ATLAS_UPLO  Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT  N,
   ATL_CSZT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   const TYPE *B,
   ATL_CSZT ldb,
   #ifdef Conj_
   const TYPE rbeta,
      #define IFLG 1
   #else
      #define IFLG 0
   const SCALAR beta,
   #endif
   TYPE *C,
   ATL_CSZT ldc
)
{
   void *vp=NULL;
   ablk2cmat_t blk2C_b0, blk2C_b1;
   ATL_SZT sz, szA, szB, nfnblks, nnblks;
   unsigned int szC;
   unsigned int mb, nb, KB;
   TYPE *aC, *bC;  /* access- and block-major C */
   TYPE *aA, *aB;  /* access-major A & B */
   void (*syput)(ATL_CSZT N, const TYPE *D, const SCALAR beta0, TYPE *A,
                 ATL_CSZT lda);
   void (*geputT)(ATL_CSZT M, ATL_CSZT N, const TYPE *A, ATL_CSZT lda,
                  const SCALAR beta, TYPE *C, ATL_CSZT ldc);
   #ifdef TCPLX
      TYPE *aCr;
      const TYPE ZERO[2] = {ATL_rzero, ATL_rzero}, *ONE;
      #ifdef Conj_
         const TYPE beta[2]={rbeta, ATL_rzero};
      #endif
   #else
      #define ZERO ATL_rzero
      #define ONE ATL_rone
      #define aCr aC
   #endif
   opinfo_t op;

   if (!N)
      return(0);
   if (SCALAR_IS_ZERO(alpha) || !K)
   {
      if (SCALAR_IS_ONE(beta))
         return(0);
      Mjoin(PATL,trscal)(Uplo, N, N, beta, C, ldc);
      #ifdef Conj_
         Mjoin(PATLU,zero)(N, C+1, (ldc+1)<<1);
      #endif
      return(0);
   }
   if (K < 3)
   {
      #if !defined(TCPLX) || defined(Conj_)
         #ifdef Conj_
         if (K == 1 && rbeta == ATL_rone && TA == AtlasNoTrans)
         #else
         if (K == 1 && SCALAR_IS_ONE(beta))
         #endif
         {
            ATL_SZT incA, incB;
            if (TA == AtlasTrans || TA == AtlasConjTrans)
            {
               incA = lda;
               incB = ldb;
            }
            else
            {
               incA = 1;
               incB = 1;
            }
            #ifdef Conj_
               Mjoin(PATL,her2)(Uplo, N, alpha, A, incA, B, incB, C, ldc);
            #else
               Mjoin(PATL,syr2)(Uplo, N, alpha, A, incA, B, incB, C, ldc);
            #endif
            return(0);
         }
      #endif
      #ifdef Conj_
         Mjoin(PATL,refher2k)(Uplo, TA, N, K, alpha, A, lda, B, ldb,
                              rbeta, C, ldc);
      #else
         Mjoin(PATL,refsyr2k)(Uplo, TA, N, K, alpha, A, lda, B, ldb,
                              beta, C, ldc);
      #endif
      return(0);
   }
   blk2C_b0 = Mjoin(PATL,opsyr2kInfo)(&op, IFLG, TA, N, K, lda, ldb, ldc,
                                      alpha, beta);
   if (!blk2C_b0)
      return(2);
   blk2C_b1 = op.blk2C;
   KB = op.KB;
   #ifdef TCPLX
      ONE = op.ONE;
      #ifdef Conj_
         if (rbeta == ATL_rone)
         {
            syput = ((Uplo == AtlasUpper) ?
                     Mjoin(PATL,her2k_putU_b1):Mjoin(PATL,her2k_putL_b1));
            geputT = Mjoin(PATL,geput1H_b1);
         }
         else if (rbeta == ATL_rzero)
         {
            syput = ((Uplo == AtlasUpper) ?
                     Mjoin(PATL,her2k_putU_b0):Mjoin(PATL,her2k_putL_b0));
            geputT = Mjoin(PATL,geput1H_b0);
         }
         else
         {
            syput = ((Uplo == AtlasUpper) ?
                     Mjoin(PATL,her2k_putU_bXi0):Mjoin(PATL,her2k_putL_bXi0));
            geputT = (rbeta == ATL_rnone) ?
               Mjoin(PATL,geput1H_bN) : Mjoin(PATL,geput1H_bX);
         }
      #else
         if (beta[1] == ATL_rzero) /* real beta */
         {
            const TYPE rb=(*beta);
            if (rb == ATL_rone)
            {
               syput = ((Uplo == AtlasUpper) ?
                        Mjoin(PATL,syr2k_putU_b1):Mjoin(PATL,syr2k_putL_b1));
               geputT = Mjoin(PATL,geput1T_b1);
            }
            else if (rb == ATL_rzero)
            {
               syput = ((Uplo == AtlasUpper) ?
                        Mjoin(PATL,syr2k_putU_b0):Mjoin(PATL,syr2k_putL_b0));
               geputT = Mjoin(PATL,geput1T_b0);
            }
            else if (rb == ATL_rnone)
            {
               syput = ((Uplo == AtlasUpper) ?
                        Mjoin(PATL,syr2k_putU_bn1):Mjoin(PATL,syr2k_putL_bn1));
               geputT = Mjoin(PATL,geput1T_bN);
            }
            else
            {
               syput = ((Uplo == AtlasUpper) ?
                       Mjoin(PATL,syr2k_putU_bXi0):Mjoin(PATL,syr2k_putL_bXi0));
               geputT = Mjoin(PATL,geput1T_bX);
            }
         }
         else
         {
            syput = ((Uplo == AtlasUpper) ?
                     Mjoin(PATL,syr2k_putU_bX):Mjoin(PATL,syr2k_putL_bX));
            geputT = Mjoin(PATL,geput1T_bX);
         }
      #endif
   #else
      if (beta == ATL_rone)
      {
         syput = ((Uplo == AtlasUpper) ?
                  Mjoin(PATL,syr2k_putU_b1):Mjoin(PATL,syr2k_putL_b1));
         geputT = Mjoin(PATL,geput1T_b1);
      }
      else if (beta == ATL_rzero)
      {
         syput = ((Uplo == AtlasUpper) ?
                  Mjoin(PATL,syr2k_putU_b0):Mjoin(PATL,syr2k_putL_b0));
         geputT = Mjoin(PATL,geput1T_b0);
      }
      else
      {
         syput = ((Uplo == AtlasUpper) ?
                  Mjoin(PATL,syr2k_putU_bX):Mjoin(PATL,syr2k_putL_bX));
         geputT = (beta == ATL_rnone) ?
            Mjoin(PATL,geput1T_bN) : Mjoin(PATL,geput1T_bX);
      }
   #endif
   nfnblks = op.nfnblks;
   nnblks = nfnblks + op.npnblks;
   if (nfnblks)
      mb = nb = op.nb;
   else
   {
      mb = op.pmb;
      nb = op.pnb;
   }
/*
 * This algorithm requires allocating all of A & B + 2 blks for C
 */
   szC = op.szC;
   szA = nfnblks*op.szA + op.npnblks*op.pszA;
   szB = nfnblks*op.szB + op.npnblks*op.pszB;
   sz = ATL_MulBySize(szA+szB+szC+szC) + 4*ATL_Cachelen;
   if (sz <= ATL_MaxMalloc)
      vp = malloc(sz);
   if (!vp)
      return(1);

   bC = ATL_AlignPtr(vp);
   aC = bC + (szC SHIFT);
   aC = ATL_AlignPtr(aC);
   #ifdef TCPLX
     aCr = aC + szC;
     aA = aCr + szC;
   #else
      aA = aC + szC;
   #endif
   aA = ATL_AlignPtr(aA);
   aB = aA + (szA SHIFT);
   aB = ATL_AlignPtr(aB);
   szB = (nfnblks) ? (op.szB SHIFT) : 0;  /* inc between B blks */
   szA = (nfnblks) ? (op.szA SHIFT) : 0;

   if (Uplo == AtlasLower)
   {
      ATL_SZT j;
      ATL_CUINT NB = op.nb;
      TYPE *wAj=aA, *wBj=aB;

      for (j=0; j < nnblks; j++, C += nb*((ldc+1)SHIFT), wAj += szA, wBj += szB)
      {
         TYPE *wAi=wAj+szA, *wBi=wBj+szB;
         TYPE *pC = C;
         ATL_SZT i;
         ATL_UINT nb, n;

         if (j < nfnblks)
            nb = n = NB;
         else
         {
            nb = op.pnb;
            n = op.nF;
         }
/*
 *       First do diagonal block
 */
         Mjoin(PATL,opblk)(&op, j, j, A, B, NULL, wAj, wAj,
                           wBj, wBj+szB, aCr, aC);
         #ifdef TCPLX
            blk2C_b0(n, n, ONE, aCr, aC, ZERO, bC, n);
         #else
            blk2C_b0(n, n, ONE, aC, ZERO, bC, n);
         #endif
         syput(n, bC, beta, pC, ldc);
         pC += n SHIFT;
         for (i=j+1; i < nnblks; i++, wAi += szA, wBi += szB)
         {
            ATL_CSZT in = (i < nnblks-1) ? i+1 : 0;
            ATL_UINT mb, m;
            if (i < nfnblks)
               mb = m = NB;
            else
            {
               mb = op.pmb;
               m = op.mF;
            }
            Mjoin(PATL,opblk)(&op, j, i, NULL, B, NULL, wAj, wAi,
                              wBi, wBj, aCr, aC);
            #ifdef TCPLX
               blk2C_b0(n, m, ONE, aCr, aC, ZERO, bC, n);
            #else
               blk2C_b0(n, m, ONE, aC, ZERO, bC, n);
            #endif
            if (m == n)
            {
               Mjoin(PATL,sqtrans)(n, bC, n);
               #ifdef Conj_  /* need hermitian transpose for HER2K! */
                  Mjoin(PATLU,scal)(n*n, ATL_rnone, bC+1, 2);
               #endif
            }
            else
               geputT(m, n, bC, n, beta, pC, ldc);
            Mjoin(PATL,opblk)(&op, i, j, A, NULL, NULL, wAi, wAj+szA,
                              wBj, wBi+szB, aCr, aC);
            if (m == n)
            {
               #ifdef TCPLX
                  blk2C_b1(n, n, ONE, aCr, aC, ONE, bC, n);
               #else
                  blk2C_b1(n, n, ONE, aC, ONE, bC, n);
               #endif
               Mjoin(PATL,geadd)(n, n, ONE, bC, n, beta, pC, ldc);
            }
            else
               #ifdef TCPLX
                  blk2C_b1(m, n, ONE, aCr, aC, ONE, pC, ldc);
               #else
                  blk2C_b1(m, n, ONE, aC, ONE, pC, ldc);
               #endif
            pC += m SHIFT;
         }
         A = NULL;
         B = NULL;
      }
   }
   else
   {
      ATL_SZT j;
      ATL_CUINT NB = Mmin(mb,nb);
      TYPE *wAj=aA, *wBj=aB;

      for (j=0; j < nnblks; j++, C += nb*(ldc SHIFT), wAj += szA, wBj += szB)
      {
         const TYPE *aa=A, *bb=B;
         TYPE *wAi=aA, *wBi=aB;
         TYPE *pC = C;
         ATL_SZT i;
         ATL_UINT n;

         if (j < nfnblks)
            nb = n = NB;
         else
         {
            nb = op.pnb;
            n = op.nF;
         }
         for (i=0; i < j; i++, wAi += szA, wBi += szB)
         {
            ATL_UINT m;
            if (i < nfnblks)
               mb = m = NB;
            else
            {
               mb = op.pmb;
               m = op.mF;
            }
            Mjoin(PATL,opblk)(&op, j, i, aa, NULL, NULL, wAj, wAi,
                              wBi, wBj, aCr, aC);
            aa = NULL;
            #ifdef TCPLX
               blk2C_b0(n, m, ONE, aCr, aC, ZERO, bC, n);
            #else
               blk2C_b0(n, m, ONE, aC, ZERO, bC, n);
            #endif
            if (m == n)
            {
               Mjoin(PATL,sqtrans)(n, bC, n);
               #ifdef Conj_  /* need hermitian transpose for HER2K! */
                  Mjoin(PATLU,scal)(n*n, ATL_rnone, bC+1, 2);
               #endif
            }
            else
               geputT(m, n, bC, n, beta, pC, ldc);
            Mjoin(PATL,opblk)(&op, i, j, NULL, bb, NULL, wAi, wAj+szA,
                              wBj, wBi+szB, aCr, aC);
            bb = NULL;
            if (m == n)
            {
               #ifdef TCPLX
                  blk2C_b1(n, n, ONE, aCr, aC, ONE, bC, n);
               #else
                  blk2C_b1(n, n, ONE, aC, ONE, bC, n);
               #endif
               Mjoin(PATL,geadd)(n, n, ONE, bC, n, beta, pC, ldc);
            }
            else
               #ifdef TCPLX
                  blk2C_b1(m, n, ONE, aCr, aC, ONE, pC, ldc);
               #else
                  blk2C_b1(m, n, ONE, aC, ONE, pC, ldc);
               #endif
            pC += m SHIFT;
         }
/*
 *       Do diagonal block to copy new A block
 */
         Mjoin(PATL,opblk)(&op, j, j, aa, bb, NULL, wAj, wAj,
                           wBj, wBj+szB, aCr, aC);
         #ifdef TCPLX
            blk2C_b0(n, n, ONE, aCr, aC, ZERO, bC, n);
         #else
            blk2C_b0(n, n, ONE, aC, ZERO, bC, n);
         #endif
         syput(n, bC, beta, pC, ldc);
      }
   }
   free(vp);
   return(0);
}
