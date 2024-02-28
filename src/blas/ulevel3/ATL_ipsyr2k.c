/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_level1.h"
#include "atlas_kernel3.h"
#include "atlas_reflevel3.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))

#ifdef Conj_
int Mjoin(PATL,ipher2k)
#else
int Mjoin(PATL,ipsyr2k)
#endif
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
   ATL_SZT sz, szA, szB, nfnblks, nnblks, nkblks;
   unsigned int szC, MVS=0;
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
   ipinfo_t ip;
   #ifdef Conj_
      const enum ATLAS_TRANS TB = (TA == AtlasConjTrans) ?
         AtlasNoTrans:AtlasConjTrans;
   #else
      const enum ATLAS_TRANS TB = (TA == AtlasTrans) ?
         AtlasNoTrans:AtlasTrans;
   #endif
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
   Mjoin(PATL,ipmenInfo)(&ip, TA, TB, N, N, K, lda, ldb, ldc, alpha, ZERO);
   blk2C_b0 = ip.blk2c;
   blk2C_b1 = ip.blk2c_b1;
   #ifdef TCPLX
      ONE = ip.ONE;
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
   nfnblks = ip.nfnblks;
   nnblks = nfnblks + ip.npnblks;
   nkblks = ip.nfkblks + 1;
/*
 * This algorithm requires allocating all of A & B + 2 blks for C
 */
   szC = ip.szC;
   if (nnblks > 1)
   {
      szA = (nfnblks*ip.szA + ip.npnblks*ip.pszA)*nkblks;
      szB = (nfnblks*ip.szB + ip.npnblks*ip.pszB)*nkblks;
      MVS=3;
   }
   else if (nfnblks)
   {
     szA = ip.szA;
     szB = ip.szB;
   }
   else
   {
     szA = ip.pszA;
     szB = ip.pszB;
   }
   sz = ATL_MulBySize(szA+szB+szC+szC) + 3*ATL_Cachelen;
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
/*
 * Now, set sz[A,B] for incr workspace ptrs
 */
   if (nnblks > 1)
   {
      szA = ip.szA;
      szB = ip.szB;
      szA *= nkblks SHIFT;
      szB *= nkblks SHIFT;
   }
   else   /* just reusing same A/B space for 1 blk C! */
      szA = szB = 0;

   if (Uplo == AtlasLower)
   {
      ATL_iptr_t imsk = ~(0L);
      ATL_SZT j;
      ATL_CUINT NB = (nfnblks) ? ip.nb : ip.pnb;
      TYPE *wAj=aA, *wBj=aB;

      for (j=0; j < nnblks; j++, wAj += szA, wBj += szB)
      {
         TYPE *wAi=wAj+szA, *wBi=wBj+szB;
         TYPE *pC = C;
         ATL_SZT incA, incB, i;
         ATL_UINT n;

         if (j < nfnblks)
         {
            n = NB;
            incA = ip.incAm & imsk;
            incB = ip.incBn & imsk;
         }
         else
         {
            n = ip.nF;
            incA = ip.pincAm & imsk;
            incB = ip.pincBn & imsk;
         }
/*
 *       First do diagonal block
 */
         ip.ldc = n;
         Mjoin(PATL,iploopsK)(&ip, j, j, A, B, bC, MVS, wAj, wBj, aCr, aC,
                              ZERO, blk2C_b0);
         syput(n, bC, beta, pC, ldc);
         pC += n SHIFT;
         A += incA;
         B += incB;
         for (i=j+1; i < nnblks; i++, wAi += szA, wBi += szB)
         {
            ATL_UINT m;
            if (i < nfnblks)
            {
               m = NB;
               incA = ip.incAm & imsk;
               incB = ip.incBn & imsk;
            }
            else
            {
               m = ip.mF;
               incA = ip.pincAm & imsk;
               incB = ip.incBn & imsk;
            }
            ip.ldc = n;
            Mjoin(PATL,iploopsK)(&ip, j, i, NULL, B, bC, 3, wAj, wBi, aCr, aC,
                                 ZERO, blk2C_b0);
            if (m == n)
            {
               Mjoin(PATL,sqtrans)(n, bC, n);
               #ifdef Conj_  /* need hermitian transpose for HER2K! */
                  Mjoin(PATLU,scal)(n*n, ATL_rnone, bC+1, 2);
               #endif
            }
            else
               geputT(m, n, bC, n, beta, pC, ldc);
            if (m == n)
            {
               Mjoin(PATL,iploopsK)(&ip, i, j, A, NULL, bC, 3, wAi, wBj,
                                    aCr, aC, ONE, blk2C_b1);
               Mjoin(PATL,geadd)(n, n, ONE, bC, n, beta, pC, ldc);
            }
            else
            {
               ip.ldc = ldc;
               Mjoin(PATL,iploopsK)(&ip, i, j, A, NULL, pC, 3, wAi, wBj,
                                    aCr, aC, ONE, blk2C_b1);
            }
            pC += m SHIFT;
            A += incA;
            B += incB;
         }
         A = NULL;
         B = NULL;
         imsk = 0L;
         C += NB*((ldc+1)SHIFT);
      }
   }
   else   /* C is symmetric matrix stored in Upper format */
   {
      ATL_CSZT incB = (nfnblks) ? ip.incBn : ip.pincBn;
      ATL_CSZT incA = (nfnblks) ? ip.incAm : ip.pincAm;
      ATL_SZT j;
      ATL_CUINT NB = (nfnblks) ? ip.nb : ip.pnb;
      TYPE *wAj=aA, *wBj=aB;

      for (j=0; j < nnblks; j++, C += NB*(ldc SHIFT), wAj += szA, wBj += szB)
      {
         const TYPE *aa=A, *bb=B;
         TYPE *wAi=aA, *wBi=aB;
         TYPE *pC = C;
         ATL_SZT i, incA;
         ATL_UINT n;

         if (j < nfnblks)
         {
            n = NB;
            incA = ip.incAm;
         }
         else
         {
            n = ip.nF;
            incA = ip.pincAm;
         }
         for (i=0; i < j; i++, wAi += szA, wBi += szB)
         {
            ATL_UINT m;
            ip.ldc = n;
            Mjoin(PATL,iploopsK)(&ip, j, i, aa, NULL, bC, 3, wAj, wBi, aCr, aC,
                                 ZERO, blk2C_b0);
            aa = NULL;
            if (i < nfnblks)
               m = NB;
            else
               m = ip.mF;
            if (m == n)
            {
               Mjoin(PATL,sqtrans)(n, bC, n);
               #ifdef Conj_  /* need hermitian transpose for HER2K! */
                  Mjoin(PATLU,scal)(n*n, ATL_rnone, bC+1, 2);
               #endif
            }
            else
               geputT(m, n, bC, n, beta, pC, ldc);

            if (m == n)
            {
               Mjoin(PATL,iploopsK)(&ip, i, j, NULL, bb, bC, 3, wAi, wBj,
                                    aCr, aC, ONE, blk2C_b1);
               Mjoin(PATL,geadd)(n, n, ONE, bC, n, beta, pC, ldc);
            }
            else
            {
               ip.ldc = ldc;
               Mjoin(PATL,iploopsK)(&ip, i, j, NULL, bb, pC, 3, wAi, wBj,
                                    aCr, aC, ONE, blk2C_b1);
            }
            bb = NULL;
            pC += m SHIFT;
         }
/*
 *       Do diagonal block to copy new A block
 */
         ip.ldc = n;
         Mjoin(PATL,iploopsK)(&ip, j, j, aa, bb, bC, MVS, wAj, wBj, aCr, aC,
                              ZERO, blk2C_b0);
         syput(n, bC, beta, pC, ldc);
         B += incB;
         A += incA;
      }
   }
   free(vp);
   return(0);
}
