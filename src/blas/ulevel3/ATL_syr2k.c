#include "atlas_amm.h"
#include "atlas_level1.h"
#include "atlas_level3.h"
#if 0
   #define USEREF 1
   #include "atlas_reflevel3.h"
#endif
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#ifdef Conj_
   #define SYR2K Mjoin(PATL,her2k)
   #define IFLG 4
#else
   #define SYR2K Mjoin(PATL,syr2k)
   #define IFLG 0
#endif
#ifdef Conj_
int Mjoin(PATL,opher2k_wrap)
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
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   return(Mjoin(PATL,opher2k)(Uplo, TA, N, K, alpha, A, lda, B, ldb,
                              *beta, C, ldc));
}
#else
typedef int (*syr2k_t)(const enum ATLAS_UPLO, const enum ATLAS_TRANS,
                       ATL_CSZT, ATL_CSZT, const SCALAR,
                       const TYPE*, ATL_CSZT, const TYPE*, ATL_CSZT,
                       const SCALAR, TYPE*, ATL_CSZT);
void Mjoin(PATL,syr2k_rec)
(
   ATL_CUINT bv, /* 0:Upper?, 1:Transpose?, 2:HER2K? */
   ATL_CSZT  N,
   ATL_CSZT K,  /* if NoTrans, ncols of A&B, else nrows A&B */
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   const TYPE *B,
   ATL_CSZT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc,
   void *syr2kOP
)
{
   const enum ATLAS_UPLO UL=(bv&1) ? AtlasUpper:AtlasLower;
   #ifdef TCPLX
      const TYPE ONE[2] = {ATL_rone, ATL_rzero};
      enum ATLAS_TRANS TA;
      if (bv&4)
         TA = (bv&2) ? AtlasConjTrans : AtlasNoTrans;
      else
         TA = (bv&2) ? AtlasTrans : AtlasNoTrans;
   #else
      const enum ATLAS_TRANS TA = (bv&2) ? AtlasTrans : AtlasNoTrans;
      #define ONE ATL_rone
   #endif
   if (N == 1)  /* can use dot product! */
   {
      ATL_SZT incA, incB;
      if (bv&2) /* K is rows */
      {
         incA = 1;
         incB = 1;
      }
      else
      {
         incA = lda;
         incB = ldb;
      }
      #ifdef TCPLX
         if (bv&4) /* dotc */
         {
            TYPE dot[2];
            register TYPE rd;
            const register TYPE ra=(*alpha), ia=alpha[1], rb=(*beta);
            if (incA == 1 && incB == 1 && ia == ATL_rzero)
            {
               rd = Mjoin(PATLU,dot)(K+K, A, incA, B, incB);
               rd *= ra;
            }
            else
            {
               register TYPE id;
               Mjoin(PATL,dotc_sub)(K, A, incA, B, incB, dot);
               rd = *dot;
               id = dot[1];
               if (bv&2)
                  rd = rd*ra - id*ia;
               else
                  rd = rd*ra + id*ia;
            }
            rd += rd;
            if (rb == ATL_rzero)
               *C = rd;
            else
               *C = *C * rb + rd;
            C[1] = ATL_rzero;
         }
         else /* dotu */
         {
            TYPE dot[2];
            register TYPE rc, ic, rd, id;
            const register TYPE ra=(*alpha),ia=alpha[1], rb=(*beta),ib=beta[1];
            Mjoin(PATL,dotu_sub)(K, A, incA, B, incB, dot);
            rd = *dot;
            id = dot[1];
            rc = rd*ra - id*ia;  /* multiply by alpha */
            ic = rd*ia + id*ra;
            rc += rc;            /* double it */
            ic += ic;
            if (rb != ATL_rzero || ib != ATL_rzero)  /* apply beta to C */
            {                                        /* and add into answer */
               rd = *C;
               id = C[1];
               rc += rb*rd - ib*id;
               ic += rb*id + ib*rd;
            }
            *C = rc;
            C[1] = ic;
         }
      #else
         TYPE dot;
         dot = Mjoin(PATL,dot)(K, A, incA, B, incB);
         dot *= alpha;
         dot += dot;
         if (beta == ATL_rzero)
            *C = dot;
         else
            *C = *C * beta + dot;
      #endif
       return;
   }
/*
 * If K small enough, need the outer-product version of syrk
 */
   if (K <= ATL_VWopgen_MAX_KB)
   {
      syr2k_t opsyr2k = syr2kOP;
      if (!opsyr2k(UL, TA, N, K, alpha, A, lda, B, ldb, beta, C, ldc))
         return;
   }
/*
 * Try the inner product version for any case outer-product fails on
 */
   {
      #ifdef TCPLX
      if (bv&4)
      {
         if (!Mjoin(PATL,ipher2k)(UL, TA, N, K, alpha, A, lda, B, ldb,
                                  *beta, C, ldc));
            return;
      }
      else
      #endif
      {
         if (!Mjoin(PATL,ipsyr2k)(UL, TA, N, K, alpha, A, lda, B, ldb,
                                  beta, C, ldc));
            return;
      }
   }
/*
 *    Splitting K or N reduces mem req by same amount, so split largest.
 *    When we split along N, we wind up with two syr2k calls, and two calls
 *    to GEMM.
 */
      if (N >= K)
      {
         ATL_SZT NL=N>>1, NR=N-NL;
         const TYPE *an, *bn;
         enum ATLAS_TRANS TA, TB;
         if (bv&2) /* A & B are KxN */
         {
            an = A + (NL SHIFT)*lda;
            bn = B + (NL SHIFT)*ldb;
            #ifdef TCPLX
               TA = (bv&4) ? AtlasConjTrans:AtlasTrans;
            #else
               TA = AtlasTrans;
            #endif
               TB = AtlasNoTrans;
         }
         else      /* A & B are NxK */
         {
            TA = AtlasNoTrans;
            #ifdef TCPLX
               TB = (bv&4) ? AtlasConjTrans:AtlasTrans;
            #else
               TB = AtlasTrans;
            #endif
            an = A + (NL SHIFT);
            bn = B + (NL SHIFT);
         }
         Mjoin(PATL,syr2k_rec)(bv, NL, K, alpha, A, lda, B, ldb,
                               beta, C, ldc, syr2kOP);
         Mjoin(PATL,syr2k_rec)(bv, NR, K, alpha, an, lda, bn, ldb,
                               beta, C+((ldc+1)SHIFT)*NL, ldc, syr2kOP);
         if (bv&1) /* Upper */
         {
            C += ldc*(NL SHIFT);
            Mjoin(PATL,ammm)(TA, TB, NL, NR, K, alpha, A, lda, bn, ldb,
                             beta, C, ldc);
            Mjoin(PATL,ammm)(TA, TB, NL, NR, K, alpha, B, ldb, an, lda,
                             ONE, C, ldc);
         }
         else
         {
            C += (NL SHIFT);
            Mjoin(PATL,ammm)(TA, TB, NR, NL, K, alpha, an, lda, B, ldb,
                             beta, C, ldc);
            Mjoin(PATL,ammm)(TA, TB, NR, NL, K, alpha, bn, ldb, A, lda,
                             ONE, C, ldc);
         }
         return;
      }
/*
 * If we reach here, split K, resulting in 2 syr2k, 2nd with BETA=ONE
 */
   {
      ATL_SZT KL=(K>>1), KR = K - KL;
      Mjoin(PATL,syr2k_rec)(bv, N, KL, alpha, A, lda, B, ldb,
                            beta, C, ldc, syr2kOP);
      #ifdef TCPLX
         KL += KL;
      #endif
      if (bv&2)  /* K controls rows of A&B */
      {
         A += KL;
         B += KL;
      }
      else       /* K controls cols of A&B */
      {
         A += KL * lda;
         B += KL * ldb;
      }
      Mjoin(PATL,syr2k_rec)(bv, N, KR, alpha, A, lda, B, ldb, ONE, C, ldc,
                            syr2kOP);
   }
}
   #ifndef TCPLX
      #undef ONE
   #endif
#endif
void SYR2K
(
   const enum ATLAS_UPLO Uplo,
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
   #else
      const SCALAR beta,
   #endif
   TYPE *C,
   ATL_CSZT ldc
)
{
   void *vp;
   #ifdef Conj_
      const TYPE beta[2] = {rbeta, ATL_rzero};
   #endif
   ATL_UINT bv;

   #ifdef USEREF
      #ifdef Conj_
         Mjoin(PATL,refher2k)(Uplo, TA, N, K, alpha, A,lda, B,ldb,rbeta,C,ldc);
      #else
         Mjoin(PATL,refsyr2k)(Uplo, TA, N, K, alpha, A,lda, B,ldb, beta,C,ldc);
      #endif
      return;
   #endif
   if (!N)
      return;
   if (SCALAR_IS_ZERO(alpha) || !K)
   {
      #ifdef Conj_
         if (rbeta == ATL_rone)
            return;
      #else
         if (SCALAR_IS_ONE(beta))
            return;
      #endif
      Mjoin(PATL,trscal)(Uplo, N, N, beta, C, ldc);
      #ifdef Conj_
         Mjoin(PATLU,zero)(N, C+1, (ldc+1)<<1);
      #endif
      return;
   }
/*
 * Outer-product syr2k handles degenerate K, along with any 1-block K
 */
   if (K < 3) /* ATL_VWopgen_MAX_KB) */
   #ifdef Conj_
      if (!Mjoin(PATL,opher2k)(Uplo, TA, N, K, alpha, A, lda, B, ldb,
                               rbeta, C, ldc))
         return;
      vp = Mjoin(PATL,opher2k_wrap);
      bv = 4;
      if (TA == AtlasConjTrans)
         bv |= 2;
   #else
      if (!Mjoin(PATL,opsyr2k)(Uplo, TA, N, K, alpha, A, lda, B, ldb,
                               beta, C, ldc))
         return;
      vp = Mjoin(PATL,opsyr2k);
      bv = 0;
      #ifdef TCPLX
         if (TA == AtlasTrans)
      #else
         if (TA == AtlasTrans || TA == AtlasConjTrans)
      #endif
            bv |= 2;
   #endif
   bv |= (Uplo == AtlasUpper);
   Mjoin(PATL,syr2k_rec)(bv, N, K, alpha, A, lda, B, ldb, beta, C, ldc, vp);
}
