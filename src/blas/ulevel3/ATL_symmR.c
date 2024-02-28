#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_level1.h"
#include "atlas_level2.h"
#include "atlas_level3.h"
#include "atlas_reflevel3.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
/*
 * On the Right, symmetrix matrix is GEMM's B, with K=N
 */
#ifdef Conj_
int Mjoin(PATL,ophemmR)
#else
int Mjoin(PATL,opsymmR)
#endif
(
   const enum ATLAS_UPLO  Uplo,
   ATL_CSZT  M,
   ATL_CSZT N,
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
   ATL_SZT sz, szW;
   void *vp=NULL;
   TYPE *aS, *aB, *aC, *bS;
   #ifdef TCPLX
      TYPE *rS, *rC;
   #endif
   opinfo_t oi;
   int Mjoin(PATL,opsymmInfo)
      (opinfo_t *op, ATL_UINT bv,
       ATL_CSZT M, ATL_CSZT N, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
       const SCALAR alpha, const SCALAR beta);

   if (!M || !N)
      return(0);
   if (N == 1) /* symmetric matrix is dense scalar! */
   {
      #ifdef TCPLX
         const register TYPE ra=(*alpha), ia=alpha[1], rA = *A;
         #ifdef Conj_
            const TYPE X[2] = {ra*rA, ia*rA};
         #else
            const register TYPE  iA=A[1];
            TYPE X[2];
            *X = ra*rA - ia*iA;
            X[1] = ra*iA + ia*rA;
         #endif
         Mjoin(PATL,axpby)(M, X, B, 1, beta, C, 1);
      #else
         Mjoin(PATL,axpby)(M, alpha*(*A), B, 1, beta, C, 1);
      #endif
      return(0);
   }
   #ifdef Conj_
   if (Mjoin(PATL,opsymmInfo)(&oi, (Uplo == AtlasUpper)?6:4, M, N, lda, ldb,
                              ldc, alpha, beta))
   #else
   if (Mjoin(PATL,opsymmInfo)(&oi, (Uplo == AtlasUpper)?2:0, M, N, lda, ldb,
                              ldc, alpha, beta))
   #endif
      return(1);

   szW = oi.szB;
   szW = Mmax(szW, oi.szA);
   sz = oi.szB + szW + oi.szC;
   sz = ATL_MulBySize(sz) + 3*ATL_Cachelen + oi.exsz;
   vp = malloc(sz);
   ATL_assert(vp);
   aS = ATL_AlignPtr(vp);
   #ifdef TCPLX
      rS = aS + oi.szB;
      aB = rS + oi.szB;
      aB = ATL_AlignPtr(aB);
      aC = aB + (szW<<1);
      aC = ATL_AlignPtr(aC);
      rC = aC + oi.szC;
   #else
      aB = aS + oi.szB;
      aB = ATL_AlignPtr(aB);
      aC = aB + szW;
      aC = ATL_AlignPtr(aC);
   #endif
   if (Uplo == AtlasUpper)
   {
      #ifdef TCPLX
         #ifdef Conj_
            #define SYCPY Mjoin(PATL,hecpyUNB)
         #else
            #define SYCPY Mjoin(PATL,sycpyUNB)
         #endif
         if (alpha[1] == ATL_rzero)
         {
            const register TYPE ral=(*alpha);
            if (ral == ATL_rone)
               Mjoin(SYCPY,_a1)(N, alpha, A, lda, aB, N);
            else if (ral == ATL_rnone)
               Mjoin(SYCPY,_an)(N, alpha, A, lda, aB, N);
            else
               Mjoin(SYCPY,_ar)(N, alpha, A, lda, aB, N);
         }
         else
         #ifdef Conj_
            Mjoin(SYCPY,_a1)(N, oi.ONE, A, lda, aB, N);
         #else
            Mjoin(SYCPY,_aX)(N, alpha, A, lda, aB, N);
         #endif
         #undef SYCPY
      #else
         if (alpha == ATL_rone)
            Mjoin(PATL,sycpyUNB_a1)(N, alpha, A, lda, aB, N);
         else if (alpha == ATL_rnone)
            Mjoin(PATL,sycpyUNB_an)(N, alpha, A, lda, aB, N);
         else
            Mjoin(PATL,sycpyUNB_aX)(N, alpha, A, lda, aB, N);
      #endif
      #if 0
         Mjoin(PATL,geprint)("aG", N, N, A, lda);
         Mjoin(PATL,geprint)("aS", N, N, aB, N);
      #endif
   }
   else /* Uplo == AtlasLower */
   {
      #ifdef TCPLX
         #ifdef Conj_
            #define SYCPY Mjoin(PATL,hecpyLNB)
         #else
            #define SYCPY Mjoin(PATL,sycpyLNB)
         #endif
         if (alpha[1] == ATL_rzero)
         {
            const register TYPE ral=(*alpha);
            if (ral == ATL_rone)
               Mjoin(SYCPY,_a1)(N, alpha, A, lda, aB, N);
            else if (ral == ATL_rnone)
               Mjoin(SYCPY,_an)(N, alpha, A, lda, aB, N);
            else
               Mjoin(SYCPY,_ar)(N, alpha, A, lda, aB, N);
         }
         else
         #ifdef Conj_
            Mjoin(SYCPY,_a1)(N, oi.ONE, A, lda, aB, N);
         #else
            Mjoin(SYCPY,_aX)(N, alpha, A, lda, aB, N);
         #endif
         #undef SYCPY
      #else
         if (alpha == ATL_rone)
            Mjoin(PATL,sycpyLNB_a1)(N, alpha, A, lda, aB, N);
         else if (alpha == ATL_rnone)
            Mjoin(PATL,sycpyLNB_an)(N, alpha, A, lda, aB, N);
         else
            Mjoin(PATL,sycpyLNB_aX)(N, alpha, A, lda, aB, N);
      #endif
      #if 0
         Mjoin(PATL,geprint)("aG", N, N, A, lda);
         Mjoin(PATL,geprint)("aS", N, N, aB, N);
      #endif
   }
   #ifdef TCPLX
      oi.b2blk(N, N, oi.ONE, aB, N, rS, aS);
      Mjoin(PATL,oploopsM)(&oi, 0, 0, B, NULL, C, 0, aB, aS, rC, aC);
   #else
      oi.b2blk(N, N, ATL_rone, aB, N, aS);
      Mjoin(PATL,oploopsM)(&oi, 0, 0, B, NULL, C, 0, aB, aS, aC, aC);
   #endif
   free(vp);
   return(0);
}

#ifndef Conj_
   #ifndef TCPLX
      #define ONE ATL_rone
      #define TA AtlasTrans
   #endif
void Mjoin(PATL,symmR_rec)
(
   ATL_CUINT flg, /* bitvec 0:Upper? 1:HEMM? */
   ATL_CSZT  M,
   ATL_CSZT N,
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
/*
 * If the symmetric matrix is smaller than our rank-K GEMM support, use
 * outer product gemm to solve this problem.
 */
   if (N <= ATL_VWopgen_LAST_KB)
   {
      const enum ATLAS_UPLO Uplo = (flg&1) ? AtlasUpper:AtlasLower;
      #ifdef TCPLX
         if (flg&2) /* really a HEMM */
         {
            if (!Mjoin(PATL,ophemmR)(Uplo, M, N, alpha, A, lda, B, ldb,
                                     beta, C, ldc))
               return;
         }
         else
            if (!Mjoin(PATL,opsymmR)(Uplo, M, N, alpha, A, lda, B, ldb,
                                     beta, C, ldc))
               return;
      #else
         if (!Mjoin(PATL,opsymmR)(Uplo, M, N, alpha, A, lda, B, ldb,
                                  beta, C, ldc))
            return;
      #endif
   }
/*
 * To avoid recurring and redundantly copying matrices as we average various
 * size GEMM's performance, attempt to halt recursion on symmetric dims by
 * allocating enough space to copy entire symmetric matrix to dense, and then
 * call GEMM only one time.  Do not do this if the general matrix consists
 * only of a single K-panel: in this case better to recur down and use opsymm.
 * Reason is copy of A a big deal, and copying up front brings it through cache
 * twice, wheareas recursion will provide automatic blocking, and avoid this
 * for small problems where it matters.
 */
   else if (M > ATL_VIEW_BEST_MB)
   {
      if (!Mjoin(PATL,ipsymmR)(flg, M, N, alpha,A,lda,B,ldb,beta,C,ldc))
         return;
   }
/*
 * Otherwise, recur on symmetric matrix until we can stop one of above using:
 * C0 = B0 * A00  + beta*C0  (symm)
 * C0 = B1 * A10 + C0        (gemm)
 * C1 = B0 * A10^T + beta*C1 (gemm)
 * C1 = B1 * A11 + C1        (symm)
 */
   {
      #ifdef TCPLX
         const TYPE ONE[2] = {ATL_rone, ATL_rzero};
         const enum ATLAS_TRANS TA = (flg&2) ? AtlasConjTrans : AtlasTrans;
      #endif
      const int NL=(N>>1), NR = N-NL;
      const TYPE *B1 = B + ldb*(NL SHIFT);
      TYPE *C1 = C + ldc*(NL SHIFT);

      ATL_assert(NL);  /* debugging, remove later */
      Mjoin(PATL,symmR_rec)(flg, M, NL, alpha, A, lda, B, ldb,
                            beta, C, ldc);
      if (flg&1)  /* Upper matrix */
      {
         const TYPE *A01 = A + NL*(lda SHIFT);
         Mjoin(PATL,ammm)(AtlasNoTrans, TA, M, NL, NR, alpha,
                          B1, ldb, A01, lda, ONE, C, ldc);
         Mjoin(PATL,ammm)(AtlasNoTrans, AtlasNoTrans, M, NR, NL, alpha,
                          B, ldb, A01, lda, beta, C1, ldc);
      }
      else         /* Lower matrix */
      {
         const TYPE *A10 = A + (NL SHIFT);
         Mjoin(PATL,ammm)(AtlasNoTrans, AtlasNoTrans, M, NL, NR, alpha,
                          B1, ldb, A10, lda, ONE, C, ldc);
         Mjoin(PATL,ammm)(AtlasNoTrans, TA, M, NR, NL, alpha,
                          B, ldb, A10, lda, beta, C1, ldc);
      }
      Mjoin(PATL,symmR_rec)(flg, M, NR, alpha, A+(NL SHIFT)*(lda+1), lda,
                            B1, ldb, ONE, C1, ldc);
   }
}
   #ifdef TCPLX
   #else
      #undef TA
   #endif
#else
   void Mjoin(PATL,symmL_rec)
      (ATL_CUINT, ATL_CSZT, ATL_CSZT, const SCALAR, const TYPE*, ATL_CSZT,
       const TYPE*, ATL_CSZT, const SCALAR, TYPE*, ATL_CSZT);
#endif
