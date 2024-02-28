/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_level1.h"
#include "atlas_level2.h"
#include "atlas_level3.h"
#include "atlas_reflevel3.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#if 0
   #define USEREF 1
   #include "atlas_reflevel3.h"
#endif

/*
 * On the left, symmetrix matrix is GEMM's A, with K=M
 */
#ifdef Conj_
int Mjoin(PATL,ophemmL)
#else
int Mjoin(PATL,opsymmL)
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
   TYPE *aA, *aB, *aC, *bA;
   #ifdef TCPLX
      TYPE *rA, *rC;
   #endif
   opinfo_t oi;
   int Mjoin(PATL,opsymmInfo)
      (opinfo_t *op, ATL_UINT bv,
       ATL_CSZT M, ATL_CSZT N, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
        const SCALAR alpha, const SCALAR beta);

   if (!M || !N)
      return(0);
   if (M == 1) /* symmetric matrix is dense scalar! */
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
         Mjoin(PATL,axpby)(N, X, B, ldb, beta, C, ldc);
      #else
         Mjoin(PATL,axpby)(N, alpha*(*A), B, ldb, beta, C, ldc);
      #endif
      return(0);
   }
   #ifdef Conj_
   if (Mjoin(PATL,opsymmInfo)(&oi, (Uplo == AtlasUpper)?7:5, M, N, lda, ldb,
                              ldc, alpha, beta))
   #else
   if (Mjoin(PATL,opsymmInfo)(&oi, (Uplo == AtlasUpper)?3:1, M, N, lda, ldb,
                              ldc, alpha, beta))
   #endif
      return(1);

   szW = oi.szB;
   szW = Mmax(szW, oi.szA);
   sz = oi.szA + szW + oi.szC;
   sz = ATL_MulBySize(sz) + ATL_Cachelen + oi.exsz;
   vp = malloc(sz);
   ATL_assert(vp);
   aA = ATL_AlignPtr(vp);
   #ifdef TCPLX
      rA = aA + oi.szA;
      aB = rA + oi.szA;
      aC = aB + (szW<<1);
      rC = aC + oi.szC;
   #else
      aB = aA + oi.szA;
      aC = aB + szW;
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
               Mjoin(SYCPY,_a1)(M, alpha, A, lda, aB, M);
            else if (ral == ATL_rnone)
               Mjoin(SYCPY,_an)(M, alpha, A, lda, aB, M);
            else
               Mjoin(SYCPY,_ar)(M, alpha, A, lda, aB, M);
         }
         else
         #ifdef Conj_
            Mjoin(SYCPY,_a1)(M, oi.ONE, A, lda, aB, M);
         #else
            Mjoin(SYCPY,_aX)(M, alpha, A, lda, aB, M);
         #endif
         #if 0
            Mjoin(PATL,geprint)("aG", M, M, A, lda);
            Mjoin(PATL,geprint)("aS", M, M, aB, M);
         #endif
         #undef SYCPY
      #else
         if (alpha == ATL_rone)
            Mjoin(PATL,sycpyUNB_a1)(M, alpha, A, lda, aB, M);
         else if (alpha == ATL_rnone)
            Mjoin(PATL,sycpyUNB_an)(M, alpha, A, lda, aB, M);
         else
            Mjoin(PATL,sycpyUNB_aX)(M, alpha, A, lda, aB, M);
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
               Mjoin(SYCPY,_a1)(M, alpha, A, lda, aB, M);
            else if (ral == ATL_rnone)
               Mjoin(SYCPY,_an)(M, alpha, A, lda, aB, M);
            else
               Mjoin(SYCPY,_ar)(M, alpha, A, lda, aB, M);
         }
         else
         #ifdef Conj_
            Mjoin(SYCPY,_a1)(M, oi.ONE, A, lda, aB, M);
         #else
            Mjoin(SYCPY,_aX)(M, alpha, A, lda, aB, M);
         #endif
         #if 0
            Mjoin(PATL,geprint)("aG", M, M, A, lda);
            Mjoin(PATL,geprint)("aS", M, M, aB, M);
         #endif
         #undef SYCPY
      #else
         if (alpha == ATL_rone)
            Mjoin(PATL,sycpyLNB_a1)(M, alpha, A, lda, aB, M);
         else if (alpha == ATL_rnone)
            Mjoin(PATL,sycpyLNB_an)(M, alpha, A, lda, aB, M);
         else
            Mjoin(PATL,sycpyLNB_aX)(M, alpha, A, lda, aB, M);
      #endif
   }
   #ifdef TCPLX
      oi.a2blk(M, M, oi.ONE, aB, M, rA, aA);
      Mjoin(PATL,oploopsN)(&oi, 0, 0, NULL, B, C, 0, aA, aB, rC, aC);
   #else
      oi.a2blk(M, M, ATL_rone, aB, M, aA);
      Mjoin(PATL,oploopsN)(&oi, 0, 0, NULL, B, C, 0, aA, aB, aC, aC);
   #endif
   free(vp);
   return(0);
}

#ifndef Conj_
   #ifdef TCPLX
   #else
      #define ONE ATL_rone
      #define TA AtlasTrans
   #endif
void Mjoin(PATL,symmL_rec)
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
 * outer product gemm to solve this problem by reflecting block matrix to dense.
 */
   if (M <= ATL_VWopgen_LAST_KB)
   {
      const enum ATLAS_UPLO Uplo = (flg&1) ? AtlasUpper : AtlasLower;
      #ifdef TCPLX
         if (flg&2) /* really HEMM */
         {
            if (!Mjoin(PATL,ophemmL)(Uplo, M, N, alpha, A, lda, B, ldb,
                                     beta, C, ldc))
               return;
         }
         else
            if (!Mjoin(PATL,opsymmL)(Uplo, M, N, alpha, A, lda, B, ldb,
                                     beta, C, ldc))
               return;
      #else
         if (!Mjoin(PATL,opsymmL)(Uplo, M, N, alpha, A, lda, B, ldb,
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
 * Reason is copy of symat a big cost, and copying up front brings it through
 * cache twice, wheareas recursion will provide automatic blocking, and avoid
 * this for small problems where it matters.  Unlike more square cases,
 * due to lack of symmetric block reuse, ipsymm has no performance benefit
 * over continuing recursion.
 */
   else if (N > ATL_VIEW_BEST_NB)
   {
      if (!Mjoin(PATL,ipsymmL)(flg, M, N, alpha,A,lda,B,ldb,beta,C,ldc))
         return;
   }
/*
 * If neither of those work, split along M using recursive steps:
 *    C0 = A00 * B0 + beta*C0 (symm)
 *    C0 = A10^T * B1 + C0    (gemm)
 *    C1 = A10 * B0 + beta*C1 (gemm)
 *    C1 = A11 * B1 + C1      (symm)
 */
   {
      #ifdef TCPLX
         const TYPE ONE[2] = {ATL_rone, ATL_rzero};
         const enum ATLAS_TRANS TA = (flg&2) ? AtlasConjTrans : AtlasTrans;
      #endif
      const int ML=(M>>1), MR = M-ML;
      const TYPE *B1 = B + (ML SHIFT);
      TYPE *C1 = C + (ML SHIFT);

      ATL_assert(ML);  /* debugging, remove later */
      Mjoin(PATL,symmL_rec)(flg, ML, N, alpha, A, lda, B, ldb,
                            beta, C, ldc);
      if ((flg&1))  /* Upper Matrix */
      {
         const TYPE *A01 = A + ML*(lda SHIFT);
         Mjoin(PATL,ammm)(AtlasNoTrans, AtlasNoTrans, ML, N, MR, alpha,
                          A01, lda, B1, ldb, ONE, C, ldc);
         Mjoin(PATL,ammm)(TA, AtlasNoTrans, MR, N, ML, alpha,
                          A01, lda, B, ldb, beta, C1, ldc);
      }
      else          /* Lower matrix */
      {
         const TYPE *A10 = A + (ML SHIFT);
         Mjoin(PATL,ammm)(AtlasNoTrans, AtlasNoTrans, MR, N, ML, alpha,
                          A10, lda, B, ldb, beta, C1, ldc);
         Mjoin(PATL,ammm)(TA, AtlasNoTrans, ML, N, MR, alpha, A10, lda,
                          B1, ldb, ONE, C, ldc);
      }
      Mjoin(PATL,symmL_rec)(flg, MR, N, alpha, A+(ML SHIFT)*(lda+1),
                            lda, B1, ldb, ONE, C1, ldc);
   }
}
   #ifndef TCPLX
      #undef TA
   #endif
#else
   void Mjoin(PATL,symmL_rec)
      (ATL_CUINT, ATL_CSZT, ATL_CSZT, const SCALAR, const TYPE*, ATL_CSZT,
       const TYPE*, ATL_CSZT, const SCALAR, TYPE*, ATL_CSZT);
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
);
#ifdef Conj_
void Mjoin(PATL,hemm)
#else
void Mjoin(PATL,symm)
#endif
(
   const enum ATLAS_SIDE  Side,
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
   TYPE *aA, *aB, *aC, *bA;
   opinfo_t oi;
   #ifdef Conj_
      ATL_UINT bv = (Uplo == AtlasUpper)|2;
   #else
      ATL_UINT bv = (Uplo == AtlasUpper);
   #endif

   #ifdef USEREF
   {
      #ifdef Conj_
         Mjoin(PATL,refhemm)(Side, Uplo, M, N, alpha, A, lda, B, ldb,
                             beta, C, ldc);
      #else
         Mjoin(PATL,refsymm)(Side, Uplo, M, N, alpha, A, lda, B, ldb,
                             beta, C, ldc);
      #endif
      return;
   }
   #endif
   if (!M || !N)
      return;
   if (SCALAR_IS_ZERO(alpha))  /* doing nothing but scaling C */
   {
      if (!SCALAR_IS_ONE(beta))  /* if there is scaling of C to do */
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);

      return;
   }

   if (Side == AtlasLeft) /* K=M, A is GEMM's A */
   {
      #if defined(Conj_) || !defined(TCPLX)
      if (N == 1) /* really HEMV/SYMV */
      {
         #ifdef Conj_
            Mjoin(PATL,hemv)(Uplo, M, alpha, A, lda, B, 1, beta, C, 1);
         #else
            Mjoin(PATL,symv)(Uplo, M, alpha, A, lda, B, 1, beta, C, 1);
         #endif
         return;
      }
      #endif
      Mjoin(PATL,symmL_rec)(bv, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   }
   else                   /* Right: K=N, A is GEMM's B */
   {
      #if !defined(TCPLX)
      if (M == 1) /* really HEMV/SYMV */
      {
/*
 *       Form C^H = S^H B^H (getting symm on left for hemv),
 *       then convert C by conjugating.  This code works, but the
 *       extra work makes usually it slower than refblas or gemm.
 *       So, above outside cpp if doesn't allow this case.
 */
         #ifdef Conj_
         {
            void *vp;
            const TYPE CALP[2]={*alpha,-alpha[1]},
                       ONE[3]={ATL_rone, ATL_rzero, ATL_rzero};
            const ATL_iptr_t N2 = N+N;
            TYPE *x, *y;
            vp = malloc(ATL_MulBySize(N2)+ATL_Cachelen*2);
            ATL_assert(vp);
            x = ATL_AlignPtr(vp);
            y = x + N2;
            y = ATL_AlignPtr(y);
            Mjoin(PATL,moveConj)(N, CALP, B, ldb, x, 1);
            Mjoin(PATL,hemv)(Uplo, N, ONE, A, lda, x, 1, ONE+1, y, 1);
            Mjoin(PATL,axpbyConj)(N, ONE, y, 1, beta, C, ldc);
            free(vp);
            return;
         }
         #else
            Mjoin(PATL,symv)(Uplo, N, alpha, A, lda, B, ldb, beta, C, ldc);
            return;
         #endif
      }
      #endif
      Mjoin(PATL,symmR_rec)(bv, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
   }
}
