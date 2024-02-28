#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_lvl2.h"
#include "atlas_level1.h"
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))
#if 0
   #define USEREF 1
   #include "atlas_reflevel3.h"
#endif

static int ATL_ammm
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
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
/*
 * Just do a scale and return
 */
   if (SCALAR_IS_ZERO(alpha) || !K)
   {
      if (SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gezero)(M, N, C, ldc);
      else if (!SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);
      return(0);
   }
/*
 * Scope for degenerate cases that should call Level-2 BLAS; these
 * routines assert they work, since their workspace is O(N) and so are
 * not allowed to fail.
 */
   if (K == 1)  /* really a GER */
   {
   #ifdef TCPLX
      if (!SCALAR_IS_ONE(beta))  /* can't use GER for beta != 1 */
      {
         int i;
         const register TYPE ral=alpha[0], ial=alpha[1];
         const size_t ldc2 = ldc+ldc;
         ATL_CSZT incB = ((TB == AtlasNoTrans) ? ldb : 1)SHIFT;
         const register TYPE
            cjm = (TB==AtlasConj || TB==AtlasConjTrans) ? ATL_rnone:ATL_rone;
         TYPE *X=(TYPE*)A;
         void *vp=NULL;
/*
 *       Copy A if it's a row or if it must be conjugated
 */
         if (TA == AtlasTrans || TA == AtlasConjTrans || TA == AtlasConj)
         {
            vp = malloc(ATL_MulBySize(M)+ATL_Cachelen);
            ATL_assert(vp);
            X = ATL_AlignPtr(vp);
            if (TA == AtlasTrans)
               Mjoin(PATL,copy)(M, A, lda, X, 1);
            else if (TA == AtlasConjTrans)
               Mjoin(PATL,copyConj)(M, A, lda, X, 1);
            else
               Mjoin(PATL,copyConj)(M, A, 1, X, 1);
         }
         for (i=0; i < N; i++, B += incB, C += ldc2)
         {
            TYPE scal[2];
            register TYPE rb=(*B), ib=cjm*B[1];
            scal[0] = rb*ral - ib*ial;
            scal[1] = rb*ial + ib*ral;
            Mjoin(PATL,axpby)(M, scal,  X, 1, beta, C, 1);
         }
         if (vp)
            free(vp);
      }
      else  /* BETA=1, can use GERU/GERC */
      {
         if (TA == AtlasConjTrans || TA == AtlasConj)  /* must copyConj A */
         {
            void *vp;
            TYPE *X;
            const TYPE ONE[2] = {ATL_rone, ATL_rzero};
            vp = malloc(ATL_MulBySize(M)+ATL_Cachelen);
            ATL_assert(vp);
            X = ATL_AlignPtr(vp);
            Mjoin(PATL,moveConj)(M, alpha, A, (TA == AtlasConj) ? 1:lda, X, 1);
            if (TB == AtlasConjTrans || TB == AtlasConj)  /* use gerc */
               Mjoin(PATL,gerc)(M, N, ONE, X, 1, B,
                                (TB == AtlasConj) ? ldb : 1, C, ldc);
            else
               Mjoin(PATL,geru)(M, N, ONE, X, 1, B,
                                (TB == AtlasNoTrans) ? ldb : 1, C, ldc);
            free(vp);
         }
         else if (TB == AtlasConjTrans || TB == AtlasConj)  /* use gerc */
            Mjoin(PATL,gerc)(M, N, alpha, A, (TA == AtlasNoTrans) ? 1 : lda,
                             B, (TB == AtlasConj) ? ldb : 1, C, ldc);
         else /* use geru */
            Mjoin(PATL,geru)(M, N, alpha, A, (TA == AtlasNoTrans) ? 1 : lda,
                            B, (TB == AtlasNoTrans) ? ldb : 1, C, ldc);
      }
   #else
      if (!SCALAR_IS_ONE(beta)) /* can't use GER for beta != 1 */
      {
         int i;
         ATL_CSZT incA = ((TA == AtlasNoTrans) ? 1 : lda);
         ATL_CSZT incB = ((TB == AtlasNoTrans) ? ldb : 1);
         for (i=0; i < N; i++, B += incB, C += ldc)
            Mjoin(PATL,axpby)(M, alpha * *B,  A, incA, beta, C, 1);
      }
      else
         Mjoin(PATL,ger)(M, N, alpha, A, (TA == AtlasNoTrans) ? 1 : lda,
                         B, (TB == AtlasNoTrans) ? ldb : 1, C, ldc);
   #endif
      return(0);
   }
   if (K == 2)
      return(Mjoin(PATL,ammm_rk2)(TA, TB, M, N, alpha, A, lda, B, ldb,
                                  beta, C, ldc));
   if (N == 1)  /* Really GEMV with A as matrix, A & C as vectors */
   {
   #ifdef TCPLX
      TYPE *X = (TYPE*)B;
      void *vp = NULL;
      int incX = 1;
/*
 *    Copy B if it we need to conjugate it
 */
      if (TB == AtlasConj || TB == AtlasConjTrans)
      {
         vp = malloc(ATL_MulBySize(K)+ATL_Cachelen);
         ATL_assert(vp);
         X = ATL_AlignPtr(vp);
         Mjoin(PATL,copyConj)(K, B, (TB == AtlasConj) ? 1:ldb, X, 1);
      }
      else
         incX = (TB == AtlasNoTrans) ? 1:ldb;
      if (TA == AtlasNoTrans || TA == AtlasConj)
         Mjoin(PATL,gemv)(TA, M, K, alpha, A, lda, X, incX, beta, C, 1);
      else /* if (TA == AtlasTrans || TA == AtlasConjTrans) */
         Mjoin(PATL,gemv)(TA, K, M, alpha, A, lda, X, incX, beta, C, 1);
      if (vp)
         free(vp);
   #else
      if (TA == AtlasNoTrans)
         Mjoin(PATL,gemv)(AtlasNoTrans, M, K, alpha, A, lda, B,
                          (TB == AtlasNoTrans) ? 1:ldb, beta, C, 1);
      else
         Mjoin(PATL,gemv)(AtlasTrans, K, M, alpha, A, lda, B,
                          (TB == AtlasNoTrans) ? 1:ldb, beta, C, 1);
   #endif
      return(0);
   }
   if (M == 1)  /* Really GEMV with B as matrix, A & C as vectors */
   {
   #ifdef TCPLX
      TYPE *X = (TYPE*)A;
      void *vp = NULL;
      int incX = 1;
/*
 *    Copy A if it we need to conjugate it
 */
      if (TA == AtlasConj || TA == AtlasConjTrans)
      {
         vp = malloc(ATL_MulBySize(K)+ATL_Cachelen);
         ATL_assert(vp);
         X = ATL_AlignPtr(vp);
         Mjoin(PATL,copyConj)(K, A, (TA == AtlasConj) ? lda:1, X, 1);
      }
      else
         incX = (TA == AtlasNoTrans) ? lda:1;
      if (TB == AtlasNoTrans)
         Mjoin(PATL,gemv)(AtlasTrans, K, N, alpha, B, ldb, X, incX,
                          beta, C, ldc);
      else if (TB == AtlasConj)
         Mjoin(PATL,gemv)(AtlasConjTrans, K, N, alpha, B, ldb, X, incX,
                          beta, C, ldc);
      else if (TB == AtlasTrans)
         Mjoin(PATL,gemv)(AtlasNoTrans, N, K, alpha, B, ldb, X, incX,
                          beta, C, ldc);
      else if (TB == AtlasConjTrans)
         Mjoin(PATL,gemv)(AtlasConj, N, K, alpha, B, ldb, X, incX,
                          beta, C, ldc);
      if (vp)
         free(vp);
   #else
      if (TB == AtlasNoTrans)
         Mjoin(PATL,gemv)(AtlasTrans, K, N, alpha, B, ldb, A,
                          (TA == AtlasNoTrans) ? lda:1, beta, C, ldc);
      else
         Mjoin(PATL,gemv)(AtlasNoTrans, N, K, alpha, B, ldb, A,
                          (TA == AtlasNoTrans) ? lda:1, beta, C, ldc);
   #endif
      return(0);
   }
/*
 * Special case mainly for real LU, where K==4, N==4 beta=1.0, TA==AtlasNoTrans;
 * Can do a no-copy update with only one loop.
 */
   #ifndef TCPLX
   if (K == 4 && N == 4 && TA==AtlasNoTrans && SCALAR_IS_ONE(beta))
   {
      int Mjoin(PATL,rk4n4)(enum ATLAS_TRANS,ATL_CSZT,const SCALAR,
          const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,TYPE*,ATL_CSZT);
      if (!Mjoin(PATL,rk4n4)(TB, M, alpha, A, lda, B, ldb, C, ldc))
         return(0);
   }
   #endif
/*
 * Rank-K can fail to allocate space, so return success/failure
 */
   if (K <= ATL_VWopgen_MAX_KB)
      return(Mjoin(PATL,ammm_rkK)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                  beta, C, ldc));
/*
 * Handle case that is really an inner product shape (M<=MB, N<=NB, large K)
 * This case not allowed to fail since it requires only 3*NB^2 workspace.
 */
   if (M <= ATL_VWipgen_MAX_MB && N <= ATL_VWipgen_MAX_NB)
      if (!Mjoin(PATL,ammm_IP)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                               beta, C, ldc))
         return(0);
/*
 * If B/C have only one column panel, call special low-workspace (3NB^3)
 * code for additional performance.  This shape occurs in left-looking LU.
 */
   if (N <= ATL_VWipgen_BEST_NB)
      return(Mjoin(PATL,ammm_tN)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                 beta, C, ldc));
/*
 * Next two loop orderings are general case, so use whichever uses least
 * workspace
 */
   if (M > N)
      return(Mjoin(PATL,ammmKMNK)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                  beta, C, ldc));
   return(Mjoin(PATL,ammmKNMK)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                               beta, C, ldc));
}
/*
 * Recur to get K below this value; this puts a ceiling on workspace and
 * usually improves performance (in huge problems, reduces TLB pressure)
 */
#define ATL_MAX_RK 3000


/*
 * This routine uses recursion to cut the dimensions of the matrices until
 * workspace requirements are low enough that a call to ATL_ammm succeeds
 */
void Mjoin(PATL,ammm)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
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
   #ifdef USEREF
   Mjoin(PATL,refgemm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   return;
   #endif
   if (!M || !N)
      return;
/*
 * Cases where all we must do is possibly scale and return
 */
   if (SCALAR_IS_ZERO(alpha) || !K)
   {
      if (SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gezero)(M, N, C, ldc);
      else if (!SCALAR_IS_ONE(beta))
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);
      return;
   }
/*
 * Our stopping criteria is if ATL_ammm signals success in mallocing mem
 */
   if (K <= ATL_MAX_RK)
   {
      if(!ATL_ammm(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc))
         return;
   }
/*
 * =========================================================================
 * Otherwise, problem too large, so we'll recursively divide its largest dim
 * =========================================================================
 */
/*
 * if K is tied for largest, cut it, since it reduces size of A & B
 * NOTE: C always uses only NB^2 workspace, so only A/B matters.
 */
   if (K > ATL_MAX_RK || (K >= N && K >= M))
   {
      const size_t kL=(K>>4)<<3, kR=K-kL;
      #ifdef TCPLX
         const TYPE ONE[2] = {ATL_rone, ATL_rzero};
      #else
         #define ONE ATL_rone
      #endif

      Mjoin(PATL,ammm)(TA, TB, M, N, kL, alpha, A, lda, B, ldb, beta, C, ldc);
      if (TA == AtlasNoTrans || TA == AtlasConj)
         A += (lda*kL)SHIFT;
      else
         A += kL SHIFT;
      if (TB == AtlasNoTrans)
         B += kL SHIFT;
      else
         B += (ldb*kL) SHIFT;
      Mjoin(PATL,ammm)(TA, TB, M, N, kR, alpha, A, lda, B, ldb,
                       ONE, C, ldc);
      #ifdef ONE
         #undef ONE
      #endif
   }
   else if (N >= M)  /* cutting N */
   {
      const size_t nL = (N>>1), nR = N-nL;
      Mjoin(PATL,ammm)(TA, TB, M, nL, K, alpha, A, lda, B, ldb, beta, C, ldc);
      if (TB == AtlasNoTrans || TB == AtlasConj)
         B += (ldb*nL)SHIFT;
      else
         B += nL SHIFT;
      C += (ldc*nL)SHIFT;
      Mjoin(PATL,ammm)(TA, TB, M, nR, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
   else  /* cutting M */
   {
      const size_t mL = (M>>1), mR = M-mL;
      Mjoin(PATL,ammm)(TA, TB, mL, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      if (TA == AtlasNoTrans || TA == AtlasConj)
         A += mL SHIFT;
      else
         A += (mL*lda)SHIFT;
      C += mL SHIFT;
      Mjoin(PATL,ammm)(TA, TB, mR, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
}

void Mjoin(PATL,gemm)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                      ATL_CSZT M, ATL_CSZT N, ATL_CSZT K, const SCALAR alpha,
                      const TYPE *A, ATL_CSZT lda, const TYPE *B,
                      ATL_CSZT ldb, const SCALAR beta, TYPE *C, ATL_CSZT ldc)
{
   Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
