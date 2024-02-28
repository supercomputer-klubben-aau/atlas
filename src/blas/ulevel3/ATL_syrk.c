/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_level2.h"
#include "atlas_level3.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_sqsyrk.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_umsyrk.h))
#ifndef Conj_
   #define ATL_DECL_ 1
#endif
#include Mstr(Mjoin(ATLAS_PRE,syrk_view.h))
#ifndef Conj_
   #undef ATL_DECL_
#endif
#if 0
   #define USEREF 1
   #include "atlas_reflevel3.h"
#endif
/*
 * Service routine, particularly for parallel.  Takes its blocking from ip
 * (assuming that is what is being used below diagonal blocks
 * flag bits, meaning if set (opposite if unset):
 * 0/1: C is upper
 * 1/2: TA == AtlasNoTrans
 * 2/4: use beta=0 syrk kernel (else use beta=1)
 *
 * If (Uplo==Upper && sy2blk)
 *    (1) flag&1 == 1; (2) wU non-NULL; (3) sy2blk is beta=0.
 * ==> wU can be aliased wt wS if you don't need correct wS on output
 */
#ifdef Conj_
   #define syrk_amm  Mjoin(PATL,herk_amm)
   #define opsyrk  Mjoin(PATL,opherk)
   #define ipsyrk  Mjoin(PATL,ipherk)
   #define SYRK Mjoin(PATL,herk)
   #define herk_FLG 1
   #define umsyrk  Mjoin(PATL,umherk)
   #define sqsyrk  Mjoin(PATL,sqherk)
#else
   #define syrk_amm  Mjoin(PATL,syrk_amm)
   #define opsyrk  Mjoin(PATL,opsyrk)
   #define ipsyrk  Mjoin(PATL,ipsyrk)
   #define SYRK Mjoin(PATL,syrk)
   #define herk_FLG 0
   #define umsyrk  Mjoin(PATL,umsyrk)
   #define sqsyrk  Mjoin(PATL,sqsyrk)
#endif

int syrk_amm
(
   const enum ATLAS_UPLO  Uplo,
   const enum ATLAS_TRANS TA,
   ATL_iptr_t  N,
   ATL_iptr_t K,
   #ifdef Conj_
   const TYPE ralpha,
   #else
   const SCALAR alpha,
   #endif
   const TYPE *A,
   ATL_iptr_t lda,
   #ifdef Conj_
   const TYPE rbeta,
   #else
   const SCALAR beta,
   #endif
   TYPE *C,
   ATL_iptr_t ldc
)
{
   size_t sz, szA, szB, szC, szS, nnblks, extra;
   TYPE *wA, *wB, *wC, *wS, *wCs, *rC, *rCs, *rS;
   double timG;
   int nb, nbS, flg, idx;
   int i, k, ierr;
   #ifdef Conj_
      TYPE alpha[2]={ralpha, ATL_rzero}, beta[2]={rbeta, ATL_rzero};
      const enum ATLAS_TRANS TB=(TA==AtlasNoTrans)?AtlasConjTrans:AtlasNoTrans;
   #else
      const enum ATLAS_TRANS TB = (TA==AtlasNoTrans) ? AtlasTrans:AtlasNoTrans;
   #endif
   ipinfo_t ip, ipmen;
   cm2am_t sy2blk;
   ablk2cmat_t blk2sy, blk2c;

   #ifdef USEREF
   {
      #ifdef Conj_
         Mjoin(PATL,refherk)(Uplo, TA, N, K, ralpha, A, lda, rbeta, C, ldc);
      #else
         Mjoin(PATL,refsyrk)(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
      #endif
      return(0);
   }
   #endif
   if (N == 0)     /* no output */
      return(0);   /* so return with no-op */
   if (K < 3 || SCALAR_IS_ZERO(alpha))      /* degenerate case handled by*/
   {                                        /* ipsyrk */
      #ifdef Conj_
         ipsyrk(Uplo, TA, N, K, ralpha, A, lda, rbeta, C, ldc);
      #else
         ipsyrk(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
      #endif
      return(0);
   }

/*
 * Inner-product version calls nothing by SYRK kernels, and is called in
 * LAPACK's default left-looking Cholesky.  _IP handles L1/L2BLAS cases too.
 */
   timG = Mjoin(PATL,ipsyrkInfo)(&ip,herk_FLG, TA, N, K, lda, ldc, alpha, beta);
   #if 1
   if (timG < 0.0)
   {
      #ifdef Conj_
         ipsyrk(Uplo, TA, N, K, ralpha, A, lda, rbeta, C, ldc);
      #else
         ipsyrk(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
      #endif
      return(0);
   }
   #endif
   nb = (ip.nfnblks) ? ip.nb : ip.pnb;
/*
 * Outer product version will call outer-product-optimized amm, and since
 * N is too large to call ipsyrk, we expect amm perf to dominate this case.
 * Outer product SYRK is used in right-looking Cholesky.
 */
   #if 1
   if (K <= ATL_VWopgen_BEST_KB)
   {
      int ierr;
   #ifdef Conj_
      ierr = opsyrk(Uplo, TA, N, K, ralpha, A, lda, rbeta, C, ldc);
   #else
      ierr = opsyrk(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
   #endif
      if (ierr)
      {
         if (ierr == 1)  /* if we failed to malloc */
            return(1);   /* tell recursion to continue */
         if (nb >= N)    /* if we can do it with ipsyrk or opsyrk */
         {               /* use ipsyrk if opsyrk fails */
            #ifdef Conj_
               ipsyrk(Uplo, TA, N, K, ralpha, A, lda, rbeta, C, ldc);
            #else
               ipsyrk(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
            #endif
            return(0);
         }
      }
      else
         return(0);
   }
   #endif
/*
 * Will eventually need syrk timed for all square blocks to select best case.
 * For now, just pretend syrk time doesn't matter
 */
/*
 * We will decide between sqsyrk and umsyrk based on the pref data from
 * amm_syrkPerf.h.
 */
#if 1  /* apply umsyrk*/
   #ifdef Conj_
      ierr = umsyrk(Uplo, TA, N, K, ralpha, A, lda, rbeta, C, ldc);
   #else
      ierr = umsyrk(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
   #endif
#else /* sqsyrk */
   #ifdef Conj_
      ierr = sqsyrk(Uplo, TA, N, K, ralpha, A, lda, rbeta, C, ldc);
   #else
      ierr = sqsyrk(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
   #endif
#endif
   return(ierr);
}

static void syrk_rec
(
   const enum ATLAS_UPLO  Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT N,
   ATL_CSZT K,
   #ifdef Conj_
   const TYPE alpha,
   #else
   const SCALAR alpha,
   #endif
   const TYPE *A,
   ATL_iptr_t lda,
   #ifdef Conj_
   const TYPE beta,
   #else
   const SCALAR beta,
   #endif
   TYPE *C,
   ATL_CSZT ldc
)
{
   #if defined(Conj_) || defined(TREAL)
      #define ONE ATL_rone
   #else
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
   #endif
/*
 * If we cannot solve problem at this size, recursively cut dims until we can
 */
   if (syrk_amm(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc))
   {
      if (N > K)
      {
         ATL_iptr_t nL=N>>1, nR=N-nL;
         const TYPE *Ar = A+((TA==AtlasNoTrans)?(nL SHIFT):lda*(nL SHIFT));
         TYPE *Ct = C + (ldc+1)*(nL SHIFT);
         const enum ATLAS_TRANS TB =
         #ifdef Conj_
               (TA == AtlasNoTrans) ? AtlasConjTrans:AtlasNoTrans;
            const TYPE galp[2]={alpha,ATL_rzero}, gbet[2]={beta,ATL_rzero};
         #else
               (TA == AtlasNoTrans) ? AtlasTrans:AtlasNoTrans;
            #define galp alpha
            #define gbet beta
         #endif
         syrk_rec(Uplo, TA, nL, K, alpha, A, lda, beta, C, ldc);
         syrk_rec(Uplo, TA, nR, K, alpha, Ar, lda, beta, Ct, ldc);
         if (Uplo == AtlasUpper)
            Mjoin(PATL,gemm)(TA, TB, nL, nR, K, galp, A, lda, Ar, lda,
                             gbet, Ct-(nL SHIFT), ldc);
         else
            Mjoin(PATL,gemm)(TA, TB, nR, nL, K, galp, Ar, lda, A, lda,
                             gbet, C+(nL SHIFT), ldc);
         #ifndef Conj_
            #undef galp
            #undef gbet
         #endif
      }
      else
      {
         ATL_iptr_t kL=K>>1, kR=K-kL, incA=((TA == AtlasNoTrans)?lda:1)SHIFT;
         syrk_rec(Uplo, TA, N, kL, alpha, A, lda, beta, C, ldc);
         syrk_rec(Uplo, TA, N, kR, alpha, A+incA*kL, lda, ONE, C, ldc);
      }
   }
}
#if defined(Conj_) || defined(TREAL)
   #undef ONE
#endif

void SYRK
(
   const enum ATLAS_UPLO Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT N,
   ATL_CSZT K,
   #ifdef Conj_
   const TYPE alpha,
   #else
   const SCALAR alpha,
   #endif
   const TYPE *A,
   ATL_CSZT lda,
   #ifdef Conj_
   const TYPE beta,
   #else
   const SCALAR beta,
   #endif
   TYPE *C,
   ATL_CSZT ldc
)
{
#if 0
   #ifdef Conj_
      Mjoin(PATL,herk_APR)(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
   #else
      Mjoin(PATL,syrk_APR)(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
   #endif
   return;
#endif
   syrk_rec(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);

}
