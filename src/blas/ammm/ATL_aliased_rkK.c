#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_level1.h"
#include "atlas_lvl2.h"

static void *FixVector(enum ATLAS_TRANS TX, ATL_CSZT N, const SCALAR alpha,
                       const TYPE *X, ATL_CSZT incX)
{
   void *vx;
   TYPE *x;
   vx = malloc(ATL_MulBySize(N)+ATL_Cachelen);
   ATL_assert(vx);
   x = ATL_AlignPtr(vx);
   #ifndef TCPLX
      if (SCALAR_IS_ONE(alpha))
         Mjoin(PATL,copy)(N, X, incX, x, 1);
      else
         Mjoin(PATL,cpsc)(N, alpha, X, incX, x, 1);
   #else
      if (SCALAR_IS_ONE(alpha))
      {
         if (TX == AtlasTrans || TX == AtlasNoTrans)
            Mjoin(PATL,copy)(N, X, incX, x, 1);
         else
            Mjoin(PATL,copyConj)(N, X, incX, x, 1);
      }
      else
      {
         if (TX == AtlasTrans || TX == AtlasNoTrans)
            Mjoin(PATL,cpsc)(N, alpha, X, incX, x, 1);
         else
            Mjoin(PATL,moveConj)(N, alpha, X, incX, x, 1);
      }
   #endif
   return(vx);
}

/*
 * This entry makes rkK safe for L3kernel aliased calls.  It handles
 * only the aliasing required by the L3kernels, namely square blocks
 * less than ATLAS's largest blocking factor for the square dimensions,
 * with one of A/B aliased with C, and aliased by having
 * either A == C or B == C (i.e., not a partial overlap).  When A==C,
 * M=K < ATL_rkAMM_LASTKB; when B==C, N=K < ATL_rkAMM_LASTKB.
 */
int Mjoin(PATL,ammm_aliased_rkK)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alp,
   const TYPE *A,
   ATL_CSZT lda,
   const TYPE *B,
   ATL_CSZT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   #ifdef TCPLX
      const TYPE ONE[2] = {ATL_rone, ATL_rzero};
   #else
      #define ONE ATL_rone
   #endif
   void *vp=NULL;

   if (K == 0 || SCALAR_IS_ZERO(alp))
   {
      if (SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gezero)(M, N, C, ldc);
      else if (!SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);
      return(0);
   }
   if (K == 1)
   {
      #ifdef TCPLX
         const SCALAR alpha = alp;
      #else
         TYPE alpha = alp;
      #endif
      void *xp=NULL, *yp=NULL;
      TYPE *a = (TYPE*)A, *b = (TYPE*)B;
      size_t ldA=lda, ldB=ldb;

      if (A == C)
      {
         xp = FixVector(TA, M, alpha, A,
                        (TA == AtlasTrans || TA == AtlasConjTrans) ? lda:1);
         a = ATL_AlignPtr(xp);
         alpha = ONE;
         if (TA == AtlasConjTrans || TA == AtlasTrans)
         {
            TA = AtlasTrans;
            ldA = 1;
         }
         else if (TA == AtlasConj)
            TA = AtlasNoTrans;
      }
      if (B == C)
      {
         yp = FixVector(TB, N, alpha, B,
                        (TB == AtlasTrans || TB == AtlasConjTrans) ? 1:ldb);
         b = ATL_AlignPtr(yp);
         alpha = ONE;
         if (TB == AtlasConjTrans)
            TB = AtlasTrans;
         else if (TB == AtlasConj || TB == AtlasNoTrans)
         {
            TA = AtlasNoTrans;
            ldB = 1;
         }
      }
      Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, a, ldA, b, ldB, beta, C, ldc);
      if (xp)
         free(xp);
      if (yp)
         free(yp);
      return(0);
   }

   if (K == 2)
   {
      #ifdef TCPLX
         const SCALAR alpha = alp;
      #else
         TYPE alpha = alp;
      #endif
/*
 *    If BETA != 1, ammm_rk2 will copy all inputs and thus aliasing safe
 */
      if (!SCALAR_IS_ONE(beta))
         Mjoin(PATL,ammm_rk2)(TA, TB, M, N, alpha, A, lda, B, ldb, beta,
                              C, ldc);
/*
 *    For beta = 1, copy aliased input array(s) and then call GER2
 */
      else
      {
         void *wp=NULL, *xp=NULL, *yp=NULL, *zp=NULL;
         TYPE *w, *x, *y, *z;
         ATL_SZT incX, incY;
      #ifndef TCPLX
         if (A == C)
         {
            ATL_CSZT incx = (TA == AtlasTrans) ? lda:1;
      #else
         if (A == C || TA == AtlasConjTrans || TA == AtlasConj)
         {
            ATL_CSZT incx = (TA == AtlasTrans || TA == AtlasConjTrans) ? lda:1;
      #endif
            wp = FixVector(TA, M, alpha, A, incx);
            xp = FixVector(TA, M, alpha, A+((incx==1 ? lda:1)SHIFT), incx);
            alpha = ONE;
            w = ATL_AlignPtr(wp);
            incX = 1;
         }
         else  /* don't need to copy A */
         {
            w = (TYPE*)A;
            if (TA == AtlasNoTrans)
            {
               x = (TYPE*)(A + (lda SHIFT));
               incX = 1;
            }
            else /* if (TA == AtlasTrans) */
            {
               x = (TYPE*)(A + (1 SHIFT));
               incX = lda;
            }
         }
         if (B == C)
         {
            ATL_CSZT incy = (TB == AtlasTrans || TB == AtlasConjTrans) ? 1:ldb;
            yp = FixVector(TB, N, alpha, A, incy);
            zp = FixVector(TB, N, alpha, A+(((incy==1)?ldb:1)SHIFT), incy);
            y = ATL_AlignPtr(yp);
            z = ATL_AlignPtr(zp);
         #ifdef TCPLX
            if (TB == AtlasConj)
               TB = AtlasNoTrans;
            else if (TB == AtlasConjTrans)
               TB = AtlasTrans;
         #endif
            incY = 1;
            alpha = ONE;

         }
         else  /* no need to copy B */
         {
            y = (TYPE*)B;
         #ifdef TCPLX
            if (TB == AtlasNoTrans || TB == AtlasConj)
         #else
            if (TB == AtlasNoTrans)
         #endif
            {
               incY = ldb;
               z = (TYPE*)(B + (1 SHIFT));
            }
            else
            {
               incY = 1;
               z = (TYPE*)(B + (ldb SHIFT));
            }
         }
         #ifndef TCPLX
            Mjoin(PATL,ger2)(M, N, alpha, w, incX, y, incY, ONE,
                             x, incX, z, incY, C, ldc);
         #else
            if (TB == AtlasNoTrans || TB == AtlasTrans)
               Mjoin(PATL,ger2u)(M, N, alpha, w, incX, y, incY, ONE,
                                 x, incX, z, incY, C, ldc);
            else
               Mjoin(PATL,ger2c)(M, N, alpha, w, incX, y, incY, ONE,
                                 x, incX, z, incY, C, ldc);
         #endif
         if (wp)
            free(wp);
         if (xp)
            free(xp);
         if (yp)
            free(yp);
         if (zp)
            free(zp);
         return(0);
      }
      return(0);
   }
/*
 * For K > 3, ATL_ammm_rkK is safe for these precise aliasing conditions
 */
   ATL_assert(!Mjoin(PATL,ammm_rkK)(TA, TB, M, N, K, alp, A, lda, B, ldb,
                                    beta, C, ldc));
   return(0);
}
