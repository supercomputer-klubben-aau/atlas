/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine Petitet
 */
#include "atlas_misc.h"
#include "atlas_level3.h"
#if 0
   #define USEREF 1
   #include "atlas_reflevel3.h"
#endif

static unsigned int trsm_recR
/*
 * RETURNS: recursion depth (starting from 0)
 */
(
   const int             bv,    /* 0:Right,1:Upper,2:TransA,3:Conj, 4:NonUnit */
   ATL_CSZT               N,    /* size of triangle, will recur on this dim */
   ATL_CSZT               NRHS, /* don't recur on RHS */
   const SCALAR           alpha,
   const TYPE             *T,
   ATL_CSZT               ldt,
   ATL_CSZT               rincT, /* rowinc: 1/2 for NoTrans, else ldt2 */
   ATL_CSZT               cincT, /* colinc: ldt2 for NoTrans, else 1/2 */
   TYPE                   *Z,
   ATL_CSZT               ldz
)
{
   unsigned int iret=0;
   /*if (Mjoin(PATL,trsmR_IP)(bv, N, NRHS, alpha, T, ldt, Z, ldz))*/
   if (Mjoin(PATL,iptrsmR)(bv, N, NRHS, alpha, T, ldt, Z, ldz))
   {
      #ifdef TCPLX
         enum ATLAS_TRANS TB;
         const TYPE none[2] = {ATL_rnone, ATL_rzero};
         const TYPE one[2]  = {ATL_rone,  ATL_rzero};
      #else
         #define none ATL_rnone
         #define one ATL_rone
         const enum ATLAS_TRANS TB = (bv&4) ? AtlasTrans : AtlasNoTrans;
      #endif
      ATL_CSZT NL = (N>>1), NR = N - NL;
      const TYPE *T0, *T1;
      TYPE *C, *Z0;
      const TYPE *B;
      ATL_UINT LNUT;        /* is this LN or UT case? */

      #ifdef TCPLX
         if (bv&4)
            TB = (bv&8) ? AtlasConjTrans : AtlasTrans;
         else
            TB = (bv&8) ? AtlasConj : AtlasNoTrans;
      #endif
      LNUT = (bv&6);             /* get Upper(2) and Trans(4) bits */
      LNUT = (!LNUT)|(LNUT==6);  /* is this LN or UT case? */
      if (LNUT)                  /* LN or UT case */
      {
         T0 = T + NR*((ldt+1)SHIFT);
         T1 = T;
         B = T + NR*rincT;
         Z0 = Z + NR*(ldz SHIFT);
         C = Z;
      }
      else                       /* LT or UN case */
      {
         Z0 = Z;
         T0 = T;
         B = T + NL*rincT;
         C = Z + NL*(ldz SHIFT);
         T1 = B + NL*cincT;
      }
      iret = 2 + trsm_recR(bv, NL, NRHS, alpha, T0, ldt, rincT, cincT, Z0,ldz);
      Mjoin(PATL,gemm)(AtlasNoTrans, TB, NRHS, NR, NL, none, Z0, ldz, B, ldt,
                       alpha, C, ldz);
      iret += trsm_recR(bv, NR, NRHS, one, T1, ldt, rincT, cincT, C, ldz);
      return(iret);
   }
   return(iret);
}

static unsigned int trsm_recL
/*
 * RETURNS: recursion depth (starting from 0)
 */
(
   const int             bv,    /* 0:Right,1:Upper,2:TransA,3:Conj, 4:NonUnit */
   ATL_CSZT               N,    /* size of triangle, will recur on this dim */
   ATL_CSZT               NRHS, /* don't recur on RHS */
   const SCALAR           alpha,
   const TYPE             *T,
   ATL_CSZT               ldt,
   ATL_CSZT               rincT, /* rowinc: 1/2 for NoTrans, else ldt2 */
   ATL_CSZT               cincT, /* colinc: ldt2 for NoTrans, else 1/2 */
   TYPE                   *Z,
   ATL_CSZT               ldz
)
{
   unsigned int iret=0;
   /*if (Mjoin(PATL,trsmL_IP)(bv, N, NRHS, alpha, T, ldt, Z, ldz))*/
   if (Mjoin(PATL,iptrsmL)(bv, N, NRHS, alpha, T, ldt, Z, ldz))
   {
      #ifdef TCPLX
         enum ATLAS_TRANS TA;
         const TYPE none[2] = {ATL_rnone, ATL_rzero};
         const TYPE one[2]  = {ATL_rone,  ATL_rzero};
      #else
         #define none ATL_rnone
         #define one ATL_rone
         const enum ATLAS_TRANS TA = (bv&4) ? AtlasTrans : AtlasNoTrans;
      #endif
      ATL_CSZT NL = (N>>1), NR = N - NL;
      const TYPE *T0, *T1;
      TYPE *C, *Z0;
      const TYPE *A, *B;
      ATL_UINT LNUT;        /* is this LN or UT case? */

      #ifdef TCPLX
         if (bv&4)
            TA = (bv&8) ? AtlasConjTrans : AtlasTrans;
         else
            TA = (bv&8) ? AtlasConj : AtlasNoTrans;
      #endif
      LNUT = (bv&6);             /* get Upper(2) and Trans(4) bits */
      LNUT = (!LNUT)|(LNUT==6);  /* is this LN or UT case? */
      if (LNUT)                  /* LN or UT case */
      {
         T0 = T;
         A = T + NL*rincT;
         Z0 = Z;
         C = Z + (NL SHIFT);
         T1 = A + NL*cincT;
      }
      else                       /* LT or UN case */
      {
         T0 = T + NR*((ldt+1)SHIFT);
         T1 = T;
         A = T + NR*rincT;
         C = Z;
         Z0 = Z + (NR SHIFT);
      }
      iret = 2 + trsm_recL(bv, NL, NRHS, alpha, T0, ldt, rincT, cincT, Z0,ldz);
      if (bv&1)  /* Right case */
         Mjoin(PATL,gemm)(AtlasNoTrans, TA, NR, NRHS, NL, none, Z0, ldz, A, ldt,
                          alpha, C, ldz);
      else       /* Left case */
         Mjoin(PATL,gemm)(TA, AtlasNoTrans, NR, NRHS, NL, none, A, ldt, Z0, ldz,
                          alpha, C, ldz);
      iret += trsm_recL(bv, NR, NRHS, one, T1, ldt, rincT, cincT, C, ldz);
      return(iret);
   }
   return(iret);
}
void Mjoin(PATL,trsm)
(
   const enum ATLAS_SIDE  SD,
   const enum ATLAS_UPLO  UL,
   const enum ATLAS_TRANS TA,
   const enum ATLAS_DIAG  DI,
   ATL_CSZT               M,
   ATL_CSZT               N,
   const SCALAR           alpha,
   const TYPE             *A,
   ATL_CSZT               lda,
   TYPE                   *B,
   ATL_CSZT               ldb
)
/*
 * Purpose
 * =======
 * ATL_trsm solves one of the matrix equations
 *    op( A ) * X = alpha * B,   or  X * op( A ) = alpha * B,
 *
 * where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 * non-unit, upper or lower triangular matrix and op( A ) is one of
 *    op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
 *
 * The matrix X is overwritten on B.
 */
{
   ATL_SZT incR, incC;
   ATL_UINT r, bv = (SD == AtlasRight) ? 1 : 0;
            /* 0:Right, 1:Upper, 2:TransA, 3: Conj, 4:NonUnit */

   #ifdef USEREF
   {
      Mjoin(PATL,reftrsm)(SD, UL, TA, DI, M, N, alpha, A, lda, B, ldb);
      return;
   }
   #endif
   if (!M || !N)  /* if either array is of size 0 */
      return;     /* return as there is nothing to do */
   if (SCALAR_IS_ZERO(alpha))                   /* alpha == 0 means */
   {
      Mjoin(PATL,gescal)(M, N, alpha, B, ldb);  /* only need to scale B */
      return;
   }
   if (UL == AtlasUpper)
   {
      bv |= 2;
      incR = lda SHIFT;
      incC = 1 SHIFT;
   }
   else
   {
      incR = 1 SHIFT;
      incC = lda SHIFT;
   }
   #ifdef TCPLX
      if (TA == AtlasConjTrans)
         bv |= (8|4);
      else
         bv |= TA == (AtlasTrans) ? 4 : 0;
   #else
      bv |= (TA == AtlasNoTrans) ? 0:4;
   #endif
   bv |= (DI == AtlasNonUnit) ? 16 : 0;
/*
 * LEFT: op(A) * X = alpha * B, A is MxM, X is MxN, B is MxN
 */
   if (SD == AtlasLeft)
      r = trsm_recL(bv, M, N, alpha, A, lda, incR, incC, B, ldb);
   else
      r = trsm_recR(bv, N, M, alpha, A, lda, incR, incC, B, ldb);
#ifdef DEBUG
   printf("Solved with %u recursive calls.\n", r);
#endif
}
