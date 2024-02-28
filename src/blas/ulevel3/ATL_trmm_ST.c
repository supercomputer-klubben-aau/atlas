/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Majedul Sujon
 */
#include "atlas_misc.h"
#include "atlas_lvl3.h"
#ifdef Right_
   #define ATL_trmm_ST Mjoin(PATL,trmmR_ST)
   #define ATL_utrmm_LN Mjoin(PATL,utrmmR_LN)
   #define ATL_utrmm_LT Mjoin(PATL,utrmmR_LT)
   #define ATL_utrmm_UN Mjoin(PATL,utrmmR_UN)
   #define ATL_utrmm_UT Mjoin(PATL,utrmmR_UT)
   #ifdef TCPLX
      #define ATL_utrmm_LC Mjoin(PATL,utrmmR_LC)
      #define ATL_utrmm_UC Mjoin(PATL,utrmmR_UC)
   #endif
#else
    #define ATL_trmm_ST Mjoin(PATL,trmmL_ST)
   #define ATL_utrmm_LN Mjoin(PATL,utrmmL_LN)
   #define ATL_utrmm_LT Mjoin(PATL,utrmmL_LT)
   #define ATL_utrmm_UN Mjoin(PATL,utrmmL_UN)
   #define ATL_utrmm_UT Mjoin(PATL,utrmmL_UT)
   #ifdef TCPLX
      #define ATL_utrmm_LC Mjoin(PATL,utrmmL_LC)
      #define ATL_utrmm_UC Mjoin(PATL,utrmmL_UC)
   #endif
#endif
/*#define DEBUG 1*/
/*
 * This is basically a wrapper function to TRMM of small triangle
 */
int ATL_trmm_ST
(
   ATL_UINT bv,  /* 0:Right, 1:Upper, 2:TransA, 3: Conj, 4:NonUnit */
   ATL_CSZT               NT,  /* number of diagonals in triangle */
   ATL_CSZT               NR,  /* number of right hand sides */
   const SCALAR           alpha,
   const TYPE             *A,
   ATL_CSZT               lda,
   TYPE                   *B,
   ATL_CSZT               ldb
)
{
   unsigned int ctu; /* bits for conj-trans-upper */
   int res;
/*
 * prototypes for trmm micro-kernels
 */
   int (*utrmm)(const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx);

   int ATL_utrmm_LN(const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx);
   int ATL_utrmm_UT(const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx);
   int ATL_utrmm_LT(const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx);
   int ATL_utrmm_UN(const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx);
   #ifdef TCPLX
   int ATL_utrmm_UC(const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx);
   int ATL_utrmm_LC(const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx);
   #endif
/*
 *       C T U
 * LNUT  0 0 0 -- LN = 0
 *       0 1 1 -- UT = 3
 *       1 1 1 -- UC = 7
 *
 * LTUN  0 1 0 -- LT = 2
 *       1 1 0 -- LC = 6
 *       0 0 1 -- UN = 1
 *
 */
      ctu = (bv >> 1) & 7;

      if (!ctu)
         utrmm = ATL_utrmm_LN;
      else if (ctu == 2)
         utrmm = ATL_utrmm_LT;
      #ifdef TCPLX
         else if (ctu == 6)
            utrmm = ATL_utrmm_LC;
      #endif
      else if (ctu == 1)
         utrmm = ATL_utrmm_UN;
      else if (ctu == 3)
         utrmm = ATL_utrmm_UT;
      #ifdef TCPLX
      else if (ctu == 7)
         utrmm = ATL_utrmm_UC;
      #endif
      else
      {
           fprintf(stderr,"**** Unknown TRMM case!!\n");
           /*assert(0);*/
           return(1);
      }
   if (bv & 1)
      res = utrmm((bv&16)?AtlasNonUnit:AtlasUnit, NR, NT, alpha, A, lda, B, ldb);
   else
      res = utrmm((bv&16)?AtlasNonUnit:AtlasUnit, NT, NR, alpha, A, lda, B, ldb);
   return(res);
}
