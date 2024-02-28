#include "atlas_misc.h"
#include "atlas_lvl3.h"
/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 */
#ifdef Right_
   #define ATL_trsm_IP Mjoin(PATL,trsmR_IP)
    #define ATL_sminfoN Mjoin(PATL,sminfoR_LN)
    #define ATL_sminfoT Mjoin(PATL,sminfoR_LT)
    #define ATL_utrsm_alloc Mjoin(PATL,utrsmR_alloc)
    #define INCR 1
#else
    #define ATL_trsm_IP Mjoin(PATL,trsmL_IP)
    #define ATL_sminfoN Mjoin(PATL,sminfoL_LN)
    #define ATL_sminfoT Mjoin(PATL,sminfoL_LT)
    #define ATL_utrsm_alloc Mjoin(PATL,utrsmL_alloc)
    #define INCR ldb
#endif
int ATL_trsm_IP
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
   int RECUR;
   RECUR=0;
   #ifdef DEBUG
   printf("D=(%u,%u), SD=%c, UP=%c TA=%c DI=%c, alpha=%e\n", NT, NR,
          (bv&1)?'R':'L', (bv&2)?'U':'L', (bv&4)?'T':'N', (bv&16)?'N':'U',
          SVAL alpha);
   #endif
   if (NT == 1)
   {
      if (bv&16)   /* if non-unit, work to do, else nothing to do at M=1! */
      #ifdef TCPLX
      {
         TYPE inv[2], tmp;
         *inv = *A;
         tmp = A[1];
         inv[1] = (bv&8) ? -tmp : tmp;
         if (SCALAR_IS_ONE(alpha))
         {  CXINV(inv, inv); }
         else
            Mjoin(PATL,cplxdivide)(1, inv, (TYPE*)alpha, 1, inv, 1);
         Mjoin(PATL,scal)(NR, inv, B, INCR);
      }
      #else
         Mjoin(PATL,scal)(NR, alpha / *A, B, INCR);
      #endif
      else if (!SCALAR_IS_ONE(alpha))  /* unit may still need to scale! */
         Mjoin(PATL,scal)(NR, alpha, B, INCR);
   }
   else if (RECUR)  /* if told to */
      return(1);    /* do pure recursion */
#if 1     /* default case uses utrsm */
   else
   {
      ATL_UINT LN;
      void *vp;
      TYPE *diag, *T, *R, *w;
      #ifdef TCPLX
         enum ATLAS_TRANS TRANSA= ((bv&12)==12) ? AtlasConjTrans :  AtlasTrans;
      #else
         #define TRANSA AtlasTrans
      #endif
      void (*utrsm)(sminfo_t*,const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx,
         TYPE *diag, TYPE *L, TYPE *RW, TYPE *w);
      sminfo_t si;
      int ATL_sminfoN
         (sminfo_t *ip, ATL_CUINT bv, ATL_CSZT N, ATL_CSZT R,
          const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
      int ATL_sminfoT
         (sminfo_t *ip, ATL_CUINT bv, ATL_CSZT N, ATL_CSZT R,
          const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);

      LN = bv & 6;
      LN = (!LN) |(LN == 6);
      if (LN) /* Lower,NoTrans & Upper,Trans use LN case */
         ATL_sminfoN(&si, bv, NT, NR, alpha, lda, ldb);
      else              /* Upper,NoTrans & Lower,Trans use LT case */
         ATL_sminfoT(&si, bv, NT, NR, alpha, lda, ldb);
      if (NT > si.kb)
         return(2);
      #ifdef DEBUG
         fprintf(stderr, "NT=%u, kb=%u\n", NT, si.kb);
      #endif
      vp = ATL_utrsm_alloc(&si, NT, &diag, &T, &R, &w);
      if (!vp)
         return(1);
      utrsm = si.utrsm;
      if (bv & 1)
         utrsm(&si, (bv&16)?AtlasNonUnit:AtlasUnit, NR, NT, alpha, A, lda,
               B, ldb, diag, T, R, w);
      else
         utrsm(&si, (bv&16)?AtlasNonUnit:AtlasUnit, NT, NR, alpha, A, lda,
               B, ldb, diag, T, R, w);
      free(vp);
   }
#else
   else
   {
      enum ATLAS_TRANS TA;
      #ifdef TCPLX
         if ((bv&12)==12)
            TA == AtlasConjTrans;
         else
      #endif
      TA = (bv&4) ? AtlasTrans:AtlasNoTrans;
      #ifdef DEBUG
      printf("D=(%u,%u), UP=%c TA=%c DI=%c, alpha=%e\n", M, N, (bv&2)?'U':'L',
             TA==AtlasTrans?'T':'N', (bv&16)?'N':'U', alpha);
      #endif
      Mjoin(PATL,trsm_APR)(AtlasLeft, (bv&2)?AtlasUpper:AtlasLower, TA,
         (bv&16)?AtlasNonUnit:AtlasUnit, M, N, alpha, A, lda, B, ldb);
   }
#endif
   return(0);
}
