/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Majedul Sujon
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_lvl3.h"
#ifdef Right_
   #define ATL_iptrsm Mjoin(PATL,iptrsmR)
   #define ATL_sminfoN Mjoin(PATL,sminfoR_LN)
   #define ATL_sminfoT Mjoin(PATL,sminfoR_LT)
   #define ATL_utrsm_alloc Mjoin(PATL,utrsmR_alloc)
   #define INCR 1
   #define ATL_trsm_alloc Mjoin(PATL,trsmR_alloc)
   #define ATL_trsm_LNUT Mjoin(PATL,trsmR_LNUT)
   #define ATL_trsm_LTUN Mjoin(PATL,trsmR_LTUN)
#else
   #define ATL_iptrsm Mjoin(PATL,iptrsmL)
   #define ATL_sminfoN Mjoin(PATL,sminfoL_LN)
   #define ATL_sminfoT Mjoin(PATL,sminfoL_LT)
   #define ATL_utrsm_alloc Mjoin(PATL,utrsmL_alloc)
   #define INCR ldb
   #define ATL_trsm_alloc Mjoin(PATL,trsmL_alloc)
   #define ATL_trsm_LNUT Mjoin(PATL,trsmL_LNUT)
   #define ATL_trsm_LTUN Mjoin(PATL,trsmL_LTUN)
#endif

/*#define DEBUG 1 */

int ATL_iptrsm
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
/*
 * Handle degenrated cases
 */
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
/*
 *    general case
 */
      ATL_UINT LN, ctu;
      ATL_SZT Tsz, Rsz;
      void *vp;
      ipinfo_t ip;
      sminfo_t si;
      TYPE *diag, *T, *L, *R, *w;
      #ifdef TCPLX
         TYPE ONE[2] = {ATL_rone, ATL_rzero};
         TYPE ZERO[2] = {ATL_rzero, ATL_rzero};
         TYPE NONE[2] = {ATL_rnone, ATL_rzero};
      #else
         #define ONE ATL_rone
         #define ZERO ATL_rzero
         #define NONE ATL_rnone
      #endif
      enum ATLAS_TRANS TA;
      #ifdef TCPLX
         if ((bv&12)==12)
            TA = AtlasConjTrans;
         else
      #endif
      TA = (bv&4) ? AtlasTrans:AtlasNoTrans;
      void (*utrsm)(sminfo_t*,const enum ATLAS_DIAG,ATL_CINT N, ATL_CINT R,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda,TYPE *X, ATL_CSZT ldx,
         TYPE *diag, TYPE *L, TYPE *RW, TYPE *w);
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
/*
 *    ST case
 */
      if (NT <= si.kb)
      {
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
         return(0);
      }
/*
 *    outer loop
 */
      else
      {
         if (bv&1)
            Mjoin(PATL,ipnekInfo)(&ip, AtlasNoTrans, TA, NR, NT, NT, ldb, lda,
                                  ldb, NONE, alpha ); /* alpha=-1, beta=alpha*/
         else  /* Left */
            Mjoin(PATL,ipmekInfo)(&ip, TA, AtlasNoTrans, NT, NR, NT, lda, ldb,
                                  ldb, NONE, alpha); /* alpha=-1, beta=alpha*/
         #ifdef DEBUG
            fprintf(stderr, "MB=%d, NB=%d, KB=%d\n", ip.mb, ip.nb, ip.kb);
         #endif

         vp = ATL_trsm_alloc(&ip, &si, NT, &Tsz, &Rsz, &diag, &L, &R, &w);
         if (!vp)
            return(1); /* not enough memory, reucrse again  */
         int (*ammtrsm)(ipinfo_t *ip, sminfo_t *si,const enum ATLAS_DIAG Diag,
               ATL_CSZT R,
               ATL_CSZT N, const SCALAR alpha, const TYPE *A, ATL_CSZT lda,
               TYPE *X, ATL_CSZT ldx, ATL_CSZT Tsz, ATL_CSZT Rsz, TYPE *diag,
               TYPE *L, TYPE *RW, TYPE *w);

         ctu = (bv>>1)&7; /*masked three bits responsible for Conj Trans Upper*/
         if (!ctu || ctu==3 || ctu == 7)
            ammtrsm = ATL_trsm_LNUT;
         else
            ammtrsm = ATL_trsm_LTUN;
/*
 *       apply outer-loop trsm
 */
         if (bv & 1)
            ammtrsm(&ip, &si, (bv&16)?AtlasNonUnit:AtlasUnit, NR, NT, alpha, A,
                  lda, B, ldb, Tsz, Rsz, diag, L, R, w);
         else
            ammtrsm(&ip, &si, (bv&16)?AtlasNonUnit:AtlasUnit, NT, NR, alpha, A,
                  lda, B, ldb, Tsz, Rsz, diag, L, R, w);
         free(vp);
      }
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
