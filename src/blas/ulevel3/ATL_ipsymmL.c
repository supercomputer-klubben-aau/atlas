/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
/*
 * SYMM involves a symmetric matrix (A) and dense matrix (B).
 * This routine uses normal GEMM to do a SYMM by reflecting the symmetric matrix
 * into dense access-major storage.  The recursive symm divides A
 * which means eventually the shape is dominated by B, and
 * so the double copy of B inherent in recursion becomes a dominant cost.
 * This routine avoids this problem by copying A up front, and then putting
 * the dense A as inner matrix, so that B is copied only once, 1 panel at
 * a time.
 */
#include "atlas_misc.h"
#define ATL_GLOBIDX 1
#include "atlas_amm.h"
#undef ATL_GLOBIDX
#define ATL_WANT_ILCM 1
#include "atlas_iopt.h"
#undef ATL_WANT_ILCM
#include "atlas_level1.h"
#include "atlas_level2.h"
#include "atlas_level3.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))
/*
 * inline copy never really wins, and recursion beats symm_tN, so this should
 * always be 1
 */
#define SEP_CPY 1
/*
 * These required for testing only
 */
#define DEBUG 1
#ifdef DEBUG
   #include "atlas_bitvec.h"
   #include "atlas_tst.h"
#endif
static void ATL_cpsy2amm
(
   ipinfo_t *ip,
   ATL_CUINT flg, /* bitvec: 0:Upper?; 1:HEMM? */
   ATL_iptr_t N,                         /* order of A */
   const TYPE *S, ATL_iptr_t lds,        /* symmetric in lower storage */
   TYPE *W, ATL_iptr_t ldW,              /* ldWxkb workspace */
   TYPE *w                               /* spc for entire S in block storage */
)
/*
 * b = alpha * S, b is in access-major storage, S is upper or lower col-major
 * NOTE: this is most general routine, where the LxL workspace W is used
 *       to copy square superblocks from the original matrix.  This allows
 *       us to reflect symmetric matrices to make full, and follow the diagonal
 *       easily even if MB & NB are not multiples.  For the case where MB & NB
 *       are multiples (or equal), we can use less workspace & copies.
 */
{
   ATL_iptr_t j;
   cm2am_t cpN, cpT;
   const unsigned int mb=ip->mb, kb = ip->kb;
   const ATL_iptr_t nmblks = ip->nfmblks + ip->npmblks, nkblks = ip->nfkblks+1;
   ATL_iptr_t iB, kB; /* block counters */
   void (*cpyBlk)
      (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpU, ATL_iptr_t IBLK,
       ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D, ATL_iptr_t ldd);

   cpN = Mjoin(PATL,ipsyGetCopyA)(ip, flg, &cpT);
   #ifdef TCPLX
      if (flg&1)
         cpyBlk = (flg&2) ? Mjoin(PATL,heLU2ipBlk) : Mjoin(PATL,syLU2ipBlk);
      else
         cpyBlk = (flg&2) ? Mjoin(PATL,heLL2ipBlk) : Mjoin(PATL,syLL2ipBlk);
   #else
      cpyBlk = (flg&1) ? Mjoin(PATL,syLU2ipBlk) : Mjoin(PATL,syLL2ipBlk);
   #endif
   for (kB=0; kB < nkblks; kB++)
      for (iB=0; iB < nmblks; iB++)
         cpyBlk(ip, 1, cpN, cpT, iB, kB, S, w, W, ldW);
}

int Mjoin(PATL,ipsymmL)
(
   ATL_CUINT flg, /* bitvec 0:Upper? 1:HEMM? */
   ATL_CSZT  M,
   ATL_CSZT N,
   const SCALAR alpha,
   const TYPE *S,
   ATL_CSZT lds,
   const TYPE *G,
   ATL_CSZT ldg,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   #ifdef TCPLX
      TYPE *rC;
   #else
      #define rC wC
   #endif
   ipinfo_t ip;
   ATL_iptr_t sz, szA, szB, ldd, tkblks;
   unsigned int mb, extra;
   void *vp=NULL;
   TYPE *wd, *wA; /* work for diagonal blks, work for copied A */
   TYPE *wB, *wC; /* work for B & C matrices; both aliased with wd */
   #ifdef DEBUG
      ATL_BV_t *bv;
   #endif
   Mjoin(PATL,ipgenInfo)(&ip, 0, AtlasNoTrans, AtlasNoTrans, M, N, M, lds, ldg,
                         ldc, alpha, beta);
/*
 * This routine attempts to malloc space for 1 block of C, 1 column panel of B,
 * and entire matrix A.  Also needs space for kb*mb workspace, which is
 * cannot be overlapped with any workspace when doing A copy inline
 */
   ATL_assert(ip.npmblks < 2);
   mb = (ip.nfmblks) ? ip.mb : ip.pmb;
   sz = (ip.nfkblks) ? ip.kb : ip.kb0;
   if (mb <= M)
      ldd = mb;
   else
   {
      ldd = ATL_MulBySize(M) + ATL_Cachelen - 1;
      ldd = ATL_MulByCachelen(ATL_DivByCachelen(ldd));
      ldd = ATL_DivBySize(ldd);
   }
   tkblks = ip.nfkblks+1;
   szA = tkblks*(ip.nfmblks*ip.szA + ip.npmblks*ip.pszA);
   szB = (ip.nfnblks) ? ip.szB : ip.pszB;
   szB *= tkblks;
   sz *= ldd;
   #if SEP_CPY
      extra = ip.exsz;
      if (sz > extra+szB+ip.szC)
         extra += sz - extra - szA - ip.szC;
   #else
      extra = Mmax(sz, ip.exsz);
   #endif

   sz = ATL_MulBySize(extra + szA + szB + ip.szC) + 3*ATL_Cachelen;
   if (sz <= ATL_MaxMalloc)
      vp = malloc(sz);
   if (!vp)
      return(1);
    wA = ATL_AlignPtr(vp);
    wB = wA + (szA SHIFT);
    wB = ATL_AlignPtr(wB);
    wC = wB + (szB SHIFT);
    wC = ATL_AlignPtr(wC);
    #ifdef TCPLX
       rC = wC + ip.szC;
       #if !SEP_CPY
          wd = rC + ip.szC
       #endif
    #elif !SEP_CPY
       wd = wC + ip.szC;
    #endif
    #if SEP_CPY
       wd = wB;
    #endif
   #if SEP_CPY
      ATL_cpsy2amm(&ip, flg, N, S, lds, wd, ldd, wA);
   #else
      Mjoin(PATL,ipsymmL_tN_wrk)(&ip, flg+256, S, G, beta, C, wd, wA,wB,rC,wC);
   #endif
   #if defined(DEBUG) && 0
      #ifdef TCPLX
         bv = Mjoin(PATL,ipsycmpBV)(2, 0.0, &ip, 1|((flg&1)<<2)|((flg&2)<<2),
                                    wA, M, S);
      #else
         bv = Mjoin(PATL,ipsycmpBV)(2, 0.0, &ip, 1|((flg&1)<<2), wA, M, S);
      #endif
      if (bv)
      {
         ATL_print2dBV(M, M, bv);
         free(bv);
         printf("\nDONE PRINt2DBV\n\n");
      }
      else
         printf("\nNO ERRORS FOUND IN SYCP!\n\n");
   #endif
   #if 0
      Mjoin(PATL,ipprint)(stdout, &ip, 1, "SY", M, M, wA);
   #endif
   #if SEP_CPY
      Mjoin(PATL,iploopsNMK)(&ip, 0, 0, NULL, G, C, 1|2|4, wA, wB, rC, wC, beta,
                             ip.blk2c);
   #else
      Mjoin(PATL,iploopsNMK)(&ip, 0, 1, NULL, G, C, 1|2|4, wA, wB, rC, wC, beta,
                             ip.blk2c);
   #endif
   free(vp);
   return(0);
}
#undef ONE
