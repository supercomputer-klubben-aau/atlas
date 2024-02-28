/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
/*
 * SYMM involves a dense matrix (A) and a symmetric matrix (B).
 * This routine uses normal GEMM to do a SYMM by reflecting the symmetric matrix
 * into dense access-major storage.  The recursive symm divides B
 * which means eventually the shape is dominated by A, and
 * so the double copy of A inherent in recursion becomes a dominant cost.
 * This routine avoids this problem by copying B up front, and then putting
 * the reflected B as inner matrix, so that A is copied only once, 1 panel at
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
 * These required for testing only
 */
/* #define DEBUG 1 */
#ifdef DEBUG
   #include "atlas_bitvec.h"
   #include "atlas_tst.h"
#endif
static void ATL_cpsy2amm
(
   ipinfo_t *ip,
   ATL_CUINT flg, /* bitvec: 0:Upper?; 1:HEMM? */
   const TYPE *S, ATL_iptr_t lds,        /* symmetric in lower storage */
   TYPE *W,                              /* nbxkb workspace */
   TYPE *w                               /* spc for entire S in block storage */
)
/*
 * w = alpha * S, w is in access-major storage, S is upper or lower col-major
 * NOTE: this is most general routine, where the nbxkb W allows us to reflect
 *       any block containing the diagonal.  For the case where MB & NB
 *       are multiples, we could optimize this by not doing the intermediate
 *       copy to W, but since SYMM is not that important, we have not bothered
 *       to write this special case.
 */
{
   ATL_iptr_t j;
   cm2am_t cpN, cpT;
   const ATL_iptr_t nnblks = ip->nfnblks + ip->npnblks, nkblks = ip->nfkblks+1;
   ATL_iptr_t jB, kB; /* block counters */
   void (*cpyBlk)
      (ipinfo_t *ip, ATL_CUINT, cm2am_t cpN, cm2am_t cpU, ATL_iptr_t KBLK,
       ATL_iptr_t JBLK, const TYPE *B, TYPE *W, TYPE *D);

   cpN = Mjoin(PATL,ipsyGetCopyB)(ip, flg, &cpT);
   #ifdef TCPLX
      if (flg&1)
         cpyBlk = (flg&2) ? Mjoin(PATL,heRU2ipBlk) : Mjoin(PATL,syRU2ipBlk);
      else
         cpyBlk = (flg&2) ? Mjoin(PATL,heRL2ipBlk) : Mjoin(PATL,syRL2ipBlk);
   #else
      cpyBlk = (flg&1) ? Mjoin(PATL,syRU2ipBlk) : Mjoin(PATL,syRL2ipBlk);
   #endif
   for (jB=0; jB < nnblks; jB++)
      for (kB=0; kB < nkblks; kB++)
         cpyBlk(ip, 1, cpN, cpT, kB, jB, S, w, W);
}

int Mjoin(PATL,ipsymmR)
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
   unsigned int nb, extra;
   void *vp=NULL;
   TYPE *wd, *wA; /* work for diagonal blks, work for copied A */
   TYPE *wB, *wC; /* work for B & C matrices; both aliased with wd */
   #ifdef DEBUG
      ATL_BV_t *bv;
   #endif
   Mjoin(PATL,ipgenInfo)(&ip, 0, AtlasNoTrans, AtlasNoTrans, M, N, N, ldg, lds,
                         ldc, alpha, beta);
/*
 * This routine attempts to malloc space for 1 block of C, 1 column panel of A,
 * and entire matrix B.  Also needs space for kb*nb workspace, which is
 * overlapped with all non-B workspace.
 */
   ATL_assert(ip.npmblks < 2);
   nb = (ip.nfnblks) ? ip.nb : ip.pnb;
   sz = (ip.nfkblks) ? ip.kb : ip.kb0;
   tkblks = ip.nfkblks+1;
   szB = tkblks*(ip.nfnblks*ip.szB + ip.npnblks*ip.pszB);
   szA = (ip.nfmblks) ? ip.szA : ip.pszA;
   szA *= tkblks;
   extra = ip.exsz;
   sz *= nb;
   if (sz > extra+szA+ip.szC)
      extra += sz - extra - szA - ip.szC;
   sz = ATL_MulBySize(extra + szA + szB + ip.szC) + 3*ATL_Cachelen;
   if (sz <= ATL_MaxMalloc)
      vp = malloc(sz);
   if (!vp)
      return(1);
    wB = ATL_AlignPtr(vp);
    wA = wB + (szB SHIFT);
    wd = wA = ATL_AlignPtr(wA);
    wC = wA + (szA SHIFT);
    wC = ATL_AlignPtr(wC);
    #ifdef TCPLX
       rC = wC + ip.szC;
    #endif
   ATL_cpsy2amm(&ip, flg, S, lds, wd, wB);
   #if defined(DEBUG) && 1
      #ifdef TCPLX
         bv = Mjoin(PATL,ipsycmpBV)(2, 0.0, &ip, 2|((flg&1)<<2)|((flg&2)<<2),
                                    wB, N, S);
      #else
         bv = Mjoin(PATL,ipsycmpBV)(2, 0.0, &ip, 2|((flg&1)<<2), wB, N, S);
      #endif
      if (bv)
      {
         ATL_print2dBV(N, N, bv);
         free(bv);
         printf("\nDONE PRINt2DBV\n\n");
      }
      else
         printf("\nNO ERRORS FOUND IN SYCP!\n\n");
   #endif
   #if 0
      Mjoin(PATL,ipprint)(stdout, &ip, 2, "SY", N, N, wB);
   #endif
   Mjoin(PATL,iploopsMNK)(&ip, 0, 0, G, NULL, C, 1|2|8, wA, wB, rC, wC, beta,
                          ip.blk2c);
   free(vp);
   return(0);
}
#undef ONE
