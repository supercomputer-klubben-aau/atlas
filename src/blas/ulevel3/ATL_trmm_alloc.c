/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Majedul Sujon
 * Code contributers : Majedul Sujon, Rakib Hasan
 */
#include "atlas_misc.h"
#include "atlas_amm.h"

#ifdef Right_
   #define ATL_trmm_alloc Mjoin(PATL,trmmR_alloc)
#else
   #define ATL_trmm_alloc Mjoin(PATL,trmmL_alloc)
#endif

/* #define MEM_DEBUG 1*/
void* ATL_trmm_alloc
(
   const ipinfo_t *ip,  /* gemm's ipinfo */
   tminfo_t *si,       /* trmm's ipinfo */
   ATL_CSZT N,
   ATL_SZT *Tsz,
   ATL_SZT *Rsz,
   ATL_SZT *Csz,
   TYPE **bw,  /* workspace for TRMM: single B-block and result (C) in TRMM*/
   TYPE **L,   /* whole A-matrix */
   TYPE **R,   /* row/col panel of B (shared with TRMM if possible) */
   TYPE **w    /* workspace for the result of single blk (C in gemm)*/
)
{
   void *vp;
   size_t szL, szR, szW;
   const int mu = ip->mu, nu = ip->nu, ku = ip->ku;
   const int mb = ip->mb;
   const int nb = ip->nb;
   const int vlen = ip->vlen;
   const size_t nmb = (N + mb - 1) / mb;
   const size_t nnb = (N + nb-1) / nb;
   const unsigned int szA = ip->szA, szB = ip->szB, szC = ip->szC;

   const int tvlen = si->vlen;
   const int tmu = si->mu;
   const int tnu = si->nu;
   const int tku = si->ku;
   const int tnmu = (mb+tmu-1)/tmu;
   const int tnnu = (nb+tnu-1)/tnu;

   int exsz;

   #ifdef Right_
   const int tnku = (nb+tku-1)/tku;
/*
 *    storage calculation for TRMM
 */
      *Tsz = nb*nb; /* space for triangle A block in TRMM */
/*
 *    NOTE: Tsz needs to be multiple of gemm's vlen
 */
      *Tsz = ((*Tsz + vlen-1)/vlen)*vlen; /* make it mult of vlen */
      *Rsz = tmu*tnmu*tku*tnku; /* space for a B-block in TRMM */
      if (si->flg&1) /* no space needed for c in trmm if shared with gemm */
         *Csz = 0;
      else
         *Csz = (((tmu*tnu+tvlen-1)/tvlen)*tvlen)*tnmu*tnnu;
/*
 *    bw : TRMM: B-copy +  result of C
 *    L:   Whole A-matrix minus diagonal blocks + space for diagoal blks in TRMM
 *         = nb*nb*(nnb-1)*(nnb>>1) + (*Tsz)*nnb
 *    R:   B col-panel = mbnb*nnb
 *    w:   workspace for C-result of single blk in GEMM = mbnb
 */
      szL = szB*(((nnb)*(nnb-1))>>1) + (*Tsz)*nnb;
      szR = szA*nnb;
      szW = Mmax(szC,mb*nb);
      exsz = tnu*tku + Mmax(nu*ku,tnu*tku) + (mu+mu)*nu;
   #else
   const int tnku = (mb+tku-1)/tku;
      *Tsz = mb*mb; /* workspace for triangle A block in TRMM */
      *Tsz = ((*Tsz + vlen-1)/vlen)*vlen; /* make it mult of vlen */
      *Rsz = tnu*tnnu*tku*tnku; /* workspace for a B-block in TRMM*/
      if (si->flg&1) /* no space needed for c in trmm if shared with gemm */
         *Csz = 0;
      else
         *Csz = (((tmu*tnu+tvlen-1)/tvlen)*tvlen)*tnmu*tnnu; /* size for C in TRMM*/
/*
 *    bw : TRMM: B-copy +  result of C
 *    L:   Whole A-matrix minus diagonal blocks + space for diagoal blks in TRMM
 *         = mb*mb*(nmb-1)*(nmb>>1) + (*Tsz)*nmb
 *    R:   B col-panel = mbnb*nmb
 *    w:   workspace for C-result of single blk in GEMM = mbnb
 */
      szL = (szA*(((nmb-1)*nmb)>>1)) + (*Tsz)*nmb;
      szR = szB*nmb;
      szW = Mmax(szC,mb*nb);
      exsz = tmu*tku + Mmax(mu*ku,tmu*tku) + (mu+mu)*nu;
   #endif
   vp = malloc( ATL_MulBySize( (*Rsz) + (*Csz) + szL + szR + szW + exsz )
                  + 4*ATL_Cachelen );
   if (!vp) return(NULL); /* not enough memory */

   *bw = ATL_AlignPtr(vp);
   *L = (*bw) + (((*Rsz)+(*Csz)) SHIFT);
   *L = ATL_AlignPtr(*L);
   *R = (*L) + (szL SHIFT);
   *R = ATL_AlignPtr(*R);
   *w = (*R) + (szR SHIFT);
   *w = ATL_AlignPtr(*w);
   return(vp);
}
