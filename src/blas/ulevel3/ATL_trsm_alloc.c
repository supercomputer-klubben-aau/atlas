/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Majedul Sujon
 * Code contributers : Majedul Sujon, Rakib Hasan
 */
#include "atlas_misc.h"
#include "atlas_amm.h"

#ifdef Right_
   #define ATL_trsm_alloc Mjoin(PATL,trsmR_alloc)
#else
   #define ATL_trsm_alloc Mjoin(PATL,trsmL_alloc)
#endif

/* #define DEBUG 1*/

void* ATL_trsm_alloc
(
   const ipinfo_t *ip,  /* ipinfo for gemm */
   const sminfo_t *si,  /* ipinfo for trsm */
   ATL_CSZT N,          /* number of triangles */
   ATL_SZT *Tsz,        /* size of triangle blk used in utrsm */
   ATL_SZT *Rsz,        /* R size using in utrsm */
   TYPE **Diag,         /* workspace for diagonals */
   TYPE **L,            /* workspace for A: both gemm and trsm */
   TYPE **R,            /* workspace for B rowspan or colspan */
   TYPE **w             /* workspace for C block in gemm, utrsm: R & w*/
)
{
   void *vp = NULL;
   ATL_SZT sz;
   ATL_SZT szL, szR, szW, exsz;
   ATL_CSZT szD = N;

   const int tmu = si->mu;
   const int tnu = si->nu;
   const int tku = si->ku;
   const int mb = ip->mb;
   const int nb = ip->nb;
   const int mu = ip->mu;
   const int nu = ip->nu;
   const int ku = ip->ku;
   const unsigned int szA = ip->szA;
   const unsigned int szB = ip->szB;
   const unsigned int szC = ip->szC;

#ifdef Right_
   /*const size_t nbnb = ((nb * nb + ip->vlen-1)/ip->vlen)*ip->vlen;*/
   const int nnb = (N + nb - 1) / nb;
   const int tnnu = (nb + tnu - 1) / tnu;
/*
 * size calculated for utrsm
 * NOTE: even though Tsz is used for TRSM, we need to make it multiple of
 * gemm's vlen so that the starting address of the workspace of A in gemm
 * can be aligned to vector length
 */
   *Tsz = tnu*tnu*tnnu*(tnnu+1)>>1;
   *Tsz = ((*Tsz + ip->vlen-1)/ip->vlen)*ip->vlen; /* make it mult of vlen */
   *Rsz = tmu*tnu*tnnu;
   szL = (szB*nnb*(nnb-1)>>1) + (*Tsz)*nnb;
   szR = szA*nnb;
   szW = (Mmax(*Rsz,szC))<<1;
   exsz = tnu*tku + Mmax(mu*ku,tmu*tku) + (mu+mu)*nu;
#else
   const int nmb = (N + mb - 1) / mb;
   const int tnmu = (mb + tmu - 1) / tmu;

   *Tsz = tmu*tmu * (tnmu*(tnmu+1))>>1;
   *Tsz = ((*Tsz + ip->vlen-1)/ip->vlen)*ip->vlen; /* make it mult of vlen */
   *Rsz = tmu*tnu * tnmu;
   szL = ((szA*nmb*(nmb-1))>>1) + (*Tsz)*nmb;
   szR = szB*nmb;
   szW = (Mmax(*Rsz,szC))<<1;
   exsz = tmu*tku + Mmax(mu*ku,tmu*tku) + (mu+mu)*nu;
#endif
   sz = ATL_MulBySize( szD + szL + szR + szW + exsz) + 5*ATL_Cachelen;
   if (sz <= ATL_MaxMalloc)
      vp = malloc(sz);
   if (!vp) return(NULL); /* not enough memory */
   *Diag = vp;
   *L = (*Diag) + (szD SHIFT);
   *L = ATL_AlignPtr(*L);
   *R = (*L) + (szL SHIFT);
   *R = ATL_AlignPtr(*R);
   *w = (*R) + (szR SHIFT);
   *w = ATL_AlignPtr(*w);
   return(vp);
}
