/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lapack.h"
#include "atlas_amm.h"

#ifndef UAMMID
   #define UAMMID 0
#endif
#define PRID Mjoin(Mjoin(atlas_,PRE),Mjoin(Mjoin(u,UAMMID),amm_))
#include Mstr(Mjoin(PRID,blk.h))
void Mjoin(PATL,laswp_amm)
(
   ATL_CINT N,   /* row-length to swap */
   TYPE *A,      /* column-major matrix holding top block */
   ATL_CINT lda, /* stride between row elts in A */
   ATL_CINT K1,  /* First elt of ipiv for which a swap will be done */
   ATL_CINT K2,  /* Last elt of ipiv for which a row interchange will be done */
   int *ipiv,    /* Vector of pivot indices.  Only K1-K2 are accessed. */
   ATL_CINT inci,/* piv stride; If inci<0, pivs applied in reverse order */
   ATL_CUINT nb, /* nb */
   ATL_CUINT mu, /* M unrolling used by amm kernel */
   ATL_CUINT nu, /* N unrolling used by amm kernel */
   ATL_CUINT bb, /* # of ipiv bits encoding major block number (mbn) */
   ATL_CUINT sb, /* # of ipiv bits encoding subblock number (sbn) */
   ATL_CUINT rb, /* # of ipiv bits encoding subblock row (sbr) */
   TYPE *b,      /* block-major panel wt trailing blocks */
   TYPE *bi      /* cplx only: imaginary panel */
)
{
   int i, i1, i2, KeepOn;
   ATL_CUINT bmsk=(1<<bb)-1, smsk=(1<<sb)-1, bb_sb=bb+sb;
   ATL_CUINT nbnb=nb*nb, munu=mu*nu;
   #ifdef TCPLX
      ATL_CUINT lda2 = lda+lda;
   #endif
   if (K2 < K1)
      return;
   if (inci >= 0)
   {
      ipiv += K1*inci;
      i1 = K1;
      i2 = K2-1;
   }
   else
   {
      ipiv -= (K2-1) * inci;
      i1 = K2-1;
      i2 = K1;
   }
   i = i1;
   do
   {
     ATL_UINT mbn, sbn, sbr, k, j;
     k = *ipiv;
     mbn = k & bmsk;
     sbn = (k>>bb) & smsk;
     sbr = (k>>bb_sb);
     k = sbr + sbn*mu + mbn*nb;
     if (k < K2)
        Mjoin(PATL,swap)(N, A+(i SHIFT), lda, A+(k SHIFT), lda);
     else
     {
        #ifdef TCPLX
           k = mbn*nbnb+sbn*munu + sbr;
           k = (mbn*nbnb)<<1;
           j = sbn*munu + sbr;
           Mjoin(PATLU,swap)(N, A+(i SHIFT), lda2, A+k+j, lda);
           Mjoin(PATLU,swap)(N, A+(i SHIFT)+1, lda2, A+k+j+nbnb, lda);
        #else
           k = mbn*nbnb+sbn*munu + sbr;
           Mjoin(PATL,swap)(N, A+i, lda, b+k, mu);
        #endif
     }
     ipiv += inci;
     if (inci >=0)
        KeepOn = (++i <= i2);
     else
        KeepOn = (--i >= i2);
   }
   while(KeepOn);
}
