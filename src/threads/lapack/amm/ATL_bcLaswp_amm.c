/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_pca.h"
#include "atlas_bcamm.h"
void Mjoin(PATL,bcLaswp_amm)
(
   ATL_bcpiv_t *bp,
   ATL_CINT N,     /* row-length to swap */
   ATL_CINT coff,  /* column offset to enable partial swap */
   TYPE *A,        /* column-major matrix holding top block */
   ATL_CINT lda,   /* stride between row elts in A */
   ATL_CINT K1,    /* First elt of ipiv for which a swap will be done */
   ATL_CINT K2,    /* Last elt of ipiv for which a row intrchg will be done */
   ATL_CINT lpj,   /* local panel no. to apply pivots on */
   ATL_CINT nnu,   /* no. of nus in current panel */
   ATL_CINT ZOFF   /* offset for complex block from real */
)
{
   int i, i1, i2, KeepOn, ii;
   void **larrs = bp->larrs;
   int* lldps = bp->lldps;
   ATL_UINT *ipiv = bp->ipiv;
   ATL_CINT inci = bp->inci;
   ATL_CUINT P=bp->R, B=bp->B, MU=bp->MU;
   ATL_CUINT nbRNK=bp->nbRNK, nbSBN=bp->nbSBN, nbSBR=bp->nbSBR;
   ATL_CUINT BB=bp->neMB, PB = P*B, NN=nnu*bp->neSB;
   ATL_CUINT mskSBR=(1<<nbSBR)-1, mskRNK=(1<<nbRNK)-1, mskSBN=(1<<nbSBN)-1;
   ATL_CUINT cmaj_off = coff*MU;
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
   A += coff*(lda SHIFT);  /* adjust for partial block swap */
   i = i1;
   ii = 0;
   do
   {
      ATL_UINT rank, sbn, sbr, lbn = ipiv[0];
      ATL_SZT k;
/*
 *    Decode ipiv entry
 */
      rank = lbn & mskRNK;   /* rank of row owner */
      lbn >>= nbRNK;
      sbr = lbn & mskSBR;    /* subblock row */
      lbn >>= nbSBR;
      sbn = lbn & mskSBN;    /* subblock number */
      lbn >>= nbSBN;         /* local block number */
      k = lbn*PB + rank*B + sbn*MU+sbr;  /* indxl2g wt psrc=0 */
      if (k < K2)
         Mjoin(PATL,swap)(N, A+(ii SHIFT), lda, A+((k-K1) SHIFT), lda);
      else  /* (rank,lbn,sbn,sbr) */
      {
         TYPE *w;
         k = ((lbn*BB) SHIFT)+sbn*NN+sbr;
         w = larrs[rank];
         w += ((lpj*lldps[rank]) SHIFT);  /* move to the correct panel */
         w += cmaj_off;                   /* adjust for partial block */
         w += k;
         #ifdef TCPLX
            Mjoin(PATL,swap_cplx2real)(N, A+ii+ii, lda, w, MU, w+ZOFF, MU);
         #else
            Mjoin(PATL,swap)(N, A+ii, lda, w, MU);
         #endif
      }
      ipiv += inci;
      ii++;
      if (inci >=0)
         KeepOn = (++i <= i2);
      else
         KeepOn = (--i >= i2);
   }
   while(KeepOn);
}
