/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lapack.h"
/*
 * Takes ipiv split in (mbn,sbn,sbr) bitpattern, and translates it back
 * to standar getrf row coordinates.  See ATL_blkIPiv_amm for more detail
 */
void ATL_unblkIpiv_amm
(
   int n,      /* number of pivot entries to unblock */
   int *ipiv,  /* pivot array to encode wt block bit patterns */
   const int nb,
   const int mu,
   const int bb,
   const int sb,
   const int rb
)
{
   unsigned int k;
   const unsigned int bmsk=(1<<bb)-1, smsk=(1<<sb)-1, bb_sb=bb+sb;
   for (k=0; k < n; k++)
   {
      const unsigned int i = ipiv[k];
      unsigned int mbn, sbn, sbr;
      mbn = i & bmsk;
      sbn = (i>>bb) & smsk;
      sbr = (i>>bb_sb);
      ipiv[k] = sbr + sbn*mu + mbn*nb;
/*printf("decode: (%d,%d,%d) = %d\n", mbn, sbn, sbr, ipiv[k]);*/
   }
}

