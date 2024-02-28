/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lapack.h"
static void GetPivBitCount
(
   int nb,   /* block factor to use */
   int mu,   /* mu used by amm kernel */
   int *BB,  /* number of bits to encode block number */
   int *SB,  /* number of bits to encode subblock number */
   int *RB   /* number of bits to encode row index in sub-block */
)
{
   int bb, sb, rb, nmu;


   if (mu > 1)
      for (rb=0; (1<<rb) < mu; rb++);
   else
      rb = 1;

   nmu = nb / mu;
   if (nmu > 1)
      for (sb=0; (1<<sb) < nmu; sb++);
   else
      sb = 1;

   bb = (sizeof(int)<<3) - rb - sb;
   *RB = rb;
   *SB = sb;
   *BB = bb;
}

/*
 * Split row index (I) into 3 bit patterns: (sbr,sbn,mbn) where:
 * + (subblock row, subblock number, main block number)
 * + For column major we compute:
 *     I = r + sbn*mu + bn*nb:
 * + For C-format we compute:
 *     panoff = bn * nbnb;
 *     blkoff = sbn*mu*nu;
 *     rowoff = r;
 * Bad news, is now both indices must be computed using * and +, but in
 * col-major encoding, C-format computation requires in-the-loop division,
 * so this seems better.
 *
 *       mbn = i / nb
 *       sbn = (i-mbn*nb) / mu;
 *       sbr = (i - mbn*nb - sbn*mu);
 *
 *    We will lose at most 2 bits of range, which should never be a problem
 *    in practice, since this algorithm will fail malloc long before we are
 *    using all but two bits of the integer range
 *
 *    Always adjusts nb entries starting from ipiv.
 */
void ATL_blkIpiv_amm
(
   int N,      /* # of ipiv entries to encode */
   int nb,
   int mu,
   int *ipiv,  /* INPUT/OUTPUT: pivot array to encode wt block bit patterns */
   int iadj,   /* offset to add to entries prior to encoding */
   int *BB,
   int *SB,
   int *RB
)
{
   unsigned int k, bb_sb;
   int bb, sb, rb;   /* bits to encode block #, subblk#, row in subblk */

   GetPivBitCount(nb, mu, BB, SB, RB);
   bb = *BB; sb = *SB; rb = *RB;
   bb_sb = bb + sb;

   for (k=0; k < N; k++)
   {
      unsigned int i, mbn, sbn, sbr;
      i = ipiv[k] + iadj;
      mbn = i / nb;
      i -= mbn*nb;
      sbn = i / mu;
      sbr = i - sbn*mu;
/*printf("encode: i=%d (%d+%d) = (%d,%d,%d)\n", i, ipiv[k], iadj, mbn, sbn, sbr);*/
      ipiv[k] = (sbr << bb_sb) | (sbn << bb) | mbn;
   }
}

