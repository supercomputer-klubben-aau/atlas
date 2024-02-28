/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
#include <errno.h>

/*
 * wrapper routine for all funcs that affect only one bit
 * flg: 1:IsBitSet, 2:SetBit, 4:UnsetBit, 512: unsafe (no lock/unlock)
 */
long ATL_tScopeBitBV(void *vp, unsigned long bit, unsigned int flg)
/*     0       1     2       3      4       5         6          7
 * 1. <P> <gnbits> <sumSz> <lckSz> <bvSz> <nLrgBlks> <nSmlBlks> <b>
 * 2. local lock storage area : *each lock* rounded up to 128 bytes
 * 3. locBV area: <nunset> [serial BV] -> <nleft> <nbits> <BV> => 128 rounded
 */
{
   if (vp)
   {
      unsigned long *lp = ATL_AlignSafeLS(vp), *lck, *bv;
      const unsigned long P=lp[0], ngbits=lp[1], sumSz=lp[2], lckSz=lp[3],
                          bvSz=lp[4], nlrg=lp[5], nsml=lp[6], B=lp[7];
      const unsigned long bigN=nlrg*(B+1);
      unsigned long idx;
      long obit, lbit;

      ATL_assert(bit < ngbits);
      if (bit > bigN)    /* have both large & small blocks */
      {
         lbit = (bit - bigN);
         idx = lbit / B;
         lbit -= idx*B;
         idx = nlrg + idx;
      }
      else if (bit == bigN)  /* have all large blocks */
      {
         lbit = 0;
         idx = nlrg;
      }
      else                  /* have only large blocks */
      {
         idx = bit/(B+1);
         lbit = bit - idx*(B+1);
      }
      lck = lp + sumSz + idx*lckSz;
      bv = lp + sumSz + lckSz*P + idx*bvSz;
      if (!(flg&512))                       /* lock unless unsafe bit (512) */
         ATL_assert(!ATL_lock(lck));        /* is set */

      if (flg&1)
         obit = ATL_IsBitSetBV(bv+1, lbit);
      else if (flg&2)
      {
         obit = ATL_SetBitBV(bv+1, lbit);
         if (!obit)
           bv[0]--;
      }
      else /* if (flg&4) */
      {
         obit = ATL_UnsetBitBV(bv+1, lbit);
         if (obit)
           bv[0]++;
      }

      if (!(flg&512))                     /* if we locked for safety */
         ATL_assert(!ATL_unlock(lck));    /* must unlock here */
      return(obit);                       /* return original value */
   }
   return(0);
}
