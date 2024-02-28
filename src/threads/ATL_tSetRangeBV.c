/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
unsigned long ATL_tSetRangeBV(void *bv, unsigned int *NBITS, unsigned long pos,
                              unsigned long setmsk)
/*
 * On input *nbits is the max number of bits in bv to set to the bit pattern
 * in the least significant bits of mask.  On exit, the number of bits actually
 * set (may be less if some bits owned by another mutex; only one lock will
 * be made).
 * RETURNED: original bitpattern, *nbits set to the number of bits
 */
{
   if (bv)
   {
      void *lck;
      ATL_BV_t *lbv;
      unsigned long *gp=ATL_AlignSafeLS(bv);
      unsigned long nbits=(*NBITS), b=gp[7], ornk, elt, lret, nn, i, v, msk, n2;
      long nch;
      const unsigned long P=gp[0], gnbits=gp[1], sumSz=gp[2], lckSz=gp[3],
         bvSz=gp[4], nlrg=gp[5], nsm=gp[6], lrgN=(b+1)*nlrg;

      if (pos >= gnbits)
      {
         *NBITS = 0;
         return(0);
      }
      ornk  = gnbits - pos;
      nbits = Mmin(nbits, ornk);  /* replace with local comp? */
      if (lrgN > pos) /* contained in a large block */
      {
         ornk = pos/(++b);
         pos -= ornk*b;      /* local position */
      }
      else
      {
         pos -= lrgN;
         ornk = pos/b;
         pos -= ornk*b;  /* local position */
         ornk += nlrg;
      }
      v = b - pos;           /* # of bits from pos prot by this lock */
      nbits = Mmin(nbits, v);/* max bits: number prot by this lock */
      lck = gp + sumSz + ornk*lckSz;
      lbv = gp + sumSz + P*lckSz + ornk*bvSz;
      elt = pos >> shBV;     /* local elt */
      pos &= modmskBV;       /* pos within elt */
      nn = bpiBV-pos;        /* bits left in this elt */
      nn = Mmin(nn, nbits);  /* # of bits from this elt */
      msk = (1L<<nn)-1;      /* only lower nn bits set */

      ATL_lock(lck);
      v = lbv[2+elt];
      lret = (v>>pos)&msk;   /* lret holds original low nn bits */
      msk = ~(msk<<pos);     /* all set except for target range */
      lbv[2+elt] = (v&msk)|(setmsk<<pos);
      n2 = nbits - nn;       /* # of bits in next element */
      if (n2)
      {
         msk = (1L<<n2)-1;  /* low n2 bits set */
         v = lbv[3+elt];
         lret |= (v&msk)<<nn;
         lbv[3+elt] = (v&(~msk))|(setmsk>>nn);
      }
/*
 *    Compute change in unset
 */
      for (nch=i=0; i < nbits; i++)
         nch += ((lret>>i)&1) - ((setmsk>>i)&1);
      lbv[0] += nch;
      ATL_unlock(lck);

      *NBITS = nbits;
      return(lret);
   }
   *NBITS = 0;
   return(0);
}
