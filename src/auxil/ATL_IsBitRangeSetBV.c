/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
/*
 * RETURNS: 1 if all bits in range [b0,bn] are set, 0 otherwise
 */
int ATL_IsBitRangeSetBV(ATL_BV_t *bv, unsigned long b0, unsigned long bN)
{
/*
 * False if range includes bits that don't exist, or bad range
 */
   if (bv && bN < bv[0] && b0 <= bN)
   {
      const ATL_BV_t n=bN>>shBV, nr=(bN+1)&modmskBV, fskip=b0&modmskBV;
      const ATL_BV_t allset=allsetBV;
      b0 >>= shBV;
      bv++;
      if (fskip)  /* partial first result */
      {
         ATL_BV_t v;
         v = bv[b0] | ((1L<<fskip)-1);  /* set 1st fskip bits */
         if (b0 == n && nr) /* 1st is also last with remainder */
         {
            v |= ~(allset>>(bpiBV-nr));
            return(v == allset);
         }
         if (v != allset)
            return(0);
         b0++;
      }
      while (b0 < n)  /* handle full blocks */
      {
         if (bv[b0] != allset)
            return(0);
         b0++;
      }
      if (nr)  /* partial result */
      {
         ATL_BV_t v;
         v = bv[b0] | (~(allset>>(bpiBV-nr)));
         return(v == allset);
      }
   }
   return(0);
}
