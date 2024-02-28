/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
void ATL_InitBV
(
   unsigned long nbits,             /* max # of bits BV needs to store */
   ATL_BV_t *bv,
   ATL_BV_t msk
)
/*
 * Initializes bv, with each full element of bv getting msk.  Last partial
 * block gets lower part of msk.
 */
{
   if (nbits)
   {
      const ATL_BV_t n=nbits>>shBV, nr=nbits&modmskBV;
      ATL_BV_t i;

      ATL_assert(bv);
      bv[0] = nbits;
      for(bv++,i=0; i < n; i++)
         bv[i] = msk;
      if (nr)
         bv[i] = ((1L<<nr)-1) & msk;
   }
   else if (bv)
      bv[0] = 0;
}
