/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
void ATL_UnsetAllBitsBV(ATL_BV_t *bv)
{
   if (bv)
   {
      const ATL_BV_t nbits=bv[0], n=nbits>>shBV, nr=nbits-(n<<shBV);
      ATL_BV_t i;
      for(bv++,i=0; i < n; i++)
         bv[i] = 0;
      if (nr)
         bv[i] &= ~((1L<<nr)-1);
   }
}
