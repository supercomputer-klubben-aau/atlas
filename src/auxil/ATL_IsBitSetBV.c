/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
int ATL_IsBitSetBV(ATL_BV_t *bv, unsigned long pos)
{
   if (bv && bv[0] > pos)
   {
      const unsigned int iv=(pos>>shBV)+1, il=pos&modmskBV;
      return((bv[iv]>>il)&1L);
   }
   return(0);
}
