/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
int ATL_SetBitBV(ATL_BV_t *bv, unsigned long pos)
{
   if (bv && bv[0] > pos)
   {
      ATL_BV_t v;
      const unsigned long iv=(pos>>shBV)+1, il=pos&modmskBV;
      v = bv[iv];
      bv[iv] = v | (1L<<il);
      return((v>>il)&1);
   }
   return(-1);
}
