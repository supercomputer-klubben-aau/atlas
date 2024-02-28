/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
int ATL_UnsetBitBV(ATL_BV_t *bv, unsigned long pos)
{
   if (bv && bv[0] > pos)
   {
      ATL_BV_t v;
      const unsigned int iv=pos>>shBV, il=pos&modmskBV;
      v = bv[iv+1];
      bv[iv+1] = v & ~(1L<<il);
      return((v>>il)&1L);
   }
   return(-1);
}
