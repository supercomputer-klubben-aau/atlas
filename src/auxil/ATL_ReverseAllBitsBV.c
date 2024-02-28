/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
void ATL_ReverseAllBitsBV(ATL_BV_t *bv)
{
  ATL_BV_t n=bv[0], nl=n&modmskBV, i;
  bv++;
  n >>= shBV;
  for (i=0; i < n; i++)
     bv[i] = ~bv[i];
  if (nl)
  {
     ATL_BV_t msk = allsetBV>>(bpiBV-nl), v=bv[i];
     bv[i] = ((~v) & msk) | (v & ~msk);
  }
}
