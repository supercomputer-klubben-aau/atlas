/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
ATL_BV_t ATL_GetEltBV(ATL_BV_t *bv, unsigned long ielt)
{
   if (bv)
      if ((ielt<<shBV) < bv[0])
         return(bv[ielt+1]);
   return(0);
}
