/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
ATL_BV_t ATL_SetEltBV
   (ATL_BV_t *bv, unsigned long ielt, ATL_BV_t val)
{
   if (bv)
   {
      if ((ielt<<shBV) < bv[0])
      {
         ATL_BV_t ret=bv[ielt+1];
         bv[ielt+1] = val;
         return(ret);
      }
   }
   return(0);
}
