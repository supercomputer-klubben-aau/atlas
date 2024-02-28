/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
unsigned long ATL_GetTotBitsBV(ATL_BV_t *bv)
{
   if (bv)
      return(bv[0]);
   return(0);
}
