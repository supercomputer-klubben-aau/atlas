/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
unsigned long ATL_tGetTotBitsBV(void *bv)
{
   if (bv)
   {
      unsigned long *lp = ATL_AlignSafeLS(bv);
      return(lp[1]);
   }
   return(0);
}
