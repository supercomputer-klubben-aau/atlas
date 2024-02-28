/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
void ATL_tFreeBV(void *bv)
{
   if (bv)
   {
      unsigned long *lp = ATL_AlignSafeLS(bv);
      const unsigned long P=lp[0], sumSz=lp[2], lckSz=lp[3];
      unsigned int i;
      lp += sumSz;
      for (i=0; i < P; i++, lp += lckSz)
         ATL_lock_destroy(lp);
      free(bv);
   }
}
