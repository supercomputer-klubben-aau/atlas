/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
long ATL_tGlb2locBV(ATL_BV_t *lbv, void *gbv, unsigned long pos)
/*
 * Combines all P individual global bitvecs into one local bitvec.
 * Mutexes are not locked, so info may be wrong.
 * RETURNS: number of unset bits in new local bitvec
 */
{
   long nunset = 0;
   if (gbv)
   {
      long *gp=ATL_AlignSafeLS(gbv), *tp;
      const long P=gp[0], nbits=gp[1], inc=gp[4];
      long i, nb;

      ATL_assert(lbv);
      ATL_assert(lbv[0] >= nbits+pos);
      nb = pos;
      tp = gp + gp[2] + P*gp[3];
      for (i=0; i < P; i++)
      {
         nunset += tp[0];
         nb += ATL_IncorpBV(lbv, tp+1, nb);
         tp += inc;
      }
   }
   return(nunset);
}
