/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
long ATL_tSetUnsetBitBV(void *vp, unsigned int rank)
{
   if (vp)
   {
      long *gp=ATL_AlignSafeLS(vp), *bv, *bv0, *lck0, *lck;
      unsigned int vrk, i;
      const unsigned long P=gp[ATL_TBV_P], sumsz=gp[ATL_TBV_SUMSZ],
         lcksz=gp[ATL_TBV_LCKSZ], bvsz=gp[ATL_TBV_BVSZ];

      if (rank >= P)
         rank = rank - (rank / P)*P;
/*
 *    Make an unsafe pass looking for unset bits.  Note that the return of -1
 *    only means no bits possibly unset if this is used on a vector that
 *    starts all unset, and sets bits during execution (the usual case).
 *    If other threads can unset bits, we may return -1 with unset bits due
 *    to another thread unsetting a rank's bit after we check it.
 */
      lck0 = gp + sumsz;
      bv0 = lck0 + lcksz*P;
      lck = lck0 + rank*lcksz;
      bv  = bv0  + rank*bvsz;
      vrk = rank;
      for (i=0; i < P; i++)
      {
         if (bv[0])  /* found a bv with unset bits */
         {
            ATL_lock(lck);
            if (bv[0])
            {
               long ret, nlrg, nsm, b;
               (*bv)--;  /* subtract bit I'm setting from nunset */
               bv++;
               ret = ATL_FindFirstUnsetBitBV(bv, 0);
               #ifdef DEBUG
                  ATL_assert(ret >= 0);
               #endif
               ATL_SetBitBV(bv, ret);
               ATL_unlock(lck);

               b = gp[ATL_TBV_B];
               nlrg = gp[ATL_TBV_NLGB];
               nlrg = Mmin(nlrg, vrk);
               nsm = vrk - nlrg;
               ret += nlrg*(b+1) + nsm*b;
               return(ret);
            }
            ATL_unlock(lck);
         }
         if (++vrk != P)
         {
            bv  += bvsz;
            lck += lcksz;
         }
         else
         {
            bv = bv0;
            lck = lck0;
            vrk = 0;
         }
      }
   }
   return(-1);
}
