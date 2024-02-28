/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
long ATL_tFindUnsetBitBV(void *vp, unsigned int rank)
/*
 * RETURNS: first (unsafe) unset bit, starting from rank's portion
 */
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
      bv0 = gp + sumsz + lcksz*P;
      bv  = bv0  + rank*bvsz;
      vrk = rank;
      for (i=0; i < P; i++)
      {
         if (bv[0])  /* found a bv with unset bits */
         {
            long ret;
            bv++;
            ret = ATL_FindFirstUnsetBitBV(bv, 0);
            #ifdef DEBUG
               ATL_assert(ret >= 0);
            #endif
            if (ret >= 0)
            {
               long nlrg, nsm, b;
               b = gp[ATL_TBV_B];
               nlrg = gp[ATL_TBV_NLGB];
               nlrg = Mmin(nlrg, vrk);
               nsm = vrk - nlrg;
               ret += nlrg*(b+1) + nsm*b;
               return(ret);
            }
         }
         if (++vrk != P)
            bv  += bvsz;
         else
         {
            bv = bv0;
            vrk = 0;
         }
      }
   }
   return(-1);
}
