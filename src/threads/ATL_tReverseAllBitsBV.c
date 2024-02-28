/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
void ATL_tReverseAllBitsBV(void *vp, unsigned int flg)
{
   if (vp)
   {
      long *gp=ATL_AlignSafeLS(vp), *bv;
      void *lk;
      const unsigned long P=gp[ATL_TBV_P], sumsz=gp[ATL_TBV_SUMSZ],
         lcksz=gp[ATL_TBV_LCKSZ], bvsz=gp[ATL_TBV_BVSZ],
         nlg=gp[ATL_TBV_NLGB], B=gp[ATL_TBV_B];
      long b;
      int i;
      lk = gp + sumsz;
      bv = lk + lcksz*P;
      for (i=0; i < P; i++)
      {
         long b = (i < nlg) ? b+1 : b;
         if (flg&1)
            ATL_lock(lk);
         bv[0] = b - bv[0];
         ATL_ReverseAllBitsBV(bv+1);
         if (flg&1)
            ATL_unlock(lk);
         bv += bvsz;
         lk += lcksz;
      }
   }
}
