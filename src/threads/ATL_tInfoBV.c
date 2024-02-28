/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
long ATL_tInfoBV(const void *vp, long what)
/*
 * Provides unsafe answers to global questions
 */
{
   if (vp)
   {
      const long *gp=ATL_AlignSafeLS(vp);
      if (what <= ATL_TBV_B)
         return(gp[what]);
      else if (what == ATL_TBV_NUNSET)
      {
         const unsigned long P=gp[ATL_TBV_P], sumsz=gp[ATL_TBV_SUMSZ],
               lcksz=gp[ATL_TBV_LCKSZ], bvsz=gp[ATL_TBV_BVSZ];
         const unsigned long *tp = gp + sumsz + lcksz*P;
         long nun=0;
         int i;
         for (i=0; i < P; i++, tp += bvsz)
            nun += *tp;
         return(nun);
      }
   }
   return(-1);
}
