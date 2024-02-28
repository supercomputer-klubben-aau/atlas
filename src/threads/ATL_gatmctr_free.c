/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_gatmctr.h"
void ATL_gatmctr_free(void *ac)
{
   if (ac)
   {
      long *lp = ac;
      const long flg = lp[ATL_GAC_FLG];

      if (GAC_FLAG_SET(flg, ATL_GAC_NOPRV))
      {
         const long inc=lp[3];
         const unsigned int P=lp[ATL_GAC_P];
         unsigned int i;
         ac = ATL_AddBytesPtr(ac, ATL_SAFELS+lp[2]);
         for (i=0; i < P+1; i++, ac = ATL_AddBytesPtr(ac, inc))
            ATL_atmctr_destroy(ac);
         ac = ATL_SubBytesPtr(lp, lp[4]);
      }
      else if (GAC_FLAG_SET(flg, ATL_GAC_NOPUB))
         ac = ATL_SubBytesPtr(ac, lp[2]);
      else
      {
         const unsigned int P=lp[ATL_GAC_P];
         unsigned int i;
         const long inc = lp[ATL_GAC_INC];
         void *lck;

         lck = ATL_AddBytesPtr(lp, (ATL_SAFELS<<2));
         for (i=0; i < P; i++, lck = ATL_AddBytesPtr(lck, inc))
            ATL_lock_destroy(lck);
         ac = ATL_SubBytesPtr(ac, lp[4]);
      }
      free(ac);
   }
}

