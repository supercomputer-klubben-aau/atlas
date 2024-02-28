#include "atlas_atmctr.h"
void *ATL_atmctr_init(void *vp, long cnt)
{
   if (vp)
   {
      long *lp = ATL_AlignSafeLS(vp);
      void *lck;

      lp[0] = cnt;
      lp[1] = ((size_t)vp) - ((size_t)lp); /* will be negative or 0 */
      #if !ATL_ATM_ASM
         lck = ATL_IncBySafeLS(lp);
         ATL_lock_init(lck);
      #endif
      return(lp);
   }
   return(NULL);
}

void *ATL_atmctr_new(long cnt)
{
   void *vp, *lck;
   long *lp;
   #if ATL_ATM_ASM
      vp = malloc(ATL_SAFELS+ATL_SAFELS);
   #else
      vp = malloc(ATL_SAFELS+ATL_SAFELS+sizeof(ATL_lock_t));
   #endif
   ATL_assert(vp);
   return(ATL_atmctr_init(vp, cnt));
}
