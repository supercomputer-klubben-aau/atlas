#include "atlas_atmctr.h"
#if !ATL_ATM_ASM
long ATL_atmctr_dec(void *ac)
{  /* RETURNS: old value of count */
   void *lck=ATL_IncBySafeLS(ac);
   long ret, *lp=ac;
   ATL_lock(lck);
   ret = *lp;
   if (ret > 0)
      *lp = ret-1;
   ATL_unlock(lck);
   return(ret);
}
#endif
