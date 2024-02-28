#include "atlas_atmctr.h"
#if !ATL_ATM_ASM
long ATL_atmctr_set(void *ac, long val)
{  /* RETURNS: old value of count */
   void *lck=ATL_IncBySafeLS(ac);
   long ret, *lp=ac;
   ATL_lock(lck);
   ret = *lp;
   *lp = val;
   ATL_unlock(lck);
   return(ret);
}
#endif
