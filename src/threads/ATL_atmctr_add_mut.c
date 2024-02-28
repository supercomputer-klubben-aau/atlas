#include "atlas_atmctr.h"
#if !ATL_ATM_ASM
long ATL_atmctr_add(void *ac, unsigned long val)
{  /* RETURNS: old value of count */
   void *lck=ATL_IncBySafeLS(ac);
   long ret, *lp=ac;
   ATL_lock(lck);
   ret = *lp;
   *lp = ret + val;
   ATL_unlock(lck);
   return(ret);
}
#endif
