#include "atlas_atmctr.h"
void  ATL_atmctr_free(void *ac)
{
   long *lp=ac, off=lp[1];
   #if !ATL_ATM_ASM
      ATL_lock_destroy(ATL_IncBySafeLS(ac));
   #endif
   ac = ATL_AddBytesPtr(ac, lp[1]);
   free(ac);
}
