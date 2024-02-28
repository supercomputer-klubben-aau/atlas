#include "atlas_tsched_dly.h"
#ifdef ATL_AMM_NSCHED
unsigned int ATL_GetSchedDelayIdx(float dly)
/*
 * Searches ATL_AMM_NSCHED-entry array ATL_AMM_SCHEDdly for last entry <= dly.
 * Uses binary search with assumption time is strictly increasing, which may
 * not always be true due to timing variance & kernel choice, but it should
 * be good enough for estimation.  Don't want the overhead of linear search
 * everytime we want to estimate parallel schedule overhead!
 */
{
   unsigned int i0=0, dist=ATL_AMM_NSCHED>>1, iret;
   float mid;
/*
 * Early exit for out-of-range delay
 */
   if (ATL_AMM_NSCHED<2 || dly <= *ATL_AMM_SCHEDdly)
      return(0);
   else if (dly >= ATL_AMM_SCHEDdly[ATL_AMM_NSCHED-1])
      return(ATL_AMM_NSCHED-1);
/*
 * DELAY[0] < dly < DELAY[N-1], do binary search
 */
   do
   {
      iret = i0 + dist;
      mid=ATL_AMM_SCHEDdly[iret];
      if (mid < dly) /* change i0 */
         i0 = iret;
      else if (mid == dly)
         return(iret);
      dist >>= 1;     /* eliminated half of entries, endpt dist halved too */
   }
   while (dist > 1);
   iret = (mid > dly) ? iret-1 : iret;
   return(iret);
}
#endif
