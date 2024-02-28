#include "atlas_cbc.h"
#include "atlas_threads.h"
/*
 * Used to post done to boolean sync array, by reversing present condition.
 * should only be called by the thread given by rank.
 * RETURNS: new value of post variable
 */
#ifdef DOFENCE
   #undef ATL_cbc_post
   char ATL_cbc_post(ATL_CUINT rank, void *vchk)
#else
   char ATL_post(ATL_CUINT rank, void *vchk)
#endif
{
   volatile char *bchk = (vchk) ? (volatile char*)vchk : ATL_TP_PTR->bchkin;
   ATL_CUINT II = rank<<ATL_chksh;
   const char newv = !bchk[II];

   #if defined(DOFENCE) && ATL_CBC_WEAK
      #if !ATL_CBC_NOBAR
         ATL_wmembar;
      #else
         ATL_mutex_lock(ATL_TP_PTR->cbcmut);   /* make my writes visible */
      #endif
   #endif
   bchk[II] = newv;
   #if defined(DOFENCE) && ATL_CBC_WEAK && ATL_CBC_NOBAR
      ATL_mutex_unlock(ATL_TP_PTR->cbcmut);
   #endif
   return(newv);
}
