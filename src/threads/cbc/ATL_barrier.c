#include "atlas_threads.h"
#include "atlas_cbc.h"
#include "atlas_misc.h"
/*
 * Use cache-based communication to perform a barrier for P threads.
 * This code works on any system with coherent caches (weakly-ordered OK).
 * For weakly-ordered caches, if this barrier is protecting memory, you
 * need to use the DOFENCE (cbc_) variant.
 */
#ifdef DOFENCE
   #ifdef ATL_cbc_barrier
      #undef ATL_cbc_barrier
   #endif
   #ifdef DOPOST
      #define mybarr ATL_cbc_barrier
   #else
      #define mybarr ATL_cbc_barrier_nopost0
   #endif
#else
   #ifdef DOPOST
      #define mybarr ATL_barrier
   #else
      #define mybarr ATL_barrier_nopost0
   #endif
#endif
/*
 * The DOFENCE varient must do a memory fence on systems wt weakly-ordered
 * caches.  On such systems w/o membar support (eg., POWER), we need
 * something for correctness, so I do lock & unlock a mutex to force
 * flushing the write buffers.  The code is ugly, but allows me to use
 * the same model of comminication w/o writing special code everywhere.
 *
 * It hurt me more to write this piece of crap than it hurts you to read it,
 * so stuff you and your judgement!
 */
void mybarr
(
   ATL_CUINT P,     /* # of threads to barrier */
   ATL_CUINT iam,   /* rank of calling thread in barrier */
   void *vchk
)
{
   volatile char *bchk = (vchk) ? (volatile char*)vchk : ATL_TP_PTR->bchkin;
   ATL_CUINT II = iam<<ATL_chksh;
   const char newv = !bchk[II];

   if (iam)
   {
      #if defined(DOFENCE) && ATL_CBC_WEAK /* need to membar */
         #if ATL_CBC_NOBAR
            ATL_mutex_lock(ATL_TP_PTR->cbcmut);
         #else
            ATL_wmembar;
         #endif
      #endif
      bchk[II] = newv;
      #if defined(DOFENCE) && ATL_CBC_WEAK && ATL_CBC_NOBAR
         ATL_mutex_unlock(ATL_TP_PTR->cbcmut);
      #endif
      while (*bchk != newv);
      #if defined(DOFENCE) && ATL_CBC_WEAK
         #if !ATL_CBC_NOBAR
            ATL_rmembar;
         #else
            ATL_mutex_lock(ATL_TP_PTR->cbcmut);   /* crappy ATL_rmembar */
            ATL_mutex_unlock(ATL_TP_PTR->cbcmut); /* with huge overhead! */
         #endif
      #endif
   }
   else
   {
      int i;
      for (i=1; i < P; i++)
      {
         ATL_CUINT d = i<<ATL_chksh;
         while (bchk[d] != newv);
      }
      #ifdef DOPOST
         #if defined(DOFENCE) && ATL_CBC_WEAK
            #if !ATL_CBC_NOBAR
               ATL_wmembar;
            #else
               ATL_mutex_lock(ATL_TP_PTR->cbcmut);
            #endif
         #endif
         *bchk = newv;
      #endif
      #if defined(DOFENCE) && ATL_CBC_WEAK
         #if !ATL_CBC_NOBAR
            ATL_rmembar;
         #else
            ATL_mutex_unlock(ATL_TP_PTR->cbcmut);
         #endif
      #endif
   }
}
