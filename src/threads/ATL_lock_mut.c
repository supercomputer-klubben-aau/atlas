#include "atlas_tprim.h"
int ATL_lock(void *vp)
{
#ifdef ATL_WINTHREADS  /* if not using pthreads, use AtomicCount to sim mut */
   #error "should not compile this file under Windows!"
#elif 0 & defined(ATL_OS_OSX) && defined(ATL_SSE1)  /* assume OSX better? */
   #error "should not compile this file under OS X in x86!"
#elif defined(ATL_OMP_THREADS)
   omp_set_lock(vp);
   return(0);
#else
   return(pthread_mutex_lock(vp));
#endif
}
