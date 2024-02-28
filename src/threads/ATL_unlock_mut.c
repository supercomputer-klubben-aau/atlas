#include "atlas_tprim.h"
int ATL_unlock(void *vp)
{
#if defined(ATL_WINTHREADS) /* || (defined(ATL_OS_OSX) && defined(ATL_SSE1)) */
   #error "should not compile this file under Windows!"
#elif defined(ATL_OMP_THREADS)
   omp_unset_lock(vp);
   return(0);
#else
   return(pthread_mutex_unlock(vp));
#endif
}
