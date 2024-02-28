#include "atlas_threads.h"
void ATL_lock_destroy(void *vp)
{
#ifdef ATL_USE_mut
   #if defined(ATL_WINTHREADS) || (defined(ATL_OS_OSX) && defined(ATL_SSE1))
      #error "should not compile this file under Windows!"
   #elif defined(ATL_OMP_THREADS)
      omp_destroy_lock(vp);
   #else
      ATL_assert(!pthread_mutex_destroy(vp));
   #endif
#else
   *((unsigned char)vp) = 0;
#endif
}

