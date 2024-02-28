#include "atlas_threads.h"
void ATL_lock_init(void *vp)
{
/*
 * On Windows, use known-good x86 code.  OS X's mutex have horrible scaling,
 * so use homebrewed code instead
 */
#ifdef ATL_USE_mut
   #if defined(ATL_WINTHREADS) /* || (defined(ATL_OS_OSX)&&defined(ATL_SSE1)) */
   #error "should not compile this file under Windows!"
   #elif defined(ATL_OMP_THREADS)
      omp_init_lock(vp);
   #else
      ATL_assert(!pthread_mutex_init(vp, NULL));
   #endif
#else
   *((unsigned char)vp) = 0;
#endif
}
