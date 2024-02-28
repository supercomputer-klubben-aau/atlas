#include "atlas_tprim.h"
int ATL_trylock(void *vp)
/*
 * return 0 if lock not required, else return non-zero
 */
{
#if defined(ATL_WINTHREADS) /*  || (defined(ATL_OS_OSX) && defined(ATL_SSE1)) */
   #error "should not compile this file under Windows!"
   return(ATL_DecAtomicCount(vp));
#elif defined(ATL_OMP_THREADS)
   return(omp_test_lock(vp));
#else
   return(!pthread_mutex_trylock(vp));
#endif
}
