#define Mstr2(m) # m
#define Mstr(m) Mstr2(m)
void ATL_xerbla(int p, char *rout, char *form, ...);
#define ATL_assert(n_) \
{ \
   if (!(n_)) \
   { \
      ATL_xerbla(0, __FILE__, "assertion %s failed, line %d of file %s\n", \
                 Mstr(n_), __LINE__, __FILE__); \
   } \
}
#include <stdlib.h>
#ifndef ATL_NOAFFINITY
#include "atlas_taffinity.h"
#endif
#ifdef ATL_WINTHREADS
   #include <windows.h>
#endif
#if defined(ATL_NOAFFINITY)
   #ifndef ATL_WINTHREADS
      #include <pthread.h>
   #endif
   #ifdef ATL_TUNING
      #define ATL_thread_start ATL_thread_start_noaff
   #endif
#elif defined(ATL_PAFF_SETAFFNP)
   #define _GNU_SOURCE 1 /* what manpage says you need to get CPU_SET */
   #define __USE_GNU   1 /* what actually works on linuxes I've seen */
   #define __USE_XOPEN2K8/* needed to avoid undef locale_t on some systems */
   #include <sched.h>    /* must be inced wt above defs before pthread.h */
   #include <pthread.h>
#elif defined(ATL_PAFF_SETPROCNP)
   #include <pthread.h>
#else
   #define ATL_PAFF_SELF
#endif
#include "atlas_tprim.h"
int ATL_thread_start(ATL_thread_t *thr, int proc, int JOINABLE,
                     void *(*rout)(void*), void *arg)
/*
 * Creates a thread that will run only on processor proc.
 * RETURNS: 0 on success, non-zero on error
 * NOTE: present implementation dies on error, so 0 is always returned.
 */
{
#ifdef ATL_WINTHREADS
   #ifdef ATL_WIN32THREADS
      DWORD thrID;
   #else
      unsigned thrID;
   #endif

   #ifdef ATL_NOAFFINITY
      #ifdef ATL_WIN32THREADS
         thr->thrH = CreateThread(NULL, 0, rout, arg, 0, &thrID);
      #else
         thr->thrH = (HANDLE)_beginthreadex(NULL, 0, rout, arg, 0, &thrID);
      #endif
      ATL_assert(thr->thrH);
   #else
      thr->rank = proc;
      #ifdef ATL_WIN32THREADS
         thr->thrH = CreateThread(NULL, 0, rout, arg, CREATE_SUSPENDED, &thrID);
      #else
         thr->thrH = (HANDLE)_beginthreadex(NULL, 0, rout, arg,
                                            CREATE_SUSPENDED, &thrID);
      #endif
      ATL_assert(thr->thrH);
      #ifdef ATL_RANK_IS_PROCESSORID
         ATL_assert(SetThreadAffinityMask(thr->thrH, ((long long)1)<<proc)));
         thr->affID = proc;
      #else
         thr->affID = ATL_affinityIDs[proc%ATL_AFF_NUMID];
         ATL_assert(SetThreadAffinityMask(thr->thrH,
                    (((long long)1)<<(thr->affID))));
      #endif
      ATL_assert(ResumeThread(thr->thrH) == 1);
   #endif
#elif defined(ATL_OMP_THREADS)
   fprintf(stderr, "Should not call thread_start when using OpenMP!");
   ATL_assert(0);
#elif 0 && defined(ATL_OS_OSX)  /* unchecked special OSX code */
/* http://developer.apple.com/library/mac/#releasenotes/Performance/RN-AffinityAPI/_index.html */
   pthread_attr_t attr;
   #define ATL_OSX_AFF_SETS 2       /* should be probed for */
   thread_affinity_policy ap;

   ap.affinity_tag = proc % ATL_OSX_AFF_SETS;
   ATL_assert(!pthread_attr_init(&attr));
   if (JOINABLE)
      ATL_assert(!pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE));
   else
      ATL_assert(!pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED));
   pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); /* no chk, OK to fail */

   ATL_assert(!pthread_create(&thr->thrH, &attr, rout, arg));
   ATL_assert(!thread_policy_set(thr->thrH, THREAD_AFFINITY_POLICY,
                                 (integer_t*)&ap,
                                 THREAD_AFFINITY_POLICY_COUNT));
   ATL_assert(!pthread_attr_destroy(&attr));
#else
   pthread_attr_t attr;
   #ifndef ATL_NOAFFINITY
      #if defined(ATL_PAFF_SETAFFNP)
         cpu_set_t *cpuset=NULL;
      #elif defined(ATL_PAFF_PLPA)
         plpa_cpu_set_t cpuset;
      #elif defined(ATL_PAFF_CPUSET) /* untried FreeBSD code */
         cpuset_t mycpuset;
      #endif
      #ifdef ATL_RANK_IS_PROCESSORID
         const int affID = proc;
      #else
         const int affID = ATL_affinityIDs[proc%ATL_AFF_NUMID];
      #endif
      #ifdef ATL_PAFF_SELF
         thr->affID = -affID-1; /* affinity must be set by created thread */
      #endif
   #endif
   thr->rank = proc;
   ATL_assert(!pthread_attr_init(&attr));
   if (JOINABLE)
   {
      #ifdef IBM_PT_ERROR
         ATL_assert(!pthread_attr_setdetachstate(&attr,
                                                 PTHREAD_CREATE_UNDETACHED));
      #else
         ATL_assert(!pthread_attr_setdetachstate(&attr,
                                                 PTHREAD_CREATE_JOINABLE));
      #endif
   }
   else
      ATL_assert(!pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED));
   pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); /* no chk, OK to fail */
   #ifdef ATL_PAFF_SETAFFNP
/*
 *    On PowerPC/Linux, pthread_attr_setaffinity_np sometimes reallocs() the
 *    cpuset variable, thus it must be malloced and not taken from stack!
 */
      cpuset = malloc(sizeof(cpu_set_t));
      CPU_ZERO(cpuset);
      CPU_SET(affID, cpuset);
      ATL_assert(!pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t),cpuset));
      free(cpuset);
   #elif defined(ATL_PAFF_SETPROCNP)
      ATL_assert(!pthread_attr_setprocessor_np(&attr, (pthread_spu_t)affID,
                                               PTHREAD_BIND_FORCED_NP));
   #endif
   ATL_assert(!pthread_create(&thr->thrH, &attr, rout, arg));
   ATL_assert(!pthread_attr_destroy(&attr));
#endif
   return(0);
}
