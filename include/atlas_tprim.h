#ifndef ATLAS_TPRIM_H
   #define ATLAS_TPRIM_H
/*
 * This file defines just a few primitive threading routines necessary
 * for us to perform the parallel runs necessary to test/time ATLAS.
 * Therefore, all of the files prototyped here can be built prior to any
 * ATLAS tuning, other than the creation of the two files in BLDdir/include:
 * 1. atlas_pthreads.h  : defines ATL_NTHREAD/NTHRPOW2
 * 2. atlas_taffinity.h : defines how to support affinity
 */
#define ATL_MEMBAR_ONLY 1
#include "atlas_cbc.h"
#include "atlas_pthreads.h"
#include "atlas_taffinity.h"
#include "atlas_aux.h"

#if defined(ATL_OS_Win64) || defined(ATL_OS_WinNT)
   #ifdef ATL_USE64BITS
      #define ATL_WIN64THREADS 1
      #define ATL_WINTHREADS 1
   #else
      #define ATL_WIN32THREADS 1
      #define ATL_WINTHREADS 1
   #endif
#endif
#include "atlas_pthreads.h" /* gened file defs ATL_NTHREADS & ATL_NTHRPOW2 */
#ifdef ATL_WINTHREADS
   #include <windows.h>
   typedef struct
   {
      void *vp;      /* ptr to extra info */
      HANDLE thrH;   /* handle to thread */
      int rank;      /* my rank */
      int P;         /* # of processors in this call */
      int affID;     /* affinity id for this core */
   } ATL_thread_t;
   #ifndef CREATE_SUSPENDED
      #define CREATE_SUSPENDED 0x00000004
   #endif
   #ifndef WAIT_FAILED
      #define WAIT_FAILED (~0)
   #endif
#elif defined(ATL_OMP_THREADS)
   #include <omp.h>
   typedef struct
   {
      void *vp;      /* ptr to extra info */
      int rank;      /* my rank */
      int P;         /* # of processors in this call */
      int affID;     /* < 0: not set, affID=1-affID, else affID */
   } ATL_thread_t;
  #ifndef ATL_lock_t
     #define ATL_lock_t omp_lock_t
  #endif
#else
   #define ATL_USE_POSIXTHREADS 1
   #include <pthread.h>
   typedef struct
   {
      pthread_t thrH;/* handle of thread */
      void *vp;      /* ptr to extra info */
      int rank;      /* my rank */
      int P;         /* # of processors in this call */
      int affID;     /* < 0: not set, affID=1-affID, else affID */
   } ATL_thread_t;
  #ifndef ATL_lock_t
     #define ATL_lock_t pthread_mutex_t
  #endif
#endif

int ATL_thread_start
   (ATL_thread_t *thr, int proc, int JOINABLE, void *(*rout)(void*), void*);
int ATL_thread_join(ATL_thread_t*); /* waits on completion of thread */

int ATL_lock(void *lck);
int ATL_trylock(void *lck);
/*
 * On weakly-ordered systems, unlock must flush store buffs before or at
 * same time as lock is written, and must prevent local reads from moving
 * above the unlock
 */
int ATL_unlock(void *lck);
void ATL_lock_init(void *lck);
void ATL_lock_destroy(void *lck);
#if 0
void *ATL_newLocks(unsigned int nlocks);
void *ATL_freeLocks(unsigned int nlocks, void *vp);
#endif

#ifdef ATL_SAFELS
   #if ATL_SAFELS > 128
      #if ATL_SAFELS == 256
         #define ATL_SAFELS 256L
         #define ATL_SAFELS_SH 8
      #elif ATL_SAFELS == 512
         #define ATL_SAFELS 512L
         #define ATL_SAFELS_SH 9
      #else
         #error "CRAZY ALERT: ATL_SAFELS has grown beyond 512!"
      #endif
   #else
      #undef ATL_SAFELS
   #endif
#endif
#ifndef ATL_SAFELS
   #define ATL_SAFELS 128L
   #define ATL_SAFELS_SH 7
#endif
#define ATL_AddBytesPtr(vp_, by_) (void*) ((size_t)(vp_) + (size_t)by_)
#define ATL_SubBytesPtr(vp_, by_) (void*) ((size_t)(vp_) - (size_t)by_)
#define ATL_AlignSafeLS(vp_) (void*) \
   ((((size_t)(vp_))+ATL_SAFELS-1L)&(~((size_t)ATL_SAFELS-1L)))
#define ATL_IncBySafeLS(vp_) ATL_AddBytesPtr(vp_, ATL_SAFELS)
#define ATL_DecBySafeLS(vp_) ATL_AddBytesPtr(vp_, -ATL_SAFELS)
/*
 * Helper function that starts P threads, waits for them to complete,
 * and returns.  Used in initial tuning before we've set up full thread
 * framework.
 */
#ifdef ATL_DEF_RUNTHR
#include <assert.h>
static void ATL_runThreads(int P, void *(*rout)(void*), void *arg)
{
   int i;
   ATL_thread_t tp[ATL_NTHREADS];
#ifdef ATL_OMP_THREADS
   omp_set_num_threads(P);
   #pragma omp parallel
   {
      assert(omp_get_num_threads() == P);
      i = omp_get_thread_num();
      #ifdef ATL_PAFF_SELF
         ATL_setmyaffinity(i);
      #endif
      tp[i].vp = arg;
      tp[i].rank = i;
      tp[i].P = P;
      rout(tp+i);
   }
#else
   for (i=0; i < P; i++)
   {
      tp[i].vp = arg;
      assert(!ATL_thread_start(tp+i, i, 1, rout, tp+i));
   }
   for (i=0; i < P; i++)
      assert(!ATL_thread_join(tp+i));
#endif
}
#endif

#endif
