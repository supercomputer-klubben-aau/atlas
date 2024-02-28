#ifndef ATLAS_CBC_H
#define ATLAS_CBC_H

#ifndef ATL_MEMBAR_ONLY
   #include "atlas_misc.h"
#endif
/*
 ******************************************************************************
 * This file prototypes cache-based communication (CBC) routines.
 * The idea is exploit cache coherence mechanisms for thread synchronization
 * (barrier operations) and communication (combine/reduction operations).
 * In the best case, this improves performance from software levels to hardware.
 * If the call being replaced doesn't require real OS intervention (eg., mutex)
 * and is implemented by a reasonable intelligent OS (eg., Linux, not OS X),
 * then the speed improvement is still large due to the difference between
 * cache coherence speed and the speed of explicit mutex calls.
 *
 * The routines prototyped in this file come into two classes:
 * (1) Those with _cbc_ in the name enforce strong memory ordering in order
 *     to allow the primitive to guarantee that at the end of the sync, all
 *     cores have seen all the writes that preceeded the sync.
 * (2) Those without the _cbc_ just perform the sync/combine, but don't
 *     use a memory fence.
 *
 * Imagine you have a set of threads copying some data into a new format,
 * and another (possibly overlapping) set of threads that need to wait
 * until this copy is complete before they use the produced/copied data.
 * In this case, the sync will need to guarantee that the writes done by
 * the producing (copying) threads can be seen by consuming threads, so
 * the _cbc_ version of the routine must be used.
 *
 * Imagine a producing core copies a large amount of data, and then posts
 * to a boolean sync variable that he has produced the required data.  If
 * the data and sync variable are on the same cacheline, then as far as I
 * know this works correctly on all systems (this is how the scalar combines
 * work without using memory fences).
 *
 * In the more general case where the produced data is on independent cache
 * lines, then the sync will not work on systems with weakly-ordered caches.
 * Coherence protocols can be split into two types:
 * (1) STRONGLY ORDERED CACHES: guarantee that if a thread writes to cachelines
 *     A and then B, no other core can "see" the write of B before they see
 *     the write of A, since that's the order the producing thread wrote them.
 * (2) WEAKLY ORDERED CACHES: if a producing thread writes A and then B,
 *     other cores may see these writes in any order.  This means the fact
 *     that you have seen the change of B (eg. our sync variable saying we are
 *     are done producing data) does not mean it is safe to begin reading A!
 *
 * Weak ordering allows caches to retain data in private write buffers until
 * another core *writes* them (who cares if they read the stale value?), which
 * allows for less delays on write buffer flushes.
 *
 * In order to make CBC work when the data and sync variable are on separate
 * cachelines, weakly-ordered systems must use memory barriers to ensure that
 * seeing the change in B ensures that all prior writes can be seen.
 *
 * AMD64/IA32 have strongly-ordered caches, as long as you don't use the
 * special weakly-ordered stores.  I believe ATLAS does not use weakly-ordered
 * stores on x86 (I think this instructions were brought in with the PHI),
 * so presently the _cbc_ routines are simply aliased to their non-cbc
 * brethren on x86
 *
 * POWERPC/ARM64 feature weakly-ordered cache systems, and so we need _cbc_
 * routs with the correct memory and OOE fences to guarantee correct behavior.
 *
 * This file defines the following integer macros
 * ATL_CBC_WEAK : defined to 0 on systems known to be strongly ordered, else 1
 * ATL_CBC_STRONG: !ATL_CBC_WEAK
 * ATL_CBC_RBAR : 1 if ATL_rmembar is defined, else 0
 * ATL_CBC_WBAR : 1 if ATL_wmembar is defined, else 0
 * ATL_CBC_RWBAR: 1 if ATL_membar is defined, else 0
 * ATL_CBC_NOBAR: 1 if there is no safe memory barrier method (must use mutex)
 *
 * In addition, only strongly ordered caches, the macro ATL_ALIAS_CBC will
 * be defined, which just means that we don't compile the cbc_ variants of
 * some functions.  Remembering that ATL_barrier is a barrier using CBC
 * communication, and ATL_cbc_barrier is a barrier that *also* does a
 * memory fence, then these two funcs are the same thing on a strongly-ordered
 * system, so rather than compiling two routs, we'll just alias them.
 *
 * This file also provides macros for barriers that can be used:
 * ATL_membar : force all prior local stores globally visible before executing
 *              any following local stores, and do not allow any local loads
 *              to be hoisted above membar.
 * ATL_rmembar: do not allow any local loads to be hoisted above rmembar
 * ATL_wmembar: force all prior local stores globally visible before executing
 *              any following local stores
 *
 *
 ******************************************************************************
 */
/*
 * By default, we assume weakly-ordered with mutexes required.  This will
 * be overwridden for known systems below
 */
#define ATL_CBC_WEAK  1  /* assume weakly-ordered (worse case) */
#define ATL_CBC_RBAR  0
#define ATL_CBC_WBAR  0
#define ATL_CBC_RWBAR 0
#define ATL_CBC_NOBAR 1 /* assume no memory barrier support */

/*
 * All known x86 use strongly-ordered caches, so CBC is safe.
 * There are some instructions on the PHI that seem to violate strong ordering,
 * but ATLAS does not presently use them.  If these inst become important,
 * we may need to define the barriers, which I've outlined below (untested)
 */
#if defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)
   #undef ATL_CBC_WEAK
   #define ATL_CBC_WEAK  0
   #define ATL_ALIAS_CBC 1
   #ifdef ATL_ALIAS_CBC
      #define ATL_membar
      #define ATL_rmembar
      #define ATL_wmembar
   #else
      #define ATL_membar __asm__ __volatile__ ("mfence" : : : "memory")
      #define ATL_rmembar __asm__ __volatile__ ("lfence" : : : "memory")
      #define ATL_wmembar __asm__ __volatile__ ("sfence" : : : "memory")
   #endif
/*
 * ARM has weakly-ordered cache, so CBC must use explicit membarrier to work.
 * This memory barrier code provided by David Nuechterlein, who has gotten
 * CBC-based codes to work based on it.  Only defined for GNUC because I
 * need a way to do inline assembly, can support other compilers if given info.
 */
#elif defined(ATL_ARCH_ARM64) || defined(ATL_ARCH_ARMv7)
   #if __GNUC__
      #undef  ATL_CBC_RBAR
      #define ATL_CBC_RBAR  1
      #undef  ATL_CBC_WBAR
      #define ATL_CBC_WBAR  1
      #undef  ATL_CBC_RWBAR
      #define ATL_CBC_RWBAR 1
      #undef  ATL_CBC_NOBAR
      #define ATL_CBC_NOBAR 0
      #define ATL_membar __asm__ __volatile__ ("dmb sy" : : : "memory")
      #define ATL_rmembar __asm__ __volatile__ ("isb sy" : : : "memory")
      #define ATL_wmembar __asm__ __volatile__ ("dmb st" : : : "memory")
   #endif
/*
 * CBC came about after Itanium was essentially dead, so it has never been
 * tested there, so don't enable it.  We have this code here in case Itanium
 * becomes important enough to test.
 */
#elif defined(ATL_ARCH_IA64Itan) || defined(ATL_ARCH_IA64Itan2)
   #ifdef __GNUC__
      #define ATL_membar __asm__ __volatile__ ("mf")
   #endif
/*
 * On PowerPC and POWER I've never succeeded in getting any memory barrier
 * to work correctly.  IBM docs I've seen essentially say "here's how it should
 * work, but it doesn't, and we aren't going to tell you how it does work".
 * So, this code is just as a starting point if we find some docs or want
 * to do some experimentation later.
 */
#elif defined(ATL_ARCH_PPCG4) || defined(ATL_ARCH_PPCG5)
   #ifdef __GNUC__
      #define ATL_membar __asm__ __volatile__ ("sync")
   #endif
#elif defined(ATL_ARCH_POWER3) || defined(ATL_ARCH_POWER4) || \
      defined(ATL_ARCH_POWER5) || defined(ATL_ARCH_POWER6) || \
      defined(ATL_ARCH_POWER7) || defined(ATL_ARCH_POWER8)
   #ifdef __GNUC__
      #define ATL_membar __asm__ __volatile__ ("dcs")
   #endif
#endif
#ifndef ATL_CBC_WEAK
   #error "Malformed atlas_cbc.h!"
#endif
#if ATL_CBC_WEAK
   #define ATL_CBC_STRONG 0
#else
   #define ATL_CBC_STRONG 1
#endif
#if !ATL_CBC_RBAR && defined(ATL_rmembar)
   #undef ATL_rmembar
#endif
#if !ATL_CBC_WBAR && defined(ATL_wmembar)
   #undef ATL_wmembar
#endif
#if !ATL_CBC_RWBAR && defined(ATL_membar)
   #undef ATL_membar
#endif
/*
 * If read/write membarrier not defined, defined them as global membarrier
 */
#if !ATL_CBC_RBAR && ATL_CBC_RWBAR
   #define ATL_rmembar ATL_membar
#endif
#if !ATL_CBC_WBAR && ATL_CBC_RWBAR
   #define ATL_wmembar ATL_membar
#endif
#if !ATL_CBC_RWBAR && ATL_CBC_WBAR && ATL_CBC_RBAR
   #define ATL_membar { ATL_wmembar ; ATL_rmembar; }
#endif

#ifndef ATL_MEMBAR_ONLY
   void ATL_barrier(ATL_CUINT P, ATL_CUINT IAM, void*);
   void ATL_barrier_nopost0(ATL_CUINT P, ATL_CUINT IAM, void*);
   #ifdef ATL_ALIAS_CBC
      #define ATL_cbc_barrier ATL_barrier
      #define ATL_cbc_barrier_npost0 ATL_barrier_nopost0
   #else
      void ATL_cbc_barrier(ATL_CUINT P, ATL_CUINT IAM, void*);
      void ATL_cbc_barrier_nopost0(ATL_CUINT P, ATL_CUINT IAM, void*);
   #endif
   int ATL_icomb_sum(ATL_CUINT P, ATL_CUINT IAM, int val, void*);
   int ATL_icomb_max(ATL_CUINT P, ATL_CUINT IAM, int val, void*);
   int ATL_icomb_min(ATL_CUINT P, ATL_CUINT IAM, int val, void*);

   float ATL_scomb_sum(ATL_CUINT P, ATL_CUINT IAM, float val, void*);
   float ATL_scomb_max(ATL_CUINT P, ATL_CUINT IAM, float val, void*);
   float ATL_scomb_min(ATL_CUINT P, ATL_CUINT IAM, float val, void*);

   double ATL_dcomb_sum(ATL_CUINT P, ATL_CUINT IAM, double val, void*);
   double ATL_dcomb_max(ATL_CUINT P, ATL_CUINT IAM, double val, void*);
   double ATL_dcomb_min(ATL_CUINT P, ATL_CUINT IAM, double val, void*);

   void ATL_ccomb_sum(ATL_CUINT P, ATL_CUINT IAM, float *val, void*);
   void ATL_ccomb_max(ATL_CUINT P, ATL_CUINT IAM, float *val, void*);
   void ATL_ccomb_min(ATL_CUINT P, ATL_CUINT IAM, float *val, void*);
   void ATL_zcomb_sum(ATL_CUINT P, ATL_CUINT IAM, double *val, void*);
   void ATL_zcomb_max(ATL_CUINT P, ATL_CUINT IAM, double *val, void*);
   void ATL_zcomb_min(ATL_CUINT P, ATL_CUINT IAM, double *val, void*);

   int ATL_scomb_iamax(ATL_CUINT P, ATL_CUINT iam, int idx,
                       float *val, void *vchk);
   int ATL_scomb_iamax_nopost0(ATL_CUINT P, ATL_CUINT iam, int *idx,
                               float *val, void *vchk);
   #ifdef ATL_ALIAS_CBC
      #define ATL_scbc_comb_iamax_nopost0 ATL_scomb_iamax_nopost0
   #else
      int ATL_scbc_comb_iamax_nopost0(ATL_CUINT P, ATL_CUINT iam, int *idx,
                                      float *val, void *vchk);
   #endif
   int ATL_dcomb_iamax(ATL_CUINT P, ATL_CUINT iam, int idx,
                       double *val, void *vchk);
   int ATL_dcomb_iamax_nopost0(ATL_CUINT P, ATL_CUINT iam, int *idx,
                               double *val, void *vchk);
   #ifdef ATL_ALIAS_CBC
      #define ATL_dcbc_comb_iamax_nopost0 ATL_dcomb_iamax_nopost0
   #else
      int ATL_dcbc_comb_iamax_nopost0(ATL_CUINT P, ATL_CUINT iam, int *idx,
                                      double *val, void *vchk);
   #endif
   int ATL_ccomb_iamax(ATL_CUINT P, ATL_CUINT iam, int idx,
                       float *val, void *vchk);
   int ATL_ccomb_iamax_nopost0(ATL_CUINT P, ATL_CUINT iam, int *idx,
                               float *val, void *vchk);
   #ifdef ATL_ALIAS_CBC
      #define ATL_ccbc_comb_iamax_nopost0 ATL_ccomb_iamax_nopost0
   #else
      int ATL_ccbc_comb_iamax_nopost0(ATL_CUINT P, ATL_CUINT iam, int *idx,
                                      float *val, void *vchk);
   #endif
   int ATL_zcomb_iamax(ATL_CUINT P, ATL_CUINT iam, int idx,
                       double *val, void *vchk);
   int ATL_zcomb_iamax_nopost0(ATL_CUINT P, ATL_CUINT iam, int *idx,
                               double *val, void *vchk);
   #ifdef ATL_ALIAS_CBC
      #define ATL_zcbc_comb_iamax_nopost0 ATL_zcomb_iamax_nopost0
   #else
      int ATL_zcbc_comb_iamax_nopost0(ATL_CUINT P, ATL_CUINT iam, int *idx,
                                      double *val, void *vchk);
   #endif

/*
 * Reverse my entry in boolean sync array: RETURNS: new boolean value
 */
   char ATL_post(ATL_CUINT rank, void *vchk);
   #ifdef ATL_ALIAS_CBC
      #define ATL_cbc_post ATL_post
   #else
      char ATL_cbc_post(ATL_CUINT rank, void *vchk);
   #endif
#endif /* end ifndef ATL_MEMBAR_ONLY */

#endif /* end ifndef ATL_CBC_H */
