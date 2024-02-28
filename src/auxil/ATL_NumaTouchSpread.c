/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

#include "atlas_dsysinfo.h"
#include "atlas_threads.h"
typedef struct
{
   size_t N;    /* total size in bytes */
   char *buff;  /* buffer to be spread over cores */
} ATL_TNUMA_t;

void ATL_NumaCoreTouch(ATL_LAUNCHSTRUCT_t *lp, void *vp)
{
   ATL_thread_t *tp=vp;
   ATL_TNUMA_t *np=(ATL_TNUMA_t*)lp->opstruct;
   const size_t inc = ((size_t)ATL_pgsz)*tp->P, N = np->N;
   char *cp = np->buff;
   size_t J;

   for (J=((size_t)ATL_pgsz)*tp->rank; J < N; J += inc)
      cp[J] = 0;
}

/*
 * This routine takes recently allocated buffer b, splits it into pgsz
 * chunks, and spreads the pages over all cores evenly assuming a first-touch
 * allocation strategy.  This is used to avoid having all memory owned by
 * one core, which is a disaster on a system that doesn't use a distributed
 * directory (particularly AMD HTassist).
 */
void ATL_NumaTouchSpread(size_t N, void *buff)
{
   ATL_TNUMA_t num;
   num.N = N;
   num.buff = buff;
   ATL_goparallel(ATL_NTHREADS, ATL_NumaCoreTouch, &num, NULL);
}
