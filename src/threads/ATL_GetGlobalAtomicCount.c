/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2014, 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_threads.h"
int ATL_GetGlobalAtomicCount(void *vp, int rank)
/*
 * This routine returns a global counter that has been distributed over
 * P local counters
 */
{
   int i, j, P, b, icnt, extra, nL, *ip=vp, *iloc;
   void **acnts;

   P = ip[0];
   b = ip[1];
   extra = ip[2];
   nL = ip[3];
   iloc = ip+4;
/*
 * See if I can get the index from purely local information
 */
   if (rank < P && rank >= 0 && nL)
   {
      j = iloc[rank];
      if (j)
      {
         j += b * rank + Mmin(rank, extra);
/*fprintf(stderr, "%d: j=%d, LRET\n", rank, j);*/
         return(j);
      }
   }
   acnts = (void**) (ip+4+(((P+3)>>2)<<2));
/*
 * Otherwise, find an atomic counter that still has count
 */
   for (i=0; i < P; i++)
   {
/*
 *    If I got a counter value, convert it from local to global
 */
      icnt = (rank+i)%P;
      if (j = ATL_GetAtomicCount(acnts[icnt]))
      {
         j += nL + b*icnt + Mmin(icnt,extra);
         break;
      }
   }
/*fprintf(stderr, "%d: j=%d, icnt=%d, b=%d P=%d, e=%d\n", rank, j, icnt, b, P, extra);*/
   return(j);
}
