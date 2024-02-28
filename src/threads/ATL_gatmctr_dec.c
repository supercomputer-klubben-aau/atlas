/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_gatmctr.h"  /* see to understand usage & data struct */

static INLINE long allPub(unsigned int rank, long *lp)
{
   const unsigned long P=lp[ATL_GAC_P];
   unsigned long i, PP;
   const long incCnt = lp[3];
   long *lowp = ATL_IncBySafeLS(lp), *acp = ATL_AddBytesPtr(lowp, lp[2]), *ac;
   unsigned long vrk;

   ac = ATL_AddBytesPtr(acp, incCnt*P);
   if (!(*ac))
      return(0);
/*
 * Use owned counter until exhausted
 */
   if (rank < P)
   {
      ac = ATL_AddBytesPtr(acp, incCnt*rank);
      if (*ac)
      {
         long lret;
         lret = ATL_atmctr_dec(ac);
         if (lret)
            return(lret+lowp[rank]);
      }
      vrk = (rank == P-1) ? 0 : rank+1;
      PP = P-1;
   }
   else
   {
      vrk = rank - (rank/P)*P;
      PP = P;
   }
/*
 * See if anybody else has count I can steal
 */
   ac = ATL_AddBytesPtr(acp, incCnt*vrk);
   for (i=0; i < PP; i++)
   {
      if (*ac)
      {
         long lret;
         lret = ATL_atmctr_dec(ac);
         if (lret)
            return(lret+lowp[vrk]);
      }
      if (++vrk != P)
         ac = ATL_AddBytesPtr(ac, incCnt);
      else
      {
         vrk = 0;
         ac = acp;
      }
   }
/*
 * Finally, check 1/0 counter
 */
   ac = ATL_AddBytesPtr(acp, incCnt*P);
   if (*ac)
      return(ATL_atmctr_dec(ac));
   return(0);
}
static INLINE long allPrv(unsigned int rank, long *lp)
{
   const unsigned long P=lp[ATL_GAC_P];
/*
 * Private counters are all that are available
 */
   if (rank < P)  /* so only P owners can hope for jobs */
   {
      lp = ATL_AddBytesPtr(lp, (rank+1)*ATL_SAFELS);
      if (*lp)
      {
         long lret;
         lret = *lp;
         if (lret)
         {
            *lp = lret - 1;
            return(lret + lp[1]);
         }
      }
   }
   return(0);
}

/*
 * Helper function that steals work from another thread.  Will make only
 * one pass through public queues in search of work, and return 0 if none
 * available at this time.
 */
static INLINE long tryPub
   (unsigned int rank, volatile long *gp, long *lp)
{
   const long P=gp[ATL_GAC_P], inc=gp[ATL_GAC_INC];
   unsigned int i;
/*
 * If I have my own public counter, prefentially use that up first
 */
   if (rank < P)
   {
      volatile long *lowPubP=ATL_AddBytesPtr(lp, inc*rank+ATL_SAFELS),
                            *cntPubP= ATL_IncBySafeLS(lowPubP);
      if (*cntPubP >= *lowPubP)  /* any public count left? */
      {
         void *lck = ATL_IncBySafeLS(cntPubP);
         long cnt;

         ATL_lock(lck);
         cnt = *cntPubP;
         if (cnt >= *lowPubP)
         {
            *cntPubP = cnt - 1;
            ATL_unlock(lck);
            return(cnt);
         }
         else
            ATL_unlock(lck);
      }
   }
/*
 * Make one pass through all P public counts, then give up
 */
   {
      unsigned int vrk = (rank >= P) ?  rank-P : rank;
      volatile long *cntPubP;
      cntPubP = ATL_AddBytesPtr(lp, (ATL_SAFELS<<1)+inc*vrk);
      for (i=0; i < P; i++)
      {
         volatile long *lowPubP = ATL_DecBySafeLS(cntPubP);

         if (rank >= P || i != rank)
         {
            if (*cntPubP >= *lowPubP)  /* any public count left? */
            {
               char *lck=ATL_IncBySafeLS(cntPubP);
               long lowPub, cnt;

               ATL_lock(lck);
               lowPub = *lowPubP;
               cnt = *cntPubP;
               if (cnt >= lowPub)
               {
                  *cntPubP = cnt - 1;
                  ATL_unlock(lck);
                  return(cnt);
               }
               else
                  ATL_unlock(lck);
            }
         }
         if (++vrk < P)
            cntPubP = ATL_AddBytesPtr(cntPubP, inc);
         else
         {
            vrk = 0;
            cntPubP = ATL_AddBytesPtr(lp, (ATL_SAFELS<<1));
         }
      }
/*
 *    If I'm quitting set it so we stop moving data from private to public.
 *    This avoids moving private to public when we have no workers left.
 *    QUIT should be set when the quitters can do another task, so ending
 *    load balance shouldn't matter much
 */
      gp[ATL_GAC_FLG] &= ~ATL_GAC_MOV;
   }
   return(0);
}

/*
 * Gets a single count from private range, returns 0 if none left
 */
static INLINE long
   prvCnt(unsigned int rank, volatile long *gp, long *tp)
{
   long *cntPrvP, *lowPubP;
   const long inc=gp[ATL_GAC_INC];
   long cnt, lowPub;

   if (rank > *gp)
      return(0);
   cntPrvP = ATL_AddBytesPtr(tp, inc*rank);
   lowPubP = ATL_IncBySafeLS(cntPrvP);

   cnt = *cntPrvP;
   if (cnt >= *lowPubP)
      return(0);
   *cntPrvP = cnt + 1;
   return(cnt);
}

static INLINE long prvPubCnt
   (unsigned int rank, volatile long *gp, long *lp)
/*
 * Get count from either private or public count.  Always prefers private
 * to public, and will move private work to public if the number of my
 * public jobs has changed (indicating work stealing has begun).
 * ASSUMES: rank < P.
 */
{
   long *prvP, *lowPubP;
   const long P=gp[ATL_GAC_P], inc=gp[ATL_GAC_INC], flg=gp[ATL_GAC_FLG];
   long lowPub, npriv;

   prvP = ATL_AddBytesPtr(lp, rank*inc);
   lowPubP = ATL_IncBySafeLS(prvP);
   lowPub = *lowPubP;
   npriv = *prvP;
   npriv = (npriv < lowPub) ? lowPub-npriv : 0;
   if (npriv)  /* can get job from private */
   {
      long *cntPubP = ATL_IncBySafeLS(lowPubP);
      const long lastPub=prvP[1];
/*
 *    Should we move some work from private to public queue? Only do so if
 *    this ctr allows moving work, we have extra private work,, and the number
 *    of public jobs has changes since the last time we checked (ctrs starting
 *    with 0 public work result in static scheduling on virtual rank).
 */
      if (GAC_FLAG_SET(flg, ATL_GAC_MOV) && npriv > 1 && (*cntPubP != lastPub))
      {
         long nmov, cntPub;
         void *lck = ATL_IncBySafeLS(cntPubP);
/*
 *       Lock my public ctr, and move some private jobs to it
 */
         ATL_lock(lck);
         lowPub = *lowPubP;
         npriv = *prvP;
         npriv = (npriv < lowPub) ? lowPub-npriv : 0;
         cntPub = *cntPubP;
         ATL_assert(lastPub > cntPub);
         nmov = lastPub - cntPub;  /* default: keep #of public jobs constant */
         nmov = (nmov < npriv) ? nmov : npriv-1;
         *lowPubP = lowPub - nmov;
         ATL_unlock(lck);
         prvP[1] = cntPub + nmov;  /* set lastPubCnt to new value */
         return((*prvP)++);
      }
      else
         return(prvCnt(rank, gp, lp));
   }
   return(tryPub(rank, gp, lp));
}

long ATL_gatmctr_dec(void *ac, unsigned int rank)
{
   if (ac)
   {
      long *gp=ac;
      const long P=gp[ATL_GAC_P], flg=gp[ATL_GAC_FLG];
      long ret;

      if (GAC_FLAG_SET(flg, ATL_GAC_NOPRV))
         ret = allPub(rank, ac);
      else if (GAC_FLAG_SET(flg, ATL_GAC_NOPUB))
         ret = allPrv(rank, ac);
      else
      {
         const long inc=gp[ATL_GAC_INC];
         long *lp=ATL_IncBySafeLS(gp);

         if (rank < P)
            ret = prvPubCnt(rank, gp, lp);
         else
         {
            rank = P + rank%P;
            ret = tryPub(rank, gp, lp);
         }
      }
      #ifdef DEBUG
         fprintf(stderr, "gDec=%ld\n", ret);
      #endif
      return(ret);
   }
   return(0);
}
