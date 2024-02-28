/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_gatmctr.h" /* scope to understand data struct */
void *ATL_gatmctr_init(void *ac, unsigned int P, long cnt, long flg)
{
   long *gp = ATL_AlignSafeLS(ac);
   long cntPrv, lowPub, cntPub;
   long lcntPub, lcntPrv, remPub, remPrv;
   unsigned int i;

   ATL_assert(cnt > 0);
   P = Mmax(P,1);
   P = Mmin(P,cnt);

   gp[ATL_GAC_P] = P;
   gp[ATL_GAC_FLG] = flg;

   if (GAC_FLAG_SET(flg, ATL_GAC_NOPUB))   /* private count only */
   {
      cntPub = 0;
      cntPrv = cnt;
      ATL_assert(!(flg&ATL_GAC_NOPRV));
      gp[2] = ((size_t)gp) - ((size_t)ac);
   }
   else if (GAC_FLAG_SET(flg, ATL_GAC_NOPRV))  /* public count only */
   {
      long inc, incLow;

      cntPrv = 0;
      cntPub = cnt;
      ATL_assert(!(flg&ATL_GAC_NOPUB));

      #if ATL_ATM_ASM
         inc = ATL_SAFELS;
      #else
         inc = ATL_SAFELS + (long) ATL_AlignSafeLS(sizeof(ATL_lock_t));
      #endif
      incLow = sizeof(long)*P;
      incLow = (long) ATL_AlignSafeLS(incLow);
      gp[2] = incLow;   /* inc to skip over low array to atmctrs */
      gp[3] = inc;      /* inc between atmctrs */
      gp[4] = ((size_t)gp) - ((size_t)ac);
   }
   else /* work stealing divides count between private/public */
   {
      long inc;
      if (cnt < (P<<2))         /* if we can't have both pub & prv */
         cntPub = cnt;         /* just make all jobs public */
      else
      {
         cntPub = P*(P-1);       /* ideal scheduling case */
         cntPrv = cnt>>1;
         cntPub = Mmin(cntPrv, cntPub);
      }
      cntPrv = cnt - cntPub;

      inc = sizeof(ATL_lock_t);           /* sz of native lock struct */
      inc = (long) ATL_AlignSafeLS(inc);  /* make mul of SAFELS */
      inc += ATL_SAFELS + (ATL_SAFELS<<1);
      gp[ATL_GAC_P] = P;
      gp[ATL_GAC_INC] = inc;
      gp[ATL_GAC_N] = cnt;
      gp[ATL_GAC_FLG] = flg;
      gp[4] = ((size_t)gp) - ((size_t)ac);
   }

   lcntPub = cntPub / P;
   remPub = cntPub - lcntPub*P;
   lcntPrv = cntPrv / P;
   remPrv = cntPrv - lcntPrv*P;

   cnt=1;
   for (i=0; i < P; i++)
   {
      long mypub, myprv;

      myprv = (i < remPrv) ? lcntPrv+1 : lcntPrv;
      mypub = (i < remPub) ? lcntPub+1 : lcntPub;
      ATL_gatmctr_initByRank(gp, i, cnt, myprv, mypub);
      cnt += myprv + mypub;
   }
   return(gp);
}

size_t ATL_gatmctr_sizeof(unsigned int P, long flg)
/*
 * RETURNS: size of specified global atomic counter, not including any
 *          alignment guard.
 */
{
   size_t sz;
   if (GAC_FLAG_SET(flg, ATL_GAC_NOPRV))   /* all-public */
   {
      long incCnt=ATL_SAFELS, incLow;
      #if !ATL_ATM_ASM
         incCnt += (long) ATL_AlignSafeLS(sizeof(ATL_lock_t));
      #endif
      incLow = sizeof(long)*P;
      incLow = (long) ATL_AlignSafeLS(incLow);
      sz = ATL_SAFELS + incLow + (P+1)*incCnt;
   }
   else if (GAC_FLAG_SET(flg, ATL_GAC_NOPUB))  /* all private */
      sz = (P+1)*ATL_SAFELS;
   else
   {
      sz = sizeof(ATL_lock_t);             /* sz of native lock struct */
      sz = (size_t) ATL_AlignSafeLS(sz);   /* make mul of SAFELS */
      sz += (ATL_SAFELS<<1)+ATL_SAFELS;
      sz = sz*P + ATL_SAFELS;
   }
   return(sz);
}

void *ATL_gatmctr_alloc(unsigned int nctr, unsigned int P, long cnt, long flg,
                        long *inc)
{
   if (cnt && nctr)
   {
      size_t sz;
      void *vp;
      void *gp;
      unsigned int i;

      sz = ATL_gatmctr_sizeof(P, flg);
      if (inc)
         *inc = sz;
      vp = malloc(ATL_SAFELS+nctr*sz);
      ATL_assert(vp);
      return(vp);
   }
   else
      *inc = 0;
   return(NULL);
}

void ATL_gatmctr_initByRank(void *ac, unsigned int rank, long lowPrv,
                            long npriv, long npub)
{
   if (ac)
   {
      long *ap = ac;
      const long P=ap[ATL_GAC_P];

      if (rank < P)
      {
         const long flg=ap[ATL_GAC_FLG];

         if (GAC_FLAG_SET(flg, ATL_GAC_NOPRV))
         {
            const long incLow=ap[2], incCnt=ap[3];
            long *lp=ATL_IncBySafeLS(ac), *atc;
            npub += npriv;
            if (rank)
            {
               atc = ATL_AddBytesPtr(lp, incLow+incCnt*rank);
               lp[rank] = lowPrv-1;
               ATL_atmctr_init(atc, npub);
            }
            else /* rank == 0, reserve one count for last */
            {
               atc = ATL_AddBytesPtr(lp, incLow);
               *lp = lowPrv;
               ATL_assert(npub && lowPrv == 1);
               ATL_atmctr_init(atc, --npub);
               atc = ATL_AddBytesPtr(lp, incLow+incCnt*P);
               ATL_atmctr_init(atc, 1);
            }
         }
         else if (GAC_FLAG_SET(flg, ATL_GAC_NOPUB))
         {
            long *lp = ATL_AddBytesPtr(ac, (rank+1)*ATL_SAFELS);
            npriv += npub;
            *lp   = npriv;
            lp[1] = lowPrv-1;
         }
         else
         {
            const long inc=ap[ATL_GAC_INC];
            const long lowPub = lowPrv+npriv;
            ap = ATL_AddBytesPtr(ap, ATL_SAFELS+rank*inc);
            *ap = lowPrv;
            ap[1] = (lowPub + npub) - 1;
            ap = ATL_IncBySafeLS(ap);
            *ap = lowPub;
            ap = ATL_IncBySafeLS(ap);
            *ap = (lowPub + npub) - 1;
            ap = ATL_IncBySafeLS(ap);
            ATL_lock_init(ap);
         }
      }
   }
}

void *ATL_gatmctr_new
(
   unsigned int P,    /* # of worker threads  */
   long N,    /* Count from [N,1] (0: already done) */
   long flg
)
{
   long *lp=NULL;

   if (!N)
      return(NULL);
   P = Mmax(P,1);
   P = Mmin(P,N);
   if (N)
   {
      lp = ATL_gatmctr_alloc(1, P, N, flg, NULL);
      lp = ATL_gatmctr_init(lp, P, N, flg);
   }
   return(lp);
}
