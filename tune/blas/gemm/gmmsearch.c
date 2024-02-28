/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016, 2015 R. Clint Whaley
 */
#include "atlas_cache.h"
#include "atlas_misc.h"
#include "atlas_mmtesttime.h"
/*
 * Flag for iFKO (must be power of 2).
 */
 #define FKO_FLAG 262144
static int Mylcm(const int M, const int N)
/*
 * Returns least common multiple (LCM) of two positive integers M & N by
 * computing greatest common divisor (GCD) and using the property that
 * M*N = GCD*LCM.
 */
{
   register int tmp, max, min, gcd=0;

   if (M != N)
   {
      if (M > N) { max = M; min = N; }
      else { max = N; min = M; }
      if (min > 0)  /* undefined for negative numbers */
      {
         do  /* while (min) */
         {
            if ( !(min & 1) ) /* min is even */
            {
               if ( !(max & 1) ) /* max is also even */
               {
                  do
                  {
                     min >>= 1;
                     max >>= 1;
                     gcd++;
                     if (min & 1) goto MinIsOdd;
                  }
                  while ( !(max & 1) );
               }
               do min >>=1 ; while ( !(min & 1) );
            }
/*
 *          Once min is odd, halve max until it too is odd.  Then, use
 *          property that gcd(max, min) = gcd(max, (max-min)/2)
 *          for odd max & min
 */
MinIsOdd:
            if (min != 1)
            {
               do  /* while (max >= min */
               {
                  max -= (max & 1) ? min : 0;
                  max >>= 1;
               }
               while (max >= min);
            }
            else return( (M*N) / (1<<gcd) );
            tmp = max;
            max = min;
            min = tmp;
         }
         while(tmp);
      }
      return( (M*N) / (max<<gcd) );
   }
   else return(M);
}
double TimeMMKernel_KB
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int FORCETIME,               /* 1: ignore any prior output file */
   ATL_mmnode_t *mmp,           /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
/*
 * If kernel has property KRUNTIME, try timing it with compile- and run-time K,
 * and if compile-time is more than 2% faster, turn off KRUNTIME
 */
{
   double mf;
   const int kb0 = mmp->kbB;
   if (mmp->ID && mmp->kbmax && kb > mmp->kbmax)  /* genned codes can adapt */
      return(0.0);
   mmp->mbB = mb;
   mmp->nbB = nb;
   mmp->kbB = kb;
/*
 * If it's a generated kernel, regen in case we need KB to match KU, and
 * to specialize it to KB (assuming compile-time is faster)
 */
   if (!mmp->ID)
   {
      if (mmp->genstr)
         free(mmp->genstr);
      if (mmp->rout)
         free(mmp->rout);
      if (FLAG_IS_SET(mmp->flag, MMF_KUISKB))
         mmp->kbmax = mmp->kbmin = mmp->ku = kb;
      mmp->rout = DupString("ATL_tmp.c");
      mmp->genstr = MMGetGenString(pre, mmp);
   }
   mf = TimeMMKernel(verb, FORCETIME, mmp, pre, mb, nb, kb, beta, mflop,cflush);
   if (FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
   {
      double mfC;
      mmp->flag &= ~(1<<MMF_KRUNTIME);
      mfC = TimeMMKernel(verb, FORCETIME, mmp, pre, mb, nb, kb, beta,
                         mflop, cflush);
      if (mfC <= 1.02*mf)
         mmp->flag |= (1<<MMF_KRUNTIME);
      else
      {
         if (verb)
            printf("      Forcing K compile-time, mfC=%.2f, mfR=%.2f\n",
                   mfC, mf);
         mf = mfC;
      }
   }
   return(mf);
}
ATL_mmnode_t *MergeCases
(
   int imf,
   ATL_mmnode_t *bs0, /* queue of cases */
   ATL_mmnode_t *bs1  /* queue of cases */
)
/*
 * Merges two queues of matmul kern cases.  Cases are not winnowed, but
 * duplicates are not allowed, so if two entries have the same kbB, then
 * we take the one with best mflop[imf].  If imf < 0, then we do indeed
 * allow duplicates of kbB.
 * NOTE: does not change bs0 or bs1.
 * ASSUMES: both bs0 & bs1 are in kb-increasing order.
 * RETURNS: base ptr to merged queue
 */
{
   ATL_mmnode_t *mb=NULL, *mp;
   while (bs0 || bs1)
   {
      ATL_mmnode_t *p;
      if (bs0 && bs1)
      {
         if (bs0->kbB < bs1->kbB)
         {
            p = CloneMMNode(bs0);
            bs0 = bs0->next;
         }
         else if (bs0->kbB > bs1->kbB)
         {
            p = CloneMMNode(bs1);
            bs1 = bs1->next;
         }
         else /* they are equal, must take best performer, or both */
         {
/*
 *          If we are taking both, special case can't use general completion
 */
            if (imf < 0)
            {
               p = CloneMMNode(bs0);
               bs0 = bs0->next;
               p->next = CloneMMNode(bs1);
               bs1 = bs1->next;
               if (mb)
                  mp->next = p;
               else
                  mb = p;
               mp = p->next;
               continue;
            }
/*
 *          Taking only the best performer, but moving both base ptrs
 */
            else
            {
/*
 *             If they are equal, take the KRUN=1 case if it exists, else
 *             take the most flexible one or one requiring the least cleanup
 */
               if (bs0->mflop[imf] == bs1->mflop[imf])
               {
                  if (FLAG_IS_SET(bs0->flag, MMF_KRUNTIME))
                     p = bs0;
                  else if (FLAG_IS_SET(bs1->flag, MMF_KRUNTIME))
                     p = bs1;
                  else if (bs0->ku < bs1->ku)
                     p = bs0;
                  else if (bs1->ku < bs0->ku)
                     p = bs1;
                  else
                  {
                     const int u0=Mmax(bs0->mu, bs0->nu),
                               u1=Mmax(bs1->mu, bs1->nu);
                     p = (u0 <= u1) ? bs0 : bs1;
                  }
               }
               else
                  p = (bs0->mflop[imf] > bs1->mflop[imf]) ? bs0 : bs1;
               p = CloneMMNode(p);
               bs0 = bs0->next;
               bs1 = bs1->next;
            }
         }
      }
      else if (bs0)
      {
         p = CloneMMNode(bs0);
         bs0 = bs0->next;
      }
      else /* if (bs1) */
      {
         p = CloneMMNode(bs1);
         bs1 = bs1->next;
      }
      if (mb)
      {
         mp->next = p;
         mp = p;
      }
      else
        mp = mb = p;
   }
   return(mb);
}
ATL_mmnode_t *WinnowCases
(
   int imf,
   ATL_mmnode_t *mb  /* queue of cases */
)
/*
 * Removes any case that runs slower than a smaller case
 * RETURNS: mb with queue bad kernels deleted
 * NOTE: mb can never change, since by def nothing smaller than 1st case
 */
{
   ATL_mmnode_t *prev = mb, *mp;

   if (!mb)
      return(NULL);
   mp = mb->next;
   while (mp)
   {
      if (mp->mflop[imf] <= prev->mflop[imf])  /* kill slow KB */
         mp = prev->next = KillMMNode(mp);
      else
      {
         prev = mp;
         mp = mp->next;
      }
   }
   return(mb);
}
int KernelIsUnique(ATL_mmnode_t *mmb, ATL_mmnode_t *mmp)
/*
 * Determines if mmp is the first mention of a unique kernel in mmb, or not.
 * For user cases (ID > 0), (ID,flag) together make a unique kernel.
 * For user generated cases, if they match on : mu,nu,ku,VLEN,flag
 *
 * RETURNS: 0 if mmp appears in mmb before mmp, else 1
 */
{
   ATL_mmnode_t *mp;
   if (mmp == mmb)
      return(1);
   for (mp=mmb; mp && mp != mmp; mp = mp->next)
      if (MMKernsPerfSame(mmp, mp))
         return(0);
   return(1);  /* didn't find it, must be first time in list */
}
double CacheRatio_all3pAB(size_t CS, size_t mb, size_t nb, size_t kb,
                          size_t mu, size_t nu)
{ /* RETURNS: ratio of util cache A,B,C + next A & B */
   double dret = 2.0*kb*(mb+nb)+mb*nb;
   return(dret/CS);
}
double CacheRatio_all3pA(size_t CS, size_t mb, size_t nb, size_t kb,
                         size_t mu, size_t nu)
{ /* RETURNS: ratio of utilized cache to keep all 3 mm ops + next A in CS */
   double dret = kb*(mb+mb+nb)+mb*nb;
   return(dret/CS);
}
double CacheRatio_all3pB(size_t CS, size_t mb, size_t nb, size_t kb,
                         size_t mu, size_t nu)
{ /* RETURNS: ratio of utilized cache to keep all 3 mm ops + next A in CS */
   double dret = kb*(mb+nb+nb)+mb*nb;
   return(dret/CS);
}

double CacheRatio_all3(size_t CS, size_t mb, size_t nb, size_t kb,
                       size_t mu, size_t nu)
{ /* RETURNS: ratio of utilized cache to keep all 3 mm ops in CS */
   double dret = 1.0*kb*(mb+nb)+mb*nb;
   return(dret/CS);
}

double CacheRatio_one(size_t CS, size_t mb, size_t nb, size_t kb,
                      size_t mu, size_t nu)
{ /* RETURNS: ratio of util cache for B + working set of A/C */
   double dret = kb*nb + 2.0*(mu*kb + mu*nu);
   return(dret/CS);
}

double CacheRatio_ws(size_t CS, size_t mb, size_t nb, size_t kb,
                     size_t mu, size_t nu)
{ /* RETURNS: ratio of working set of all matmul ops to CS */
   double dret = 2.0*(mu*nu + nu*kb) + mu*kb;
   return(dret/CS);
}

typedef void (*BudgetFunc_t)(double, size_t, size_t, size_t, size_t,
                             size_t*, size_t*, size_t*);

#define MAXNB 360
void GetBlkFromBudget_allP(char extra, double thresh, size_t CS,
                           size_t mu, size_t nu, size_t ku,
                           size_t *MB, size_t *NB, size_t *KB)
{
   double (*cacheRatio)(size_t,size_t,size_t,size_t,size_t,size_t);
   size_t mb=mu, nb=nu, kb=ku;
   int MGROW, NGROW, KGROW;

   if (extra == 'A')
      cacheRatio = CacheRatio_all3pA;
   else if (extra == 'B')
      cacheRatio = CacheRatio_all3pB;
   else if (extra == '2')
      cacheRatio = CacheRatio_all3pAB;
   else  /* nothing extra */
      cacheRatio = CacheRatio_all3;
   do
   {
      size_t mn=mb+mu, nn=nb+nu, kn=kb+ku;
      MGROW = mn < MAXNB && (cacheRatio(CS, mn, nb, kb, mu, nu) <= thresh);
      NGROW = nn < MAXNB && (cacheRatio(CS, mb, nn, kb, mu, nu) <= thresh);
      KGROW = kn < MAXNB && (cacheRatio(CS, mb, nb, kn, mu, nu) <= thresh);
      if (KGROW && ((!MGROW && !NGROW) || (kn <= nn && kn <= mn)))
         kb = kn;
      else if (MGROW && (!NGROW || mb < nb))
         mb = mn;
      else if (NGROW)
         nb = nn;
   }
   while (MGROW | NGROW | KGROW);
   *MB = mb;
   *NB = nb;
   *KB = kb;
}

void GetBlkFromBudget_all3pA(double thresh, size_t CS,
                           size_t mu, size_t nu, size_t ku,
                           size_t *MB, size_t *NB, size_t *KB)
{
   GetBlkFromBudget_allP('A', thresh, CS, mu, nu, ku, MB, NB, KB);
}
void GetBlkFromBudget_all3pB(double thresh, size_t CS,
                           size_t mu, size_t nu, size_t ku,
                           size_t *MB, size_t *NB, size_t *KB)
{
   GetBlkFromBudget_allP('B', thresh, CS, mu, nu, ku, MB, NB, KB);
}
void GetBlkFromBudget_all3pAB(double thresh, size_t CS,
                           size_t mu, size_t nu, size_t ku,
                           size_t *MB, size_t *NB, size_t *KB)
{
   GetBlkFromBudget_allP('2', thresh, CS, mu, nu, ku, MB, NB, KB);
}
void GetBlkFromBudget_all3(double thresh, size_t CS,
                           size_t mu, size_t nu, size_t ku,
                           size_t *MB, size_t *NB, size_t *KB)
{
   GetBlkFromBudget_allP('N', thresh, CS, mu, nu, ku, MB, NB, KB);
}

void GetBlkFromBudget_one(double thresh, size_t CS,
                           size_t mu, size_t nu, size_t ku,
                           size_t *MB, size_t *NB, size_t *KB)
{
   size_t mb=mu, nb=nu, kb=ku;
   int NGROW, KGROW;
   do
   {
      size_t nn=nb+nu, mn=(nn/mu)*mu, mn1 = ((nn+mu-1)/mu)*mu, kn=kb+ku;
      if (mn1 - nn <= nn - mn || !mn)
         mn = mn1;
      NGROW = nn < MAXNB && (CacheRatio_one(CS, mn, nn, kb, mu, nu) <= thresh);
      KGROW = kn < MAXNB && (CacheRatio_one(CS, mb, nb, kn, mu, nu) <= thresh);
      if (NGROW && ((nn < kn && mn < kn) || !KGROW))
      {
         nb = nn;
         mb = mn;
      }
      else if (KGROW)
         kb = kn;
   }
   while (NGROW | KGROW);
   *MB = mb;
   *NB = nb;
   *KB = kb;
}

void GetBlkFromBudget_ws(double thresh, size_t CS,
                          size_t mu, size_t nu, size_t ku,
                          size_t *MB, size_t *NB, size_t *KB)
{
   size_t mb=mu, nb=nu, kb=ku;
   int KGROW;
   do
   {
      size_t kn=kb+ku, mn=(kn > mu)?(kn/mu)*mu:mu, nn=(kn>nu)?(kn/nu)*nu:nu;
      KGROW = kn < MAXNB && (CacheRatio_ws(CS, mn, nn, kn, mu, nu) <= thresh);
      if (KGROW)
      {
         mb = mn;
         nb = mn;
         kb = kn;
      }
   }
   while (KGROW);
   *MB = mb;
   *NB = nb;
   *KB = kb;
}

int Blk2Case(size_t CS, int mb, int nb, int kb, int mu, int nu)
{
   size_t sz3 = mb*nb + kb*(mb+nb);
   if (CacheRatio_all3pAB(CS, mb, nb, kb, mu, nu) < 1.0)
      return(4);
   else if (CacheRatio_all3pA(CS, mb, nb, kb, mu, nu) < 1.0 &&
            CacheRatio_all3pB(CS, mb, nb, kb, mu, nu) < 1.0)
      return(3);
   else if (CacheRatio_all3(CS, mb, nb, kb, mu, nu) < 1.0)
      return(0);
   else if (CacheRatio_one(CS, mb, nb, kb, mu, nu) < 1.0)
      return(1);
   return(2);
}
ATL_mmnode_t *FindBestCacheBudgetCase
(
   int verb,
   char pre,
   BudgetFunc_t GetBlocking,     /* func ptr to budget function */
   double thresh,                /* max ratio of cache to fill */
   size_t CS,                    /* size of cache we are optimizing for */
   int imf,                      /* entry in mflop[] to use */
   ATL_mmnode_t *mmb             /* list of cases to try */
)
/*
 * RETURNS: clone of best-peforming kernel in mmb for kb=kb, mb & nb
 *          near-square and within budget
 */
{
   ATL_mmnode_t *mmB=NULL, *mp, *p;
   double mf, mfB=0.0;

   printf("Finding best case for cache budget case=%d, CS=%.0f elts\n",
          imf, CS*thresh);
   for (mp=mmb; mp; mp = mp->next)
   {
      size_t mb, nb, kb, ku=mp->ku;
      if (!mp->ID && FLAG_IS_SET(mp->flag, MMF_KUISKB))
         ku = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen : 1;

      GetBlocking(thresh, CS, mp->mu, mp->nu, ku, &mb, &nb, &kb);
      p = CloneMMNode(mp);  /* can't use mp, since may switch KRUNTIME */
      mf = TimeMMKernel_KB(verb, 0, p, pre, mb, nb, kb, 1, 0, -1);
      printf("   RT='%s' B=(%d,%d,%d), MFLOP=%.2f\n", p->rout,
             (int)mb, (int)nb, (int)kb, mf);
      if (mf > mfB)
      {
         if (mmB)
            KillMMNode(mmB);
         p->mbB = mb;
         p->nbB = nb;
         p->kbB = kb;
         mmB = p;
         mfB = mmB->mflop[imf] = mf;
      }
      else
         KillMMNode(p);
   }
   printf("BEST CASE %s: mb=%d, nb=%d, kb=%d, RTK=%d, MFLOP=%.2f\n\n",
          mmB->rout ? mmB->rout : "GENNED",
          mmB->mbB, mmB->nbB, mmB->kbB, FLAG_IS_SET(mmB->flag, MMF_KRUNTIME),
          mmB->mflop[imf]);
   mmB->next = NULL;
   return(mmB);
}

ATL_mmnode_t *FindBestCacheBudgetCases
(
   int verb,
   char pre,
   size_t CS,                    /* size of cache we are optimizing for */
   ATL_mmnode_t *mmb             /* list of cases to try */
)
/*
 * This case attempts to find the best kernel for 3-5 cases of interest:
 * (1) All 3 matrices fit in CS -- this case is designed for when we wish
 *     to reuse at least one of the matrices *across* mmkern calls.  It is
 *     particularly good for complex arithmetic, or when CS is large enough
 *     that A&B are reused so much internally to a mmkern call that it makes
 *     sense to retain C in cache for the next mmkern call in K-loop.
 * (2) All of B fits in CS, and so does the working set of A&C.  This one
 *     reuses all internal ops from L1, but won't allow any full op reuse
 *     across multiple GEMM calls.  Usually best for small-to-medium cache
 *     sizes.
 * (3) All of working set of A/B/C fit in cache.  This case provides maximal
 *     NB, where only the mu*KB panel of A is reused from the L1 internally
 *     to the algorthm.  It is essentially an L2-blocked algorithm internally,
 *     but can be useful on those archs where the best sustained bandwidth
 *     comes from one L1 load (A) and 1 L2 load (B).
 * If 3 blocks fit with smallest dim > 16, then also add:
 * (4) All 3 blocks fit in cache, with room for the next A or B.
 * If this fits in cache with smallest dim > 16, then add:
 * (5) All 3 blocks fit in cache, with room for the next A and B.
 */
{
   ATL_mmnode_t *mm3, *mm1, *mmw;
   mm3 = FindBestCacheBudgetCase(verb, pre, GetBlkFromBudget_all3, 1.0,
                                 CS, 1, mmb);
   mm1 = FindBestCacheBudgetCase(verb, pre, GetBlkFromBudget_one, 1.0,
                                 CS, 2, mmb);
   mmw = FindBestCacheBudgetCase(verb, pre, GetBlkFromBudget_ws, .90,
                                 CS, 3, mmb);
   mm3->next = mm1;
   mm1->next = mmw;
   mmw->next = NULL;

/*
 * If cache large enough, try fitting 4 & 5 blocks in it
 */
   if (mm3->mbB > 16 && mm3->nbB > 16 && mm3->kbB > 16)
   {
      ATL_mmnode_t *mp;
      mp = FindBestCacheBudgetCase(verb, pre, (mm3->mbB > mm3->nbB) ?
              GetBlkFromBudget_all3pA:GetBlkFromBudget_all3pB,
              0.95, CS, 4, mmb);
      mmw->next = mp;
      if (mp->mbB > 16 && mp->nbB > 16 && mp->kbB > 16)
      {
         ATL_mmnode_t *mm5;
         mm5 = FindBestCacheBudgetCase(verb, pre, GetBlkFromBudget_all3pAB,
                                       0.95, CS, 5, mmb);
         mp->next = mm5;
      }
   }
   {
      int i;
      ATL_mmnode_t *mp;
      char *exp[5] = {"3BLKS", "1BLKS", "0BLKS", "4BLKS", "5BLKS"};

      for (i=0, mp=mm3; mp; i++, mp=mp->next)
         printf("%s: RT='%s' B=(%d,%d,%d), RK=%d MFLOP=%.2f\n",
                exp[i], mp->rout, mp->mbB, mp->nbB, mp->kbB,
                FLAG_IS_SET(mp->flag, MMF_KRUNTIME), mp->mflop[i+1]);
      printf("\n");
   }
   return(mm3);
}

int NumberBetaFails(FILE *fperr, char pre, int nb, ATL_mmnode_t *p)
{
   const int mu = p->mu, nu = p->nu, ku = p->ku;
   int mb = ((nb+mu-1)/mu)*mu, kb = ((nb+ku-1)/ku)*ku;
   int i, nfail = 0;

   if (kb < p->kbmin)
      kb = p->kbmin;
   if (p->kbmax && kb > p->kbmax)
      kb = p->kbmax;
   nb = ((nb+nu-1)/nu)*nu;
   for (i=(-1); i < 2; i++)
   {
      if (pre == 'z' || pre == 'c')  /* beta=0 case tests all 3 real betas */
         i = 0;
      if (MMKernelFailsTest(pre, mb, nb, kb, i, p))
      {
         if (fperr)
         {
            char *sp;
            fprintf(fperr, "FAIL: B=(%d,%d,%d), rout='%s', genstr='%s'\n",
                    mb, nb, kb, p->rout?p->rout:"NULL",
                    p->genstr?p->genstr:"NULL");
            sp = MMGetTestString(pre, mb, nb, kb, i, p);
            fprintf(fperr,"   '%s'\n", sp);
            free(sp);
         }
         nfail++;
      }
      if (pre == 'z' || pre == 'c')
         break;
   }
   return(nfail);
}

/*
 * RETURNS: 0 if bcast slower than ld/splat combo
 */
int UseBcast(int flg, int verb, char pre, int kb, int nreg, int VL)
{
   int nu=VL, mu, mb, nb, TEST=flg&1;
   ATL_mmnode_t *mpBC, *mpNO;
   double mfBC, mfNO;

   if (flg&FKO_FLAG)        /* temporary for Majedul */
      return(1);            /* replace wt analysis later */
   if (VL < 2)
      return(1);
   mu = (nreg-nu-1) / (nu+1);
   if (mu < 2)
      mu = (nreg - 2) / (nu+1);
   mu = (mu) ? mu : 1;
   mb = ((kb+mu-1)/mu)*mu;
   nb = ((kb+nu-1)/nu)*nu;
   mpBC = MMGetNodeGEN(pre, 0, nb, mu*VL, nu, 1, VL, 0, 0, 0, NULL);
   mpNO = MMGetNodeGEN(pre, 1, nb, mu*VL, nu, 1, VL, 0, 0, 0, NULL);
   printf("TIMING BCAST VS SPLAT MVEC WITH: B=(%d,%d,%d) U=(%d,%d,1)\n",
          mb, nb, kb, mu, nu);
   mfBC = TimeMMKernel(verb, 1, mpBC, pre, mb, nb, kb, 1, 0, -1);
   printf("   BCAST = %.0f MFLOP\n", mfBC);
   mfNO = TimeMMKernel(verb, 1, mpNO, pre, mb, nb, kb, 1, 0, -1);
   printf("   SPLAT = %.0f MFLOP\n", mfNO);
   if (TEST)
   {
      printf("   TESTING . . .");
      assert(!NumberBetaFails(stderr, pre, nb, mpBC));
      printf("  . . . ");
      assert(!NumberBetaFails(stderr, pre, nb, mpNO));
      printf("   PASS!\n");
   }
   KillMMNode(mpBC);
   KillMMNode(mpNO);
   if (mfNO > mfBC)
   {
      printf("VLD/VSPLAT PROVIDES %.4f SPEEDUP\n", mfNO/mfBC);
      return(0);
   }
   printf("VBCAST PROVIDES %.4f SPEEDUP\n", mfBC/mfNO);
   return(1);
}

ATL_mmnode_t *DoSyrkMUNU
   (int flg, int vrb, char pre, int nreg, int nb, int kb, int VL)
{
   ATL_mmnode_t *pM, *pB, *pK, *mp;
   double mfB = 0.0;
   int i, j, uB=VL, bcB=1, kvecB=0, kbB=0;
   const char *frm="   %c%c %4d %4d %3d %9.0f\n";
   const int bf=(flg&FKO_FLAG)?FKO_FLAG:0;

   pM = MMGetNodeGEN(pre, 0|bf, 0, VL, 1, 1, VL, 0, 0, 0,
                     DupString("ATL_tmp.c"));
   pB = MMGetNodeGEN(pre, 1|bf, 0, VL, 1, 1, VL, 0, 0, 0,
                     DupString("ATL_tmp.c"));
   pK = MMGetNodeGEN(pre, 0|bf, 0, 1, 1, VL, VL, 1, 0, 0,
                     DupString("ATL_tmp.c"));
   pM->blask = pB->blask = pK->blask = 1;
   pK->next = pM;
   pM->next = pB;

   pM->ku = pB->ku = 1;
   pM->kbB = pB->kbB = kb;
   pK->ku = VL;
   pK->kbB = ((kb+VL-1)/VL)*VL;
   free(pM->genstr); free(pB->genstr); free(pK->genstr);

   printf("Full search for SYRK kernel for nb=%d, NREG=%d, VLEN=%d\n",
          nb, nreg, VL);
   printf("   VD   NB   KB  NU      MFLOP\n");
   printf("   ==  ===  ===  ==  =========\n");
   for (i=1; i <= nreg; i++)
   {
      for (j=1; j <= nreg; j++)
      {
         if (i*j+1 > nreg || j > nb)
            continue;
         if (i*VL == j) /* legal M-vectorized SYRK kernel */
         {
            const int b = (nb > j) ? (nb/j)*j : j;
            double mf;

            pB->mu = pB->nu = pM->mu = pM->nu = j;
            pB->mbB = pB->nbB = pM->mbB = pM->nbB = b;
            if (i*j+i+j > nreg)
            {
               pB->flag |= 1<<MMF_BREG1;
               pM->flag |= 1<<MMF_BREG1;
            }
            else
            {
               pB->flag &= ~(1<<MMF_BREG1);
               pM->flag &= ~(1<<MMF_BREG1);
            }
            pM->genstr = MMGetGenString(pre, pM);
            pB->genstr = MMGetGenString(pre, pB);
            mf = TimeMMKernel(vrb, 0, pM, pre, b, b, kb, 1, 0, -1);
            assert(!MMKernelFailsTest(pre, b, b, kb, 1, pM));
            printf(frm, 'M', 'b', b, kb, j, mf);
            if (mf > mfB)
            {
               mfB = mf;
               uB = j;
               kvecB = 0;
               bcB = 0;
               kbB = kb;
            }
            mf = TimeMMKernel(vrb, 0, pB, pre, b, b, kb, 1, 0, -1);
            assert(!MMKernelFailsTest(pre, b, b, kb, 1, pB));
            printf(frm, 'M', 's', b, nb, j, mf);
            if (mf > mfB)
            {
               mfB = mf;
               uB = j;
               kvecB = 0;
               bcB = 1;
               kbB = kb;
            }
            free(pM->genstr);
            free(pB->genstr);
         }
         if (i == j)    /* legal k-vectorized SYRK kernel */
         {
            const int b = (nb > j) ? (nb/j)*j : j;
            double mf;

            pK->mbB = pK->nbB = b;
            pK->mu = pK->nu = j;
            if (i*j+i+j > nreg)
               pK->flag |= 1<<MMF_BREG1;
            else
               pK->flag &= ~(1<<MMF_BREG1);
            pK->genstr = MMGetGenString(pre, pK);
            mf = TimeMMKernel(vrb, 0, pK, pre, b, b, pK->kbB, 1, 0, -1);
            assert(!MMKernelFailsTest(pre, b, b, pK->kbB, 1, pK));
            printf(frm, 'K', ' ', b, nb, j, mf);
            if (mf > mfB)
            {
               mfB = mf;
               uB = j;
               kvecB = VL;
               bcB = 1;
               kbB = pK->kbB;
            }
            free(pK->genstr);
         }
      }
   }
   pM->genstr = pK->genstr = pB->genstr = NULL;
   KillAllMMNodes(pK);
   printf("Done.\n");
   pK = MMGetNodeGEN(pre, bcB|bf, 0, uB, uB, kvecB ? kvecB:1, VL, kvecB,
                      0, 0, NULL);
   if (uB*uB+uB+uB > nreg)
   {
      pK->flag |= 1<<MMF_BREG1;
      if (pK->genstr)
         free(pK->genstr);
      pK->genstr = MMGetGenString(pre, pK);
   }
   if (pre == 'z' || pre == 'c')
      pK->flag |= 1<<MMF_COMPLEX;
   pK->blask = 1;
   pK->mflop[0] = mfB;
   pK->mbB = pK->nbB = nb;
   pK->kbB = kbB;
   pK->blask = 1;
   return(pK);
}

void DoSyrkKU(int vrb, char pre, ATL_mmnode_t *mb)
/*
 * For KVEC, we currently force ku=vlen, so this routine does nothing.
 * For other kernels, see if 2 <= ku <= 4 provide speedup over default ku=1
 */
{
   int ku, kuB;
   double mf, mfB, mf0;
   if (FLAG_IS_SET(mb->flag, MMF_KVEC))
      return;
   mf0 = mfB = mb->mflop[0];
   printf("TUNING KU, CURRENTLY KU=%d, mf=%.2f\n", mb->ku, mfB);
   for (ku=2; ku <= 4; ku++)
   {
      char *gen0=mb->genstr;
      int ku0=mb->ku;
      mb->ku = ku;
      mb->genstr = MMGetGenString(pre, mb);
      mf = TimeMMKernel(vrb, 0, mb, pre, mb->mbB, mb->nbB, mb->kbB, 1, 0, -1);
      printf("ku=%u: mf=%.2f, SPEEDUP=%.4f\n", ku, mf, mf/mfB);
      if (mf >= 1.005*mfB) /* require further unrolling to prove .5% speedup */
      {
         free(gen0);
         mfB = mf;
      }
      else
      {
         free(mb->genstr);
         mb->genstr = gen0;
         mb->ku = ku0;
      }
   }
   if (mfB == mf0)
      printf("NO SPEEDUP, KEEPING KU=%u\n\n", mb->ku);
   else
      printf("KU=%u PROVIDES SPEEDUP=%.4f\n\n", mb->ku, mfB/mf0);
}

ATL_mmnode_t *FullSrchMUNU(int flg, int verb, char pre, int nreg, int nb,
                           int VL, int KVEC)
{
   char fn[32];
   ATL_mmnode_t *mmp;
   double mf, mfB=0.0;
   const int CHK=(flg&1), ku = (KVEC) ? VL : 1;
   int n, i, j, mbB, nbB, kbB, muB=1, nuB=1, b1B=0;
   char ch;
   const int bf=(flg&FKO_FLAG)?FKO_FLAG:0;

   assert(VL < 1000 && nreg < 1000);  /* don't overflow fn len */
   if (VL < 2)
      ch = 'U';
   else
      ch = (KVEC) ? 'K':'M';
   if (bf)
      sprintf(fn, "gAMMUR_%c%d_%d_fko.sum", ch, VL, nreg);
   else
      sprintf(fn, "gAMMUR_%c%d_%d.sum", ch, VL, nreg);
   mmp = ReadMMFileWithPath(pre, "res", fn);
   if (mmp)
   {
      MMFillInGenStrings(pre, mmp);
      TimeNegMMKernels(0, verb, 0, mmp, pre, 1, 0, -1);
      WriteMMFileWithPath(pre, "res", fn, mmp);
      return(mmp);
   }
   assert(nb%VL == 0);
   mmp = MMGetNodeGEN(pre, 0|bf, nb, 1, 1, ku, 1, KVEC, 0, 0,
                      DupString("ATL_Xamm_munu.c"));
   mmp->rout[4] = pre;
   mmp->mbB = mmp->nbB = mmp->kbB = nb;
   mmp->vlen = VL;
   if (KVEC)
      mmp->flag |= (1<<MMF_KVEC);
   else if (!UseBcast(flg, verb, pre, nb, nreg, VL))
      mmp->flag |= (1<<MMF_NOBCAST);
/*
 * Try all MU/NU unrollings
 */
   printf("Full search on MUxNU for nb=%d, NREG=%d, VLEN=%d, KVEC=%d\n",
          nb, nreg, VL, KVEC);
   for (i=1; i <= nreg; i++)
   {
      for (j=1; j <= nreg; j++)
      {
         int mbu, nbu, mu, nu;
         if (i*j+Mmin(i,j)+1 > nreg)
            continue;
         mu = (KVEC) ? i : i*VL;
         nu = j;
         if (i*j+i+j > nreg)
            mmp->flag |= 1<<MMF_BREG1;
         else
            mmp->flag &= ~(1<<MMF_BREG1);
         mmp->mu = mu;
         mmp->nu = nu;
         if (mmp->genstr)
           free(mmp->genstr);
         mbu = (nb >= mu) ? (nb/mu)*mu : mu;
         nbu = (nb >= nu) ? (nb/nu)*nu : nu;
         if (bf) mmp->flag |= 1<<MMF_FKO;
         mmp->genstr = MMGetGenString(pre, mmp);
         mf = TimeMMKernel(verb, 0, mmp, pre, mbu, nbu, nb, 1, 0, -1);
         printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
         if (CHK)
         {
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 1, mmp));
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 0, mmp));
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, -1, mmp));
         }
         if (mf > mfB)
         {
            mbB = mbu;
            nbB = nbu;
            kbB = nb;
            muB = mu;
            nuB = nu;
            mfB = mf;
            b1B = ((mmp->flag)>>MMF_BREG1)&1;
         }
      }
   }
DONE:
   assert(mfB > 0.0);
   i = FLAG_IS_SET(mmp->flag, MMF_NOBCAST);
   i |= (b1B) ? 2 : 0;
   KillMMNode(mmp);
   mmp = MMGetNodeGEN(pre, i|bf, nb, muB, nuB, ku, VL, KVEC, 0, 0, NULL);
   mmp->mbB = mbB;
   mmp->nbB = nbB;
   mmp->kbB = kbB;
   mmp->mflop[0] = mfB;
   printf("BEST FULL-SEARCH CASE IS B=(%d,%d,%d), U=(%d,%d) MFLOP=%.2f\n\n",
          mbB, nbB, kbB, muB, nuB, mfB);
   WriteRefreshedMMFileWithPath(pre, "res", fn, mmp);
   return(mmp);
}

ATL_mmnode_t *SrchNU(int flg, int verb, char pre, int nreg, int nb, int VL,
                     int I)
/*
 * M-vectorized search for with mu=I*VLEN.  It allows us to find a case
 * that can handle smaller blocks with VLEN is long.
 */
{
   char fn[32];
   ATL_mmnode_t *mmp;
   double mf, mfB=0.0;
   const int CHK=(flg&1), mu = I*VL;
   int n, i, j, mbB, nbB, kbB, nuB=1, b1B=0, mbu;
   const int bf=(flg&FKO_FLAG)?FKO_FLAG:0;

   if (bf)
      sprintf(fn, "gAMMUR_MU%d_M%d_%d_fko.sum", I, VL, nreg);
   else
      sprintf(fn, "gAMMUR_MU%d_M%d_%d.sum", I, VL, nreg);
   mmp = ReadMMFileWithPath(pre, "res", fn);
   if (mmp)
   {
      MMFillInGenStrings(pre, mmp);
      TimeNegMMKernels(0, verb, 0, mmp, pre, 1, 0, -1);
      WriteRefreshedMMFileWithPath(pre, "res", fn, mmp);
      return(mmp);
   }
   mmp = MMGetNodeGEN(pre, 0|bf, nb, 1, 1, 1, 1, 0, 0, 0,
                      DupString("ATL_Xamm_munu.c"));
   mmp->rout[4] = pre;
   mmp->mbB = mmp->nbB = mmp->kbB = nb;
   mmp->vlen = VL;
   mmp->mu = mu;
   mbu = (nb >= mu) ? (nb/mu)*mu : mu;
   if (!UseBcast(flg, verb, pre, nb, nreg, VL))
      mmp->flag |= (1<<MMF_NOBCAST);
/*
 * Try all powers of 2 MU/NU unrollings
 */
   printf("Searching M-VEC MU=%d xNU case for mb=%d, kb=%d, NREG=%d, VLEN=%d\n",
          I, mbu, nb, nreg, VL);
   for (j=1; j <= nreg; j++)
   {
      int nbu, nu;
      if (I*j+Mmin(I,j)+1 > nreg)
         continue;
      nu = j;
      mmp->nu = nu;
      if (mmp->genstr)
        free(mmp->genstr);
      nbu = (nb >= nu) ? (nb/nu)*nu : nu;
      if (I*j+I+j > nreg)
         mmp->flag |= 1<<MMF_BREG1;
      else
         mmp->flag &= ~(1<<MMF_BREG1);
      if (bf) mmp->flag |= 1<<MMF_FKO;
      mmp->genstr = MMGetGenString(pre, mmp);
      mf = TimeMMKernel(verb, 0, mmp, pre, mbu, nbu, nb, 1, 0, -1);
      printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", I, j, mf);
      if (CHK)
      {
         assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 1, mmp));
         assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 0, mmp));
         assert(!MMKernelFailsTest(pre, mbu, nbu, nb, -1, mmp));
      }
      if (mf > mfB)
      {
         mbB = mbu;
         nbB = nbu;
         kbB = nb;
         nuB = nu;
         mfB = mf;
         b1B = ((mmp->flag)>>MMF_BREG1)&1;
      }
   }
   i = FLAG_IS_SET(mmp->flag, MMF_NOBCAST);
   i |= (b1B) ? 2 : 0;
   KillMMNode(mmp);
   if (mfB == 0.0)  /* no legal kern found! */
   {
      printf("NO LEGAL KERNS FOR MU=%d!\n", I);
      return(NULL);
   }
   mmp = MMGetNodeGEN(pre, i|bf, nb, mu, nuB, 1, VL, 0, 0, 0, NULL);
   mmp->mbB = mbB;
   mmp->nbB = nbB;
   mmp->kbB = kbB;
   mmp->mflop[0] = mfB;
   assert(mfB > 0.0);
   printf("BEST MU=%d CASE IS B=(%d,%d,%d) U=(%d,%d), MFLOP=%.2f\n\n",
          I, mbB, nbB, kbB, mu, nuB, mfB);
   WriteRefreshedMMFileWithPath(pre, "res", fn, mmp);
   return(mmp);
}

ATL_mmnode_t *SrchMUNUp2(int flg, int verb, char pre, int nreg, int nb,
                         int VL, int KVEC)
{
   char fn[32];
   ATL_mmnode_t *mmp;
   double mf, mfB=0.0;
   const int CHK=(flg&1), ku = (KVEC) ? VL : 1;
   int n, i, j, mbB, nbB, kbB, muB=1, nuB=1, b1B=0;
   char ch;
   const int bf=(flg&FKO_FLAG)?FKO_FLAG:0;

   if (VL < 2)
      ch = 'U';
   else
      ch = (KVEC) ? 'K':'M';
   if (bf)
      sprintf(fn, "gAMMURP2_%c%d_%d_fko.sum", ch, VL, nreg);
   else
      sprintf(fn, "gAMMURP2_%c%d_%d.sum", ch, VL, nreg);
   mmp = ReadMMFileWithPath(pre, "res", fn);
   if (mmp)
   {
      MMFillInGenStrings(pre, mmp);
      TimeNegMMKernels(0, verb, 0, mmp, pre, 1, 0, -1);
      WriteRefreshedMMFileWithPath(pre, "res", fn, mmp);
      return(mmp);
   }
   mmp = MMGetNodeGEN(pre, 0|bf, nb, 1, 1, ku, 1, KVEC, 0, 0,
                      DupString("ATL_Xamm_munu.c"));
   mmp->rout[4] = pre;
   mmp->mbB = mmp->nbB = mmp->kbB = nb;
   mmp->vlen = VL;
   if (KVEC)
      mmp->flag |= (1<<MMF_KVEC);
   else if (!UseBcast(flg, verb, pre, nb, nreg, VL))
      mmp->flag |= (1<<MMF_NOBCAST);
/*
 * Try all powers of 2 MU/NU unrollings
 */
   printf("Searching PWR-2 MUxNU cases for nb=%d, NREG=%d, VLEN=%d, KVEC=%d\n",
          nb, nreg, VL, KVEC);
   for (i=1; i <= nreg; i += i)
   {
      for (j=1; j <= nreg; j += j)
      {
         int mbu, nbu, mu, nu;
         if (i*j+Mmin(i,j)+1 > nreg)
            continue;
         if (i*j+i+j > nreg)
            mmp->flag |= 1<<MMF_BREG1;
         else
            mmp->flag &= ~(1<<MMF_BREG1);
         mu = (KVEC) ? i : i*VL;
         nu = j;
         mmp->mu = mu;
         mmp->nu = nu;
         if (mmp->genstr)
           free(mmp->genstr);
         mbu = (nb >= mu) ? (nb/mu)*mu : mu;
         nbu = (nb >= nu) ? (nb/nu)*nu : nu;
         if (bf) mmp->flag |= 1<<MMF_FKO;
         mmp->genstr = MMGetGenString(pre, mmp);
         mf = TimeMMKernel(verb, 0, mmp, pre, mbu, nbu, nb, 1, 0, -1);
         printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
         if (CHK)
         {
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 1, mmp));
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 0, mmp));
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, -1, mmp));
         }
         if (mf > mfB)
         {
            mbB = mbu;
            nbB = nbu;
            kbB = nb;
            muB = mu;
            nuB = nu;
            mfB = mf;
            b1B = ((mmp->flag)>>MMF_BREG1)&1;
         }
      }
   }
DONE:
   assert(mfB != 0.0);
   i = FLAG_IS_SET(mmp->flag, MMF_NOBCAST);
   i |= (b1B) ? 2 : 0;
   KillMMNode(mmp);
   mmp = MMGetNodeGEN(pre, i|bf, nb, muB, nuB, ku, VL, KVEC, 0, 0, NULL);
   mmp->mbB = mbB;
   mmp->nbB = nbB;
   mmp->kbB = kbB;
   mmp->mflop[0] = mfB;
   printf("BEST POW2 CASE IS B=(%d,%d,%d) U=(%d,%d), MFLOP=%.2f\n\n",
          mbB, nbB, kbB, muB, nuB, mfB);
   WriteRefreshedMMFileWithPath(pre, "res", fn, mmp);
   return(mmp);
}

ATL_mmnode_t *SrchMUNU(int flg, int verb, char pre, int nreg, int nb,
                       int VL, int KVEC)
{
   ATL_mmnode_t *mmp, *mmp2;
   char fn[32];
   double mf, mfB=0.0;
   const int CHK=(flg&1), ku = (KVEC) ? VL : 1;
   #if (defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)) && !defined(ATL_AVX)
      int DO1D=1;
   #else
      int DO1D=(nreg < 9 || nreg < VL);
   #endif
   int n, i, j, kb, mbB, nbB, kbB, muB=1, nuB=1, b1B=0;
   char ch;
   const int bf=(flg&FKO_FLAG)?FKO_FLAG:0;

   if (flg&2)
      return(FullSrchMUNU(flg, verb, pre, nreg, nb, VL, KVEC));
   mmp2 = SrchMUNUp2(flg, verb, pre, nreg, nb, VL, KVEC);
   assert(VL < 1000 && nreg < 1000);  /* don't overflow fn len */
   if (VL < 2)
      ch = 'U';
   else
      ch = (KVEC) ? 'K':'M';
   if (bf)
      sprintf(fn, "gAMMUR_%c%d_%d_fko.sum", ch, VL, nreg);
   else
      sprintf(fn, "gAMMUR_%c%d_%d.sum", ch, VL, nreg);
   mmp = ReadMMFileWithPath(pre, "res", fn);
   if (mmp)
   {
      KillAllMMNodes(mmp2);
      MMFillInGenStrings(pre, mmp);
      TimeNegMMKernels(0, verb, 0, mmp, pre, 1, 0, -1);
      WriteMMFileWithPath(pre, "res", fn, mmp);
      return(mmp);
   }
   mmp = MMGetNodeGEN(pre, 0|bf, nb, 1, 1, ku, 1, KVEC, 0, 0,
                      DupString("ATL_Xamm_munu.c"));
   mmp->rout[4] = pre;
   mmp->mbB = mmp->nbB = mmp->kbB = nb;
   mmp->vlen = VL;
   if (KVEC)
      mmp->flag |= (1<<MMF_KVEC);
   else if (!UseBcast(flg, verb, pre, nb, nreg, VL))
      mmp->flag |= (1<<MMF_NOBCAST);
/*
 * Try all near-square register blocking cases
 */
   printf("Finding best MUxNU case for nb=%d, NREG=%d, VLEN=%d, KVEC=%d\n",
          nb, nreg, VL, KVEC);
   for (n=4; n <= nreg; n++)
   {
      int mbu, nbu, mu, nu;
      for (j=1; j*j < n; j++);
      i = n / j;
      if (nb%i || nb%j)
         continue;
      mu = (KVEC) ? i : i*VL;
      nu = j;
      if (i*j+i+j > nreg)
         mmp->flag |= 1<<MMF_BREG1;
      else
         mmp->flag &= ~(1<<MMF_BREG1);
      mmp->mu = mu;
      mmp->nu = nu;
      if (mmp->genstr)
        free(mmp->genstr);
      mbu = (nb >= mu) ? (nb/mu)*mu : mu;
      nbu = (nb >= nu) ? (nb/nu)*nu : nu;
      if (bf) mmp->flag |= 1<<MMF_FKO;
      mmp->genstr = MMGetGenString(pre, mmp);
      mf = TimeMMKernel(verb, 0, mmp, pre, mbu, nbu, nb, 1, 0, -1);
      printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
      if (CHK)
      {
         assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 1, mmp));
         assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 0, mmp));
         assert(!MMKernelFailsTest(pre, mbu, nbu, nb, -1, mmp));
      }
      if (mf > mfB)
      {
         mbB = mbu;
         nbB = nbu;
         kbB = nb;
         muB = mu;
         nuB = nu;
         mfB = mf;
         b1B = ((mmp->flag)>>MMF_BREG1)&1;
      }
   }
/*
 * For non-AVX x86, try 1-D cases since they are 2-operand assemblies; always
 * try 1-D for low registers
 */
   if (DO1D)
   {
      printf("BEST NEAR-SQUARE CASE IS MU=%d, NU=%d, MFLOP=%.2f\n\n",
             muB, nuB, mfB);
      printf("Finding best 1-D outer loop unrolling for nb=%d\n", nb);
      for (n=2; n <= nreg; n++)
      {
         int mbu, nbu, mu, nu;
         i = 1; j = n;
         if (nb % n)
            continue;
         mu = (KVEC) ? i : i*VL;
         nu = mmp->nu = j;
         if (mmp->genstr)
           free(mmp->genstr);
         if (bf) mmp->flag |= 1<<MMF_FKO;
         mmp->genstr = MMGetGenString(pre, mmp);
         mbu = (nb >= mu) ? (nb/mu)*mu : mu;
         nbu = (nb >= nu) ? (nb/nu)*nu : nu;
         mf = TimeMMKernel(verb, 0, mmp, pre, mbu, nbu, nb, 1, 0, -1);
         printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
         if (CHK)
         {
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 1, mmp));
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 0, mmp));
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, -1, mmp));
         }
         if (mf > mfB)
         {
            muB = i;
            nuB = j;
            mfB = mf;
            b1B = 0;
         }
         i = n; j = 1;
         mu = mmp->mu = (KVEC) ? i : i * VL;
         nu = mmp->nu = j;
         mbu = (nb >= mu) ? (nb/mu)*mu : mu;
         nbu = (nb >= nu) ? (nb/nu)*nu : nu;
         if (mmp->genstr)
           free(mmp->genstr);
         if (bf) mmp->flag |= 1<<MMF_FKO;
         mmp->genstr = MMGetGenString(pre, mmp);
         mf = TimeMMKernel(verb, 1, mmp, pre, mbu, nbu, nb, 1, 0, -1);
         printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
         if (CHK)
         {
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 1, mmp));
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, 0, mmp));
            assert(!MMKernelFailsTest(pre, mbu, nbu, nb, -1, mmp));
         }
         if (mf > mfB)
         {
            muB = i;
            nuB = j;
            mfB = mf;
            b1B = 0;
         }
      }
   }

   assert(mfB > 0.0);
   i = FLAG_IS_SET(mmp->flag, MMF_NOBCAST);
   i |= (b1B) ? 2 : 0;
   KillMMNode(mmp);
   mmp = MMGetNodeGEN(pre, i|bf, nb, muB, nuB, ku, VL, KVEC, 0, 0, NULL);
   mmp->mbB = mbB;
   mmp->nbB = nbB;
   mmp->kbB = kbB;
   mmp->mflop[0] = mfB;
   if (mmp->mflop[0] < mmp2->mflop[0])
   {
      printf("Taking pow2 srch (%d,%d:%.0f) over square (%d,%d:%.0f)\n",
             mmp2->mu, mmp2->nu, mmp2->mflop[0],
             mmp->mu, mmp->nu, mmp->mflop[0]);
      KillMMNode(mmp);
      mmp = mmp2;
   }
   else
      KillAllMMNodes(mmp2);
   WriteRefreshedMMFileWithPath(pre, "res", fn, mmp);
   printf("BEST CASE IS B=(%d,%d,%d), U=(%d,%d), MFLOP=%.2f\n\n",
          mmp->mbB, mmp->nbB, mmp->kbB, mmp->mu, mmp->nu, mmp->mflop[0]);
   return(mmp);
}

#if 0
/* this idea doesn't really work */
int FindKvecXover(int flg, int verb, char pre, int nreg, int VL, int nb)
/*
 * On some machines, vectorizing the M dim will win for small problems,
 * due to not needing to the summation at the end of the K-loop.  However,
 * once this cost is dominated by the K-loop, K dim vectorization can
 * start to win, possibly by reducing the C write traffic, as well as
 * avoiding vec bcast inside the K-loop.
 *
 * For finding this crossover, we use best pw2 M/K, since we can be sure
 * they can always use a pwr2 block factor for direct comparison.  If pwr2
 * cases are much different than normal, this may cause problems!
 * IDEA: Read in normal MUNU results, and don't use this test if gap is wide.
 */
{
   ATL_mmnode_t *mmM, *mmK;
   double mfM, mfK;
   int b, b0;

   mmM = SrchMUNUp2(flg, verb, pre, nreg, nb, VL, 0);
   mmK = SrchMUNUp2(flg, verb, pre, nreg, nb, VL, 1);
   printf("FINDING NB CROSSOVER FOR M- AND K-VECTORIZATION:");
   b = Mmax(mmM->mu, mmK->mu);
   b = Mmax(b, mmM->nu);
   b = Mmax(b, mmK->nu);
   b = Mmax(b,16);
   b0 = b;
   while (b < 512)
   {
      mfM = TimeMMKernel(verb, 0, mmM, pre, b, b, b, 1, 0, -1);
      mfK = TimeMMKernel(verb, 0, mmK, pre, b, b, b, 1, 0, -1);
      printf("   B=%d, mflopM=%.0f, mflopK=%.0f\n", b, mfM, mfK);
      if (mfK > mfM*1.02)
         break;
      b += b;
   }
   if (b == b0)
   {
      printf("K-VECTORIZATION BETTER FROM FIRST BLOCK!\n");
      b = 1;            /* Kvec always better */
   }
   else if (b == 512)
   {
      printf("M-VECTORIZATION ALWAYS BETTER\n");
      b = 0;            /* Kvec never better */
   }
   else
      printf("K-VEC BEGINS WINNING AROUND %d\n", b);
   return(b);
}
#endif

void FindDefMUNU(int flg, int verb, char pre, int nb, int *NREG, int *VLEN)
{
   ATL_mmnode_t *mp, *mmM, *mmK;
   int nreg=(*NREG), VL=(*VLEN), chkNR=0, chkVL=0;

   if (nreg < 1)
      nreg = GetNumVecRegs(pre);
   if (nreg < 1)
   {
      #ifdef ATL_GAS_x8632
         nreg = 8;
      #else
         nreg = 16;
      #endif
      chkNR = 1;
   }
   if (VL < 1)
      VL = GetNativeVLEN(pre);
   if (!VL)
   {
      VL = (pre == 'c' || pre == 's') ? 4:2;
      chkVL = 1;
   }
/*
 * Always do full search for low number of registers, where this is only search
 */
   if (!chkNR && nreg > 0 && nreg < 20)
      flg |= 2;
   mmM = SrchMUNU(flg, verb, pre, nreg, nb, VL, 0);
   mmK = SrchMUNU(flg, verb, pre, nreg, nb, VL, 1);
   printf("MVEC: B=(%d,%d,%d) mu=%d, nu=%d, MFLOP=%.0f\n",
          mmM->mbB,  mmM->nbB,  mmM->kbB, mmM->mu, mmM->nu, mmM->mflop[0]);
   printf("KVEC: B=(%d,%d,%d) mu=%d, nu=%d, MFLOP=%.0f\n",
          mmK->mbB,  mmK->nbB,  mmK->kbB, mmK->mu, mmK->nu, mmK->mflop[0]);
/*
 * After this, fastest code in mmM, slowest mmK
 */
   if (mmK->mflop[0] > mmM->mflop[0])
   {
      mp = mmK;
      mmK = mmM;
      mmM = mp;
   }
   KillAllMMNodes(mmK);
/*
 * If we only guessed a lower bound on # regs, try some searches with
 * increasing regs
 */
   if (chkNR)
   {
      const int KVEC = FLAG_IS_SET(mmM->flag, MMF_KVEC);
      int i, nr = nreg+nreg;
      printf("NREG=%d, U=(%d,%d): MFLOP=%.0f\n",
             nreg, mmM->mu, mmM->nu, mmM->mflop[0]);
      for (i=0; i < 4; i++)  /* sanity check for stopping */
      {
         mp = SrchMUNU(flg, verb, pre, nr, nb, VL, KVEC);
         printf("NREG=%d, U=(%d,%d): MFLOP=%.0f\n",
                nr, mp->mu, mp->nu, mp->mflop[0]);
         if (mp->mu*mp->nu + mp->mu + 1 <= (nr>>1)) /* did not use more regs */
            break;
         if (mp->mflop[0] < 1.03*mmM->mflop[0])     /* perf not better */
            break;
         KillMMNode(mmM);
         mmM = mp;
         nreg = nr;
         nr += nr;
      }
   }
/*
 * Now that we are confident in our NREG, see if we need to confirm VLEN
 */
   if (chkVL)
   {
   }
   *NREG = nreg;
   *VLEN = VL;
#if 0
/*
 * Now see if K-vec has a crossover with M-vec
 */
   FindKvecXover(flg, verb, pre, nreg, VL, nb);
#endif
   KillAllMMNodes(mmM);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set type/precision prefix (d) \n");
   fprintf(stderr, "   -r <nreg> : set max # of registers to try\n");
   fprintf(stderr, "   -V <vlen> : force vector length\n");
   fprintf(stderr, "   -b <nb>   : set initial block factor to try\n");
   fprintf(stderr, "   -v <verb> : set verbosity (1)\n");
   fprintf(stderr, "   -T 1      : test all legal kerns, and exit\n");
   fprintf(stderr,
           "   -f <flg>  : bitvec for srch control, add vals you want set:\n");
   fprintf(stderr, "        1: test all generated kernels\n");
   fprintf(stderr, "        2: do full MUxNU search\n");
   fprintf(stderr, "        4: print # of regs to res/<pre>nreg\n");
   fprintf(stderr, "        8: pref tune & time kerns in time.sum & exit\n");
   fprintf(stderr, "       16: create [d,s]flops.frc and return\n");
   fprintf(stderr, "   %d: use iFKO generator\n", FKO_FLAG);
   fprintf(stderr, "      DEFAULT: all bits unset\n");
   exit(ierr ? ierr : -1);
}

void GetFlags(int nargs, char **args, int *FLG, int *VERB, char *PRE, int *NREG,
              int *VLEN, int *NB, int *TEST)
{
   int i, flg=0, nreg=0;
   char pre = 'd';

   *VERB = 0;
   *NB = 120;
   *VLEN = 0;
   *TEST = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'p':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        pre = tolower(args[i][0]);
        assert(pre == 's' || pre == 'd' || pre == 'z' || pre == 'c');
        break;
      case 'T':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *TEST = atoi(args[i]);
         break;
      case 'V':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *VLEN = atoi(args[i]);
         break;
      case 'f':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         flg = atoi(args[i]);
         break;
      case 'r':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         nreg = atoi(args[i]);
         break;
      case 'v':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *VERB = atoi(args[i]);
         break;
      case 'b':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *NB = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!nreg)
      nreg = GetNumVecRegs(pre);
   if (pre == 's' || pre == 'c' && (*NB)%16)
      *NB = (*NB == 120) ? 240 : (((*NB+15)>>4)<<4);
   *PRE = pre;
   *FLG = flg;
   *NREG = nreg;
}

/*
 * Using best discovered kernel, figure out the largest NB < 512 that
 * gets good performance
 */
int GetMaxNB(int flag, int verb, char pre, ATL_mmnode_t *mp)
{
   int inc, i, bB=0, badrow=0;
   double mf, mfB=0.0;

   printf("FINDING RANGE OF NB FOR GENKERN MU=%d, NU=%d, %cVEC=%d:\n",
          mp->mu, mp->nu, FLAG_IS_SET(mp->flag, MMF_KVEC)?'K':'M', mp->vlen);
   inc = Mylcm(mp->mu, mp->nu);
   inc = Mylcm(inc,mp->ku);
   while (inc < 12)
      inc += inc;

   for (i=inc; i < 512; i += inc)
   {
       mf = TimeMMKernel(verb, 0, mp, pre, i, i, i, 1, 0, -1);
       printf("   NB=%d, mf=%.0f\n", i, mf);
       if (mf > mfB)
       {
          bB=i;
          mfB = mf;
       }
       else badrow++;
       if (badrow > 4)
          break;
   }
   mp->mbB = mp->kbB = mp->nbB = bB;
   mp->mflop[0] = mfB;
   printf("BEST SQUARE NB=%d (%.0f)\n", bB, mfB);
   return(bB+inc-1);
}

void FindInfo(int flag, int verb, char pre, int NB, int *NREG, int *VLEN)
{
   const int WNR=(flag&4);
   int nreg=(*NREG), vlen=(*VLEN);
   ATL_mmnode_t *mp, *mpN, *mpB;
   const int bf=(flag&FKO_FLAG)?FKO_FLAG:0;
   char *mvsum, *kvsum;
   char fn[16]={pre,'m','f','l','o','p','s','.','f','r','c','\0'};

   if (pre == 'z')
      pre = 'd';
   else if (pre == 'c')
      pre = 's';
   if (bf)
   {
      mvsum = "gmvAMMUR_fko.sum";
      kvsum = "gkvAMMUR_fko.sum";
   }
   else
   {
      mvsum = "gmvAMMUR.sum";
      kvsum = "gkvAMMUR.sum";
   }
   mp = ReadMMFileWithPath(pre, "res", mvsum);
   mpN = ReadMMFileWithPath(pre, "res", kvsum);
   if (mp && mpN)
   {
      FILE *fp;
      *VLEN = mp->vlen;
      *NREG = mp->ivar;
      MMFillInGenStrings(pre, mp);
      MMFillInGenStrings(pre, mpN);
      TimeNegMMKernels(0, verb, 0, mp, pre, 1, 0, -1);
      TimeNegMMKernels(0, verb, 0, mpN, pre, 1, 0, -1);
      WriteRefreshedMMFileWithPath(pre, "res", mvsum, mp);
      WriteRefreshedMMFileWithPath(pre, "res", kvsum, mpN);
      mpB = FindMaxMflopMMQ(mp, 0);
      mpB = CloneMMNode(mpB);
      KillAllMMNodes(mp);
      mp = FindMaxMflopMMQ(mpN, 0);
      if (mp->mflop[0] > mpB->mflop[0])
         CopyMMNode(mpB, mp);
      KillAllMMNodes(mpN);
      if (WNR)
      {
         FILE *fp;
         char fn[12];
         sprintf(fn, "res/%cnreg", pre);
         fp = fopen(fn, "w");
         fprintf(fp, "%d\n", *NREG);
         fclose(fp);
      }
      fp = fopen(fn, "r");
      if (!fp)
         goto FIND_REPS;
      KillMMNode(mpB);
      return;
   }
   else if (mp)
      KillAllMMNodes(mp);
   else if (mpN)
      KillAllMMNodes(mpN);
   FindDefMUNU(flag, verb, pre, NB, &nreg, &vlen);
   printf("\nNREG=%d, VLEN=%d\n", nreg, vlen);
   if (NB%vlen)
   {
      if ((NB+NB)%vlen == 0)
         NB += NB;
      else
         while((++NB)%vlen);
   }
/*
 * With nreg & vlen set, create output with best K- & M- vectorized code
 * in standard names with nreg in mp->ivar
 */
   mp = SrchMUNU(flag, verb, pre, nreg, NB, vlen, 0);
   mpN = SrchMUNUp2(flag, verb, pre, nreg, NB, vlen, 0);
   mpB = (mpN->mflop[0] > mp->mflop[0]) ? mpN : mp;
   mpN->next = SrchNU(flag, verb, pre, nreg, NB, vlen, 1);
   mpB = (mpB->mflop[0] >= mpN->next->mflop[0]) ? mpB : mpN->next;
   if (mpN->next)
   {
      mpN->next->next = SrchNU(flag, verb, pre, nreg, NB, vlen, 2);
      mpB = (mpB->mflop[0] >= mpN->next->next->mflop[0]) ? mpB:mpN->next->next;
      if (mpN->next->next)
          mpN->next->next->next = SrchNU(flag, verb, pre, nreg, NB, vlen, 3);
   }
   mpB = CloneMMNode(mpB);
   mp = AddUniqueMMKernCompList(mp, mpN);
   KillAllMMNodes(mpN);
   mp = ReverseMMQ(mp);
   mp->ivar = nreg;
   WriteRefreshedMMFileWithPath(pre, "res", mvsum, mp);
   KillAllMMNodes(mp);

   mp = SrchMUNU(flag, verb, pre, nreg, NB, vlen, 1);
   if (mp->mflop[0] > mpB->mflop[0])
      CopyMMNode(mpB, mp);
   mp->ivar = nreg;
   mpN = SrchMUNUp2(flag, verb, pre, nreg, NB, vlen, 1);
   if (mpN->mflop[0] > mpB->mflop[0])
      CopyMMNode(mpB, mpN);
   mp = AddUniqueMMKernCompList(mp, mpN);
   KillAllMMNodes(mpN);
   mp = ReverseMMQ(mp);
   WriteRefreshedMMFileWithPath(pre, "res", kvsum, mp);
   KillAllMMNodes(mp);
   if (WNR)
   {
      FILE *fp;
      char fn[12];
      sprintf(fn, "res/%cnreg", pre);
      fp = fopen(fn, "w");
      fprintf(fp, "%d\n", nreg);
      fclose(fp);
   }
   *VLEN = vlen;
   *NREG = nreg;
/*
 * Now, use best-performing kernel found so far to set target MFLOP timing
 * interval for all subsequent timings
 */
   FIND_REPS:
   {
      double tim, nrep=1.0, mf;
      FILE *fp;
      printf(
      "FINDING NUMBER OF FLOPS TO FORCE FOR .2 SECOND TIMING INTERVAL:\n");
      mf = (mpB->mbB*1e-6)*mpB->nbB*2.0*mpB->kbB;
      tim = mf / mpB->mflop[0];
      if (tim < 0.2)
         nrep = .22 / tim;
      fprintf(stdout, "\nINCREASING TIMING INTERVAL BY: %.2f, MFLOP=%e!\n\n",
              nrep, nrep*mf);
      fp = fopen(fn, "w");
      fprintf(fp, "%le", nrep*mf);
      fclose(fp);
   }
}

ATL_mmnode_t *GetBestKernVT(char pre, char vt, int flg)
{
   ATL_mmnode_t *mp;
   char *mvsum, *kvsum;
   const int bf=(flg&FKO_FLAG)?FKO_FLAG:0;

   if (bf)
   {
      mvsum = "gmvAMMUR_fko.sum";
      kvsum = "gkvAMMUR_fko.sum";
   }
   else
   {
      mvsum = "gmvAMMUR.sum";
      kvsum = "gkvAMMUR.sum";
   }
   if (vt == 'K')
      mp = ReadMMFileWithPath(pre, "res", kvsum);
   else
      mp = ReadMMFileWithPath(pre, "res", mvsum);
   assert(mp);
   if (mp->next)
   {
      KillAllMMNodes(mp->next);
      mp->next = NULL;
   }
   return(mp);
}

ATL_mmnode_t *GetBestKern(char pre, int flg)
{
   ATL_mmnode_t *mp, *mpB;
   mpB = GetBestKernVT(pre, 'M', flg);
   mp =  GetBestKernVT(pre, 'K', flg);
   if (mp->mflop[0] > mpB->mflop[0])
   {
      KillAllMMNodes(mpB);
      mpB = mp;
   }
   else
      KillAllMMNodes(mp);
   return(mpB);
}


void DoSquare(int flag, int verb, char pre, int nreg, int VL)
{
   ATL_mmnode_t *mp;
   int maxNB;

   mp = GetBestKern(pre, flag);
   maxNB = GetMaxNB(flag, verb, pre, mp);
}

int TestWithKU(int verb, char pre, int NB, ATL_mmnode_t *mp, FILE *fperr)
{
   int nf;
   mp->flag |= (1<<MMF_KUISKB);
   mp->kbB = mp->ku = NB;
   free(mp->rout);
   free(mp->genstr);
   mp->rout = MMGetGenName(pre, NB, mp);
   mp->genstr = MMGetGenString(pre, mp);
   nf = NumberBetaFails(fperr, pre, NB, mp);
   return(nf);
}

int CountFails(int TEST, int flg, int verb, char pre, int NB, int nreg, int VL)
{
   int i, j, ntest=0, nfail=0;
   char *frm="%8d %4d %3d %3d %3d  %2d   %c %2d  %5d\n";
   FILE *fperr;
   const int bf=(flg&FKO_FLAG)?FKO_FLAG:0;

   fperr = fopen("res/FAIL.OUT", "w");
   assert(fperr);

   assert(VL);
   assert(nreg);
   assert(NB > 0);
   printf("     NUM    B  MU  NU  KU  VL VEC BC  NPASS\n");
   printf("======== ==== === === ===  == === ==  =====\n");
   for (i=1; i <= nreg; i++)
   {
      for (j=1; j <= nreg; j++)
      {
         ATL_mmnode_t *mp;
         int nf, flg=0;
         if (i == 1 || j == 1)
         {
            if (i*j+1 > nreg)
               continue;
         }
         else if (i*j+Mmin(i,j)+1 > nreg)
               continue;
         else if (i*j+i+j > nreg)
            flg = 2;

         mp = MMGetNodeGEN(pre, flg|bf, NB, i*VL, j, 1, VL, 0, 0, 0, NULL);
         nf = NumberBetaFails(fperr, pre, NB, mp);
         printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'M', 0, 3-nf);
         nfail += nf;
         ntest += 3;
         if (!nf)  /* try fully unrolled case if rolled worked*/
         {
            nf = TestWithKU(verb, pre, NB, mp, fperr);
            printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'M', 0, 3-nf);
            nfail += nf;
            ntest += 3;
         }
         mp = KillMMNode(mp);
         if (j%VL == 0 && !bf)  /* try m-vec w/o bcast */
         {
            mp = MMGetNodeGEN(pre, 1, NB, i*VL, j, 1, VL, 0, 0, 0, NULL);
            nf = NumberBetaFails(fperr, pre, NB, mp);
            printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'M',
                   FLAG_IS_SET(mp->flag, MMF_NOBCAST), 3-nf);
            nfail += nf;
            ntest += 3;
            if (!nf)  /* try fully unrolled case if rolled worked*/
            {
               nf = TestWithKU(verb, pre, NB, mp, fperr);
               printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'M',
                      FLAG_IS_SET(mp->flag, MMF_NOBCAST), 3-nf);
               nfail += nf;
               ntest += 3;
            }
            mp = KillMMNode(mp);
         }
         if (1 /*(i*j)%VL == 0*/) /* try k-vec */
         {
            mp = MMGetNodeGEN(pre, flg|bf, NB, i, j, VL, VL, 1, 0, 0, NULL);
            nf = NumberBetaFails(fperr, pre, NB, mp);
            printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'K', 0, 3-nf);
            nfail += nf;
            ntest += 3;
            if (!nf)  /* try fully unrolled case if rolled worked*/
            {
               nf = TestWithKU(verb, pre, NB, mp, fperr);
               printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'K', 0, 3-nf);
               nfail += nf;
               ntest += 3;
            }
            mp = KillMMNode(mp);
         }
      }
   }
   if (TEST > 1)
   {
      printf("\n   .... start of prefetch tests ....\n");
/*
 *    Now test some pref on both 1- & 2-D kernels.  Prefetch of next block of
 *    A & B affect peeling, and so can cause errors (K-loop prefetch no).
 *    Fully unrolled cases don't need to be tested again, since they are
 *    totally peeled with or without prefetch.
 */
      assert(nreg > 2);
      for (i=0; i < 7; i++)  /* 1-D in both directions, and some 2-D */
      {
         const int NOBCAST = (i == 5 || i == 6); /* spec cases for no-bcast */
         int MU;
         int mu, nu, pfLS;
         if (i == 0)
         {
            mu = nreg - 2;
            nu = 1;
         }
         else if (i == 1)
         {
            mu = 1;
            nu = nreg - 2;
            if (nu > VL)
               nu = (nu / VL) * VL;
         }
         else if (i == 2)
         {
            mu = 3;
            for (nu=3; nu*mu+mu+1 < nreg; nu++);
         }
         else if (i == 4)
            mu = nu = 2;
         else if (i == 5)
         {
            nu = VL;
            mu = 3;
         }
         else if (i == 6)
         {
            nu = VL*2;
            mu = 2;
         }

         if (i <= 2)  /* 1-D case */
         {
            if (mu*nu+1 > nreg)
               continue;
         }
         else if (mu*nu+Mmin(mu,nu)+1 > nreg)
            continue;

         MU = mu*VL;

         for (pfLS=16; pfLS < 128; pfLS += pfLS)
         {
            int k, pfs[6]={1, 2, 4, 1|2, 2|4, 1|2|4};
            for (k=0; k < 5; k++)
            {
               ATL_mmnode_t *mp;
               int nf;
/*
 *             Test M-vec using broadcast
 */
               if (!NOBCAST)
               {
                  mp = MMGetNodeGEN(pre, 0|bf, NB, MU, nu, 1, VL, 0,
                                    pfs[k], pfLS, NULL);
                  nf = NumberBetaFails(fperr, pre, NB, mp);
                  printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'M',
                         FLAG_IS_SET(mp->flag, MMF_NOBCAST), 3-nf);
                  nfail += nf;
                  ntest += 3;
                  KillMMNode(mp);
               }
               if (nu%VL == 0 && !bf)  /* try m-vec using splat & w/o bcast */
               {
                  mp = MMGetNodeGEN(pre, 1, NB, MU, nu, 1, VL, 0,
                                    pfs[k], pfLS, NULL);
                  nf = NumberBetaFails(fperr, pre, NB, mp);
                  printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'M',
                         FLAG_IS_SET(mp->flag, MMF_NOBCAST), 3-nf);
                  nfail += nf;
                  nfail += nf;
                  ntest += 3;
                  KillMMNode(mp);
               }
               if (!NOBCAST && (mu*nu)%VL == 0)  /* try k-vec */
               {
                  mp = MMGetNodeGEN(pre, 0|bf, NB, mu, nu, VL, VL, 1,
                                    pfs[k], pfLS, NULL);
                  nf = NumberBetaFails(fperr, pre, NB, mp);
                  printf(frm, ntest, NB, mp->mu, mp->nu, mp->ku, VL, 'K',
                         0, 3-nf);
                  nfail += nf;
                  ntest += 3;
                  KillMMNode(mp);
               }
            }
         }
      }
   }
   if (!nfail)
   {
      printf("ALL %d TESTS PASS!\n", ntest);
      fprintf(fperr, "ALL %d TESTS PASS!\n", ntest);
   }
   else
   {
      printf("FAILED %d OF %d TESTS!!\n\n", nfail, ntest);
      fprintf(fperr, "FAILED %d OF %d TESTS!!\n\n", nfail, ntest);
   }
   fclose(fperr);
   return(nfail);
}

int FindBlockingRegions(int flag, int verb, char pre)
/*
 * This routine attempts to find the best prefetch kernel for all problem sizes
 * using timings of the fastest unprefetched kernel found in MU/NU search.
 * We can prefetch or not each of the three operands, and for each we an
 * also use ATL_pfl1 (pref to L1) or ATL_pfl2.  L2 may really be last-lvl cache.
 * We assume there are five operand size ranges of interest:
 * 1. nb <= sqrt(L1sz/5) : can fit 5 blks of A/B in L1.  May be worthwhile
 *    to prefetch next block of A & B to L1 cache
 * 2. nb <= sqrt(L1sz/4) : fit 4 blks, try fetching A or B to L1, other to L2
 * 3. nb <= sqrt(L2sz/6) : pref next A&B blocks to L2
 * 4. nb <= sqrt(L2sz/5) : pref only one of A/B to L2
 * 5. else only do inter-block prefetch
 * RETURNS: maximum NB providing speedup
 */
{
   int maxNB=0;

   return(maxNB);
}

int TryKU(char pre, int verb, ATL_mmnode_t *mp, int imf, int ku)
/*
 * Assumes mp->mflop[imf] presently holds perf of present ku unrolling.
 * Times (& tests) unrolling of ku, and retains best setting.
 * Penalizes larger unroll very slightly due to potential code size problems.
 * RETURNS: 1 if new ku is faster, else 0.
 */
{
   const double newmul = (mp->ku < ku) ? 0.99 : 1.01; /* bias big unroll */
   double mf0=mp->mflop[0], mfN;
   char *gs0=mp->genstr;
   int ku0 = mp->ku, fail;

   mp->ku = ku;
   if (mp->rout)
      free(mp->rout);
   mp->rout = MMGetGenName(pre, mp->kbB, mp);
   mp->genstr = MMGetGenString(pre, mp);
   fail = MMKernelFailsAnyBeta(pre, mp->mbB, mp->nbB, mp->kbB, mp);
   if (fail)
      printf("   ku=%d FAILS TESTS!\n", ku);
   else
   {
      mfN = TimeMMKernel(verb, 0, mp, pre, mp->mbB, mp->nbB, mp->kbB, 1, 0, -1);
      printf("      ku=%d, MFLOP=%.2f, SPEEDUP=%.3f\n", ku, mfN, mfN/mf0);
   }
   if (!fail && mfN*newmul > mf0)  /* is new timing better? */
   {
      if (gs0)
         free(gs0);
      mp->mflop[0] = mfN;
      if (ku == mp->kbB)
         mp->flag |= (1<<MMF_KUISKB);
      return(1);
   }
   else  /* original ku best */
   {
      free(mp->genstr);
      free(mp->rout);
      mp->ku = ku0;
      mp->genstr = gs0;
      mp->rout = MMGetGenName(pre, mp->kbB, mp);
   }
}

int TryPref(char pre, int verb, int vrb, ATL_mmnode_t *mp, int ibet, int imf,
            int ipf, int pfLS)
/*
 * Times mmb with present prefetch setting, and (ipf,pfLS), and returns
 * in mmb the best.  mp->mflop[imf] must contain the timing for the present
 * prefetch settings.
 * RETURNS: 1 if new settings taken, 0 if old are kept
 */
{
   char *gs0=mp->genstr, *rt0=mp->rout;
   int ipf0=mp->pref, pfLS0=mp->pfLS;
   double mf0=mp->mflop[imf], mfN;

   if (mp->pref == ipf && mp->pfLS == pfLS)
   {
      if (vrb)
         printf("      pref settings unchanged!\n");
      return(0);
   }
   mp->pref = ipf;
   mp->pfLS = pfLS;
   mp->rout = DupString("ATL_tmp.c");
   mp->genstr = MMGetGenString(pre, mp);
   mfN = TimeMMKernel(verb, 0, mp, pre, mp->mbB, mp->nbB, mp->kbB, ibet, 0, -1);
   if (vrb)
      printf("      old=(%d,%d,%.0f), new=(%d,%d,%.0f), spdup=%.3f\n",
             ipf0, pfLS0, mf0, mp->pref, mp->pfLS, mfN, mfN/mf0);
   if (mfN > mf0)  /* new settings win */
   {
      if (gs0)
         free(gs0);
      if (rt0)
         free(rt0);
      mp->mflop[imf] = mfN;
      return(1);
   }
   else  /* old settings best */
   {
      free(mp->genstr);
      mp->genstr = gs0;
      free(mp->rout);
      mp->rout = rt0;
      mp->pref = ipf0;
      mp->pfLS = pfLS0;
   }
   return(0);
}

void DoKUs(char pre, int flg, int verb, ATL_mmnode_t *mmb)
/*
 * Given routs times with KU=1 (mflop stored at imf), try various KU settings
 * and return the best
 */
{
   ATL_mmnode_t *mp;
   printf("TUNING KU FOR ALL KERNELS:\n");
   for (mp=mmb; mp; mp = mp->next)
   {
      int kumax, kumul, ku;
      const double mf0 = mp->mflop[0];

      printf("   TRYING KUs, RT='%s' B=(%d,%d,%d)\n", mp->rout,
             mp->mbB, mp->nbB, mp->kbB);
      printf("      ku=%d, MFLOP=%.2f, SPEEDUP=1.0\n", mp->ku, mf0);
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         assert(mp->ku == mp->vlen);
      else
         assert(mp->ku == 1);
      #if 0  /* present generator handles only ku=1,full */
      kumax = (mp->kbB)>>1;
      kumax = Mmin(kumax, 128);
      kumul = Mmin(mp->mu, mp->nu); /* square-friendly K unrollings */
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         kumul = Mylcm(kumul, mp->vlen);
      else
         kumul = (kumul != 1) ? kumul : 2;

      for (ku=kumul; ku <= kumax; ku += kumul)
         TryKU(pre, verb, mp, 0, ku);
      #endif
      TryKU(pre, verb, mp, 0, mp->kbB);

      printf("   DONE: KU=%d gives MFLOP=%.0f, SPEEDUP=%.3f\n",
             mp->ku, mp->mflop[0], mp->mflop[0]/mf0);
   }
   printf("DONE TUNING KU FOR ALL KERNELS.\n\n");
}

void PrintPref(FILE *fp, int ipf)
{
   if (ipf & 2)
   {
      if (ipf & 64)
         fprintf(fp, " Ab to L1");
      else if (ipf & 512)
         fprintf(fp, " Ab to LLC");
      else
         fprintf(fp, " Ab to L2");
   }
   else
      fprintf(fp, " NO Ab PREF");
   if (ipf & 4)
   {
      if (ipf & 128)
         fprintf(fp, ", Bb to L1");
      else if (ipf & 1024)
         fprintf(fp, ", Bb to LLC");
      else
         fprintf(fp, ", Bb to L2");
   }
   else
      fprintf(fp, ", NO Bb PREF");
   if (ipf & 8)
   {
      if (ipf & 2048)
         fprintf(fp, ", Ak to L1");
      else if (ipf & 8192)
         fprintf(fp, ", Ak to LLC");
      else
         fprintf(fp, ", Ak to L2");
   }
   else
      fprintf(fp, ", NO Ak PREF");
   if (ipf & 16)
   {
      if (ipf & 4096)
         fprintf(fp, ", Bk to L1");
      else if (ipf & 16384)
         fprintf(fp, ", Bk to LLC");
      else
         fprintf(fp, ", Bk to L2");
   }
   else
      fprintf(fp, ", NO Bk PREF");
   if (ipf & 1)
   {
      if (ipf & 32)
         fprintf(fp, ", C to L1");
      else if (ipf &256)
         fprintf(fp, ", C to LLC");
      else
         fprintf(fp, ", C to L2");
   }
   else
      fprintf(fp, ", NO C PREF");
}

void DoPref1(char pre, int flg, int verb, int vrb0, ATL_mmnode_t *mmb)
/*
 * Given routs times with no prefetch (mflop stored at imf), try various
 * prefetch strategies, and return the best
 */
{
   ATL_mmnode_t *mp;
                    /* Lx(A)  ; Lx(B),   Lx(A&B)       2(A),x(B),  x(A),2(B) */
   const int ipfs[11]={2|512, 4|1024, 2|4|512|1024, 2|4|1024, 2|512|4,
                  /*L2(A); L2(B);L2(A&B);1(A),2(B);2(A)1(B),    L1(A&B) */
                        2,     4, 2|4,   2|64|4,   2|4|128, 2|64|4|128};
   int k, n, vrb=(vrb0 > 1);

/*
 * First, try A/B prefetch, since they are N^3
 */
   for (n=0,mp=mmb; mp; n++,mp = mp->next);
   if (vrb0)
      printf("PREFETCH TUNING FOR ALL %d KERNELS\n", n);
   for (k=0,mp=mmb; mp; k++,mp = mp->next)
   {
      const int ipf=mp->pref;
      int npf, i;
      double mf0=mp->mflop[0];
      mp->mflop[1] = mf0;
      npf = 11;     /* try all cases */
      if (vrb)
         printf("TUNING A/B PREFETCH FOR KERN %d of %d\n", k+1, n);
/*
 *    First, try only next-block prefetch
 */
      for (i=0; i < npf; i++)
         TryPref(pre, verb, vrb, mp, 1, 0, ipf|ipfs[i], 64);
/*
 *    Now let's try K-loop prefetch, both in addition to prior setting and
 *    w/o other prefetch
 */
      for (i=0; i < npf; i++)
      {
         const int kpf = mp->pref;
         const double kmf = mp->mflop[0];
         int pf;

         mp->mflop[0] = mf0;
         mp->pref = ipf;
/*
 *       Try fetching next B working set to L1, L2, & L3
 */
         TryPref(pre, verb, vrb, mp, 1, 0, ipf|16|4096, 64);
         TryPref(pre, verb, vrb, mp, 1, 0, ipf|16, 64);
         TryPref(pre, verb, vrb, mp, 1, 0, ipf|16|16384, 64);
         pf = mp->pref;
/*
 *       Using B setting, try A working set prefetch.  Not as likely to help,
 *       since we prefetch same set MB/mu times for 1 use.
 */
         TryPref(pre, verb, vrb, mp, 1, 0, pf|16|4096, 64);
         TryPref(pre, verb, vrb, mp, 1, 0, pf|16, 64);
         TryPref(pre, verb, vrb, mp, 1, 0, pf|16|16384, 64);
/*
 *       Now try combining best-found with above
 */
         pf = mp->pref;
         if (pf != kpf && kpf)
            TryPref(pre, verb, vrb, mp, 1, 0, pf|kpf, 64);
         if (mp->mflop[0] < kmf*1.005) /* if new settings (pen expen k-pref) */
         {                              /* slower M/N fetch alone, revert */
            mp->pref = kpf;
            mp->mflop[0] = kmf;
         }
         if (vrb)
            printf("   A/B PREFETCH FOR %d SPEEDUP=%.3f\n",
                   k+1,mp->mflop[0]/mf0);
      }
      {
         const int ipf=mp->pref;
         double mf0=mp->mflop[0];
         if (vrb)
            printf("   TUNING C PREFETCH FOR %d\n", k+1);
         TryPref(pre, verb, vrb, mp, 1, 0, ipf|1, 64);
         TryPref(pre, verb, vrb, mp, 1, 0, ipf|1|32, 64);
         TryPref(pre, verb, vrb, mp, 1, 0, ipf|1|256, 64);
         if (vrb)
            printf("   C PREFETCH FOR %d SPEEDUP=%.3f\n", k+1,
                   mp->mflop[0]/mf0);
      }
      if (vrb0)
      {
         printf("   PREF:");
         PrintPref(stdout, mp->pref);
         printf("; SPD=%.2f\n", mp->mflop[0] / mp->mflop[1]);
      }
      mp->mflop[1] = 0.0;
   }
   if (vrb0)
      printf("DONE  PREFETCH TUNING FOR ALL %d KERNELS\n\n", n);
}

void DoPref(char pre, int flg, int verb, ATL_mmnode_t *mmb)
/*
 * Given routs times with no prefetch (mflop stored at imf), try various
 * prefetch strategies, and return the best
 */
{
   ATL_mmnode_t *mp;
   char *ctxt[5]={"3BLKS IN L1", "1BLKS IN L1", "0BLKS IN L1",
                  "4BLKS IN L1", "5BLKS IN L1"};
                    /* Lx(A)  ; Lx(B),   Lx(A&B)       2(A),x(B),  x(A),2(B) */
   const int ipfs[11]={2|512, 4|1024, 2|4|512|1024, 2|4|1024, 2|512|4,
                  /*L2(A); L2(B);L2(A&B);1(A),2(B);2(A)1(B),    L1(A&B) */
                        2,     4, 2|4,   2|64|4,   2|4|128, 2|64|4|128};
   int k;

   printf("START PREFETCH TUNING FOR ALL CACHE CONTEXTS\n");
/*
 * First, try A/B prefetch, since they are N^3
 */
   for (k=0,mp=mmb; mp; k++,mp = mp->next)
   {
      const int ipf=mp->pref;
      int npf, i;
      double mf0=mp->mflop[0];
      mp->mflop[1] = mf0;
      if (k <= 2) /* for cases where L1 is full */
         npf = 8; /* only try prefetches to L2 (really last level cache) */
      else if (k == 3) /* room for one block */
         npf = 9;      /* try fetching one block to L1 */
      else             /* room for 5 blocks */
         npf = 11;     /* try all cases */
      printf("   TUNING A/B PREFETCH FOR %s\n", ctxt[k]);
/*
 *    First, try only next-block prefetch
 */
      for (i=0; i < npf; i++)
         TryPref(pre, verb, 1, mp, 1, 0, ipf|ipfs[i], 64);
/*
 *    Now let's try K-loop prefetch, both in addition to prior setting and
 *    w/o other prefetch
 */
      for (i=0; i < npf; i++)
      {
         const int kpf = mp->pref;
         const double kmf = mp->mflop[0];
         int pf;

         mp->mflop[0] = mf0;
         mp->pref = ipf;
/*
 *       Try fetching next B working set to L1, L2, & L3
 */
         TryPref(pre, verb, 1, mp, 1, 0, ipf|16|4096, 64);
         TryPref(pre, verb, 1, mp, 1, 0, ipf|16, 64);
         TryPref(pre, verb, 1, mp, 1, 0, ipf|16|16384, 64);
         pf = mp->pref;
/*
 *       Using B setting, try A working set prefetch.  Not as likely to help,
 *       since we prefetch same set MB/mu times for 1 use.
 */
         TryPref(pre, verb, 1, mp, 1, 0, pf|16|4096, 64);
         TryPref(pre, verb, 1, mp, 1, 0, pf|16, 64);
         TryPref(pre, verb, 1, mp, 1, 0, pf|16|16384, 64);
/*
 *       Now try combining best-found with above
 */
         pf = mp->pref;
         if (pf != kpf && kpf)
            TryPref(pre, verb, 1, mp, 1, 0, pf|kpf, 64);
         if (mp->mflop[0] < kmf*1.005) /* if new settings (pen expen k-pref) */
         {                              /* slower M/N fetch alone, revert */
            mp->pref = kpf;
            mp->mflop[0] = kmf;
         }
      }
      printf("   A/B PREFETCH FOR %s SPEEDUP=%.3f\n", ctxt[k],mp->mflop[0]/mf0);
   }
   printf("\n");
/*
 * First, try prefC for each kernel
 */
   for (k=0,mp=mmb; mp; k++,mp = mp->next)
   {
      const int ipf=mp->pref;
      double mf0=mp->mflop[0];
      printf("   TUNING C PREFETCH FOR %s\n", ctxt[k]);
      TryPref(pre, verb, 1, mp, 1, 0, ipf|1, 64);
      TryPref(pre, verb, 1, mp, 1, 0, ipf|1|32, 64);
      TryPref(pre, verb, 1, mp, 1, 0, ipf|1|256, 64);
      printf("   C PREFETCH FOR %s SPEEDUP=%.3f\n", ctxt[k], mp->mflop[0]/mf0);
      printf("   PREF:");
      PrintPref(stdout, mp->pref);
      printf("; SPD=%.2f\n", mp->mflop[0] / mp->mflop[1]);
      mp->mflop[1] = 0.0;
   }
   printf("DONE  PREFETCH TUNING FOR ALL CACHE CONTEXTS\n\n");
}

void FixContextMflop(ATL_mmnode_t *mm3)
/*
 * Put all context mflop back to ->mflop[0]
 */
{
   ATL_mmnode_t *mp;

   mm3->mflop[0] = mm3->mflop[1]; mm3->mflop[1] = 0.0;
   mp = mm3->next;
   mp->mflop[0] = mp->mflop[2]; mp->mflop[2] = 0.0;
   mp = mp->next;
   mp->mflop[0] = mp->mflop[3]; mp->mflop[3] = 0.0;
   mp = mp->next;
   if (mp)
   {
      mp->mflop[0] = mp->mflop[4]; mp->mflop[4] = 0.0;
      mp = mp->next;
      if (mp)
      {
         mp->mflop[0] = mp->mflop[5]; mp->mflop[5] = 0.0;
      }
   }
}
void DoBlock(char pre, int flg, int verb)
{
   int L1ELTS;
   ATL_mmnode_t *mmb, *mm3, *mp;
   int b;
   char upr;
   char *ressum, *mvsum, *kvsum;

   if (flg&FKO_FLAG)
   {
      ressum = "gAMMRES_fko.sum";
      mvsum = "gmvAMMUR_fko.sum";
      kvsum = "gkvAMMUR_fko.sum";
   }
   else
   {
      ressum = "gAMMRES.sum";
      mvsum = "gmvAMMUR.sum";
      kvsum = "gkvAMMUR.sum";
   }
   mmb = TimeMMFileWithPath(pre, "res", ressum, 0, verb|1, 0, 1, 0, -1);
   if (mmb)
   {
      KillAllMMNodes(mmb);
      return;
   }
   if (pre == 'c')
      upr = 's';
   else
      upr = (pre == 'z') ? 'd' : pre;
/*
 * Compute # of L1 cache elements, and join M- & K-vec kernels found so far
 */
   L1ELTS = GetL1CacheElts(upr);
   mmb = ReadMMFileWithPath(upr, "res", mvsum);
   assert(mmb);
   mp = ReadMMFileWithPath(upr, "res", kvsum);
   ATL_LastMMNode(mmb)->next = mp;
   MMFillInGenStrings(pre, mmb);
/*
 * mm3: A,B,C in L1, mm3->next: B fits in L1,
 * mm3->next->next: mu*KB panel of A fits in L1 (L2 blocked)
 */
   mm3 = FindBestCacheBudgetCases(verb, pre, L1ELTS, mmb);
   KillAllMMNodes(mmb);
   FixContextMflop(mm3);
   mmb = mm3;
   DoKUs(pre, flg, verb, mmb);   /* find best ku */
   DoPref(pre, flg, verb, mmb);  /* find best prefetch patterns */
   WriteRefreshedMMFileWithPath(pre, "res", ressum, mmb);
   KillAllMMNodes(mmb);
}

void DoSyrkAmmUM(char pre, int flag, int verb, int nreg)
/*
 * Find best generated ammm kernel to use for SYRK with restriction that MU
 * must be a multiple of NU (or vice versa if NU > MU).
 */
{
   ATL_mmnode_t *mb, *mD, *pT;
   ATL_UINT KVEC, KB, NB, VL, i, flg, muB=1, nuB=1, flgB;
   double mfB = 0.0;

   mb = TimeMMFileWithPath(pre, "res", "gSYRKUM.sum", 0, verb|1, 0, 1, 0, -1);
   if (mb)
   {
      KillAllMMNodes(mb);
      return;
   }
/*
 * Should read in 5 kernels, we want only highest performing
 */
   mb = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
   mD = FindMaxMflopMMQ(mb, 0);
   mb = RemoveMMNodeFromQ(mb, mD);
   KillAllMMNodes(mb);
/*
 * If mu already a multiple of nu (or vice versa), no need to search!
 */
   if (!(mD->mu % mD->nu) || !(mD->nu % mD->mu))
   {
      WriteMMFileWithPath(pre, "res", "gSYRKUM.sum", mD);
      KillMMNode(mD);
      return;
   }
   VL = mD->vlen;
   NB = (mD->mbB + mD->nbB)>>1;
   KB = (mD->kbB);
   flg = mD->flag;
   KVEC = ((flg>>MMF_KVEC)&1)*VL;
   flg &= ~(1<<MMF_KUISKB);
   flg |= (1<<MMF_KRUNTIME);
   pT = MMGetNodeGEN(pre, 0, 0, VL, 1, 1, VL, 0, 0, 0,
                     DupString("ATL_tmp.c"));
   pT->flag = flg;
   printf("FULL SEARCH FOR MU%%NU == 0 (or vice versa) AMM, NREG=%u, "
          "B=(%u,%u), VLEN=%u:\n", nreg, NB, KB, VL);
   for (i=1; i <= nreg; i++)
   {
      ATL_UINT j;
      if (i > NB)
         break;
      for (j=1; j <= nreg; j++)
      {
         double mf;
         ATL_UINT mu, nu, mb, nb;
         if (j > NB || i*j+1 > nreg)
            break;
         if (i*j+i+j > nreg)
            pT->flag |= 1<<MMF_BREG1;
         else
            pT->flag &= ~(1<<MMF_BREG1);
         pT->mu = mu = (KVEC) ? i : i*VL;
         pT->nu = nu = j;
         pT->ku = (KVEC) ? KVEC : 1;
         mb = (NB/mu)*mu;
         nb = (NB/nu)*nu;
         if (mu%nu && nu%mu)
            continue;
         if (!mb || !nb)
            continue;
         if (pT->genstr)
            free(pT->genstr);
         pT->genstr = MMGetGenString(pre, pT);
         mf = TimeMMKernel(verb, 0, pT, pre, mb, nb, KB, 1, 0, -1);
         printf("   MU=%2u/%2u, NU=%2u, MFLOP=%.2f\n", i, mu, j, mf);
         if (mf > mfB)
         {
            muB = mu;
            nuB = nu;
            mfB = mf;
            flgB = pT->flag;
         }
      }
   }
   assert(mfB > 0.0);
   KillMMNode(pT);
   pT = MMGetNodeGEN(pre, flgB, mD->kbB, muB, nuB, mD->ku, VL, KVEC,
                     mD->pref, mD->pfLS, NULL);
   pT->mbB = (NB/muB)*muB;
   pT->nbB = (NB/nuB)*nuB;
   pT->mbB = mD->kbB;
   printf("BEST MU%%NU==0 CASE is B=(%u,%u,%u), U=(%u,%u), MFLOP=%.2f\n",
          pT->mbB, pT->nbB, pT->kbB, muB, nuB, mfB);
   KillMMNode(mD);
   WriteRefreshedMMFileWithPath(pre, "res", "gSYRKUM.sum", pT);
   KillMMNode(pT);
}

int getSzC(ATL_mmnode_t *mp, int nb)
{
   const int nu = mp->nu, nnu = (nb+nu-1)/nu;
   const int vlen=mp->vlen, nvec=(nu*nu+vlen-1)/vlen, blksz=nvec*vlen;
   return(nnu*nnu*blksz);
}
int getSzA(ATL_mmnode_t *mp, int nb, int kb)
{
   const int nu = mp->nu, nnu = (nb+nu-1)/nu, ku=mp->ku, nku=(kb+ku-1)/ku;
   return(nnu*nku*nu*ku);
}

void TryAllSyrkPref(int flag, int vrb, char pre, ATL_mmnode_t *mp,
                    int nb, int kb)
/*
 * This kernel is designed to do an inner-product syrk.  This means we actually
 * reuse workspace for both A & C all the time, so it doesn't make sense to
 * prefetch anything from outside the block.
 */
{
   double mf0, mfB;
   int i, ipf, ils;
   const int thsA[3] = {8|8192, 8, 8|2048};
   const int thsC[3] = {1|256, 1, 1|32};
   printf("PREFETCH TUNING NB=%d, KB=%d, U=%d, KVEC=%d, SPLAT=%d:\n", nb, kb,
          mp->nu, FLAG_IS_SET(mp->flag, MMF_KVEC)?mp->vlen:0,
          FLAG_IS_SET(mp->flag, MMF_NOBCAST));
   mp->pref = 0; mp->pfLS=64;
   mf0 = mfB = TimeMMKernel(vrb, 0, mp, pre, nb, nb, kb, 1, 0, -1);
   mp->mflop[0] = mf0;
   printf("   No prefetch, MFLOP=%.0f\n", mf0);
   printf("   TRY PREFETCH A:\n");
/*
 * Try K-loop A/At pref for LS=64 & 128, since it's N^3 cost
 */
   for (i=0; i < 3; i++)
   {
      TryPref(pre, vrb, 1, mp, 1, 0, thsA[i], 64);
      TryPref(pre, vrb, 1, mp, 1, 0, thsA[i], 128);
   }
/*
 * Now try C prefetch (N^2 cost)
 */
  printf("   TRY PREFETCH C:\n");
  ipf = mp->pref;
  ils = mp->pfLS;
  for (i=0; i < 3; i++)
     TryPref(pre, vrb, 1, mp, 1, 0, ipf|thsC[i], ils);
  printf("DONE PREFETCH=(%d,%d), FINAL MFLOP=%.0f, SPEEDUP=%.3f\n",
         mp->pref, mp->pfLS, mp->mflop[0], mp->mflop[0] / mf0);
}
double TimeSyrkKBs(int flg, int vrb, char pre, ATL_mmnode_t *mp, int nb)
{
   const unsigned int ku=mp->ku, kb=(nb+ku-1)/ku;
   unsigned int KB, sz;
   static int l1elts=0;
   double mf0, mf=0.0;

   if (l1elts == 0)
      l1elts = GetL1CacheElts((pre=='c' || pre=='s') ? 's':'d');

   mf0 = TimeMMKernel(vrb, 0, mp, pre, nb, nb, kb, 1, 0, -1);
/*
 * For problems that fit entirely in L1 (including A you're copying from),
 * find a large KB that fills L1.  This will reduce the N^2 costs, which
 * are important for these small problems.  We can do this, since this
 * kernel is intended for inner product (large K, small N).
 */
   KB = kb;
   do
   {
      KB += ku;
      sz = (getSzA(mp, nb, KB)<<1)+getSzC(mp, nb);
   }
   while (l1elts > sz);
   KB -= ku;
   if (KB > kb)
      mf = TimeMMKernel(vrb, 0, mp, pre, nb, nb, KB, 1, 0, -1);

   if (mf0 > mf)
   {
      mp->mflop[0] = mf0;
      mp->kbB = kb;
   }
   else
   {
      mp->mflop[0] = mf;
      mp->kbB = KB;
      mf0 = mf;
   }
   return(mf0);
}

int SyrkXoverUp(int flg, int vrb, char pre, ATL_mmnode_t *ms,/* syrk kerns */
                ATL_mmnode_t *mg)  /* gemm largest kb */
{
   const double mfG = mg->mflop[0];
   double mfS=ms->mflop[0], mfR; /* MFLOP for SYRK & recursion */
   const unsigned int nuS = ms->nu, kuS=ms->ku, NB=ms->nbB;
   int nb, incb;

   incb = Mylcm(nuS,kuS);
   nb = NB;
   printf("ESTIMATING SYRK RECURSIVE STOPPING POINT BY INCREASING N:\n");
   while (1)
   {
      int nbn = nb+incb, nbS, nbG;
      double mfR, mfn;            /* mflop for recursion & syrk(N/2) */
      if (nbn > 512)
         break;
      mfn = TimeSyrkKBs(flg, vrb, pre, ms, nbn);
      mfR = 0.5*(mfS + mg->mflop[0]);
      printf("   SYRK-%d = %.0f, GEMM/SYRK-%d = %.0f (S:%.0f,G:%.0f)\n",
             nbn, mfn, nb, mfR, mfS, mfG);
      if (mfR > mfn)
         break;
      mfS = mfn;
      nb = nbn;
   }
   printf("STOPPING POINT N = %d\n", nb);
   return(nb);
}
int SyrkXover(int flg, int vrb, char pre, ATL_mmnode_t *ms,/* syrk kerns */
              ATL_mmnode_t *mG)  /* gemm kernels */
{
   ATL_mmnode_t *mg=mG, *mp;
   const double mf0 = ms->mflop[0];
   double mfS=mf0, mfR; /* MFLOP for SYRK & recursion */
   const unsigned int nuS = ms->nu, kuS=ms->ku, NB=ms->nbB;
   int nb;
/*
 * We are using the SYRK kernel when we stop recurring on SYRK.  In each
 * recursion, we divide C into hi&low triangles (SYRK) and a square GEMM.
 * During recursion, each SYRK is again divided to continue to use GEMM.
 * When we stop the recursion, we call the SYRK kernel we are tuning here
 * for the high & low triangles, and then unwind all the GEMM calls.
 * We note that SYRK has roughly half the flops as GEMM, but during the
 * recursion we make two SYRK calls, so at the stopping point, we can
 * simply average the MFLOP rate of SYRK & GEMM to estimate the MFLOP rate
 * the non-recursive call would get.
 *
 * We want to estimate where we should stop our recursion based on kernel
 * timings.  At each candidate N, we can make one SYRK call, or recur.
 * If we assume this is the stopping point, then its MFLOP rate is
 * that of a single syrk call.  If we assume the N/2 (next candidate N)
 * is the stopping point, its mflop rate is (syrkMF(N/2)+gemmMF(N/2)).
 * So, we have reached the estimated crossover point when
 * syrkMF(N) >= syrkMF(N/2)+gemmMF(N/2).
 *
 * This code attempts to find this point.  Note that this is a very rough
 * estimation, since the kernel timer doesn't include cpy of A, which is
 * an important part of inner product.  SYRK has to copy only A, whereas
 * GEMM needs both A & B, so we should give SYRK the benefit of the doubt,
 * and err on the side of stopping the recursion.
 */

   nb = NB;
   mg = mG;
   printf("ESTIMATING SYRK RECURSIVE STOPPING POINT:\n");
   while (1)
   {
      int nb2 = nb>>1, nbS, nbG;
      double mfR, mfs;            /* mflop for recursion & syrk(N/2) */
      if (nb2 < nuS)
         break;
      for (; mg; mg=mg->next)     /* find gemm used by nb/2 */
         if (mg->nbB <= nb2)
            break;
      if (!mg)
         break;
      nbS = (nb2/nuS)*nuS;       /* basic syrk perf estimate */
      mfs = TimeSyrkKBs(flg, vrb, pre, ms, nbS);
      mfR = 0.5*(mfs + mg->mflop[0]);
      printf("   N=%3d (S:%d,G:%d), SYRK=%.0f, RECUR=%.0f (S:%.0f,G:%.0f)\n",
             nb, nbS, mg->nbB, mfS, mfR, mfs, mg->mflop[0]);
      if (mfR < mfS)
         break;
      mfS = mfs;
      nb = nb2;
   }
   ms->mbB = ms->nbB = ((nb+nuS-1)/nuS)*nuS;
   ms->kbB = ((nb+kuS-1)/kuS)*kuS;
   if (nb == NB)
   {
      ms->mflop[0] = mf0;
      nb = SyrkXoverUp(flg, vrb, pre, ms, mG);
   }
   else
      printf("STOPPING POINT N = %d\n", nb);
   return(nb);
}

void DoAllSyrkNB(char pre, int vrb, ATL_mmnode_t *mSQ, ATL_mmnode_t *sy, int flg)
{
   ATL_mmnode_t *mSY, *mp;
   const unsigned int nu=sy->nu, ku=sy->ku, maxD;
   if (flg&FKO_FLAG)
      mp = TimeMMFileWithPath(pre, "res", "SYRKFNL_fko.sum", 0, 1, 0, 1, 0, -1);
   else
      mp = TimeMMFileWithPath(pre, "res", "SYRKFNL.sum", 0, 1, 0, 1, 0, -1);
   if (mp)
   {
      KillAllMMNodes(mp);
      return;
   }
/*
 * Make a syrk node for all nodes in mSQ
 */
   printf("TIMING SYRK KERN FOR ALL NBs:\n");
   for (mSY=NULL,mp=mSQ; mp; mp = mp->next)
   {
      ATL_mmnode_t *syp;
      double mf, mul=1.0, mfG;

      syp = CloneMMNode(sy);
      syp->mbB = syp->nbB = ((mp->nbB+nu-1) / nu)*nu;
      syp->kbB = ((syp->kbB+ku-1)/ku)*ku;
      mf = TimeMMKernel(vrb, 0, syp, pre, syp->mbB, syp->nbB, syp->kbB, 1,0,-1);
      syp->mflop[0] = mf;
      if (syp->nbB != mp->nbB)
      {
         mul = mp->nbB;
         mul /= syp->nbB;
         mul *= mul;
      }
      syp->mflop[0] *= mul;
      printf("   NB=%4u,%4u KB=%4u,%4u: GEMM=%.0f,(%.0f), SYRK=%.0f,(%.0f)\n",
             mp->nbB, syp->nbB, mp->kbB, syp->kbB,
             mp->mflop[0]/2.0,mp->mflop[0],syp->mflop[0],mf);
      syp->next = mSY;
      mSY = syp;
   }
   printf("DONE TIMING SYRK KERN FOR ALL NBs.\n");
   if (mSY->next && mSY->nbB > mSY->next->nbB)
      mSY = ReverseMMQ(mSY);
   if(flg&FKO_FLAG)
      WriteMMFileWithPath(pre, "res", "SYRKFNLi_fko.sum", mSY);
   else
      WriteMMFileWithPath(pre, "res", "SYRKFNL.sum", mSY);
   KillAllMMNodes(mSY);
}

void DoSyrkTim(char pre, int verb, int flg)
{
   ATL_mmnode_t *tb, *tp, *sb, *syp, *syb=NULL;
   int mu, nu;
   if (flg&FKO_FLAG)
      tb = ReadMMFileWithPath(pre, "res", "syrkKi_fko.sum");
   else
      tb = ReadMMFileWithPath(pre, "res", "syrkK.sum");
   if (tb)  /* already ran! */
   {
      for (tp=tb; tp; tp = tp->next)
      {
         if (tp->mflop[0] <= 0.0)
         {
            tp->mflop[0] = TimeMMKernel(verb, 0, tp, pre, tp->mbB, tp->nbB, 0,
                                        1, 0, -1);
            if (verb)
               printf("   syrkK %dx%d: %.2f\n", tp->mbB, tp->nbB, tp->mflop[0]);
         }
      }
      if (flg&FKO_FLAG)
         WriteRefreshedMMFileWithPath(pre, "res", "syrkK_fko.sum", tb);
      else
         WriteRefreshedMMFileWithPath(pre, "res", "syrkK.sum", tb);
      KillAllMMNodes(tb);
      return;
   }
   if (flg&FKO_FLAG)
   {
      tb = ReadMMFileWithPath(pre, "res", "sqAMMRES_fko.sum");
      syp = ReadMMFileWithPath(pre, "res", "gAMSYRK_fko.sum");
   }
   else
   {
      tb = ReadMMFileWithPath(pre, "res", "sqAMMRES.sum");
      syp = ReadMMFileWithPath(pre, "res", "gAMSYRK.sum");
   }
   assert(syp);
   mu=syp->mu;
   nu=syp->nu;
   printf("\nFINDING PERFORMANCE OF SYRK KERNELS FOR MB=NB AMM:\n");
   for (tp=tb; tp; tp = tp->next)
   {
      ATL_mmnode_t *new;
      double ratio, mf;
      int mb=tp->mbB, nb=tp->nbB;

      new = CloneMMNode(syp);
      mb = ((mb+mu-1)/mu)*mu;
      nb = ((nb+nu-1)/nu)*nu;
      new->mbB = mb;
      new->nbB = nb;
      new->kbB = tp->kbB;
      if (mb == tp->mbB && nb == tp->nbB)
         ratio = 1.0;
      else
         ratio = (((double)tp->mbB)*tp->nbB) / (((double)mb)*nb);
      mf = TimeMMKernel(verb, 0, new, pre, mb, nb, tp->kbB, 1, 0, -1);
      mf *= ratio;  /* don't count useless flops */
      new->mflop[0] = mf;
      printf("   syrkK %dx%d: %.2f (%.2f)\n", tp->mbB, tp->nbB, mf, ratio);
      new->next = syb;
      syb = new;
   }
   KillAllMMNodes(tb);
   KillAllMMNodes(syp);
   syb = ReverseMMQ(syb);
   if (flg&FKO_FLAG)
      WriteRefreshedMMFileWithPath(pre, "res", "syrkK_fko.sum", syb);
   else
      WriteRefreshedMMFileWithPath(pre, "res", "syrkK.sum", syb);
   printf("DONE SYRK TIMING.\n");
   KillAllMMNodes(syb);
}

void DoSyrk(int flg, int vrb, char pre, int nreg, int nb, int VL)
{
   ATL_mmnode_t *mb, *mSQ, *mp, *mSY;
   int L1ELTS, i, maxNB, nu, ku, kb;
   char upr = (pre == 'c' || pre == 's') ? 's' : 'd';
   double mf0, mfB;
   char *ipsum, *syrksum;

   if (flg&FKO_FLAG)
   {
      ipsum = "ipmen_fko.sum";
      syrksum = "gAMSYRK_fko.sum";
   }
   else
   {
      ipsum = "ipmen.sum";
      syrksum = "gAMSYRK.sum";
   }
   mSQ = TimeMMFileWithPath(pre, "res", ipsum, 0, vrb|1, 0, 1, 0, -1);
   if (!mSQ)
      return;
   mb = TimeMMFileWithPath(pre, "res", syrksum, 0, vrb|1, 0, 1, 0, -1);
   if (mb)
   {
      DoAllSyrkNB(pre, vrb, mSQ, mb, flg);
      KillAllMMNodes(mb);
      KillAllMMNodes(mSQ);
      return;
   }
   mSQ = ReverseMMQ(mSQ);  /* sorted from large to small NB */
   mb = DoSyrkMUNU(flg, vrb, pre, nreg, mSQ->nbB, mSQ->kbB, VL);
   DoSyrkKU(vrb, pre, mb);
   mp = mb;
   nu = mp->nu;
   ku = mp->ku;
   mp->nbB = mp->mbB = nb = (mSQ->nbB / nu)*nu;
   mp->kbB = kb = (mSQ->kbB / ku)*ku;
   assert(nb && kb);
   TryAllSyrkPref(flg, vrb, pre, mp, nb, kb);
   WriteMMFileWithPath(pre, "res", syrksum, mp);
/*
 * Make a syrk node for all nodes in mSQ
 */
   DoAllSyrkNB(pre, vrb, mSQ, mb, flg);
   KillAllMMNodes(mSQ);
   KillMMNode(mb);
}

void DoPrefOnList(int flg, int vrb, char pre, int maxB)
{
   ATL_mmnode_t *mb, *mp, *pb;
   int b, binc;
   FILE *fp;

   if (flg&FKO_FLAG)
      mb = ReadMMFile("time.sum"); /* may use different name later */
   else
      mb = ReadMMFile("time.sum");
   if (!mb)
   {
      fprintf(stderr, "No files in time.sum, exiting!\n");
      remove("NB-pref.txt");
      return;
   }
   binc = Mylcm(mb->mu, mb->nu);
   for (mp=mb->next; mp; mp = mp->next)
   {
      binc = Mylcm(binc, mp->mu);
      binc = Mylcm(binc, mp->nu);
   }
   fp = fopen("NB-pref.txt", "w");
   assert(fp);
   fprintf(fp, "      ");
   for (mp=mb; mp->next; mp = mp->next)
      fprintf(fp, "-----------");
   fprintf(fp, "---------  ");
   for (mp=mb; mp->next; mp = mp->next)
      fprintf(fp, "-----------");
   fprintf(fp, "---------\n");
   fprintf(fp, "  NB");
   for (mp=mb; mp; mp = mp->next)
      fprintf(fp, "   %c%c:%2dx%2d", FLAG_IS_SET(mp->flag, MMF_KVEC) ? 'K':'M',
              FLAG_IS_SET(mp->flag,MMF_NOBCAST)?'s':'b', mp->mu, mp->nu);
   for (mp=mb; mp; mp = mp->next)
      fprintf(fp, "  PREF %c%2dx%2d",
              FLAG_IS_SET(mp->flag, MMF_KVEC) ? 'K':'M', mp->mu, mp->nu);
   fprintf(fp, "\n====");
   for (mp=mb; mp; mp = mp->next)
      fprintf(fp, "  =========");
   for (mp=mb; mp; mp = mp->next)
      fprintf(fp, "  ==== ======");
   fprintf(fp, "\n");
   fflush(fp);
   for (b=binc; b <= maxB; b += binc)
   {
      fprintf(fp, "%4d", b);
      for (mp=mb; mp; mp = mp->next)
      {
         double mf;
         mp->pref = 0;
         mp->nbB = mp->mbB = mp->kbB = b;
         if (mp->genstr)
            free(mp->genstr);
         mp->genstr = NULL;
         mf = TimeMMKernel(vrb, 0, mp, pre, b, b, b, 1, 0, -1);
         mp->mflop[0] = mf;
         fprintf(fp, "     %6.0f", mf);
      }
      DoPref1(pre, flg, vrb, 1, mb);
      for (mp=mb; mp; mp = mp->next)
         fprintf(fp, "  %4.4x %6.0f", mp->pref, mp->mflop[0]);
      fprintf(fp, "\n");
      fflush(fp);
   }
   fclose(fp);
   KillAllMMNodes(mb);  /* done with these! */
}

void SetFlops(int flg, int verb, char pre, int VL)
/*
 * This search finds a decent generated kernel, and uses it to find the
 * number of flops that must be forced to get
 * 1:  each problem above .25 sec?
 * 2:  10 timings withing 2% of each other?
 * which should be enough to avoid huge timing error bars
 */
{
   ATL_mmnode_t *mp;
   FILE *fp;
   char *rt;
   char fn[32];
   int b, i;
   int muB, nuB, kvecB, bB, bcB;
   #define TOL .015
   const double mbig=(1.0+TOL), msml=(1.0-TOL);
   double mfB, mf, mul, nmf, tim;
   const int bf=(flg&FKO_FLAG)?FKO_FLAG:0;

/*
 * Get square blocking factor that fills cache with one block
 */
   i = GetL1CacheElts(pre);
   for (b=4; b*b <= i; b += 4);
   b -= 4;
/*
 * If we have no vector length set, see if we know the platform, else
 * attempt to determine it empirically (can use gnuvec to vectorize
 * archs where ATLAS has no explicit support).
 */
   if (VL < 1)
      VL = GetNativeVLEN(pre);
   if (VL < 1)  /* no good guess, will have to probe for it */
   {
      int vl, vlB=1;
      int B = (b>>1)<<1;
      mp = MMGetNodeGEN(pre, 0, 0, 2, 2, 1, 1, 1, 0, 0, DupString("ATL_tmp.c"));
      mp->mbB = mp->nbB = 2;
      mp->kbB = 480;
      mp->flag &= ~( (1<<MMF_MVA)|(1<<MMF_MVB)|(1<<MMF_MVC) );
      mfB = 0.0;
      printf("ATTEMPTING TO DETECT VECTOR LENGTH:\n");
      for (vl=1; vl < 16; vl += vl)
      {
         mp->ku = mp->vlen = vl;
         if (vl == 2)
            mp->flag |= (1<<MMF_KVEC);
         if (mp->genstr)
            free(mp->genstr);
         mp->genstr = NULL;
         mf = TimeMMKernel(verb, 1, mp, pre, 2, 2, 480, 1, mul, 0);
         printf("   VL=%d, KVEC=1, B=(2,2,480) U=2: %.0f\n", vl, mf);
         if (mf > mfB)
         {
            mfB = mf;
            vlB=vl;
         }
      }
      printf("VLEN SET TO %d\n\n", vlB);
      VL = vlB;
      KillMMNode(mp);
   }

   mp = MMGetNodeGEN(pre, 0|bf, 0, 1, 1, VL, VL, 1, 0, 0, DupString("ATL_tmp.c"));
   mp->ku = 1;
   mp->flag &= ~( (1<<MMF_MVA)|(1<<MMF_MVB)|(1<<MMF_MVC) ); /* no cache flush */
   printf("SCOPING FOR A DECENT SQUARE REGISTER BLOCKING, VLEN=%d\n", VL);
   printf("     B   U  Mb MFLOP  Ms MFLOP  K  MFLOP\n");
   printf("   ===  ==  ========  ========  ========\n");
   bB = Mylcm(VL, 6);
   bB = ((b+bB-1)/bB)*bB;
   i = MMTimeAllVecGen(verb, 1, mp, pre, bB, bB, bB, 0, 0, 0);
   printf("%6d %3d %9.0f %9.0f %9.0f\n", i, 1,
          mp->mflop[0],mp->mflop[1],mp->mflop[2]);
   fflush(stdout);
   mfB = mp->mflop[i];
   muB = nuB = 1;
   kvecB = (i == 2) ? 1 : 0;
   bcB = (i != 1);
   for (i=2; i*i <= 64; i++)
   {
      int B = (b > i) ? (b/i)*i : i, iB;
      mp->mu = mp->nu = i;
      iB = MMTimeAllVecGen(verb, 1, mp, pre, B, B, B, 0, 0, 0);
      printf("%6d %3d %9.0f %9.0f %9.0f\n", B, i,
             mp->mflop[0],mp->mflop[1],mp->mflop[2]);
      if (mp->mflop[iB] > mfB)
      {
         mfB = mp->mflop[iB];
         muB = nuB = i;
         bB = B;
         bcB = (iB != 1);
         kvecB = (iB == 2) ? 1 : 0;
      }
      fflush(stdout);
   }
   if (mp->genstr)
      free(mp->genstr);
   mp->genstr = NULL;
/*
 * Now, pump up mflop until the run takes 3 seconds or we can get 8 timings
 * in a row that are within 1.5% of each other
 */
   mp->mflop[0] = mfB;
   mp->mflop[1] = mp->mflop[2] = 0.0;
   mp->nu = nuB;
   mp->mbB = mp->nbB = mp->kbB = bB;
   if (kvecB)
   {
      mp->ku = VL;
      mp->mu = muB;
      mp->flag |= (1<<MMF_KVEC);
      mp->flag &= ~(1<<MMF_NOBCAST);
   }
   else
   {
      mp->ku = 1;
      mp->mu = muB * VL;
      mp->flag &= ~(1<<MMF_KVEC);
      if (bcB)
         mp->flag &= ~(1<<MMF_NOBCAST);
      else
         mp->flag |= (1<<MMF_NOBCAST);
   }
/*
 * aim for a run that should take .08 sec if above timing was good
 *    (seconds / mflop) * nmfl = .08 -> nmflop = .08*mfB;
 */
   nmf = .08*mfB;
   if (nmf < 1.0)
      nmf = 1.0;
   while(1)
   {
      double mf0;
      int inmf = (int) nmf;
      tim = nmf/mfB;
      fprintf(stderr, "nmf=%d, pred time = %e seconds\n", inmf,
              1.0/(nmf*mfB));
      if (tim >= 2.0)
         break;
      mf0 = TimeMMKernel(verb, 1, mp, pre, bB, bB, bB, 1, nmf, 0);
      printf("mf=%.0f:", mf0);
      for (i=0; i < 7; i++)
      {
         double spdup;
         mf = TimeMMKernel(verb, 1, mp, pre, bB, bB, bB, 1, mul, 0);
         spdup = mf*mf0;
         printf(" %.4f%c", mf/mf0, i==6?'.':',');
         if (mf > mf0*mbig || mf < mf0*msml)
            break;
      }
      printf("\n");
      fflush(stdout);
      if (i == 7)
         break;
      nmf *= 2.0;
   }
   printf("Final timing interval: %e seconds, nmflops=%.0f\n", tim, nmf);

   if(bf)
      sprintf(fn, "%cmflopsi_fko.frc", pre);
   else
      sprintf(fn, "%cmflops.frc", pre);
   fp = fopen(fn, "w");
   fprintf(fp, "%e", nmf);
   fclose(fp);
   KillMMNode(mp);
}
int DoTrmmLeft(char pre, int flag, int verb, int nrg)
{
   ATL_mmnode_t *mb, *mD, *pT;
   ATL_mmnode_t *mbn, *mbt;
   ATL_UINT KVEC, KB, NB, VL, i, flg, muB=1, nuB=1, flgB;
   double mfB = 0.0;
/*
 * FIXME: need to update timer to have correct timer data
 */
   mbn = TimeMMFileWithPath(pre, "res", "trmmL_LN.sum", 0, verb|1, 0, 1, 0, -1);
   mbt = TimeMMFileWithPath(pre, "res", "trmmL_LT.sum", 0, verb|1, 0, 1, 0, -1);
   if (mbn && mbt) /* generated either both or none for now */
   {
      KillAllMMNodes(mbn);
      KillAllMMNodes(mbt);
      return(0);
   }
/*
 * Should read in 5 kernels, we want only highest performing
 */
   mb = ReadMMFileWithPath(pre, "res", "ipmek.sum");
   if (!mb) return(0); /* no ipmek yet! */
   mD = FindMaxMflopMMQ(mb, 0);
   mb = RemoveMMNodeFromQ(mb, mD);
   KillAllMMNodes(mb);
/*
 * If mu already a multiple of nu (or vice versa), no need to search!
 */
   if (!(mD->mu % mD->ku) || !(mD->ku % mD->mu))
   {
      /*mD->blask = ATL_KTRMM; */
      WriteMMFileWithPath(pre, "res", "trmmL_LN.sum", mD);
      WriteMMFileWithPath(pre, "res", "trmmL_LT.sum", mD);
      KillMMNode(mD);
      return(0);
   }
/*
 * FIXME: need to write search for trmm... returning error now,
 * will implement it later
 */
   return(1);
}

int DoTrmmRight(char pre, int flag, int verb, int nrg)
{
   ATL_mmnode_t *mb, *mD, *pT;
   ATL_mmnode_t *mbn, *mbt;
   ATL_UINT KVEC, KB, NB, VL, i, flg, muB=1, nuB=1, flgB;
   double mfB = 0.0;
/*
 * FIXME: need to update timer to time
 */
   mbn = TimeMMFileWithPath(pre, "res", "trmmR_LN.sum", 0, verb|1, 0, 1, 0, -1);
   mbt = TimeMMFileWithPath(pre, "res", "trmmR_LT.sum", 0, verb|1, 0, 1, 0, -1);
   if (mbn && mbt) /* generated either both or none for now */
   {
      KillAllMMNodes(mbn);
      KillAllMMNodes(mbt);
      return(0);
   }
/*
 * Should read in 5 kernels, we want only highest performing
 */
   mb = ReadMMFileWithPath(pre, "res", "ipnek.sum");
   if (!mb) return(0); /* no ipnek yet! */
   mD = FindMaxMflopMMQ(mb, 0);
   mb = RemoveMMNodeFromQ(mb, mD);
   KillAllMMNodes(mb);
/*
 * If mu already a multiple of nu (or vice versa), no need to search!
 */
   if (!(mD->mu % mD->ku) || !(mD->ku % mD->mu))
   {
      /*mD->blask = ATL_KTRMM;*/
      WriteMMFileWithPath(pre, "res", "trmmR_LN.sum", mD);
      WriteMMFileWithPath(pre, "res", "trmmR_LT.sum", mD);
      KillMMNode(mD);
      return(0);
   }
/*
 * FIXME: need to write search for trmm... returning error now,
 * will implement it later
 */
   return(1);
}

static int AdjustForTRMM(char pre, int flag, ATL_mmnode_t *mp, ATL_mmnode_t *gd)
/*
 * Takes an amm kernel, and adjusts it for use as a TRMM.  We will do following:
 * Adjust ku: KVEC must use VLEN, otherwise number between 1-4 where U%ku==0
 * or ku%U == 0.  KRUNTIME is also forced, and it doesn't need KCLEAN.
 * If the kern isn't generated, get rid of special compiler/flags.
 * We also delete any flags or special compilers, in case its not a C code.
 * if gd non-null, we also change mp's block factors to roughly match its.
 * RETURNS: old U (mu if flag&1, else nu)
 */
{
   const double mfR = mp->mflop[0];
   const unsigned int UU=(flag&1)?mp->mu:mp->nu;
   unsigned int U=UU, B;
/*
 * Make K-loop runtime
 */
   if (mp->ID)  /* if present code is hand-tuned, make it generated */
   {
      mp->ID = 0;
      if (mp->rout)
         free(mp->rout);
      mp->rout = DupString("ATL_tmp.c");
      if (mp->comp)
      {
         free(mp->comp);
         mp->comp = NULL;
      }
      if (mp->cflags)
      {
         free(mp->cflags);
         mp->cflags = NULL;
      }
   }
   mp->flag &= ~((1<<MMF_KUISKB)|(1<<MMF_KCLN));
   mp->flag |= 1<<MMF_KRUNTIME;
/*
 * For kvec kernel, KU is fixed, so we must adjust U
 */
   if (FLAG_IS_SET(mp->flag, MMF_KVEC))
   {
      const unsigned int ku = mp->vlen;
      unsigned int u, nreg=mp->mu*(mp->nu+1);

      mp->ku = ku;
      nreg += (mp->flag & (1<<MMF_BREG1)) ? 1 : mp->nu;
/*
 *    If U not multiple of ku, or vice versa, must adjust U
 */
      if ((U/ku)*ku != U && (ku/U)*U != ku)
      {
         if (U > ku)
            U = (U/ku)*ku;
         else
            U = ku;
         if (flag&1) /* U == mu */
            mp->mu = U;
         else        /* U == nu */
            mp->nu = U;
      }
   }
   else  /* MVEC kernel can adjust both U and ku, but we just reduce ku */
   {
      unsigned int ku = Mmax(mp->ku,4);
      while ((U/ku)*ku != U)
         ku--;
      mp->ku = ku;
   }
   if (flag & 1)  /* Left, MB == KB */
   {
      unsigned int sz = mp->mbB * mp->kbB;
      U = Mylcm(mp->mu, mp->ku);
      for (B=U; B*B < sz; B += U) ;
      if (B*B > sz && B > U)
         B -= U;
      mp->mbB = mp->kbB = B;
      U = mp->nu;
      B = mp->nbB;
      mp->nbB = (B >= U) ? ((B/U)*U) : U;
   }
   else           /* Right, NB == KB */
   {
      unsigned int sz = mp->nbB * mp->kbB;
      U = Mylcm(mp->nu, mp->ku);
      for (B=U; B*B < sz; B += U) ;
      if (B*B > sz && B > U)
         B -= U;
      mp->nbB = mp->kbB = B;
      U = mp->mu;
      B = mp->mbB;
      mp->mbB = (B >= U) ? ((B/U)*U) : U;
   }
   if (mp->genstr)
      free(mp->genstr);
   mp->genstr = MMGetGenString(pre, mp);

   return(UU);
}

ATL_mmnode_t *srchTrmm(char pre, int flag, int verb, int nreg,
                       int MB, int NB, ATL_mmnode_t *mmB)
/*
 * flag bitvec: 0:Left
 * Assumes mmB is best-performing case found so far, MB, NB block factors of
 * original (unadapted) TRMM case.  Searches for all U that match ku.
 * RETURNS: (possibly new) mmB case (kills original mmB in this case).
 */
{
   ATL_mmnode_t *pM, *pK;
   const char *frm="    %c %4d %4d %3d %3d %9.0f\n";
   unsigned int MU, KU=mmB->ku, bf=mmB->flag;

   pM = MMGetNodeGEN(pre, bf&(~(1<<MMF_KVEC)), 0, 0, 0, 1, mmB->vlen,
                     0, 0, 0, DupString("ATL_tmp.c"));
   pM->ku = 1;
   pK = MMGetNodeGEN(pre, bf|(1<<MMF_KVEC), 0, 0, 0, mmB->vlen, mmB->vlen,
                     mmB->vlen, 0, 0, DupString("ATL_tmp.c"));
   free(pM->genstr);
   free(pK->genstr);
   pK->ku = mmB->vlen;
   printf("Search for TRMM kernel for SD=%c, B=(%u,%u), NREG=%d, VLEN=%d\n",
          (flag&1) ? 'R':'L', MB, NB, nreg, mmB->vlen);
   printf("   VD   MB   NB  MU  NU     MFLOP\n");
   printf("   ==  ===  ===  ==  ==  =========\n");
   for (MU=1; MU < nreg; MU++)
   {
      unsigned int NU, mb, nb, kb, DOK;
      double mf;
      for (NU=1; NU < nreg; NU++)
      {
         if (MU*(NU+1)+1 > nreg)
            continue;
         if (MU*(NU+1)+NU > nreg)
         {
            pM->flag |= 1<<MMF_BREG1;
            pK->flag |= 1<<MMF_BREG1;
         }
         else
         {
            pM->flag &= ~(1<<MMF_BREG1);
            pK->flag &= ~(1<<MMF_BREG1);
         }
         pK->mu = MU;
         pM->mu = MU * pM->vlen;
         pK->nu = pM->nu = NU;
         pM->genstr = MMGetGenString(pre, pM);
         mb = (MB > MU) ? (MB/MU)*MU : MU;
         nb = (NB > NU) ? (NB/NU)*NU : NU;
         kb = (flag&1) ? mb:nb;
         /* assert(!MMKernelFailsTest(pre, mb, nb, kb, 1, pM)); */
         mf = TimeMMKernel(0, 0, pM, pre, mb, nb, kb, 1, 0, -1);
         printf(frm, 'M', mb, nb, MU, NU, mf);
         if (mf > mmB->mflop[0])
         {
            CopyMMNode(mmB, pM);
            mmB->mbB = mb;
            mmB->nbB = nb;
            mmB->kbB = kb;
            /* assert(!MMKernelFailsTest(pre, mb, nb, kb, 1, mmB)); */
            mmB->mflop[0] = mf;
         }
         free(pM->genstr);
         if (flag&1)
            DOK = (KU%MU == 0 || MU%KU == 0);
         else
            DOK = (KU%NU == 0 || NU%KU == 0);
         if (DOK)
         {
            unsigned int U;
            if (flag&1)
            {
               U = Mylcm(MU, KU);
               mb = (MB > U) ? (MB/U)*U : U;
            }
            else
            {
               U = Mylcm(NU, KU);
               nb = (NB > U) ? (NB/U)*U : U;
            }
            kb = (flag&1) ? mb:nb;
            pK->genstr = MMGetGenString(pre, pK);
            mf = TimeMMKernel(0, 0, pK, pre, mb, nb, kb, 1, 0, -1);
            printf(frm, 'K', mb, nb, MU, NU, mf);
            if (mf > mmB->mflop[0])
            {
               CopyMMNode(mmB, pK);
               mmB->mbB = mb;
               mmB->nbB = nb;
               mmB->kbB = kb;
               /* assert(!MMKernelFailsTest(pre, mb, nb, kb, 1, mmB)); */
               mmB->mflop[0] = mf;
            }
            free(pK->genstr);
         }
      }
   }
   printf("BEST TRMM CASE:\n");
   PrintMMLine(stdout, mmB);
   assert(!MMKernelFailsTest(pre, mmB->mbB, mmB->nbB, mmB->kbB, 1, mmB));
/*
 * If best case was M-vectorized, try ku=2,3,4
 */
   if (!(mmB->flag & (1<<MMF_KVEC)))
   {
      const unsigned U = (flag&1) ? mmB->mu : mmB->nu;
      const unsigned B = (flag&1) ? mmB->mbB : mmB->nbB;
      unsigned int mb, nb;
      unsigned int ku;

      pK->genstr = NULL;
      printf("\nTUNING TRMM's KU:\n");
      CopyMMNode(pK, mmB);
      free(pK->genstr);
      for (ku=2; ku < 4; ku++)
      {
         double mf;
         if (U%ku || B%ku)
            continue;
         pK->ku = ku;
         pK->genstr = MMGetGenString(pre, pK);
         mf = TimeMMKernel(0, 0, pK, pre, pK->mbB, pK->nbB, pK->kbB, 1, 0, -1);
         printf("   KU=%d, mf=%.2f\n", ku, mf);
         if (mf > mmB->mflop[0])
         {
            CopyMMNode(mmB, pK);
            mmB->mflop[0] = mf;
         }
         free(pK->genstr);
      }
   }
   pK->genstr = pM->genstr = NULL;
   KillMMNode(pK);
   KillMMNode(pM);
   printf("BEST TRMM CASE:\n");
   PrintMMLine(stdout, mmB);
   printf("\n");
   return(mmB);
}

void DoTrmmS(char pre, int flag, int verb, int nreg)
{
   ATL_mmnode_t *bp, *mmB;
   unsigned int U, ku, MB, NB, KB;
   double mf0, mfB, mf;

   bp = TimeMMFileWithPath(pre,"res",(flag&1)?"trmmKLL.sum":"trmmKRL.sum",
                           0, 0, 0, 1, 0, -1);
   if (bp)
   {
      WriteMMFileWithPath(pre,"res",(flag&1)?"trmmKLU.sum":"trmmKRU.sum", bp);
      WriteMMFileWithPath(pre,"res",(flag&1)?"trmmKLL.sum":"trmmKRL.sum", bp);
      KillAllMMNodes(bp);
      return;
   }
   bp = ReadMMFileWithPath(pre, "res", (flag&1)?"ipmek.sum":"ipnek.sum");
   if (!bp)
      return;
   printf("\nSEARCHING FOR TRMM KERNEL, SIDE=%c\n", (flag&1) ? 'L':'R');
   mmB = FindMaxMflopMMQ(bp, 0);
   KillAllMMNodes(RemoveMMNodeFromQ(bp, mmB));
   MB = mmB->mbB;
   NB = mmB->nbB;
   KB = mmB->kbB;
   U = (flag&1) ? mmB->mu : mmB->nu;
   ku = mmB->ku;
   mf0 = mmB->mflop[0];
   AdjustForTRMM(pre, flag, mmB, NULL);
/*
 * Adjust given kernel to allow it to use with TRMM, and see what perf we get
 */
   mf = TimeMMKernel(0, 0, mmB, pre, mmB->mbB, mmB->nbB, mmB->kbB, 1, 0,-1);
   printf("   AMM=%.2f, TRM=%.2f\n", mf0, mf);
   mfB = mmB->mflop[0] = mf;
/*
 * If performance drops significantly, try tuning a kernel specifically for TRMM
 */
   if (mf*1.035 < mf0)
   {
      ATL_mmnode_t *mp;
      unsigned int B, Ug;

      mmB->mflop[0] = mfB;
      bp = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
      assert(bp);
      mp = FindMaxMflopMMQ(bp, 0);
      KillAllMMNodes(RemoveMMNodeFromQ(bp, mp));
      Ug = (flag&1) ? mp->mu : mp->nu;
/*
 *    Set blocking factor to one found by [m,n]eqk
 */
      if (flag&1)  /* MB = KB */
      {
         Ug = Mylcm(mp->mu, mp->ku);
         assert(MB == KB);  /* sanity check for expected case */
         mp->kbB = (KB > Ug) ? ((KB)/Ug)*Ug : Ug;
         mp->mbB = mp->kbB;
         mp->nbB = (NB/mp->nu)*mp->nu;
      }
      else         /* NB = KB */
      {
         Ug = Mylcm(mp->nu, mp->ku);
         assert(NB == KB);  /* sanity check for expected case */
         mp->kbB = (KB > Ug) ? ((KB)/Ug)*Ug : Ug;
         mp->nbB = mp->kbB;
         mp->mbB = (MB > mp->mu) ? (MB/mp->mu)*mp->mu : mp->mu;
      }
/*
 *    If best generated kernel works w/o modification, simply use it unless
 *    its performance worse than adapted main kernel
 */
      if ((mp->flag&(1<<MMF_KRUNTIME)) &&
          ((Ug/ku)*ku == Ug || (ku/Ug)*Ug == ku))
      {
         mf = TimeMMKernel(0, 0, mp, pre, mp->mbB, mp->nbB, mp->kbB, 1, 0,-1);
         printf("   Using best generated case,\n      AMM=(%u,%u,%u) %.2f, "
                "TRM=(%u,%u,%u) %.2f\n",
                MB, NB, KB, mf0, mp->mbB, mp->nbB, mp->kbB, mf);
         if (mf > mmB->mflop[0])
         {
            KillMMNode(mmB);
            mmB = mp;
            mp->mflop[0] = mfB = mf;
         }
         else
            KillMMNode(mp);
      }
      else
      {
         ATL_mmnode_t *gp;
         unsigned int NOSRCH;

         gp = CloneMMNode(mp);
         AdjustForTRMM(pre, flag, gp, NULL);
         mf = TimeMMKernel(0, 0, gp, pre, gp->mbB, gp->nbB, gp->kbB, 1, 0,-1);
         gp->mflop[0] = mf;
         if (mmB->mflop[0] >= mf)
            KillMMNode(gp);
         else
         {
            KillMMNode(mmB);
            mmB = gp;
         }
         printf("   Using adapted generated case,\n      AMM=(%u,%u,%u) %.2f, "
                "TRM=(%u,%u,%u) %.2f\n",
                MB, NB, KB, mf0, mp->mbB, mp->nbB, mp->kbB, mf);
/*
 *       If we didn't change register blocking, just go with gAMMRES
 */
         NOSRCH = (flag&1) ? (mp->mu == gp->mu) : (mp->nu == gp->nu);
         KillMMNode(mp);
         if (!NOSRCH || mf*1.035 < mf0)
            mmB = srchTrmm(pre, flag, verb, nreg, mmB->mbB, mmB->nbB, mmB);
      }
   }
   if (!(flag&1))
      mmB->flag |= 1<<MMF_RIGHT;
   WriteMMFileWithPath(pre,"res",(flag&1)?"trmmKLL.sum":"trmmKRL.sum",mmB);
   WriteMMFileWithPath(pre,"res",(flag&1)?"trmmKLU.sum":"trmmKRU.sum",mmB);
   KillMMNode(mmB);
}

void DoTrmm(char pre, int flag, int verb, int nreg)
{
#if 1
   DoTrmmS(pre, 0, verb, nreg);
   DoTrmmS(pre, 1, verb, nreg);
#else
/*
 * FIXME: add make target to time trmm to record
 * correct perf data
 */
   assert(!DoTrmmLeft(pre,flag,verb,nreg));
   assert(!DoTrmmRight(pre,flag,verb,nreg));
#endif
}
int main(int nargs, char **args)
{
   int flg, verb, nreg, VLEN, NB, TEST;
   char pre, upr;
   GetFlags(nargs, args, &flg, &verb, &pre, &nreg, &VLEN, &NB, &TEST);
   upr = pre;
   if (pre == 'z')
      upr = 'd';
   else if (pre == 'c')
      upr = 's';
   if (TEST)
      return(CountFails(TEST, flg, verb, pre, NB, nreg, VLEN));
   if (flg&8)
   {
      DoPrefOnList(flg, verb, pre, NB);
      return(0);
   }
   if (flg&16)
   {
      SetFlops(flg, verb, upr, VLEN);
      return(0);
   }
   FindInfo(flg, verb, upr, NB, &nreg, &VLEN);
   if (flg&4)     /* If I just wanted nreg calc */
      return(0);  /* then I'm done */
   DoBlock(pre, flg, verb);
   DoSyrkAmmUM(pre, flg, verb, nreg);
   DoSyrk(flg, verb, pre, nreg, NB, VLEN);
   DoTrmm(pre, flg, verb, nreg);
   return(0);
}
