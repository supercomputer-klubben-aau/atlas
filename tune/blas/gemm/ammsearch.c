/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015, 2013, 2012 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_mmtesttime.h"

static int VLEN=0;
static int TSIZE=8;
static int IMVS=3;     /* move ptrs in timing encoded in last 3 bits: CBA */
#define KRUNMUL 1.02   /* KRUNTIME speedup increase over K-compile time */

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

ATL_mmnode_t *GetGenCases(char pre)
{
   ATL_mmnode_t *mb, *mp;
   mb = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
   if (!mb)
   {
      char ln[32];
      sprintf(ln, "make res/%cgAMMRES.sum", pre);
      assert(!system(ln));
      mb = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
   }
   assert(mb);
   MMFillInGenStrings(pre, mb);
   return(mb);
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
/*
 * Finds best blocking factors for kernel mmp trying all legal values
 * between [b0, bN]
 */
ATL_mmnode_t *BestBlocking_BFI
(
   int verb,
   char pre,
   ATL_mmnode_t *mmp,
   int b0,
   int bN,
   int minInc,  /* minimum increment to use */
   int FORCE
)
/*
 * Times all legal block factors in all dims between [b0,bN].
 * RETURNS: ptr to best performing kernel, NULL if no legal block factors
 */
{
   ATL_mmnode_t *mp;
   int mbB=0, nbB=0, kbB=0;
   int mbS=0, nbS=0, kbS=0;
   int mu = mmp->mu, nu = mmp->nu, ku = mmp->ku;
   int k0, kn, m0, mn, n0, nn, m, n, k;
   double mfB=0.0, mfS=0.0;

   if (!mmp)
      return(NULL);
   if (minInc > mu)
      mu = ((minInc+mu-1)/mu)*mu;
   if (minInc > nu)
      nu = ((minInc+nu-1)/nu)*nu;
   if (minInc > ku)
      ku = ((minInc+ku-1)/ku)*ku;
   m0 = ((b0+mu-1)/mu)*mu;
   n0 = ((b0+nu-1)/nu)*nu;
   k0 = ((b0+ku-1)/ku)*ku;
   mn = ((bN+mu-1)/mu)*mu;
   nn = ((bN+nu-1)/nu)*nu;
   kn = ((bN+ku-1)/ku)*ku;
   mp = CloneMMNode(mmp);
   if (mp->kbmax && mp->kbmax < kn)
      kn = mp->kbmax;
   if (mp->kbmin && mp->kbmin > k0)
      k0 = mp->kbmin;


   printf("SEARCH BLKING [%d - %d] for %d.%s:\n\n", b0, bN, mp->ID, mp->rout);
   printf("  MB    NB    KB        MFLOP    mbB  nbB  kbB      mflopB\n");
   printf("====  ====  ====  ===========   ==== ==== ==== ===========\n");
   for (m=m0; m <= mn; m += mu)
   {
      for (n=m0; n <= nn; n += nu)
      {
         for (k=k0; k <= kn; k += ku)
         {
            double mf;
            mf = TimeMMKernel(verb, FORCE, mp, pre, m, n, k, 1, 0, -1);
            printf("%4d %5d %5d %11.1f %4d %4d %4d %11.1f\n",
                   m, n, k, mf, mbB, nbB, kbB, mfB);
            if (mf > mfB)
            {
               mfB = mf;
               mbB = m;
               nbB = n;
               kbB = k;
            }
            if (m == n && m == k)
            {
               if (mf > mfS)
               {
                  mfS = mf;
                  mbS = m;
                  nbS = n;
                  kbS = k;
               }
            }
         }
      }
   }
   if (mfB == 0)
   {
      printf("NO KERNEL POSSIBLE FOR RANGE=[%d,%d]\n", b0, bN);
      KillMMNode(mp);
      return(NULL);
   }
   mp->mbB = mbB;
   mp->nbB = nbB;
   mp->kbB = kbB;
   mp->mflop[0] = mfB;
   printf("FOR %d.'%s': MB=%d, NB=%d, KB=%d, MFLOPS=%.1f\n",
          mp->ID, mp->rout, mbB, nbB, kbB, mfB);
   k = MMKernelFailsTest(pre, mbB, nbB, kbB, 0, mp);
   if (!k)
      k = MMKernelFailsTest(pre, mbB, nbB, kbB, 1, mp);
   if (!k)
      k = MMKernelFailsTest(pre, mbB, nbB, kbB, -1, mp);
   if (k)
   {
      printf("KERNEL FAILS TESTER FOR [M,N,K]B=%d,%d,%d\n", mbB, nbB, kbB);
      exit(k);
   }
   if (mbS == 0)
      mp->next = NULL;
   else
   {
      k = MMKernelFailsTest(pre, mbS, nbS, kbS, 0, mp);
      if (!k)
         k = MMKernelFailsTest(pre, mbS, nbS, kbS, 1, mp);
      if (!k)
         k = MMKernelFailsTest(pre, mbS, nbS, kbS, -1, mp);
      if (k)
         mp->next = NULL;
      else
      {
         mp->next = CloneMMNode(mp);
         mp->next->mbB = mbS;
         mp->next->nbB = nbS;
         mp->next->kbB = kbS;
         mp->next->mflop[0] = mfS;
      }
   }
   WriteRefreshedMMFileWithPath(pre, "res", "AMMEXBLKS.sum", mp);
   return(mp);
}

ATL_mmnode_t *TimeExtraBlockings(char pre, int verb)
{
   ATL_mmnode_t *eb;
   eb = ReadMMFileWithPath(pre, "res", "AMMEXBLKS.sum");
   if (!eb)
      return(eb);
   if (eb->mflop[0] < 0)
   {
      ATL_mmnode_t *mp;
      printf("EXTRA BLOCKING FACTOR TIMINGS:\n\n");
      if (verb)
      {
         printf("  MB    NB    KB        MFLOP\n");
         printf("====  ====  ====  ===========\n");
      }
      for (mp=eb; mp; mp = mp->next)
      {
         mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mp->mbB, mp->nbB,
                                     mp->kbB, 1, 0, -1);
         if (verb)
            printf("%4d %5d %5d %11.1f\n",
                   mp->mbB, mp->nbB, mp->kbB, mp->mflop[0]);
      }
      WriteRefreshedMMFileWithPath(pre, "res", "AMMEXBLKS.sum", eb);
   }
   return(eb);
}

ATL_mmnode_t *GetGenKernForNB(char pre, int nb)
{
   static ATL_mmnode_t *mmb=NULL;
   ATL_mmnode_t *mp=NULL, *pM=NULL, *pK=NULL;
   char upr = pre;

   if (upr == 'z')
      upr = 'd';
   if (upr == 'c')
      upr = 's';

   if (!nb)
   {
      if (mmb)
         KillAllMMNodes(mmb);
      return(NULL);
   }
   if (!mmb)
   {
      mmb = ReadMMFileWithPath(upr, "res", "gAMMRES.sum");
      assert(mmb);
      for (mp=mmb; mp; mp = mp->next)
      {
         if (mp->genstr)
            free(mp->genstr);
         if (mp->rout)
            free(mp->rout);
         mp->genstr = mp->rout=NULL;
      }
   }
/*
 * See if any existing kernel can handle this case, modulo ku/kb
 */
   for (mp=mmb; mp; mp = mp->next)
   {
      if ((nb%(mp->mu) == 0) && (nb%(mp->nu) == 0))
      {
         if (!FLAG_IS_SET(mp->flag, MMF_KVEC))
         {
            if (!pM)
               pM = mp;
         }
         else if (!pK && (nb%mp->vlen == 0))  /* K-vec must match on vlen too */
            pK = mp;
      }
      if (pK && pM)
         break;
   }
/*
 * If we've got either pM or pK, all we need to do is possibly adjust kb/ku
 * With same vector length, we'll just take the first one that works.
 */
   if (pM || pK)
   {
      if (pM)
         mp = CloneMMNode(pM);
      else
         mp = CloneMMNode(pK);
      if (FLAG_IS_SET(mp->flag,MMF_KUISKB))
         mp->kbmax = mp->kbmin = mp->ku = mp->kbB = nb;
      else if (nb%(mp->ku))
         mp->ku = (pM) ? 1 : mp->vlen;
      mp->mbB = mp->nbB = mp->kbB = nb;
      mp->rout = MMGetGenName(pre, nb, mp);
      mp->genstr = MMGetGenString(pre, mp);
      return(mp);
   }
/*
 * If changing K info isn't enough, we'll have to make a new kernel that
 * changes possibly a bunch of stuff, including vlen
 *
 * Try to find both a M- & K-vec kernel as a candidate.  Since no kernel
 * working tends to happen more with small NB, prioritize end of queue.
 */
   for (mp=mmb; mp; mp = mp->next)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         pK = mp;
      else
         pM = mp;
   }

/*
 * Find a legal M-vec (or unvectorized) kernel based on best M kernel
 */
   if (pM)
   {
      int vl=pM->vlen, u;
      mp = CloneMMNode(pM);
/*
 *    Find a vlen compatible with this MB; may become 1 -> unvectorized
 */
      while (vl > 1 && nb % vl)  /* make vlen evenly divide nb */
         vl >>= 1;
      mp->vlen = vl;
      mp->mu = (pM->mu / pM->vlen) * vl;  /* use same # of regs */
/*
 *    Make mu evenly divide MB (vl already does)
 */
      for (u=mp->mu; u > vl && nb%u; u -= vl);
      assert(nb%u == 0);
      mp->mu = u;
/*
 *    Make nu evenly divide NB
 */
      for (u=mp->nu; nb%u; u--);
      if (FLAG_IS_SET(mp->flag, MMF_NOBCAST))
      {
         if (u > vl)
            u = (u/vl)*vl;
         else
            mp->flag &= (1<<MMF_NOBCAST);
      }
      mp->nu = u;
/*
 *    Make ku evenly divide KB
 */
      if (nb%(mp->ku) != 0)
         mp->ku = 1;
      pM = mp;
   }
/*
 * Find a legal K-vec kernel if possible
 */
   if (pK)
   {
      int mu=0, nu=0, vl;
      mp = NULL;
/*
 *    Make vlen evenly divide KB
 */
      for(vl=pK->vlen; vl > 1 && nb%vl; vl >>= 1);
      if (vl > 1)
      {
         int i, j;
         float minrat=0.0;
/*
 *       Find mu & nu such that mu*nu % vl == 0
 */
         for (i=pK->mu; i > 0; i--)
         {
            for (j=pK->nu; j > 0; j--)
            {
               if ((i*j)%vl == 0 && i%nb == 0 && j%nb == 0)
               {
                  float ratio = (i+j)/(i*j);
                  if (minrat == 0.0 || ratio < minrat)
                  {
                     minrat = ratio;
                     mu = i;
                     nu = j;
                  }
               }
            }
         }
/*
 *       Do we have a legal k-vectorized kernel?
 */
         if (vl > 1 & mu > 1 && nu > 1)
         {
/*
 *          Make ku divide KB & vl
 */
            i = pK->ku;
            if (i%vl)
            {
               if (i > vl)
                  i = (i/vl)*vl;
               else
                  i = vl;
            }
            for (; i > vl && nb%i; i -= vl);
            if (nb%i)
            {
               mp = pK = CloneMMNode(pK);
               pK->vlen = vl;
               pK->mu = mu;
               pK->nu = nu;
               pK->ku = i;
            }
         }
      }
      pK = mp;
   }
/*
 * If we've found no legal kernel, create a scalar kernel that will work
 */
   if (!pM && !pK)
   {
      mp = pM = GetMMNode();
      mp->mu = (nb&3) ? 1 : 4;
      if (nb%6 == 0)
         mp->mu = 6;
      else if ((nb & 3) == 0)
         mp->mu = 4;
      if ((nb & 1) == 0)
         mp->mu = 2;
      else
         mp->mu = 1;
      mp->nu = mp->ku = mp->vlen = 1;
   }
/*
 * put best kernel in pM, free other
 */
   if (pM && pK)  /* wt both K & M, guess best */
   {
      if (pK->vlen > pM->vlen)
      {
         KillMMNode(pM);
         pM = pK;
      }
      else if (pM->vlen > pK->vlen)
         KillMMNode(pK);
      else /* vlen same, break tie on relative speed & load/use ratio */
      {
         double krat, mrat;
         krat = (pK->mu + pK->nu);
         krat /= (pK->mu * pK->nu);
         mrat = pK->mu / pK->vlen;
         mrat = (mrat+pK->nu) / (mrat*pK->nu);
         if (pM->mflop[0] != 0.0 && pK->mflop[0] != 0)
            krat *= pM->mflop[0] / pK->mflop[0];
         if (krat < 0.0)
            krat = -krat;
         if (krat < mrat)
         {
            KillMMNode(pM);
            pM = pK;
         }
         else
            KillMMNode(pK);
      }
   }
   else if (!pM)
      pM = pK;
   assert(pM);
   if (!FLAG_IS_SET(pM->flag,MMF_KRUNTIME))
      pM->kbB = nb;
   pM->rout = MMGetGenName(pre, nb, pM);
   pM->genstr = MMGetGenString(pre, pM);
   pM->mbB = pM->nbB = pM->kbB = nb;
   return(pM);
}

void PrintGen0(FILE *fp, ATL_mmnode_t *mp, int mb, int nb, int kb)
{
   fprintf(fp, "B=(%d,%d,%d), U=(%d,%d,%d), pf=(%x,%d), flg=%x",
           mb?mb:mp->mbB, nb?nb:mp->nbB, kb?kb:mp->kbB, mp->mu, mp->nu,
           FLAG_IS_SET(mp->flag, MMF_KUISKB) ? -1:mp->ku,
           mp->pref, mp->pfLS, mp->flag);
}

void PrintGen(FILE *fp, ATL_mmnode_t *mp, int mb, int nb, int kb, double mf)
/*
 * Prints single line description of mp to fp
 */
{
   fprintf(fp, "   0.");
   PrintGen0(fp, mp, mb, nb, kb);
   fprintf(fp, ": MFLOP=%.0f\n", mf);
}

ATL_mmnode_t *BestForThisNB
(
   int verb,
   char pre,
   ATL_mmnode_t *mmb,
   int nb,
   int pnb,  /* previous nb */
   int nnb   /* next nb */
)
/*
 * Times all kernels in mmb
 * RETURNS: ptr to best performing kernel, empty gen node if no user case wrks
 */
{
   ATL_mmnode_t *mmp, *mmB=NULL;
   double mf, mf0, mfB=0.0;

   printf("SCOPING FOR BEST PERFORMING KERNEL FOR NB=%d\n", nb);
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      const int kb0 = mmp->kbB, ku0 = mmp->ku;
      char *gs0=mmp->genstr;
      int kb, kbOK;
/*
 *    Choose kb, if forced only kb will do, so skip if kernel can't do it
 *    Genned kerns can be adapted, so are checked differently frm user kerns.
 */
      if (!mmp->ID)  /* kvec OK wt any mul of vlen, mvec OK wt any K */
         kbOK = FLAG_IS_SET(mmp->flag, MMF_KVEC) ? (nb%mmp->vlen == 0):1;
      else
      {
         kbOK = (mmp->kbmin) ? (nb >= mmp->kbmin) : 1;
         if (kbOK && mmp->kbmax)
            kbOK = nb <= mmp->kbmax;
      }
      if (!kbOK || ((nb/mmp->mu)*mmp->mu != nb) || ((nb/mmp->nu)*mmp->nu != nb)
          || ((nb/mmp->ku)*mmp->ku != nb) || (nb == pnb) || (nb == nnb))
      {

         printf("   %d. %s: SKIPPED, bad NB\n", mmp->ID, mmp->rout);
         continue;
      }
/*
 *    Generated files may need to get genstr and related info corrected
 */
      if (mmp->ID == 0)
      {
         if (FLAG_IS_SET(mmp->flag, MMF_KUISKB))
             mmp->kbmax = mmp->kbmin = mmp->ku = mmp->kbB = nb;
         if (!FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
             mmp->kbB = nb;
         mmp->genstr = MMGetGenString(pre, mmp);
      }
      mf0 = TimeMMKernel(verb, 0, mmp, pre, nb, nb, nb, 1, 0, -1);
      if (mmp->ID == 0)  /* put original info back in queue */
      {
         free(mmp->genstr);
         mmp->genstr = gs0;
         mmp->kbB = kb0;
         mmp->ku = ku0;
      }
/*
 *    Give bonus to K-runtime variable over K-compile time; K-runtime kernels
 *    can be used for some K-cleanup, and they can be used for any required KB
 *    as well as being typically much smaller instruction load, so they are
 *    strongly preferred
 */
      mf = FLAG_IS_SET(mmp->flag, MMF_KRUNTIME) ? mf0*KRUNMUL : mf0;
      if (mf > mfB)
      {
         mfB = mf;
         mmB = mmp;
         mmB->mbB = mmB->nbB = mmB->kbB = nb;
      }
      if (mmp->ID)
         printf("   %d. %s: kb=%d, flg=%x, MFLOP=%.2f\n",
                mmp->ID, mmp->rout, nb, mmp->flag, mf0);
      else
         PrintGen(stdout, mmp, nb, nb, nb, mf0);
   }
   if (mmB)
   {
      mmB = CloneMMNode(mmB);
      mmB->mflop[0] = mfB;
      if (!mmB->ID)  /* specialize generated code for this KB */
      {
         if (mmB->genstr)
            free(mmB->genstr);
         if (mmB->rout)
            free(mmB->rout);
         mmB->mbB = mmB->nbB = mmB->kbB = nb;
         if (FLAG_IS_SET(mmB->flag, MMF_KUISKB))
             mmB->kbmin = mmB->kbmax = mmB->ku = nb;
         mmB->rout = MMGetGenName(pre, nb, mmB);
         mmB->genstr = MMGetGenString(pre, mmB);
      }
   }
   else
   {
      mmB = GetGenKernForNB(pre, nb);
      assert(mmB);
   }
   if (MMKernelFailsAnyBeta(pre, nb, nb, nb, mmB))
   {
      printf("BEST KERNEL FAILS TESTER FOR NB=%d\n", nb);
      exit(1);
   }
   mmB->mflop[0] = TimeMMKernel(verb, 0, mmB, pre, nb, nb, nb, 1, 0, -1);
   printf("BEST KERNEL FOUND FOR NB=%d: ID#%d '%s' %.2f MFLOPS\n\n",
          nb, mmB->ID, mmB->rout, mmB->mflop[0]);
   return(mmB);
}

int DeleteBadBigNBs(ATL_mmnode_t *mmb, int *nbs)
{
   ATL_mmnode_t *best=NULL, *mmp;
   double mfB=0.0;
   int n=0;
/*
 * Find the best-performing kernel
 */
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      double mf;
      mf = mmp->mflop[0];
      if (mf > mfB)
      {
         mfB = mf;
         best = mmp;
      }
   }
/*
 * Delete all NBs larger than best
 */
   while (best->next)
   {
      best->next = KillMMNode(best->next);
      n++;
   }
   if (n)
   {
      int N = *nbs;
      N = (N >= 0) ? N : -N;
      printf("Deleted %d large, slow kernels starting at NB=%d\n",
             n, nbs[N-n+1]);
   }
   return(n);
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
   fprintf(stderr, "   -n # nb1 ... nb# : NBs to try for\n");
   fprintf(stderr, "   -N # nb1 ... nb# : force exact NBs in search\n");
   fprintf(stderr, "   -r <nreg> : set max # of registers to try\n");
   fprintf(stderr, "   -b <nb>   : set initial block factor to try\n");
   fprintf(stderr, "   -v <verb> : set verbosity (1)\n");
   exit(ierr ? ierr : -1);
}

void GetFlags(int nargs, char **args, char *PRE, int *verb, int *NREG,
              int *NB, int *CS)
{
   ATL_mmnode_t *mmb=NULL, *mp;
   int B0, BN;
   int i, j=0, n, k;
   char pre='d';
   int *nbs=NULL;
   *NREG = *NB = 0;
   *verb = 1;
   *CS = 0;
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
      case 'r':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *NREG = atoi(args[i]);
         break;
      case 'v':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *verb = atoi(args[i]);
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
   *PRE = pre;
/*
 * NREG has been stored by search in ivar.  Read it, and then zero ivar so
 * it won't propogate, making .sum files confusing
 */
   if (*NREG == 0)
   {
      ATL_mmnode_t *mp;
      char upr=pre;
      if (pre == 'z')
         upr = 'd';
      else if (pre == 'c')
         upr = 's';
      mp = ReadMMFileWithPath(upr, "res", "gmvAMMUR.sum");
      if (mp)
      {
         *NREG = mp->ivar;
         KillAllMMNodes(mp);
      }
   }
   if (*CS == 0)
      *CS = GetL1CacheElts(pre);
}
static void ApplyMoves2Flag
(
   ATL_mmnode_t *mmp,  /* kernel to set MMF_MV[A,B,C] flag bits */
   int mvBits          /* last 3 bits: MOVE_[CBA] */
)
{
   int flag = mmp->flag & (~MMF_MVSET);         /* zero existing move bits */
   mmp->flag = flag | ((mvBits & 7)<<MMF_MVA); /* put new move bits in */
}
static void ApplyMoves2Flags
(
   ATL_mmnode_t *mmb,  /* kernel to set MMF_MV[A,B,C] flag bits */
   int mvBits          /* last 3 bits: MOVE_[CBA] */
)
{
   const unsigned int mvMSK = ~MMF_MVSET, mvSET = (mvBits&7)<<MMF_MVA;
   ATL_mmnode_t *mmp;
   for (mmp=mmb; mmp; mmp = mmp->next)
      mmp->flag = ((mmp->flag) & mvMSK) | mvSET;
}

ATL_mmnode_t *GetNewKCleanGenNode
(
   char pre,
   ATL_mmnode_t *kp,  /* kernel we are generating K-cleanup for */
   int mb,
   int nb,
   int kb
)
{
   ATL_mmnode_t *p;
   const int mu=kp->mu, nu=kp->nu;
   int kvec=0, ku=1, vl, vmu;
   int kmaj = FLAG_IS_SET(kp->flag, MMF_KVEC) ? kp->vlen:0;

   if (kmaj > 1)
   {
      kvec = 1;
      ku = vl = kmaj;
   }
   else
   {
      if (kp->vlen)
      {
         vl = kp->vlen;
         if (mu%vl)
         {
            fprintf(stderr, "kern=%d '%s', mu=%d, vlen=%d\n", kp->ID,
                    kp->rout?kp->rout:"NULL", kp->mu, kp->vlen);
         }
         assert(mu%vl == 0);
      }
      else
         vl = 1;
   }
   printf("  TRY: mu=%d, nu=%d, ku=%d, vl=%d, kvec=%d\n", mu, nu, ku, vl, kvec);
   p = MMGetNodeGEN(pre, 0, 0, mu, nu, ku, vl, kvec, 0, 0, NULL);
   p->mbB = mb;
   p->nbB = nb;
   p->kbB = kb;
   return(p);
}

ATL_mmnode_t *FindDefMUNU(int verb, char pre, int nreg, int lat, int nb, int ku,
                          int *MU, int *NU)
{
   ATL_mmnode_t *mmp;
   double mf, mfB=0.0;
   int n, i, j, kb, muB=1, nuB=1, VL, chkVL=0;

/*   mmp = ReadMMFileWithPath(pre, "res", "gAMMMUNU.sum"); */
   mmp = NULL;
   if (mmp)
   {
      MMFillInGenStrings(pre, mmp);
      nb = mmp->kbB;
      if (mmp->mflop[0] < 0.0)
         mmp->mflop[0] = TimeMMKernel(verb, 1, mmp, pre, nb, nb, nb, 1, 0, -1);
      printf("READ IN BEST GENNED MU=%d, NU=%d, MFLOP=%.2f\n\n",
             mmp->mu, mmp->nu, mmp->mflop[0]);
#if 0
/*
 *    See if there is a mismatch between vector settings
 */
      if (mmp->vlen != VLEN[VECi])
      {
         printf("\n\n!!! WARNING: TURNING OFF VECTORIZATION DUE TO MISMATCHED VLEN IN 'res/%cAMMMUNU.sum!!!!\n\n", pre);
         VECi = VTSC;
      }
      *MU = mmp->mu / VLEN[VECi];
      assert(*MU);
      *NU = mmp->nu;
#endif
      return (mmp);
   }
   VL = GetNativeVLEN(pre);
   if (!VL)
      chkVL = VL = (pre == 's' || pre == 'c') ? 4 : 2;  /* temp kludge */
   mmp = MMGetNodeGEN(pre, 0, nb, 1, 1, ku, 1, 0, 0, 0,
                      DupString("ATL_Xamm_munu.c"));
   if (pre == 'z')
      mmp->rout[4] = 'd';
   else if (pre == 'c')
      mmp->rout[4] = 's';
   else
      mmp->rout[4] = pre;
   mmp->mbB = mmp->nbB = mmp->kbB = nb;
   mmp->vlen = VL;
/*
 * Try all near-square register blocking cases
 */
   printf("Finding best MUxNU case for nb=%d\n", nb);
   for (n=4; n < nreg; n++)
   {
      int mbu, nbu, mu, nu;
      for (j=1; j*j < n; j++);
      i = n / j;
      if (nb%i || nb%j)
         continue;
      mu = mmp->mu = i * VL;
      nu = mmp->nu = j;
      if (mmp->genstr)
        free(mmp->genstr);
      mbu = (nb >= mu) ? (nb/mu)*mu : mu;
      nbu = (nb >= nu) ? (nb/nu)*nu : nu;
      mmp->genstr = MMGetGenString(pre, mmp);
      mf = TimeMMKernel(verb, 1, mmp, pre, mbu, nbu, nb, 1, 0, -1);
      printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
      if (mf > mfB)
      {
         muB = i;
         nuB = j;
         mfB = mf;
      }
   }
/*
 * For non-AVX x86, try 1-D cases since they are 2-operand assemblies
 */
   #if (defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)) && !defined(ATL_AVX)
      printf("BEST NEAR-SQUARE CASE IS MU=%d, NU=%d, MFLOP=%.2f\n\n",
             muB, nuB, mfB);
      printf("Finding best 1-D outer loop unrolling for nb=%d\n", nb);
      for (n=2; n < nreg; n++)
      {
         int mbu, nbu, mu, nu;
         i = 1; j = n;
         if (nb % n)
            continue;
         mu = mmp->mu = i*VL;
         nu = mmp->nu = j;
         if (mmp->genstr)
           free(mmp->genstr);
         mmp->genstr = MMGetGenString(pre, mmp);
         mbu = (nb >= mu) ? (nb/mu)*mu : mu;
         nbu = (nb >= nu) ? (nb/nu)*nu : nu;
         mf = TimeMMKernel(verb, 1, mmp, pre, mbu, nbu, nb, 1, 0, -1);
         printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
         if (mf > mfB)
         {
            muB = i;
            nuB = j;
            mfB = mf;
         }
         i = n; j = 1;
         mu = mmp->mu = i * VL;
         nu = mmp->nu = j;
         mbu = (nb >= mu) ? (nb/mu)*mu : mu;
         nbu = (nb >= nu) ? (nb/nu)*nu : nu;
         if (mmp->genstr)
           free(mmp->genstr);
         mmp->genstr = MMGetGenString(pre, mmp);
         mf = TimeMMKernel(verb, 1, mmp, pre, mbu, nbu, nb, 1, 0, -1);
         printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
         if (mf > mfB)
         {
            muB = i;
            nuB = j;
            mfB = mf;
         }
      }
   #endif

   i = FLAG_IS_SET(mmp->flag, MMF_KVEC);
   KillMMNode(mmp);
   mmp = MMGetNodeGEN(pre, 0, nb, (i)?muB:muB*VL, nuB, ku, VL, i, 0, 0, NULL);
   WriteRefreshedMMFileWithPath(pre, "res", "gAMMMUNU.sum", mmp);
   printf("BEST CASE IS MU=%d, NU=%d, MFLOP=%.2f (%.2f)\n\n",
          muB, nuB, mf, mfB);
   *MU = muB;
   *NU = nuB;
   return(mmp);
}

#if 0
void GetMUNUbyNB(int nb, int nreg, int *MU, int *NU)
{
   int mu=(*MU), nu=(*NU), vmu=mu*VLEN;

   assert(mu && nu && !(nb%VLEN[VECi]));
   if (!(nb%vmu) && !(nb%nu))
      return;
   if (nu == 1) /* handle MUx1 by decreasing by VLEN */
   {
      int u = vmu, vlen = VLEN[VECi];
      while (u+u+1 <= nreg && nb%u)
         u += vlen;
      if (u+u+1 > nreg)
         u -= vlen;
      while (nb%u)
         u -= vlen;
      assert(u);
      *MU = u / vlen;
      return;
   }
   if (mu == 1 || nu == 1) /* handle 1-D cases by just inc/dec U */
   {
      int u = (mu == 1) ? nu : mu;
      while (u+u+1 <= nreg && nb%u)
         u++;
      if (u+u+1 > nreg)
         u--;
      while (nb%u)
         u--;
      if (mu == 1)
         *NU = u;
      else
         *MU = u;
      return;
   }
   if (nb%vmu)  /* mu can't handle NB */
   {
      int i;
/*
 *    try increasing mu until we run out of registers
 */
      for (i=mu+1; i*nu+i+nu <= nreg; i++)
         if (!(nb%(i*VLEN[VECi])))
            break;
/*
 *    Try decreasing mu until it divides
 */
      if (nb%(i*VLEN[VECi]) || i*nu+i+nu > nreg)
      {
         for (mu--; mu; mu--)
            if (!(nb%(mu*VLEN[VECi])))
               break;
      }
      else
         mu = i;
   }
   if (nb%nu)  /* nu can't handle NB */
   {
      int i;
/*
 *    try increasing nu until we run out of registers
 */
      for (i=nu+1; i*mu+i+mu <= nreg; i++)
         if (!(nb%i))
            break;
/*
 *    Try decreasing nu until it divides
 */
      if (nb%i || i*mu+i+mu > nreg)
      {
         for (nu--; nu; nu--)
            if (!(nb%nu))
               break;
      }
      else
         nu = i;
   }
   *MU = mu;
   *NU = nu;
}
#endif

int FindNBInArray(int nb, int *nbs)
/*
 * RETURNS: location+1, or 0 if not found
 */
{
   int i, n = (nbs[0] > 0) ? nbs[0] : -nbs[0];
   for (i=1; i <= n; i++)
       if (nbs[i] == nb)
          return(i);
   return(0);
}
#if 0
ATL_mmnode_t *CreateGenCasesFromNBs
(
   ATL_mmnode_t *mmb,   /* best user-contributed kernels */
   char pre,            /* precision: s/d */
   int *nbs,            /* list of desired NBs */
   int nreg,            /* upper bound on register use */
   int MU, int NU,      /* default M/N unrolling */
   int KU               /* -1 for fully unrolled, else unrolling factor */
)
/*
 * Generate a list of generated kernels, with the union of nb's in nbs
 * and mmb, and return the generated nodes for timing.
 * HERE HERE HERE: this code is crap, needs to merge both lists, not user
 * list twice.
 */
{
   ATL_mmnode_t *mp, *umb=NULL, *ap;
   int i, n = (nbs[0] > 0) ? nbs[0] : -nbs[0], ne=0, *enbs;

/*
 * Create new queue with an entry for all NBs; both lists (mmb & nbs) are
 * sorted in increasing size
 */
   if (!n && !mmb)
      return(NULL);
   n++;
   ap = mmb;  /* add ptr */
   i = 1;     /* ptr to normal block under consideration */
   do
   {
      int nb, mu=MU, nu=NU, ku;
      ATL_mmnode_t *p=NULL;
      if (ap && i < n)  /* must choose amongst blocks */
      {
         nb = ap->kbB;
         nb = Mmin(nb, nbs[i]);
         if (nb == ap->kbB)
            ap = ap->next;
         if (nb == nbs[i])
            i++;
      }
      else if (ap)
      {
         nb = ap->kbB;
         ap = ap->next;
      }
      else
         nb = nbs[i++];
      ku = (KU == -1) ? nb : KU;

/*
 *    If NB is not a multiple of VLEN, drop down to shorter ops
 */
      if (nb%VLEN[VECi])
      {
         int vl=VECi;
/*
 *       For AVX, see if dropping to SSE will fix problem
 */
         if (VECi == VTAVX && !(nb%VLEN[VTSSE]))  /* AVX can drop to SSE */
         {
            VECi = VTSSE;
            GetMUNUbyNB(nb, nreg, &mu, &nu);
            p = GetNewGenNode(pre, nb, 0, mu, nu, ku, 0);
         }
         if (!p)
         {
            VECi = VTSC;
            GetMUNUbyNB(nb, nreg, &mu, &nu);
            p = GetNewGenNode(pre, nb, 0, mu, nu, ku, 0);
         }
         VECi = vl;
      }
      else
      {
         GetMUNUbyNB(nb, nreg, &mu, &nu);
         p = GetNewGenNode(pre, nb, 0, mu, nu, ku, 0);
      }
      if (umb)
      {
         mp->next = p;
         mp = p;
      }
      else
         umb = mp = p;
   }
   while (i < n || ap);
   return(umb);
}
void SetGenVec(int verb, char pre)
/*
 * This routine uses a simple timing to be sure if vectorization helps or not
 */
{
   ATL_mmnode_t *mp;
/*
 * If vector operations are being used, make sure they work; compiler and
 * flag changes can mess them up, and in this case we'll fall back to
 * scalar generation.  Try to see if we can successfully test simplist
 * possible vector kernel, and fall back to scalar kernel if we can't
 */
   if (VLEN[VECi] < 2)
      return;
   mp = GetNewGenNode(pre, 32, 0, 1, 1, 1, 0);
   if (MMKernelFailsTest(pre, 32, 32, 32, 1, mp))
   {
      printf("ERROR: VEC='%s' FAILED, genstr='%s'!\n",VECs[VECi],mp->genstr);
      KillMMNode(mp);
/*
 *    For AVX, try falling back to SSE
 */
      if (VECi == VTAVX)
      {
         VECi = VTSSE;
         KillMMNode(mp);
         mp = GetNewGenNode(pre, 32, 0, 1, 1, 1, 0);
         if (MMKernelFailsTest(pre, 32, 32, 32, 1, mp))
            VECi = VTSC;
      }
      else
         VECi = VTSC;
   }
   KillMMNode(mp);
/*
 * For AVX, switch to SSE if AVX doesn't offer a performance advantage
 * (as on AMD Dozer), since SSE smaller code size and requires less cleanup
 */
   if (VECi == VTAVX)
   {
      double mfA, mfS, mf;
      char *sp;
      int vl;
      mp = GetNewGenNode(pre, 128, 0, 1, 4, 1, 0);
      mfA = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 1, 0, -1);
      KillMMNode(mp);
      mp = GetNewGenNode(pre, 128, 0, 2, 2, 1, 0);
      mf = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 1, 0, -1);
      KillMMNode(mp);
      if (mf > mfA)
         mfA = mf;
      vl = VECi;
      VECi = VTSSE;
      mp = GetNewGenNode(pre, 128, 0, 1, 4, 1, 0);
      mfS = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 1, 0, -1);
      KillMMNode(mp);
      mp = GetNewGenNode(pre, 128, 0, 2, 2, 1, 0);
      mf = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 1, 0, -1);
      if (mf > mfA)
         mfA = mf;
      KillMMNode(mp);
      if (mfA < 1.03*mfS)
         printf("USING SSE INSTEAD OF AVX, AVX=%.2f, SSE=%.2f\n", mfA, mfS);
      else
      {
         printf("AVX GOOD, AVX=%.2f, SSE=%.2f\n", mfA, mfS);
         VECi = vl;
      }
   }
/*
 * For any system, don't use vector instructions if they aren't faster than
 * scalar.
 */
   if (VLEN[VECi] > 1)
   {
      double mfV, mfS;
      char *sp;
      int vl;
      mp = GetNewGenNode(pre, 128, 0, 1, 4, 1, 0);
      mfV = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 1, 0, -1);
      KillMMNode(mp);
      vl = VECi;
      VECi = VTSC;
      mp = GetNewGenNode(pre, 128, 0, 1, 4, 1, 0);
      mfS = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 1, 0, -1);
      KillMMNode(mp);
      if (mfV < 1.05*mfS)
         printf("USING SCALAR INSTEAD OF VECTOR, VEC=%.2f, SCALAR=%.2f\n",
                mfV, mfS);
      else
      {
         printf("VEC GOOD, VEC=%.2f, SCALAR=%.2f\n", mfV, mfS);
         VECi = vl;
      }
   }
   printf("GENERATING WITH VEC='%s', VLEN=%d\n\n", VECs[VECi], VLEN[VECi]);
}
#endif



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

#define HUGE_NB 180
ATL_mmnode_t *WinnowHugeNB
(
   int imf,
   ATL_mmnode_t *mb  /* queue of cases */
)
/*
 * Removes any NB >= HUGE_NB that aren't at least 2% faster than smaller cases
 */
{
   ATL_mmnode_t *mp, *p, *prev=mb;
   double mfB;

   if (!mb || !mb->next)
      return(mb);
   mp = mb->next;
/*
 * Find best-performing kernel below HUGE_NB
 */
   mfB = mp->mflop[imf];
   for (mp=mb->next; mp; mp = mp->next)
   {
      if (mp->mbB < HUGE_NB && mp->nbB < HUGE_NB && mp->kbB < HUGE_NB)
         mfB = Mmax(mfB, mp->mflop[imf]);
      else
         break;
      prev = mp;
   }
/*
 * If no kernels above threshold, return original queue
 */
   if (!mp)
      return(mb);
/*
 * mp points to first NB above threshold, but there is no point in deleting
 * small NB if we leave large NB, so delete only from the end of queue
 */
  do
  {
     for (p=mp; p->next; p = p->next);
     if (p->mflop[imf] <= 1.02*mfB)
        mp = RemoveMMNodeFromQ(mp, p);
     else  /* stop removing stuff */
        break;
  }
  while (mp);
  prev->next = mp;
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

ATL_mmnode_t *MergeAndWinnowCases
(
   int verb,
   char pre,
   ATL_mmnode_t *umb, /* queue of user cases */
   ATL_mmnode_t *gmb  /* genned cases, always include NBs of umb */
)
/*
 * Merges user and gmp cases, while getting rid of cases that get worse
 * performance than their smaller blocks; FREES umb and gmb
 * RETURNS: new merged and winnowed queue
 */
{
   ATL_mmnode_t *mmb=NULL, *mmp, *gmp, *ump=umb;
   for (gmp=gmb; gmp; gmp = gmp->next)
   {
      ATL_mmnode_t *p;
      if (ump)
      {
         if (ump->kbB == gmp->kbB)
         {
            if (gmp->mflop[0] >= ump->mflop[0])
               p = CloneMMNode(gmp);
            else
               p = CloneMMNode(ump);
            ump = ump->next;
         }
         else
            p = CloneMMNode(gmp);
      }
      else
         p = CloneMMNode(gmp);
      p->next = NULL;
      if (mmb)
      {
/*
 *       If larger NB isn't faster than smaller one, kill it for nb >= 16
 */
         if (p->kbB >= 16 && mmp->mflop[0] > p->mflop[0])
            KillMMNode(p);
         else
         {
            mmp->next = p;
            mmp = p;
         }
      }
      else
         mmp = mmb = p;
   }
   KillAllMMNodes(umb);
   KillAllMMNodes(gmb);
   mmb = WinnowCases(0, mmb);
   return(mmb);
}



int FailKCleanTests(char pre, int nb, ATL_mmnode_t *kp)
/*
 *  This routine tests if a kernel is suitable for use in K-cleanup by
 *  doing testing with ku=1, kb=0, and tries all K values between 1 and nb
 *  RETURNS: 0 if kernel passes all tests, else non-zero
 */
{
   int i, beg, end, inc, kmaj = FLAG_IS_SET(kp->flag, MMF_KVEC) ? kp->vlen:0;;

   if (!FLAG_IS_SET(kp->flag, MMF_KRUNTIME) ||
       (kp->ku != 1 && kp->ku != kmaj))
      return(-1);
   printf("TESTING ID=%d, rout='%s', nb=%d, mu=%d, nu=%d for K-cleanup:\n",
          kp->ID, kp->rout, nb, kp->mu, kp->nu);

   if (kmaj > 1)
   {
      inc = beg = kmaj;
      end = ((nb+inc-1)/inc)*inc;
   }
   else
   {
      beg = inc = 1;
      end = nb;
   }
   for (i=beg; i <= end; i += inc)
   {
      int ierr;
      ierr = MMKernelFailsTest(pre, nb, nb, i, 0, kp);
      if (ierr)
      {
         printf("  K=%d: FAILED!\n", i);
         return(ierr);
      }
      else
         printf("  K=%d: PASSED!\n", i);
   }
   printf("PASSED ALL K-tests!\n\n");
   return(0);
}
ATL_mmnode_t *GetUniqueKClean(int verb, char pre, ATL_mmnode_t *mmb)
/*
 * OUTPUT: <pre>AMMKCLEAN.sum: all unique kerns to be compiled
 */
{
   ATL_mmnode_t *mp, *gmmb, *ummb, *ub, *np, **dlmm;
   int nn=0, nd=0, n=0;  /* #needed & done, total, copy of done */
   int *dl, *nl;         /* done and needed lists */
   int i;
   gmmb = ReadMMFileWithPath(pre, "res", "AMMKCLEAN.sum");
   if (gmmb)
   {
      printf("READING IN UNIQUE K-CLEANUP:\n");
      MMFillInGenStrings(pre, gmmb);
      for (mp=gmmb; mp; mp = mp->next)
      {
         int mb = (mp->nbB > mp->mu) ? (mp->nbB/mp->mu)*mp->mu : mp->mu;
         int nb = (mp->nbB > mp->nu) ? (mp->nbB/mp->nu)*mp->nu : mp->nu;
         int kb = (nb > 8) ? (nb>>2) : nb, KB = kb;
         int ku = mp->ku, kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;;
         if (kmaj > 1)
            KB = ((kb+ku-1)/ku)*ku;
         if (mp->mflop[0] < 0.0)
         {
            mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mb, nb, KB, 0, 0, -1);
            mp->mflop[0] *= (double)kb / (double)KB;
         }
         printf("   nb=%d,  kb=%d, mu=%d, nu=%d, MFLOP=%.2f\n",
                nb, kb, mp->mu, mp->nu, mp->mflop[0]);
      }
      printf("Done.\n");
      return(gmmb);
   }
/*
 * Find out how many total kernels, and how many already have their own
 * cleanup (nd, number done).  This nd may be bigger than it should, because
 * we can't guarantee they are unique
 */
   for (mp=mmb; mp; mp = mp->next, n++)
   {
      const int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
      if ((mp->ku == 1 || (kmaj == mp->ku)) &&
          FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
        nd++;
   }
   dl = malloc(8*n*sizeof(int));
   assert(dl);
   if (nd)
   {
      dlmm = malloc(nd*sizeof(ATL_mmnode_t*));
      assert(dlmm);
   }
   else
      dlmm = NULL;

   nl = dl + (n<<2);
   nd = 0;
/*
 * First, go back through kernels, and add kernels that can serve as K-cleaners
 * to the done list
 */
   for (mp=mmb; mp; mp = mp->next, n++)
   {
      const int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
      if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME) &&
          (mp->ku == 1 || kmaj == mp->ku))
      {
         const int nd4 = (nd<<2), mu=mp->mu, nu=mp->nu;
/*
 *       See if trip is already in done list if so, no new entry, just update kb
 *       and cleanup kernel entry
 */
         for (i=0; i < nd4; i += 4)
            if (mu == dl[i] && nu == dl[i+1] && kmaj == dl[i+2])
               break;
         if (i < nd4)
         {                      /* (larger NB always later) */
            dl[i+3] = mp->kbB;  /* take largest kbB that matches mu/nu */
            dlmm[i>>2] = mp;
            continue;
         }
         else
         {
            dl[nd4] = mu;
            dl[nd4+1] = nu;
            dl[nd4+2] = kmaj;
            dl[nd4+3] = mp->kbB;
            dlmm[nd++] = mp;
         }
      }
   }
/*
 * Delete any kernels from dl that fail to actually work for K cleaning
 */
   for (i=0; i < nd; i++)
   {
      if (FailKCleanTests(pre, dlmm[i]->kbB, dlmm[i]))
      {
         const int i4=(i<<2), nc=nd-i-1;
         if (nc > 0)
         {
            memcpy(dl+i4, dl+i4+4, (nc<<2)*sizeof(int));
            memcpy(dlmm[i], dlmm[i+1], nc*sizeof(ATL_mmnode_t*));
         }
         nd--;
      }
   }

/*
 * Find all unique (mu,nu,kmaj) combos that still need to to be cleaned;
 * there will be nn (# needed) of these, and we'll save (mu,nu,MAXNB) in
 * needed list (nl).
 * We use MAXNB for testing (large NB tests mosts cases of K).
 * Combos that are handled by the done list (dl) aren't added to needed list.
 */
   for (mp=mmb; mp; mp = mp->next, n++)
   {
      int mu=mp->mu, nu=mp->nu, nn4=(nn<<2), nd4=(nd<<2);
      int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
/*
 *    See if pair is already in done list or needed list, if so, no change
 */
      for (i=0; i < nd4; i += 4)
         if (mu == dl[i] && nu == dl[i+1] && kmaj == dl[i+2])
            break;
      if (i < nd4)    /* if it was found in the done list */
         continue;    /* this combo is already handled */
/*
 *    If we reach here, combo is not handled, must add to needed list
 */
      for (i=0; i < nn4; i += 4)
         if (mu == nl[i] && nu == nl[i+1] && kmaj == nl[i+2])
            break;
      if (i < nn4)            /* If already in needed list */
      {
         nl[i+3] = mp->kbB;   /* just update kb so we get largest for testing */
         continue;
      }
/*
 *    If we haven't seen this pair before, add to needed list
 */
      else
      {
         nl[nn4] = mu;
         nl[nn4+1] = nu;
         nl[nn4+2] = kmaj;
         nl[nn4+3] = mp->kbB;
         nn++;
      }
   }
/*
 * Now, create a queue of generated kernels for each needed pair, and time
 * it's maxNB performance.
 */
   gmmb = NULL;
   printf("Timing Generated K-cleanup:\n");
   for (i=0; i < nn; i++)
   {
      ATL_mmnode_t *p;
      const int i4 = (i<<2), mu=nl[i4], nu=nl[i4+1], kmaj=nl[i4+2];
      int nb = Mmax(nl[i4+3],nu), mb = (nb > mu) ? (nb/mu)*mu : mu, ku;
      const int kb = (nb > 8) ? (nb>>2) : nb;
      int vl=VLEN, vmu, KK;
      double mf;
/*
 *    HERE HERE: Improve KMAJ when generator is extended!
 */
      if (kmaj > 1)
      {
         while ((mu*nu)%vl)
            vl >>= 1;
         vmu = mu;
         ku = vl;
      }
      else
      {
         while (mu%vl)
            vl >>= 1;
         ku = 1;
      }
/* HERE HERE */
      p = MMGetNodeGEN(pre, 0, 0, mu, nu, ku, vl, kmaj, 0, 0, NULL);
      p->mbB = mb;
      p->nbB = nb;
      p->kbB = kb;
      p->flag |= (1<<MMF_KRUNTIME);
      #if 1  /* by default don't waste time testing generated code */
         assert(!FailKCleanTests(pre, nb, p));
      #endif
      KK = (kmaj < 2) ? kb : ((kb+kmaj-1)/kmaj)*kmaj;
      p->mflop[0] = TimeMMKernel(verb, 0, p, pre, mb, nb, kb, 0, 0, -1);
      if (KK != kb)
         p->mflop[0] *= (double)kb / (double)KK;
      printf("   nb=%d, kb=%d,  mu=%d, nu=%d, MFLOP=%.2f\n",
             nb, kb, mu, nu, p->mflop[0]);
      if (gmmb)
      {
         mp->next = p;
         mp = p;
      }
      else
         gmmb = mp = p;
   }
   printf("Done.\n");
/*
 * Now, add the done-list items to generated list
 */
   for (i=0; i < nd; i++)
   {
      const int i4=4*i, mu=dl[i4], nu=dl[i4+1], kmaj=dl[i4+2], nb=dl[i4+3];
      ATL_mmnode_t *prev=NULL;
/*
 *    Get a copy of done-list kern that can be added to genlist
 */
      np = CloneMMNode(dlmm[i]);
      np->next = gmmb;
      gmmb = np;
   }
   if (dlmm)
      free(dlmm);
/*
 * Now, search index file for suitable user-submitted kernels to compete
 * with existing solutions
 */
   ub = ReadMMFileWithPath(pre, "AMMCASES", "amcases.idx");
   ummb = NULL;  /* no suitable user cases to begin */
/*
 * Look through user-list for any routine with ku=1 and K-Runtime
 */
   for (mp=ub; mp; mp = mp->next)
   {
      int km = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
      if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME) &&
          (mp->ku == 1 || km == mp->ku))
      {
/*
 *       It matched our gross criteria, see if it is a required mu/nu
 */
         for (i=0; i < nn; i++)
         {
            const int i4=(i<<2), mu=nl[i4], nu=nl[i4+1], kmaj=nl[i4+2],
                      nb=nl[i4+3];
            if (mp->mu == mu && mp->nu == nu && km == kmaj)
            {
               if (!FailKCleanTests(pre, nb, mp))
               {
                  ATL_mmnode_t *p;
                  p = CloneMMNode(mp);
                  p->next = NULL;
                  p->nbB = ((nb+nu-1)/nu)*nu;
                  p->mbB = ((nb+mu-1)/mu)*mu;
                  p->kbB = nb;
                  if (ummb)
                  {
                     np->next = p;
                     np = p;
                  }
                  else
                     np = ummb = p;
                  break;
               }
            }
         }
      }
   }
   KillAllMMNodes(ub);
/*
 * If we have both user and genned code, must compare timing to select best
 */
   if (ummb)
   {
/*
 *    Now, loop over user cases and time them for comparison with genned
 */
      printf("Timing User K-cleanup:\n");
      for (mp=ummb; mp; mp = mp->next)
      {
         const int nb = mp->nbB, kb = (nb > 8) ? (nb>>2) : nb;
         const int KK = FLAG_IS_SET(mp->flag,MMF_KVEC) ?
                   kb:((kb+mp->vlen-1)/mp->vlen)*mp->vlen;
         mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mp->mbB, nb, KK, 0,0,-1);
         if (KK != kb)
            mp->mflop[0] *= (double)kb / (double)KK;
         printf("   ID=%d, nb=%d, kb=%d, mu=%d, nu=%d, MFLOP=%.2f\n",
                mp->ID, nb, kb, mp->mu, mp->nu, mp->mflop[0]);
      }
      printf("Done timing, merging lists:\n");
/*
 *    Merge generated (gmmb) and user (ummb) kerns by selecting best performing.
 *    gmmb is a superset of ummb, so what we will do is look through gmmb
 *    for matching (mu,nu,dup), time them, and if ummb is faster, replace
 *    that entry in gmmb with ummb.
 */
      while (ummb)
      {
         ATL_mmnode_t *prev=NULL;
         int mu=ummb->mu, nu=ummb->nu;
         int kmaj = FLAG_IS_SET(ummb->flag, MMF_KVEC) ? ummb->vlen:0;

         for (mp=gmmb; mp; mp = mp->next)
         {
            int km = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
            if (mp->mu == mu && mp->nu == nu && km == kmaj)
               break;
            prev = mp;
         }
         assert(mp);  /* logic error if we haven't found it */
/*
 *       If user case gets better performance, replace genned case in queue
 */
         if (ummb->mflop[0] > gmmb->mflop[0])
         {
            printf("   Replacing genned case (%.2f) with user ID %d (%.2f)\n",
                   gmmb->mflop[0], ummb->ID, ummb->mflop[0]);
            if (prev)
            {
               prev->next = ummb;
               ummb = ummb->next;
               prev->next->next = KillMMNode(mp);
            }
            else /* replace gmmb, mp pts at gmmb */
            {
               ATL_mmnode_t *up=ummb;
               ummb = ummb->next;
               up->next = KillMMNode(gmmb);
               gmmb = up;
            }
         }
         else /* user case loser, just delete it */
         {
            printf("   Preferring genned case (%.2f) over user ID %d (%.2f)\n",
                   gmmb->mflop[0], ummb->ID, ummb->mflop[0]);
            ummb = KillMMNode(ummb);
         }
      }
      printf("DONE.\n\n");
   }
   else
      printf("NO VALID USER-SUBMITTED K-CLEANUP KERNELS\n\n");
   free(dl);
   WriteRefreshedMMFileWithPath(pre, "res", "AMMKCLEAN.sum", gmmb);
   return(gmmb);
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


/*
 * Returns a non-repetitive list of user kernels (ID>0) found in rb.  Note that
 * differing compilations of the same kernel are reduced to one entry.
 * rb is left unchanged.
 */
ATL_mmnode_t *GetUniqueUserKerns(ATL_mmnode_t *rb)
{
   ATL_mmnode_t *ub=NULL, *p;

   if (!rb)
      return(NULL);
   for (p=rb; p; p = p->next)
      if (p->ID > 0)
         break;
   if (!p)
      return(NULL);
   ub = CloneMMNode(p);
   for (p=p->next; p; p = p->next)
   {
       if (p->ID > 0)
       {
          ATL_mmnode_t *np;
          int ID = p->ID;

          for (np=ub; np; np = np->next)
             if (np->ID == ID)
                break;
          if (!np)
          {
             np = CloneMMNode(p);
             np->next = ub;
             ub = np;
          }
       }
   }
   return(ub);
}

static int SelfKClean(ATL_mmnode_t *mp)
{
   if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
   {
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      {
         if (mp->ku == mp->vlen && mp->kbmin == mp->vlen)
            return(1);
      }
      else if (mp->ku == 1 && mp->kbmin < 2)
         return(1);
   }
   return(0);
}

ATL_mmnode_t *CanKClean(ATL_mmnode_t *krn, ATL_mmnode_t *cln)
/*
 * RETURNS: NULL if cln/krn cannot provide K-cleanup for krn, else ptr
 *          to krn if it can do its own K-cleanup, else ptr to cln
 */
{
/*
 * First, determine if kernel can perform its own cleaning
 */
   if (SelfKClean(krn))
      return(krn);
/*
 * Cleaner must share same mu/nu, have runtime K, and handle long enough K
 */
   if (krn->mu == cln->mu && krn->nu == cln->nu &&
       FLAG_IS_SET(cln->flag, MMF_KRUNTIME) &&
      (!cln->kbmax || cln->kbmax >= krn->kbB))
   {
      int kmaj = FLAG_IS_SET(krn->flag, MMF_KVEC) ? krn->vlen:0;
      int km = FLAG_IS_SET(cln->flag, MMF_KVEC) ? cln->vlen:0;
      if (kmaj > 1)
      {
         if (km == kmaj && cln->kbmin <= kmaj)
            return(cln);
      }
      else if (cln->ku == 1 && cln->kbmin < 2 && km < 2)
         return(cln);
   }
   return(NULL);
}

ATL_mmnode_t *FindKCleaner(ATL_mmnode_t *clnb, ATL_mmnode_t *kp)
/*
 * RETURNS: kp if kp provides its own K cleanup,
 *          else NULL if no K-cleaner for kp is found in clnb,
 *          else a ptr to the first such valid K-cleaner found in clnb
 */
{
   ATL_mmnode_t *mp, *cln;
   for (mp=clnb; mp; mp = mp->next)
   {
      ATL_mmnode_t *cln;
      cln = CanKClean(kp, mp);
      if (cln)
         return(cln);
   }
   return(NULL);
}

ATL_mmnode_t *FindAllKCleaners(ATL_mmnode_t *clnb, ATL_mmnode_t *kp)
/*
 * RETURNS: if kp provides its own K-cleaning, then kp is returned.
 *          otherwise a queue cloned nodes of all kernels in clnb that
 *          could be used to clean kp is return.
 * Cloned nodes have their blocking values set to match kp
 */
{
   ATL_mmnode_t *gdb=NULL, *mp;

   gdb = FindKCleaner(clnb, kp);
   if (gdb == kp)
      return(kp);
   else if (gdb)
   {
      mp = gdb;
      gdb = CloneMMNode(gdb);
      gdb->mbB = kp->mbB;
      gdb->nbB = kp->nbB;
      gdb->kbB = kp->kbB;
      while ((mp = FindKCleaner(mp->next, kp)))
      {
         ATL_mmnode_t *p;
         p = CloneMMNode(mp);
         p->mbB = kp->mbB;
         p->nbB = kp->nbB;
         p->kbB = kp->kbB;
         p->next = gdb;
         gdb = p;
      }
      return(gdb);
   }
   return(NULL);
}

ATL_mmnode_t *FindAllUniqueKClean(int verb, char pre, ATL_mmnode_t *mmb)
/*
 * Finds a way to clean up all kernels in mmb
 * RETURNS: list of all unique kernels required to do K-cleanup
 */
{
   ATL_mmnode_t *clnb, *mkb, *kp;
   if (verb)
      printf("FINDING K CLEANERS FOR ALL KERNELS:\n");
   clnb = ReadMMFileWithPath(pre, "res", "k1AMM.sum");
   if (clnb)
   {
      MMFillInGenStrings(pre, clnb);
      TimeNegMMKernels(0, verb, 0, clnb, pre, 1, 0, -1);
      WriteRefreshedMMFileWithPath(pre, "res", "k1AMM.sum", clnb);
      return(clnb);
   }
/*
 * mkb is the list of all candidate cleanup codes
 */
   mkb = GetWorkingUserCases(verb, pre);
   for (kp=mmb; kp; kp = kp->next)
   {
      ATL_mmnode_t *cp;
      if (FindKCleaner(clnb, kp))  /* if we've already got a K-cleaner */
         continue;                 /* for this case, skip! */
      cp = FindAllKCleaners(mkb, kp);
      if (cp)
         printf("   %s --> %s!\n", kp->rout, cp->rout);
      else
         printf("   %s --> no Kclean!\n", kp->rout);
/*
 *    For kernels that serve as their own K-cleanup, just use them wt no need
 *    to time anything else
 */
      if (cp == kp)
      {
         cp = CloneMMNode(kp);
         cp->next = clnb;
         clnb = cp;
      }
/*
 *    For kernels that must be cleaned by other kernels, we must time all
 *    candidate kernels and use the best!
 */
      else
      {
         ATL_mmnode_t *mp;
         const int mb=kp->mbB, nb=kp->nbB, kb=kp->kbB;
         int ntim;
/*
 *       Add generated case to any user cases that work
 */
         mp = GetNewKCleanGenNode(pre, kp, mb, nb, kb);
         if (!CanKClean(kp, mp))
         {
            int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
            int km = FLAG_IS_SET(kp->flag, MMF_KVEC) ? kp->vlen:0;
            fprintf(stderr, "KU=(%d,%d,%d), KV=%d, MU=(%d,%d,%d), MV=%d\n",
                    kp->mu, kp->nu, kp->ku, km,
                    mp->mu, mp->nu, mp->ku, kmaj);
            assert(CanKClean(kp,mp));
         }
         mp->next = cp;
         cp = mp;
         ntim = ATL_CountNumberOfMMNodes(cp);
         if (ntim > 1)
         {
/*
 *          Now time all kernels, and choose the fastest for cleanup.
 *          We'll use kbB as kb, even though this is larger than the code
 *          will ever be used for.  However, it will allow us to directly
 *          compare kernel and cleanup performance.
 *          A better strategy would be to time many different K cases, but
 *          I don't want to spend that amount of install time tuning and timing
 *          low-order cleanup!
 */
            printf("   CHOOSING BETWEEN %d KB=%d K-CLEANERS WITH TIMINGS:\n",
                   ntim, kb);
            for (mp=cp; mp; mp = mp->next)
            {
               double mf;
               if (mp->mbB != mb || mp->nbB != nb || mp->kbB != kb ||
                   mp->mflop[0] <= 0.0)
               {
                  mp->mbB = mb;
                  mp->nbB = nb;
                  mp->kbB = kb;
                  mf = TimeMMKernel(0, 1, mp, pre, mb, nb, kb, 0, 0, -1);
                  mp->mflop[0] = mf;
               }
               else
                  mf = mp->mflop[0];
               printf("      %d-%s: %.2f\n", mp->ID, mp->rout?mp->rout:"", mf);
            }
            mp = FindMaxMflopMMQ(cp, 0);
            printf("   USING %d-%s\n", mp->ID, mp->rout? mp->rout:"");
            cp = RemoveMMNodeFromQ(cp, mp);
            KillAllMMNodes(cp);
         }
         else
         {
            mp = cp;
            mp->mbB = mb;
            mp->nbB = nb;
            mp->kbB = kb;
         }
         mp->next = clnb;
         clnb = mp;
      }
   }
   KillAllMMNodes(mkb);
   if (verb)
      printf("DONE FINDING FULL LIST OF K-CLEANERS\n");
   WriteRefreshedMMFileWithPath(pre, "res", "k1AMM.sum", clnb);
   return(clnb);
}

ATL_mmnode_t *FindMUNU(ATL_mmnode_t *mb, int mu, int nu)
{
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
      if (mp->mu == mu && mp->nu == nu)
         return(mp);
   return(NULL);
}

ATL_mmnode_t *KCleanByNB
(
   int verb,
   char pre,
   ATL_mmnode_t *mmb, /* final kernels giving final supported NBs */
   ATL_mmnode_t *mkb  /* All necessary routs ku=1 to clean all kerns in mmb */
)
/*
 * Replicates mkb so that it includes all NBs in mmb, times K-clean,
 * **FREES** mkb, and returns by-NB list
 *
 * OUTPUT:
 *   <pre>AMMKCLEANBYNB.sum: non-unique K-clean for each NB in mmb
 *      mflop[1] contains estimated time for 1 K-it using K=MAX(kb/4,4)
 */
{
   ATL_mmnode_t *nkb=NULL, *mp, *np;
   int kb;
   double mf;

   nkb = ReadMMFileWithPath(pre, "res", "AMMKCLEANBYNB.sum");
   if (nkb)
   {
      KillAllMMNodes(mkb);
      printf("READING IN BY-NB K-CLEANUP:\n");
      MMFillInGenStrings(pre, nkb);
      for (mp=nkb; mp; mp = mp->next)
      {
         int mb=mp->mbB, nb=mp->nbB;
         kb = mp->kbB >> 2;
         kb = (kb >= 4) ? kb : 4;
         if (mp->mflop[0] < 0.0)
            mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mb, nb, kb, 0, 0, -1);
         mf = (2.0*nb)*nb;  /* flop count of gemm/kits (kb) */
         mp->mflop[1] = mf / mp->mflop[0]; /* time in microsecs for 1 k-it */
         printf("   nb=%d, kb=%d, mu=%d, nu=%d, mf=%.2f (%e Usec/Kit)\n",
                nb, kb, mp->mu, mp->nu, mp->mflop[0], mp->mflop[1]);
      }
      printf("Done.\n");
      WriteRefreshedMMFileWithPath(pre, "res", "AMMKCLEANBYNB.sum", nkb);
      return(nkb);
   }
   printf("TIMING K-CLEAN FOR ALL SUPPORTED NBs:\n");
   for (mp=mmb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      int mb = mp->mbB, nb = mp->nbB, kb;

      kb = mp->kbB >> 2;
      kb = (kb >= 4) ? kb : 4;
      if (FLAG_IS_SET(mp->flag,MMF_KVEC))
         kb = ((kb+mp->vlen-1)/mp->vlen)*mp->vlen;

      p = FindMUNU(mkb, mp->mu, mp->nu);
/*
 *    If no user cleanup exists, generate one
 */
      if (!p)
      {
         if (FLAG_IS_SET(mp->flag, MMF_KVEC))
            p = MMGetNodeGEN(pre, 0, 0, mp->mu, mp->nu, mp->vlen, mp->vlen,
                             1, 0, 0, NULL);
         else
            p = MMGetNodeGEN(pre, 0, 0, mp->mu, mp->nu, 1, mp->vlen,
                             0, 0, 0, NULL);
      }
      else
      {
         p = CloneMMNode(p);
         p->next = NULL;
      }
      p->nbB = nb; p->mbB = mb;  p->kbB = kb;
      p->mflop[0] = TimeMMKernel(verb, 0, p, pre, mb, nb, kb, 0, 0, -1);
      mf = mb*nb;
      p->mflop[1] = mf / p->mflop[0];   /* time in microseconds for 1 k it */
      printf("   mb=%d, nb=%d, kb=%d, mu=%d, nu=%d, mf=%.2f (%e Usec/Kit)\n",
             mb, nb, kb, p->mu, p->nu, p->mflop[0], mf);
      if (nkb)
      {
         np->next = p;
         np = p;
      }
      else
         nkb = np = p;
   }
   printf("DONE.\n\n");
   KillAllMMNodes(mkb);
   WriteRefreshedMMFileWithPath(pre, "res", "AMMKCLEANBYNB.sum", nkb);
   return(nkb);
}

void TimeKClean(int verb, char pre, ATL_mmnode_t *mp)
/*
 *   mp is the ku=1, KRUNTIME K-cleanup kernel for a support NB.
 *   This routine creates an output file for supported the NB, where we
 *   document the performance for all NB different KB values.  These
 *   timings can therefore precisely document how expensive K-cleanup
 *   will be for each NB.
 *   OUTPUT:
 *   <pre>AMMKCLEAN_<nb>.TIM: timing of K-clean for nb=<nb>; there are
 *   i=nb-1 timings, mflop[0] contains time to do NB-i K its.  Will use
 *   these times to get completely accurate estimate of total time for
 *   large problems (use estimated time in CLBYNB for small probs).
 */
{
   ATL_mmnode_t *mmb, *p, *np;
   char fn[32];
   int mb = mp->mbB, nb = mp->nbB, i;

   sprintf(fn, "AMMKCLEAN_%d.TIM", mp->nbB);
   mmb = ReadMMFileWithPath(pre, "res", fn);
   if (mmb)
   {
      printf("READING IN K-CLEANUP TIMINGS FOR NB=%d:\n", nb);
      MMFillInGenStrings(pre, mmb);
      for (p=mmb; p; p = p->next)
      {
         int kb = p->kbB;
         assert(nb == p->nbB && mb == p->mbB);
         if (p->mflop < 0)
            p->mflop[0] = TimeMMKernel(verb, 0, p, pre, mb, nb, kb, 0, 0, -1);
         printf("   MB=%d, NB=%d, KB=%d, mu=%d, nu=%d, MFLOP=%.2f\n",
                mb, nb, kb, p->mu, p->nu, p->mflop[0]);
      }
      printf("Done.\n\n");
      WriteRefreshedMMFileWithPath(pre, "res", fn, mmb);
      KillAllMMNodes(mmb);
      return;
   }

   printf("TIMING K-CLEANUP FOR MB=%d, NB=%d:\n", mb, nb);
   for (i=1; i <= nb; i++)  /* create queue of ascending KB */
   {
      p = CloneMMNode(mp);
      p->next = NULL;
      p->kbB = i;
      p->mflop[0] = TimeMMKernel(verb, 0, p, pre, mb, nb, i, 0, 0, -1);
      printf("   MB=%d, NB=%d, KB=%d, mu=%d, nu=%d, MFLOP=%.2f\n", mb, nb, i,
             p->mu, p->nu, p->mflop[0]);
      if (mmb)
      {
         np->next = p;
         np = p;
      }
      else
         mmb = np = p;
   }
   WriteRefreshedMMFileWithPath(pre, "res", fn, mmb);
   KillAllMMNodes(mmb);
   printf("Done.\n");
}

/*
 * Specialize the K cleanup routs in mkb to the kernels in mmb by changing
 * their block factors, and timing them.
 */
ATL_mmnode_t *SpecializeKClean
   (int verb, char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *mkb)
{
   ATL_mmnode_t *mp, *b=NULL;
   for (mp=mmb; mp; mp = mp->next)
   {
      ATL_mmnode_t *kp;
      if (SelfKClean(mp))
      {
         kp = CloneMMNode(mp);
         kp->mflop[2] = 1.0;
      }
      else
      {
         const int mb=mp->mbB, nb=mp->nbB, kb=mp->kbB;
         kp = FindKCleaner(mkb, mp);
         if (!kp)
            fprintf(stderr, "UR(%d,%d,%d), %s: NO KCLEAN",
                    mp->mu, mp->nu, mp->ku, mp->rout);
         assert(kp);
         kp = CloneMMNode(kp);
         if (kp->mbB != mb || kp->nbB != nb || kp->kbB != kb ||
             kp->mflop[0] <= 0.0)
         {
            kp->mbB = mb;
            kp->nbB = nb;
            kp->kbB = kb;
            kp->mflop[0] = TimeMMKernel(verb, 0, kp, pre, mb, nb, kb, 0, 0, -1);
            kp->mflop[2] = kp->mflop[0] / mp->mflop[0];
         }
      }
      printf("   KB=%d KCLEAN SPEEDUP: %.4f\n", kp->kbB, kp->mflop[2]);
      kp->next = b;
      b = kp;
   }
   return(b);
}
void ComputeKClean(int verb, char pre)
/*
 * This kernel finds K-cleanup for all routines present in <pre>geAMMRES.sum
 * and <pre>sqAMMRES.sum.
 *
 * OUTPUT: <pre>AMMKCLEAN.sum: all unique kerns to be compiled
 *         <pre>geAMMKCLEAN.sum:  MFLOP[2] = cleanup slowdown
 *         <pre>sqAMMKCLEAN.sum: MFLOP[2] = cleanup slowdown
 *
 * NOTE: we will time K-cleanup kernels only in BETA=0 case, and peel
 *    the first K-block rather than the last.  This will minimize the C cost,
 *    which is more appreciable for short-K.  We will actually generate all
 *    beta cases, since sometimes you need other betas (in complex, or if
 *    you can't peel first partial for some reason).
 */
{
   ATL_mmnode_t *mkb, *mmGE=NULL, *mmSQ=NULL, *mp, *mmGEk=NULL, *mmSQk=NULL;
   ATL_mmnode_t *ipb, *ipk;
   mmGEk = TimeMMFileWithPath(pre, "res", "geAMMKCLEAN.sum",
                              0, verb|1, 0, 0, 0, -1);
   mmSQk = TimeMMFileWithPath(pre, "res", "sqAMMKCLEAN.sum",
                              0, verb|1, 0, 0, 0, -1);
/*
 * Find all unique K-clean kernels
 */
   mkb = ReadMMFileWithPath(pre, "res", "k1AMM.sum");
   if (!mkb)
   {
/*
 *    Will use only one K-cleanup for any (mu,nu,kmaj) combo.  Decide between
 *    competing kernels based on timings, with larger KB more important, so we
 *    reverse the list order so that the larger block factors choose cleanup for
 *    smaller, rather than reverse.
 */
      mmGE = ReverseMMQ(ReadMMFileWithPath(pre, "res", "geAMMRES.sum"));
      mmSQ = ReverseMMQ(ReadMMFileWithPath(pre, "res", "sqAMMRES.sum"));
      assert(mmGE && mmSQ);
/*
 *    Now temporarily join rect & square lists into one, and get a list of
 *    all cleanup routines that are required.  Routines that provide their
 *    own cleanup will always be used regardless of what is in the list,
 *    and the first such kernel will appear in mkb
 */
      for (mp=mmGE; mp->next; mp = mp->next);
      mp->next = mmSQ;
      mkb = FindAllUniqueKClean(verb, pre, mmGE);
      mp->next = NULL;  /* go back to separate lists */
   }
/*
 * Use reversed lists to build lists of cleanup, which will be in correct
 * order due to the way we build them
 */
   if (!mmGEk)
   {
      if (verb)
         printf("SPECIALIZING K-CLEANERS FOR RECTANGULAR BLOCKINGS:\n");
      if (!mmGE)
         mmGE = ReverseMMQ(ReadMMFileWithPath(pre, "res", "geAMMRES.sum"));
      mmGEk = SpecializeKClean(verb, pre, mmGE, mkb);
      if (verb)
         printf("DONE SPECIALIZING K-CLEANERS FOR RECTANGULAR BLOCKINGS.\n");
      KillAllMMNodes(mmGE);
      WriteRefreshedMMFileWithPath(pre, "res", "geAMMKCLEAN.sum", mmGEk);
   }
   KillAllMMNodes(mmGEk);
/*
 * Now do same for square kernels
 */
   if (!mmSQk)
   {
      if (verb)
         printf("SPECIALIZING K-CLEANERS FOR SQUARE BLOCKINGS:\n");
      if (!mmSQ)
         mmSQ = ReverseMMQ(ReadMMFileWithPath(pre, "res", "sqAMMRES.sum"));
      mmSQk = SpecializeKClean(verb, pre, mmSQ, mkb);
      if (verb)
         printf("DONE SPECIALIZING K-CLEANERS FOR SQUARE BLOCKINGS.\n");
      KillAllMMNodes(mmSQ);
      WriteRefreshedMMFileWithPath(pre, "res", "sqAMMKCLEAN.sum", mmSQk);
   }
   KillAllMMNodes(mmSQk);
/*
 * Specialize cleaners to inner-product views.  We don't have to read this
 * for unique, since they come from geAMMRES exclusively
 */
   ipk = TimeMMFileWithPath(pre, "res", "nAMMKCLEAN.sum",0,verb|1,0,0,0,-1);
   if (!ipk)
   {
      ipb = ReverseMMQ(ReadMMFileWithPath(pre, "res", "ipnPERF.sum"));
      assert(ipb);
      if (verb)
         printf("SPECIALIZING K-CLEANERS FOR IP DEGEN N:\n");
      ipk = SpecializeKClean(verb, pre, ipb, mkb);
      if (verb)
         printf("DONE SPECIALIZING K-CLEANERS FOR IP DEGEN N.\n");
      KillAllMMNodes(ipb);
      WriteRefreshedMMFileWithPath(pre, "res", "nAMMKCLEAN.sum", ipk);
   }
   KillAllMMNodes(ipk);
   ipk = TimeMMFileWithPath(pre, "res", "mAMMKCLEAN.sum",0,verb|1,0,0,0,-1);
   if (!ipk)
   {
      ipb = ReverseMMQ(ReadMMFileWithPath(pre, "res", "ipmPERF.sum"));
      assert(ipb);
      if (verb)
         printf("SPECIALIZING K-CLEANERS FOR IP DEGEN M:\n");
      ipk = SpecializeKClean(verb, pre, ipb, mkb);
      if (verb)
         printf("DONE SPECIALIZING K-CLEANERS FOR IP DEGEN M.\n");
      KillAllMMNodes(ipb);
      WriteRefreshedMMFileWithPath(pre, "res", "mAMMKCLEAN.sum", ipk);
   }
   KillAllMMNodes(ipk);
   ipk = TimeMMFileWithPath(pre, "res", "mnAMMKCLEAN.sum",0,verb|1,0,0,0,-1);
   if (!ipk)
   {
      ipb = ReverseMMQ(ReadMMFileWithPath(pre, "res", "ipmnPERF.sum"));
      assert(ipb);
      if (verb)
         printf("SPECIALIZING K-CLEANERS FOR IP DEGEN MN:\n");
      ipk = SpecializeKClean(verb, pre, ipb, mkb);
      if (verb)
         printf("DONE SPECIALIZING K-CLEANERS FOR IP DEGEN MN.\n");
      KillAllMMNodes(ipb);
      WriteRefreshedMMFileWithPath(pre, "res", "mnAMMKCLEAN.sum", ipk);
   }
   KillAllMMNodes(ipk);

   KillAllMMNodes(mkb);
   if (verb)
      printf("\nDONE FINDING K-CLEANUP FOR EACH KB.\n");
}

void FindBestKU1
(
   int verb,
   char pre,
   int K       /* K dim, should be small, probably like 23 or 17 */
)
/*
 * Find the best possible kernel for use in low-rank update;  We only consider
 * kernels with runtime-K that handle all possible K (ku=1).  We will try
 * all legal blocking factors between 16 & 480 for this kernel, and choose
 * the one that performs best.  This kernel always used for any K not covered
 * by optimized kernels given in eAMMRES kbBs.  When we match a kbB, we
 * compare the perf of this kernel at its optimal nbB/mbB wt that of the
 * specialized kernel, and choose the best.
 *
 * OUTPUT: This routine outputs two files:
 * (1) AMMRANKK: best ku=1 kern wt best MB/NB, K=K
 * (2) AMMRANKKT: timing of this kern wt M=mbB, N=nbB, all K between 1 & maxNB
 */
{
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
      printf("   ID=%d, mb=%d, nb=%d, kb=%d, RTK=%d, MFLOP=%.2f\n", p->ID,
             (int)mb, (int)nb, (int)kb, FLAG_IS_SET(p->flag, MMF_KRUNTIME), mf);
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
         printf("%s: ID=%d, RT='%s' B=(%d,%d,%d), RK=%d MFLOP=%.2f\n",
                exp[i], mp->ID, mp->rout, mp->mbB, mp->nbB, mp->kbB,
                FLAG_IS_SET(mp->flag, MMF_KRUNTIME), mp->mflop[i+1]);
      printf("\n");
   }
   return(mm3);
}
ATL_mmnode_t *TimeKBRegion
(
   int verb,
   char pre,
   ATL_mmnode_t *mmk,            /* kernel to time throughout region */
   int kbmin,                    /* start of region */
   int kbend,                    /* largest kb in region */
   int kincD                     /* default stride between kernel timings */
)
/*
 * Returns list of timings of kernel mmk using near-square cases with KB
 * varying between kbmin - kbend.  All cases that are legal and incremented
 * by kinc are tried, as are all perfectly square cases
 */
{
   ATL_mmnode_t *mmb=NULL, *mp, *mpB=NULL;
   const int ku = mmk->ku, mu=mmk->mu, nu=mmk->nu;
   int kstart, kinc, kend, k, ksq, ksqinc;
   double mf, mfB=0.0;
/*
 * Get starting and ending point that is legal for this kernel.
 */
   kstart = Mmax(mmk->kbmin, kbmin);
   kstart = ((kstart+ku-1)/ku)*ku;
   kend = ((kbend+ku-1)/ku)*ku;
   if (mmk->kbmax)
      kend = Mmin(kend, mmk->kbmax);
   k = kstart;
/*
 * square inc always lcm(mu,nu,ku).  Normal increment is always at least
 * as big as the default stride, but must be a multiple of the kernel's ku
 */
   ksqinc = Mylcm(mu, nu);
   ksqinc = Mylcm(ksqinc, ku);
   for (kinc=ku; kinc < kincD; kinc += ku);
   if (kstart <= kend)
   {
      int kb = k;
      printf("TIMING %s mu=%d, mu=%d, For KB=[%d,%d]:\n",
             mmk->rout ? mmk->rout : "Genkern", mu, nu, kstart, kend);
      ksq = ((kstart+ksqinc-1)/ksqinc)*ksqinc;
      do
      {
         ATL_mmnode_t *p;
         const int mb=((kb+mu-1)/mu)*mu, nb=((kb+nu-1)/nu)*nu;

         p = CloneMMNode(mmk);
         mf = TimeMMKernel_KB(verb, 0, p, pre, mb, nb, kb, 1, 0, -1);
         printf("   mb=%d, nb=%d, kb=%d, KRUN=%d, MFLOP=%.2f\n",
                mb, nb, kb, FLAG_IS_SET(p->flag, MMF_KRUNTIME), mf);
         p->mflop[0] = mf;
         if (mf > mfB)
         {
            mfB = mf;
            mpB = p;
         }
         if (mmb)
         {
            mp->next = p;
            mp = p;
         }
         else
            mmb = mp = p;
         if (kb == ksq)
            ksq += ksqinc;
         if (kb == k)
            k += kinc;
         kb = Mmin(k,ksq);
      }
      while(kb <= kend);
      printf("DONE, best case %s mb=%d, nb=%d, kb=%d, MFLOP=%.2f\n\n",
             mmk->rout ? mmk->rout : "Genkern",
             mpB->mbB, mpB->nbB, mpB->kbB, mfB);
   }
   else
   {
      printf("KERNEL %s mu=%d, mu=%d, has no legal cases in KB=[%d,%d]!\n\n",
             mmk->rout ? mmk->rout : "Genkern", mu, nu, kstart, kend);
   }
   return(mmb);
}

ATL_mmnode_t *TimeAllKBRegions
(
   int verb,
   char pre,
   ATL_mmnode_t *mmk,            /* kernel to time throughout regions */
   int kb                        /* maxKB to ever try */
)
/*
 * Times mmk for KBs up to kb
 */
{
   const int ku = mmk->ku;
   ATL_mmnode_t *mmb, *mp;

   mmb = TimeKBRegion(verb, pre, mmk, 24, kb, 4);
   return(mmb);
}

ATL_mmnode_t *FindCacheBudgetCasesByKB
(
   int verb,
   char pre,
   size_t CS,                    /* size of cache we are optimizing for */
   ATL_mmnode_t *mmb             /* list of cases to try */
)
/*
 * This routine is responsible for:
 * (1) Find the best performing kernels out of mmb for our 3-5 cache budget
 *     cases
 * (2) Free mmb
 * (3) For each unique kernel, find perf of kernel for all supported KBs
 *     in the budgetary regions
 * (4) Merge these lists, and winnow underperforming cases
 * (5) RETURN: queue of all supported KBs
 */
{
   ATL_mmnode_t *mm3b, *mp;
   const char upr = (pre == 'z' || pre == 'd') ? 'd' : 's';
/*
 * We want number of real elts, not # of cplx elts!
 */
   if (pre == 'c' || pre == 'z')
      CS += CS;
/*
 * See if we just need to rerun cases
 */
   mm3b = ReadMMFileWithPath(pre, "res", "bAMMRES.sum");
   if (mm3b)
   {
      int i=0, WRT=0;
      printf("READING IN LARGE KERNEL CASES FROM res/<pre>bAMMRES:\n");
      MMFillInGenStrings(pre, mm3b);
      for (mp=mm3b; mp; mp = mp->next)
      {
         if (mp->mflop[0] <= 0.0)
         {
            mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mp->mbB, mp->nbB,
                                        mp->kbB, 1, 0, -1);
            WRT=1;
         }
         printf("   ID=%d, %s: MB=%d, NB=%d, KB=%d, KRUN=%d, MFLOP=%.2f\n",
                mp->ID, mp->rout ? mp->rout : "Gennedkern",
                mp->mbB, mp->nbB, mp->kbB, FLAG_IS_SET(mp->flag, MMF_KRUNTIME),
                mp->mflop[0]);
         i++;
      }
      if (WRT)
         WriteRefreshedMMFileWithPath(pre, "res", "bAMMRES.sum", mm3b);
      printf("DONE %d CASES.\n\n", i);
      return(mm3b);
   }
/*
 * Find best performing kernels for each of our 3-5 cache budgets
 */
   mm3b = mmb;
/*
 * Get list of performance of all-3-5 in-cache kernel in all 3-5 cache regions
 */
   mmb = TimeAllKBRegions(verb, pre, mm3b, mm3b->next->next->kbB);
   for (mp=mm3b->next; mp; mp = mp->next)
   {
      if (KernelIsUnique(mm3b, mp))
      {
         ATL_mmnode_t *p, *p2;
         p = TimeAllKBRegions(verb, pre, mp, mm3b->next->next->kbB);
         p2 = MergeCases(0, mmb, p);
         KillAllMMNodes(p);
         KillAllMMNodes(mmb);
         mmb = p2;
      }
   }
   KillAllMMNodes(mm3b);
/*
 * Now, get rid of any blocking factor that is slower than the preceeding one
 */
   mmb = WinnowCases(0, mmb);
   WriteRefreshedMMFileWithPath(pre, "res", "bAMMRES.sum", mmb);
   return(mmb);
}


/*
 * This routine finds kernels to use in low-rank-K update. For 3 <= K <= 15,
 * it tries all kernels and chooses the best performing; In this search
 * we consider only compile-time K kernels, since runtime kernels will be
 * selected by general (K>15) search.
 * (K=1 and K=2 are handled by GER and GER2).
 * We consider only user-generated kernels; for any problem sizes that are
 * not supported, we will use the normal K or K-clean routines.  This list
 * is just to allow for hand-tuning small-K special cases.
 * This routine produces output file <pre>AMMLOWK.sum
 */
ATL_mmnode_t *GetLowRankKKernels
(
   int verb,            /* verbosity */
   char pre,            /* s,d */
   int MB,              /* default mb to time with */
   int NB,              /* default nb to time with */
   ATL_mmnode_t *inb    /* all working ukerns */
)
{
   int k, ik;
   ATL_mmnode_t *rkKb=NULL, *mp;
/*
 * Get rid of all K-runtime kernels from consideration
 */
   while (inb && FLAG_IS_SET(inb->flag, MMF_KRUNTIME))
      inb = KillMMNode(inb);
   if (inb)
   {
      ATL_mmnode_t *prev=inb;
      mp = inb->next;
      while (mp)
      {
         if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
         {
            mp = KillMMNode(mp);
            prev->next = mp;
         }
         else
         {
            prev = mp;
            mp = mp->next;
         }
      }
   }
   else
      return(NULL);
   for (ik=0; ik < 2; ik++)
   {
      int kbeg, kend, kinc;
      if (!ik)
      {
         kbeg = 96;
         kend = 16;
         kinc = 16;
      }
      else
      {
         kbeg = 15;
         kend = 3;
         kinc = 1;
      }
      for (k=kbeg; k >= kend; k -= kinc)
      {
         printf("FINDING BEST USER-PROVIDED KERNEL FOR K=%d:\n", k);
         ATL_mmnode_t *best=NULL;
         for (mp=inb; mp; mp = mp->next)
         {
            const int mu = mp->mu, nu = mp->nu, ku = mp->ku;
            const int mb = (MB/mu)*mu, nb = (NB/nu)*nu;
            const int KK = (!FLAG_IS_SET(mp->flag,MMF_KVEC)) ?
                           k : ((k+ku-1)/ku)*ku;
            double mf;

            assert(mb && nb);
            if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME) || KK%ku)
            {
               printf("   skipping %d. %s, KRUN=%d, ku=%d\n", mp->ID, mp->rout,
                      FLAG_IS_SET(mp->flag, MMF_KRUNTIME), ku);
               continue;
            }
            if ((mp->kbmin && k < mp->kbmin) || (mp->kbmax && KK > mp->kbmax))
            {
               printf("   skipping %d. %s, kbmin,max=%d,%d, K=%d\n",
                      mp->ID, mp->rout, mp->kbmin, mp->kbmax, KK);
               continue;
            }
            mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, KK, 0, 0, -1);
            if (KK != k)
               mf = (mf*k) / (double)KK;
            printf("   %d. %s: mb=%d, nb=%d, MFLOP=%.2f\n", mp->ID, mp->rout,
                   mb, nb, mf);
            if (!best)
            {
               best = mp;
               mp->mflop[0] = mf;
               mp->mbB = mb;
               mp->nbB = nb;
               mp->kbB = k;
            }
            else if (best->mflop[0] < mf)
            {
               best = mp;
               mp->mflop[0] = mf;
               mp->mbB = mb;
               mp->nbB = nb;
               mp->kbB = k;
            }
         }
         if (best)
         {
            best = CloneMMNode(best);
            best->next = rkKb;
            rkKb = best;
            printf("BEST FIXED-%d KERNEL: %d. %s MFLOP=%.2f\n\n",
                   k, best->ID, best->rout, best->mflop[0]);
         }
         else
            printf("NO SPECIAL CASE for K=%d\n\n", k);
      }
   }
   KillAllMMNodes(inb);
   return(rkKb);
}

/*
 * Finds list of best run-time kernel, ranked by KU.  Higher KUs are not
 * retained unless they beat any lower-KU kernel that divides that KU evenly
 */
ATL_mmnode_t *GetRuntimeKKernels
(
   int verb,            /* verbosity */
   char pre,            /* s,d */
   int MB,              /* default mb to time with */
   int NB,              /* default nb to time with */
   ATL_mmnode_t *inb    /* all working ukerns */
)
{
   ATL_mmnode_t *mp;
   int KU;
/*
 * Get rid of all non-K-runtime kernels from consideration
 */
   while (inb && !FLAG_IS_SET(inb->flag, MMF_KRUNTIME))
      inb = KillMMNode(inb);
/*
 * Got rid of any at base, now get rid of non-K-runtime from internal nodes
 */
   if (inb)
   {
      ATL_mmnode_t *prev=inb;
      mp = inb->next;
      while (mp)
      {
         if (!FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
         {
            mp = KillMMNode(mp);
            prev->next = mp;
         }
         else
         {
            prev = mp;
            mp = mp->next;
         }
      }
   }
   else
      return(NULL);
   KU = inb->ku;
   for (mp=inb->next; mp; mp = mp->next)
      KU = Mylcm(KU, mp->ku);
   if (KU > 32)
      KU = 32;
   else
      KU = ((16+KU-1)/KU)*KU;
   printf("TRYING ALL RUNTIMEK KERNS WITH MB=%d, NB=%d, KB=%d:\n", MB, NB, KU);
   for (mp=inb; mp; mp = mp->next)
   {
      int mu=mp->mu, nu=mp->nu, ku=mp->ku;
      int mb = (MB/mu)*mu, nb = (NB/nu)*nu, kb = (KU/ku)*ku;
      double mf;
      assert(mb && nb && kb);
      mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, kb, 0, 0, -1);
      printf("   %d. %s: mb=%d, nb=%d, MFLOP=%.2f\n", mp->ID, mp->rout,
             mb, nb, mf);
      mp->mflop[0] = mf;
      mp->mbB = mb;
      mp->nbB = nb;
      mp->kbB = kb;
   }
   printf("\n");
   inb = ATL_SortMMNodesByMflop(0, inb);
   if (inb->ku == 1)
   {
      KillAllMMNodes(inb->next);
      inb->next = NULL;
   }
   else
   {
      ATL_mmnode_t *p;
/*
 *    Go thru sorted list, and kill all slower nodes that don't add new K
 */
      for (p=inb; p; p = p->next)
      {
         ATL_mmnode_t *prev=p;
         mp = p->next;
         while (mp)
         {
            if (mp->ku % p->ku == 0)
            {
               mp = KillMMNode(mp);
               prev->next = mp;
            }
            else
            {
               prev = mp;
               mp = mp->next;
            }
         }
      }
   }
   if (!inb)
      printf("NO RETAINED RUNTIME KERNELS.\n\n");
   else
   {
      printf("RETAINED RUNTIME KERNELS:\n");
      for (mp=inb; mp; mp = mp->next)
         printf("   %d. %s: ku=%d, MFLOP=%.2f\n", mp->ID, mp->rout, mp->ku,
                mp->mflop[0]);
      printf("DONE.\n");
   }
   return(inb);
}
/*
 * RETURNS: 1 if mmc is slower than any kernel in mmb
 */
int IsSlowerThanList
(
   int verb,            /* verbosity */
   char pre,            /* s,d */
   int MB,
   int NB,              /* default mb/nb to time with */
   ATL_mmnode_t *mmc,  /* candidate mmkern */
   ATL_mmnode_t *mmb   /* kernels to time candidate against */
)
{
   ATL_mmnode_t *mp;
   double mfc, mf;
   int mu, nu, ku;
   int mb, nb, kb, KB;

   if (!mmb)
      return(0);
   kb = mmc->kbB;
   mu = mmc->mu;
   nu = mmc->nu;
   mb = (MB/mu)*mu;
   nb = (NB/nu)*nu;
   assert(mb && nb && kb);
   KB = (!FLAG_IS_SET(mmc->flag,MMF_KVEC))?kb:((kb+mmc->ku-1)/mmc->ku)*mmc->ku;
   mfc = TimeMMKernel(verb, 0, mmc, pre, mb, nb, KB, 0, 0, -1);
   mfc = (kb*mfc)/(double)KB;
   mmc->mflop[1] = mfc;
   kb = mmc->kbB;
   for (mp=mmb; mp; mp = mp->next)
   {
      ku = mp->ku;
      if (mp->kbmin && kb < mp->kbmin)
         continue;
      if (mp->kbmax && kb > mp->kbmax)
         continue;
      if (kb%ku == 0 || FLAG_IS_SET(mp->flag,MMF_KVEC))
      {
         int KK = (FLAG_IS_SET(mp->flag,MMF_KVEC)) ? ((kb+ku-1)/ku)*ku : kb;
         mu = mp->mu;
         nu = mp->nu;
         mb = (MB/mu)*mu;
         nb = (NB/nu)*nu;
         assert(mb && nb);
         mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, KK, 0, 0, -1);
         mp->mflop[1] = (kb*mfc)/(double)KK;
         if (mf > mfc)
         {
            return(1);
         }
      }
   }
   return(0);
}

/*
 * Finds best-performing square cases in list of pre-existing kernels, mmb.
 * Does not modify original mmb list, and will return only kernels that are
 * faster than smaller square cases.
 * RETURNS: new list of all square cases that got best performance from
 *          original mmb.
 */
/*
 * Sets all MV[A,B,C] bits in mmb to those provided in low 3 bits of bits
 */
void ResetMoveBitsInQ(ATL_mmnode_t *mmb, int bits)
{
   while (mmb)
   {
      ATL_MMF_MVPUT(mmb->flag, bits);
      mmb = mmb->next;
   }
}
ATL_mmnode_t *FindBestSquareCases(char pre, int verb, int nregs, int maxNB,
                                  ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mp, *mmSQ=NULL, *mmGN, *prev=NULL, *mmS;
   int maxU=0, i, KM;
   int vlen;
   const char upr = (pre == 's' || pre == 'c') ? 's' : 'd';

   mmGN = GetGenCases(pre);
   assert(mmGN);
   vlen = mmGN->next->next->vlen;
   if (vlen < 1)
      vlen = 1;
   for (mp=mmb; mp; mp = mp->next)
   {
      if (mp->mu > maxU)
         maxU = mp->mu;
      if (mp->nu > maxU)
         maxU = mp->nu;
   }
   maxNB = 1.2 * maxNB;
   maxNB = ((maxNB+maxU-1) / maxU)*maxU;
/*
 * For small problems, try full generated & user kerns, since square kerns
 * are very likely to use some odd kernel that got winnowed in original file
 */
   mmS = ReadMMFileWithPath(upr, "res", "WORKING.sum");
   ATL_LastMMNode(mmS)->next = CloneMMQueue(mmGN);
   ResetMoveBitsInQ(mmS, 5);  /* mmb already set to movA=movC=1, movB=0 */
/*
 * Find best NB=4 case to start queue
 */
   mmSQ = prev = BestForThisNB(verb, pre, mmS, 4, 4-1, 4+1);
/*
 * for nb < vlen, scope every multiple of 2
 */
   KM = Mmax(vlen,16);
   for (i=6; i <= KM; i += 2)
   {
      mp = BestForThisNB(verb, pre, mmS, i, i-1, i+1);
      if (prev->mflop[0] >= mp->mflop[0])
         KillMMNode(mp);
      else
      {
         prev->next = mp;
         prev = mp;
      }
   }
   KillAllMMNodes(mmS);  /* don't need small queue anymore */
/*
 * For larger problems, only use user kernels that won against generated.
 * We add back in all generated kerns in case we need an eliminated one for
 * cleanup.  Could do same for user (use mmS for all), but mmGN is at most
 * length 5, while mmb is of any length, which is why we don't want to time
 * all those kerns.
 */
   mmb = AddUniqueMMKernCompList(mmGN, mmb);
/*
 * For square cases > VLEN, try only multiples of VLEN: since we always
 * vectorize along at least one dim, only these dims can be fully vectorized
 */
   for (i=((16+vlen)/vlen)*vlen; i <= maxNB; i += vlen)
   {
      mp = BestForThisNB(verb, pre, mmb, i, i-1, i+1);
      if (prev->mflop[0] >= mp->mflop[0])
         KillMMNode(mp);
      else
      {
         prev->next = mp;
         prev = mp;
      }
   }
   return(mmSQ);
}

ATL_mmnode_t *MergeRankKKernels
(
   int verb,            /* verbosity */
   char pre,            /* s,d */
   int MB,              /* default mb to time with */
   int NB,              /* default nb to time with */
   int maxKB,           /* largest KB to produce */
   ATL_mmnode_t *fixb,  /* rank-K fixed-K kerenls */
   ATL_mmnode_t *runb,  /* rank-K, runtime-K kernels */
   ATL_mmnode_t *sqrb   /* optimized near-square kernels */
)
{
   ATL_mmnode_t *rkb, *rkp;
   int k;
   rkp = rkb = GetMMNode();
   printf("CHOOSING BEST KERNEL FOR EACH RANK-K (3 <= K <= %d):\n", maxKB);
   for (k=3; k <= maxKB; k++)
   {
      ATL_mmnode_t *best=NULL, *p;
      double mfB=0.0, mf;
/*
 *    fixb & sqrb are in K-order, so we pop them off stack until we get to
 *    one big enough to solve the problem.  We also ignore all KRUNTIME kernels
 *    in sqrb, since they should appear in runb if they are competitive
 */
      while (fixb)
      {
         if (fixb->kbB < k || (fixb->kbmin && fixb->kbmin > k) ||
             (fixb->kbmax && fixb->kbmax < k))
            fixb = KillMMNode(fixb);
         else
            break;
      }
      while (sqrb)
      {
         if (sqrb->kbB < k || FLAG_IS_SET(sqrb->flag, MMF_KRUNTIME)
             || (sqrb->kbmin && sqrb->kbmin > k) ||
                (sqrb->kbmax && sqrb->kbmax < k))
            sqrb = KillMMNode(sqrb);
         else break;
      }
      if (fixb)
      {
         if (fixb->kbB == k)
         {
            int mu = fixb->mu, nu = fixb->nu, ku = fixb->ku;
            int mb = (NB/mu)*mu, nb = (NB/nu)*nu;
            int kb = (!FLAG_IS_SET(fixb->flag,MMF_KVEC)) ? k : ((k+ku-1)/ku)*ku;
            best = fixb;
            fixb = fixb->next;
            mfB = TimeMMKernel(verb, 0, best, pre, mb, nb, kb, 0,0, -1);
            mfB = (mfB*k)/(double)kb;
         }
      }
      if (sqrb)
      {
         if (sqrb->kbB == k)
         {
            int mu = sqrb->mu, nu = sqrb->nu, ku = sqrb->ku;
            int mb = (MB/mu)*mu, nb = (NB/nu)*nu;
            int kb = (!FLAG_IS_SET(sqrb->flag,MMF_KVEC)) ? k : ((k+ku-1)/ku)*ku;
            mf = TimeMMKernel(verb, 0, sqrb, pre, mb, nb, kb, 0, 0,-1);
            mf = (mf*k)/(double)kb;
            if (mf > mfB)
            {
               mfB = mf;
               if (best)
                  KillMMNode(best);
               best = sqrb;
               sqrb = sqrb->next;
            }
            else
               sqrb = KillMMNode(sqrb);
         }
      }
      for (p=runb; p; p = p->next)
         if ((k % p->ku == 0 || FLAG_IS_SET(p->flag,MMF_KVEC)) &&
             k > p->kbmax && (!p->kbmin || p->kbmin <= k))
            break;
      if (p)
      {
         int mu = p->mu, nu = p->nu, ku = p->ku;
         int mb = (NB/mu)*mu, nb = (NB/nu)*nu;
         int kb = (!FLAG_IS_SET(p->flag,MMF_KVEC)) ? k : ((k+ku-1)/ku)*ku;
         if (p->kbmax && p->kbmax < kb)
            mf = -1.0;
         else
            mf = TimeMMKernel(verb, 0, p, pre, mb, nb, kb, 0, 0, -1);
         mf = (mf*k)/(double)kb;
         if (mf > mfB)
         {
            mfB = mf;
            if (best)
               KillMMNode(best);
            best = CloneMMNode(p);
         }
      }
      assert(best);
      printf("   Best kernel K=%d: %d. %s (%.2f)\n",k,best->ID,best->rout,mfB);
      best->kbB = k;
      rkp->next = best;
      rkp = best;
   }
   if (sqrb)
      KillAllMMNodes(sqrb);
   if (fixb)
      KillAllMMNodes(fixb);
   rkp->next = NULL;
   printf("DONE.\n\n");
   return(KillMMNode(rkb));
}

/*
 * Complex types use the previously selected real kernels in order to
 * reduce library size (means we only have 2 precision GEMMS for 4
 * types/precisions).  The only thing that is different is we may
 * reduce max NB in order to keep complex ops in cache.
 * May want to write this as part of real tuning!
 */
int DoComplex(char pre, int verb)
{
   ATL_mmnode_t *mmb;
   char upr = (pre == 'z') ? 'd' : 's';
   exit(-1);
   mmb = ReadMMFileWithPath(pre, "res", "geAMMRES.sum");
   if (!mmb)
   {
      mmb = ReadMMFileWithPath(upr, "res", "geAMMRES.sum");
      assert(mmb);
      WriteMMFileWithPath(pre, "res", "geAMMRES.sum", mmb);
   }

}

/*
 * Creates two lists from original, which is left unchanged.
 * A list of all unique runtime kernels is RETURNED,
 * while a list of all unique KB compile kernels is provided by CBAS
 */
ATL_mmnode_t *SplitRunCompKB(ATL_mmnode_t *orig, ATL_mmnode_t **CBAS)
{
   ATL_mmnode_t *comp=NULL, *run, *bp;
   for (bp=orig; bp; bp = bp->next)
   {
      if (FLAG_IS_SET(bp->flag, MMF_KRUNTIME))
      {
         if (!MMKernCompIsPresent(run, bp))
         {
            ATL_mmnode_t *tp;
            tp = CloneMMQueue(bp);
            tp->next = run;
            run = tp;
         }
      }
      else  /* candidate for compile-time K list */
      {
         if (!MMKernCompIsPresent(comp, bp))
         {
            ATL_mmnode_t *tp;
            tp = CloneMMQueue(bp);
            tp->next = comp;
            comp = tp;
         }
      }
   }
/*
 * Put them back in original order, they were produced backwards from orig
 */
   if (comp)
      comp = ReverseMMQ(comp);
   if (run)
      run = ReverseMMQ(run);

   *CBAS = comp;
   return(run);
}
/*
 * Deletes all kernels in mmb with kbB that are not a multiple of mu
 */
ATL_mmnode_t *KillIncompatible_MK(ATL_mmnode_t *mmb)
{
   while (mmb && (mmb->kbB/mmb->mu)*mmb->mu != mmb->kbB)
      mmb = KillMMNode(mmb);
   if (mmb)
   {
      ATL_mmnode_t *mp=mmb->next, *prev=mmb;

      while (mp)
      {
         if ((mp->kbB/mp->mu)*mp->mu != mp->kbB)
             mp = prev->next = KillMMNode(mp);
         else
         {
            prev = mp;
            mp = mp->next;
         }
      }
   }
   return(mmb);
}
/*
 * Given a list of kernels being already being used by GEMM, find best
 * performing cases with MB=KB, and NB allowed to vary.  These blockings
 * of already-existing kernels will be used to build triangular & symmetric
 * amm-based routines
 */
ATL_mmnode_t *DoMKB_findNB(char pre, int verb, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mmSQ, *run, *comp, *mp;
   mmSQ = ReadMMFileWithPath(pre, "res", "mkbAMMRES.sum");
   if (mmSQ)
   {
      MMFillInGenStrings(pre, mmSQ);
      if (mmSQ->mflop[0] < 0.0)
      {
         for (mp=mmSQ; mp; mp = mp->next)
         {
            int nb=mp->kbB;
            if (mp->mflop[0] < 0.0)
               mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, nb, nb, nb,
                                           1, 0, -1);
         }
         WriteMMFileWithPath(pre, "res", "mkbAMMRES.sum", mmSQ);
      }
      return(mmSQ);
   }
   run = SplitRunCompKB(mmb, &comp);
   comp = KillIncompatible_MK(comp);
/*
 * Add all K-cleanup kernels to list of candidates
 */
   mp = ReadMMFileWithPath(pre, "res", "AMMKCLEAN.sum");
   if (mp)
   {
      ATL_mmnode_t *r, *c;
      r = SplitRunCompKB(mp, &c);
      KillAllMMNodes(mp);
      c = KillIncompatible_MK(c);
      mp = MergeCases(-1, run, r);
      KillAllMMNodes(run);
      KillAllMMNodes(r);
      run = mp;
      mp = MergeCases(-1, comp, c);
      KillAllMMNodes(comp);
      KillAllMMNodes(c);
      comp = mp;
   }
/*
 * HERE HERE HERE
 */
}
ATL_mmnode_t *DoSquare(char pre, int verb, int nregs, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mmSQ;
   const char upr = (pre == 'z' || pre == 'd') ? 'd' : 's';

   mmSQ = TimeMMFileWithPath(pre, "res", "sqAMMRES.sum", 0,verb|1, 0, 1, 0, -1);
   if (mmSQ)
      return(mmSQ);
   else
   {
      int maxNB=0;
      ATL_mmnode_t *mp;

      for (mp=mmb; mp; mp = mp->next) /* find maximum NB */
      {
          if (mp->kbB > maxNB)
             maxNB = mp->kbB;
          if (mp->mbB > maxNB)
             maxNB = mp->mbB;
          if (mp->nbB > maxNB)
             maxNB = mp->nbB;
      }
      mmSQ = FindBestSquareCases(pre, verb, nregs, maxNB, mmb);
      WriteRefreshedMMFileWithPath(pre, "res", "sqAMMRES.sum", mmSQ);
   }
   return(mmSQ);
}

void DoTRSM(char pre, int verb)
/*
 * Later, may need to put trsm 'Right' timings in mflop[1]
 */
{
   ATL_mmnode_t *tb, *tp, *sb;
   if (pre == 'z' || pre == 'c')
      return;
   tb = ReadMMFileWithPath(pre, "res", "tsAMMRES.sum");
   if (tb)  /* already ran! */
   {
      for (tp=tb; tp; tp = tp->next)
      {
         if (tp->mflop[0] <= 0.0)
         {
            tp->mflop[0] = TimeTSKernel(verb, 0, pre, tp->mbB, tp->nbB, 0);
            if (verb)
               printf("   trsmKL %dx%d: %.2f\n", tp->mbB,tp->nbB, tp->mflop[0]);
         }
      }
      WriteRefreshedMMFileWithPath(pre, "res", "tsAMMRES.sum", tb);
      KillAllMMNodes(tb);
      return;
   }

   tb = ReadMMFileWithPath(pre, "res", "sqAMMRES.sum");
   printf("\nFINDING PERFORMANCE OF TRSM KERNELS FOR SQUARE AMM:\n");
   for (tp=tb; tp; tp = tp->next)
   {
      tp->mflop[0] = TimeTSKernel(verb, 0, pre, tp->mbB, tp->nbB, 0);
      if (verb)
         printf("   trsmKL %dx%d: %.2f\n", tp->mbB, tp->nbB, tp->mflop[0]);
   }
   WriteRefreshedMMFileWithPath(pre, "res", "tsAMMRES.sum", tb);
   printf("DONE TRSM TIMING.\n");
   KillAllMMNodes(tb);
}

ATL_mmnode_t *FindBestNK_M
   (char pre, int verb, ATL_mmnode_t *mmb, unsigned int mb)
/*
 * RETURNS: Clone of node that provides best perf for a problem with MB=mb.
 *          We try finding best nb&kb by using ones near any that are near
 *          those in list of mmb.  NULL if no kern can do MB=mb.
 */
{
   ATL_mmnode_t *mp, *mpMax=NULL;
   unsigned int kbB=0, nbB=0;
   double mf, mfMax=0.0;
/*
 * Try all candidate kernels for this mb
 */
   for (mp=mmb; mp; mp = mp->next)
   {
      int kb=mp->kbB, nb=mp->nbB;
      if (mb % mp->mu)
         continue;
      mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, kb, 1, 0, -1);
/*
 *    Try NB (and if allowed KB) from larger kerns to see if they improve
 *    small mb perf
 */
      if (mp->next)  /* larger KB/NB kerns exist */
      {
         ATL_mmnode_t *mpG; /* kerns of greater size */
         int KB, lastK=kb, lastN=nb;
         const int ku=mp->ku, nu=mp->nu;
         const int CHGK=FLAG_IS_SET(mp->flag, MMF_KRUNTIME);

         for (mpG=mp->next; mpG; mpG = mpG->next)
         {
            int k = mpG->kbB, n=mpG->nbB;
            double mfG;

            if (nu <= 8)
               n = ((n+nu-1)/nu)*nu;
            else
               n = (n/nu)*nu;
            if (CHGK)
            {
               if (ku <= 8)
                  k = ((k+ku-1)/ku)*ku;
               else
                  k = (k/ku)*ku;
            }
            else
               k = kb;
            if (n <= lastN || k <= lastK)
               continue;
            mfG = TimeMMKernel(verb, 0, mp, pre, mb, n, k, 1, 0, -1);
            if (mfG > mf)
            {
               mf = mfG;
               nb = n;
               kb = k;
            }
            lastN = n;
            lastK = k;
         }
      }
      if (mf > mfMax)
      {
         nbB = nb;
         kbB = kb;
         mfMax = mf;
         mpMax = mp;
      }
   }
   if (mpMax)
   {
      int i;
      for (i=0, mp=mmb; mp != mpMax; i++, mp = mp->next);
      mp = mpMax;
      mpMax = CloneMMNode(mpMax);
      mpMax->ivar = i + 1;
      mpMax->mflop[0] = mfMax;
      mpMax->mbB = mb;
      mpMax->nbB = nbB;
      mpMax->kbB = kbB;
   }
   return(mpMax);
}

ATL_mmnode_t *FindBestMK_N
   (char pre, int verb, ATL_mmnode_t *mmb, unsigned int nb)
/*
 * RETURNS: Clone of node that provides best perf for a problem with NB=nb.
 *          We try finding best mb&kb by using ones near any that are near
 *          those in list of mmb.  NULL if no kern can do NB=nb.
 */
{
   ATL_mmnode_t *mp, *mpMax=NULL;
   unsigned int kbB=0, mbB=0;
   double mf, mfMax=0.0;
/*
 * Try all candidate kernels for this nb
 */
   for (mp=mmb; mp; mp = mp->next)
   {
      int kb=mp->kbB, mb=mp->mbB;
      if (nb % mp->nu)
         continue;
      mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, kb, 1, 0, -1);
/*
 *    Try MB (and if allowed KB) from larger kerns to see if they improve
 *    small nb perf
 */
      if (mp->next)  /* larger KB/MB kerns exist */
      {
         ATL_mmnode_t *mpG; /* kerns of greater size */
         int KB, lastK=kb, lastM=mb;
         const int ku=mp->ku, mu=mp->mu;
         const int CHGK=FLAG_IS_SET(mp->flag, MMF_KRUNTIME);

         for (mpG=mp->next; mpG; mpG = mpG->next)
         {
            int k = mpG->kbB, m=mpG->mbB;
            double mfG;

            if (mu <= 8)
               m = ((m+mu-1)/mu)*mu;
            else
               m = (m/mu)*mu;
            if (CHGK)
            {
               if (ku <= 8)
                  k = ((k+ku-1)/ku)*ku;
               else
                  k = (k/ku)*ku;
            }
            else
               k = kb;
            if (m <= lastM || k <= lastK)
               continue;
            mfG = TimeMMKernel(verb, 0, mp, pre, m, nb, k, 1, 0, -1);
            if (mfG > mf)
            {
               mf = mfG;
               mb = m;
               kb = k;
            }
            lastM = m;
            lastK = k;
         }
      }
      if (mf > mfMax)
      {
         mbB = mb;
         kbB = kb;
         mfMax = mf;
         mpMax = mp;
      }
   }
   if (mpMax)
   {
      int i;
      for (i=0, mp=mmb; mp != mpMax; i++, mp = mp->next);
      mp = mpMax;
      mpMax = CloneMMNode(mpMax);
      mpMax->ivar = i + 1;
      mpMax->mflop[0] = mfMax;
      mpMax->nbB = nb;
      mpMax->mbB = mbB;
      mpMax->kbB = kbB;
   }
   return(mpMax);
}

ATL_mmnode_t *FindBestK_MN
   (char pre, int verb, ATL_mmnode_t *mmb, unsigned int nb)
/*
 * RETURNS: Clone of node that provides best perf for a problem with MB=NB=nb,
 *          and any KB in the list of mmb.  NULL if no kern can do MB=NB=nb.
 */
{
   ATL_mmnode_t *mp, *mpMax=NULL;
   unsigned int kbB=0;
   double mf, mfMax=0.0;
   for (mp=mmb; mp; mp = mp->next)
   {
      int kb=mp->kbB;
      if (nb % mp->mu || nb % mp->nu)
         continue;
      mf = TimeMMKernel(verb, 0, mp, pre, nb, nb, kb, 1, 0, -1);
/*
 *    If we can change KB, try larger ones before giving up on this kern
 */
      if (mp->next && FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
      {
         ATL_mmnode_t *mpK;
         int KB, lastK=kb;
         const int ku=mp->ku, kbmin=mp->kbmin, kbmax=mp->kbmax;

         for (mpK=mp->next; mpK; mpK = mpK->next)
         {
            int k = mpK->kbB;
            double mfK;

            if (k < kbmin)
               k = kbmin;
            else if (kbmax && k > kbmax)
               k = kbmax;
            else if (ku <= 8 || ku > k)
               k = ((k+ku-1)/ku)*ku;
            else
            {
               k = (k/ku)*ku;
               if (k <= lastK)
                  continue;
            }
            mfK = TimeMMKernel(verb, 0, mp, pre, nb, nb, k, 1, 0, -1);
            if (mfK > mf)
            {
               mf = mfK;
               kb = k;
            }
            lastK = k;
         }
      }
      if (mf > mfMax)
      {
         kbB = kb;
         mfMax = mf;
         mpMax = mp;
      }
   }
   if (mpMax)
   {
      int i;
      for (i=0, mp=mmb; mp != mpMax; i++, mp = mp->next);
      mpMax = CloneMMNode(mpMax);
      mpMax->ivar = i + 1;
      mpMax->mflop[0] = mfMax;
      mpMax->mbB = mpMax->nbB = nb;
      mpMax->kbB = kbB;
   }
   return(mpMax);
}

void GenAllIPViews(char pre, int verb)
/*
 * Generates performance views of geAMMRES.sum
 */
{
   ATL_mmnode_t *mmb, *mp, *mpL;
   mmb = ReadMMFileWithPath(pre, "res", "geAMMRES.sum");
   assert(mmb);
   mp = TimeMMFileWithPath(pre, "res", "ipmnPERF.sum", 0, verb|1, 0, 1, 0, -1);
   if (mp)
      KillAllMMNodes(mp);
   else /* generate inner-product MB=NB, KB free view */
   {
      int i, nb, NBMAX, k;
      ATL_mmnode_t *mb=NULL;
      double mfL=0.0;

      printf("   FINDING BEST CASE FOR DEGENERATE MB=NB\n");
      i = GetOffset(mmb, &(mmb->next));
      GetIntMaxMinAtOff(mmb, i, GetOffset(mmb, &(mmb->nbB)), &NBMAX, &nb);
      GetIntMaxMinAtOff(mmb, i, GetOffset(mmb, &(mmb->mbB)), &k, &i);
      NBMAX = Mmax(NBMAX,k);
      nb = Mmin(nb, i);
      for (mpL=mmb; mpL->next; mpL=mpL->next);
      while (nb <= NBMAX)
      {
         mp = FindBestK_MN(pre, verb, mmb, nb);
         if (mp && mp->mflop[0] >= mfL)
         {
            mfL = mp->mflop[0];
            printf("      BEST CASE NB=%d: %d-%s, KB=%d  mf=%.2f.\n", mp->nbB,
                   mp->ID, mp->rout ? mp->rout:"gen", mp->kbB, mfL);
            mp->next = mb;
            mb = mp;
         }
         nb++;
      }
      mb = ReverseMMQ(mb);
      WriteMMFileWithPath(pre, "res", "ipmnPERF.sum", mb);
      KillAllMMNodes(mb);
      printf("   DONE DEGENERATE MB/NB SEARCH.\n\n");
   }

   mp = TimeMMFileWithPath(pre, "res", "ipnPERF.sum", 0, verb|1, 0, 1, 0, -1);
   if (mp)
      KillAllMMNodes(mp);
   else /* generate inner-product NB < MAXNB, MB,KB free view */
   {
      int i, nb, NBMAX, k;
      ATL_mmnode_t *mb=NULL;
      double mfL=0.0;

      printf("   FINDING BEST CASE FOR DEGENERATE NB\n");
      i = GetOffset(mmb, &(mmb->next));
      GetIntMaxMinAtOff(mmb, i, GetOffset(mmb, &(mmb->nbB)), &NBMAX, &nb);
      while (nb <= NBMAX)
      {
         mp = FindBestMK_N(pre, verb, mmb, nb);
         if (mp && mp->mflop[0] >= mfL)
         {
            mfL = mp->mflop[0];
            printf("      BEST CASE NB=%d: %d-%s, MB=%d, KB=%d  mf=%.2f.\n",
                   mp->nbB, mp->ID, mp->rout ? mp->rout:"gen",
                   mp->mbB, mp->kbB, mfL);
            mp->next = mb;
            mb = mp;
         }
         nb++;
      }
      mb = ReverseMMQ(mb);
      WriteMMFileWithPath(pre, "res", "ipnPERF.sum", mb);
      KillAllMMNodes(mb);
      printf("   DONE DEGENERATE NB SEARCH.\n\n");
   }

   mp = TimeMMFileWithPath(pre, "res", "ipmPERF.sum", 0, verb|1, 0, 1, 0, -1);
   if (mp)
      KillAllMMNodes(mp);
   else /* generate inner-product NB < MAXNB, MB,KB free view */
   {
      int i, nb, NBMAX, k;
      ATL_mmnode_t *mb=NULL;
      double mfL=0.0;

      printf("   FINDING BEST CASE FOR DEGENERATE MB\n");
      i = GetOffset(mmb, &(mmb->next));
      GetIntMaxMinAtOff(mmb, i, GetOffset(mmb, &(mmb->mbB)), &NBMAX, &nb);
      while (nb <= NBMAX)
      {
         mp = FindBestNK_M(pre, verb, mmb, nb);
         if (mp && mp->mflop[0] >= mfL)
         {
            mfL = mp->mflop[0];
            printf("      BEST CASE MB=%d: %d-%s, NB=%d, KB=%d  mf=%.2f.\n",
                   mp->mbB, mp->ID, mp->rout ? mp->rout:"gen",
                   mp->nbB, mp->kbB, mfL);
            mp->next = mb;
            mb = mp;
         }
         nb++;
      }
      mb = ReverseMMQ(mb);
      WriteMMFileWithPath(pre, "res", "ipmPERF.sum", mb);
      KillAllMMNodes(mb);
      printf("   DONE DEGENERATE MB SEARCH.\n\n");
   }

   KillAllMMNodes(mmb);
}
/*
 * Finds main amm kernels
 */
ATL_mmnode_t *DoMainMM(char pre, int verb, int nregs, int CS)
{
   ATL_mmnode_t *mmb, *sqmmb, *mp, *mmGN;
   const char upr = (pre == 'z' || pre == 'd') ? 'd' : 's';
/*
 * If we are quick returning, still must possibly retime subsearches
 */
   mmb = TimeMMFileWithPath(pre, "res", "geAMMRES.sum", 0, verb|1, 0, 1, 0, -1);
   if (mmb)
   {
      KillAllMMNodes(DoSquare(pre, verb, nregs, mmb));
      return(mmb);
   }
/*
 * Find the generated cases for each cache context
 */
   mmGN = GetGenCases(pre);
   assert(mmGN);
/*
 * Find which user-supplied kernels can compile on this platform
 */
   mmb = GetWorkingUserCases(verb, upr);
/*
 * Add gmmsearch-suggested kernels to those being considered.
 */
   mp = ATL_LastMMNode(mmGN);
   mp->next = mmb;
   mmb = mmGN;
/*
 * Find the best user-supplied cases for the three common cache blking cases
 */
   mmb = FindCacheBudgetCasesByKB(verb, pre, CS, mmb);
   #if 0  /* switching square case to rank-K */
/*
 *    Now find best square cases
 */
      sqmmb = DoSquare(pre, verb, nregs, mmb);
/*
 *     Now merge square & rect lists for best performing kernels, and write out
 */
      mp = MergeCases(0, mmb, sqmmb);
      KillAllMMNodes(sqmmb);
      KillAllMMNodes(mmb);
      mmb = WinnowCases(0, mp);
   #else
      mmb = WinnowCases(0, mmb);
   #endif
   WriteRefreshedMMFileWithPath(pre, "res", "geAMMRES.sum", mmb);
   return(mmb);
}

int KernHandlesThisKB(ATL_mmnode_t *mp, int kb)
{
   if (!mp->ID)  /* genned kernel handles all K (kvec needs padding) */
      return(1);
   if (mp->kbmax && kb > mp->kbmax)
      return(0);
/*
 * KVEC kerns can always handle problems within vlen-1 less than their kb
 * due to padding in copy routines.
 */
   if (FLAG_IS_SET(mp->flag, MMF_KVEC))
   {
      int KB;
      const int ku = mp->ku;
      if (mp->kbmin)
      {
         KB = (mp->kbmin / mp->vlen)*mp->vlen;
         if (KB < mp->kbmin)
            return(0);
      }
      KB = ((kb+mp->vlen-1)/mp->vlen)*mp->vlen;
      if ((KB/ku)*ku != kb)
         return(0);
   }
   else
   {
      const int ku = mp->ku;
      if (kb < mp->kbmin)
         return(0);
      if ((kb/ku)*ku != kb)
      {
         if (mp->ID)
            return(0);
         else if (!FLAG_IS_SET(mp->flag, MMF_KUISKB))
            return(0);  /* genned codes can adjust fully-unrolled loop */
      }
   }
   return(1);
}

/*
 * This routine tries all kernels in mmb, with K=KB, M=CEIL(KB/MU)*MU,
 * N=CEIL(KB/NU)*NU.  Kernels that can't handle KB are rejected
 * RETURNS: new mmnode ptr for best case for this KB; cannot be NULL, because
 *          1st param of mmb must be a generated kernel that works for any KB
 */
ATL_mmnode_t *BestKernForKB(int verb, char pre, ATL_mmnode_t *mmb, int KB)
{
   ATL_mmnode_t *mp, *mpB=NULL;
   double mfB=0.0;
   assert(mmb);
   printf("   FINDING BEST NEAR-SQUARE KERNEL WT KB=%d:\n", KB);
   for (mp=mmb; mp; mp = mp->next)
   {
      const int mu=mp->mu, nu=mp->nu, ku=mp->ku, ID=mp->ID;
      const int mb=((KB+mu-1)/mu)*mu, nb=((KB+nu-1)/nu)*nu;
      int kb=KB;
      const char *rt=mp->rout;
      double mf;

      if (!KernHandlesThisKB(mp, KB))
      {
         if (ID)
            printf("      %d-%s: skipped, cannot handle KB=%d\n", ID, rt, KB);
         else
         {
            printf("      0-");
            PrintGen0(stdout, mp, 0, 0, 0);
            printf(": skipped, cannot handle KB=%d\n", KB);
         }
         continue;
      }
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         kb = ((KB+mp->vlen-1)/mp->vlen)*mp->vlen;
      if (!mp->ID && FLAG_IS_SET(mp->flag, MMF_KUISKB))
      {
         mp->kbmin = kb-mp->vlen+1;
         mp->kbmax = mp->ku = kb;
         mp->kbB = KB;
         if (mp->genstr);
            free(mp->genstr);
         mp->genstr = MMGetGenString(pre, mp);
      }
      mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, KB, 0,  0, -1);
      if (mf > mfB)
      {
         mfB = mf;
         mpB = mp;
      }
      if (ID)
         printf("      %d-%s, M=%d, N=%d, K=%d: %.0f\n", ID, rt, mb,nb,KB, mf);
      else
      {
         printf("      0-");
         PrintGen0(stdout, mp, mb, nb, KB);
         printf(": %.0f\n", mf);

      }
   }
   assert(mpB);
   printf("   BEST FOR KB=%d: %d-%s (%.1f MFLOPS)\n",
          KB, mpB->ID, mpB->rout, mfB);
   mpB = CloneMMNode(mpB);
   mpB->mflop[0] = mfB;
   mpB->kbB = KB;
   mpB->mbB = ((KB+mpB->mu-1)/mpB->mu)*mpB->mu;
   mpB->nbB = ((KB+mpB->nu-1)/mpB->nu)*mpB->nu;
   return(mpB);
}

ATL_mmnode_t *DoRankK(char pre, int verb, int nregs, const ATL_mmnode_t *mainb)
{
   ATL_mmnode_t *rkb=NULL, *mp, *mmb;
   const ATL_mmnode_t *cmp;
   int maxB=0, b;
   const char upr = (pre == 'z' || pre == 'd') ? 'd' : 's';
   rkb = TimeMMFileWithPath(pre, "res", "rkAMMRES.sum", 0, verb|1, 0, 1, 0, -1);
   if (rkb)
      return(rkb);
/*
 * Find largest KB used by main kernels; we will time all near-square kernels
 * of this size and below
 */
   for (cmp=mainb; cmp; cmp = cmp->next)
      if (cmp->kbB > maxB)
         maxB = cmp->kbB;
   if (!maxB)
      maxB = 256;
/*
 * All we need main kerns for is to find maxB, so now reuse the ptr to hold
 * all user cases that work on this platform
 */
   mmb = GetGenCases(pre);
   mp = ATL_LastMMNode(mmb);
   mp->next = GetWorkingUserCases(verb, upr);
   ResetMoveBitsInQ(rkb, 5);
   printf("TUNING RANK-K, 3 <= K <= %d:\n", maxB);
   for (b = maxB; b > 2; b--)
   {
      mp = BestKernForKB(verb, pre, mmb, b);
      mp->next = rkb;
      rkb = mp;
   }
   WriteRefreshedMMFileWithPath(pre, "res", "rkAMMRES.sum", rkb);
   return(rkb);
}

int main(int nargs, char **args)
{
   char pre='d';
   int verb, nregs, nb, CS, gmu, gnu;
   char *fnout;
   ATL_mmnode_t *mmb, *rnkK, *grnkK, *emb, *mp, *mmSQ;

   GetFlags(nargs, args, &pre, &verb, &nregs, &nb, &CS);
   mmb = DoMainMM(pre, verb, nregs, CS);
   emb = TimeExtraBlockings(pre, verb);
   if (emb)
   {
      emb = SortMMQByIntVal(emb, &(emb->kbB));
      mp = MergeCases(0, mmb, emb);
      KillAllMMNodes(mmb);
      KillAllMMNodes(emb);
      mmb = WinnowCases(0, mp);
      WriteMMFileWithPath(pre, "res", "geAMMRES.sum", mmb);
   }
   rnkK = DoRankK(pre, verb, nregs, mmb);
   DoSquare(pre, verb, nregs, rnkK);
/*
 * Handle K-cleanup
 */
   ComputeKClean(verb, pre);
/*
 * Time TRSM kernels with matching block factors as square cases
 */
   DoTRSM(pre, verb);
/*
 * Join mmb & rnkK to make master list of all kernels required in this search
 */
   if (rnkK)
   {
      for (mp=mmb; mp->next; mp = mp->next);
      mp->next = rnkK;
   }
   mp = GetUniqueUserKerns(mmb);
   KillAllMMNodes(mmb);
   WriteRefreshedMMFileWithPath(pre, "res", "AMMFRCLST.sum", mp);
   KillAllMMNodes(mp);
   GetGenKernForNB(pre, 0);  /* dealloc static mmnodes */
   exit(0);
}
