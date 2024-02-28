#include "atlas_mmtesttime.h"
void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -i ammsearch.sum: repeat for multiple\n");
   fprintf(stderr, "   -o <output file>\n");
   exit(ierr ? ierr : -1);
}

ATL_mmnode_t *GetFlags(int nargs, char **args, char *PRE, char **FOUT)
{
   ATL_mmnode_t *cb=NULL, *cp;
   char *fout=NULL;
   int i;
   char pre='d';

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
      case 'i':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        if (!cb)
           cb = ReadMMFile(args[i]);
        else
        {
           cp = ATL_LastMMNode(cb);
           cp->next = ReadMMFile(args[i]);
        }
        break;
      case 'o':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        fout = DupString(args[i]);
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   assert(cb);
   for (cp=cb; cp; cp = cp->next)
   {
/*
 *    TRSM kernels are just GEMM kernels w/o K cleanup.
 *    TRMM kernels aren't, but we put a matching gemm kernel in anyway, so
 *    call them GEMM to avoid problems with copies.
 */
      if (cp->blask == ATL_KTRSM || cp->blask == ATL_KTRMM)
      {
         cp->blask = ATL_KGEMM;
         cp->ivar = 0;
      }
      else
         cp->ivar = FLAG_IS_SET(cp->flag, MMF_KCLN) ? -1 : 0;
   }
   if (!fout)
   {
      fout = malloc(18);
      assert(fout);
      strcpy(fout, "res/");
      fout[4] = pre;
      strcpy(fout+5, "AMMLIST.sum");
   }
   *FOUT = fout;
   *PRE = pre;
   return(cb);
}

void FindValidCleanKUs(ATL_mmnode_t *p, ATL_mmnode_t *mb)
/*
 * Searches master list given by mb, finding all kernel compiles that match p
 * to find all kbB that this kern comp is used with.  This allows us to compute
 * whether a given KU is legal.  Recall that we can pad kbB up to 4.  1 is
 * always a multiple of all kbB, and so can be used for cleanup always.
 * However, imagine we want to have pad to 3 for cleanup.  If any target
 * kernel is used with a kernel where kbB % 3 != 0, we would either have
 * memory overwrite on non-cleanup code, or have to have non-cleanup nodes
 * allocate the extra space.  Therefore, for the possibly incompatible KUs of
 * 4, 3, 2, we rule them out if some target kbB doesn't evenly divide them.
 * mp->flag sets: 31:4, 30:3, 29:2 -- set to 1 if OK clean mul, else 0.
 */
{
   unsigned int flg = (0x7<<29);
   ATL_mmnode_t *mp;

   if (!FLAG_IS_SET(p->flag, MMF_KCLN))
      return;
   for (mp=mb; mp; mp = mp->next)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KCLN) && MMKernCompsSame(p, mp))
      {
         const unsigned kb = mp->kbB;
         if (kb & 0x3)            /* if (mp->kbB % 4 != 0) */
            flg &= ~(1<<31);
         if (kb % 3)
            flg &= ~(1<<30);
         if (kb & 1)              /* if (mp->kbB % 2 != 0) */
            flg &= ~(1<<29);
      }
   }
   p->flag = flg | (p->flag & (~(0x7<<29)));
}

ATL_mmnode_t *GetDesiredK1(ATL_mmnode_t *mb)
/*
 * Finds all unique kernels in mb, and then takes the largest blocking factor
 * that matches that signature.
 */
{
   ATL_mmnode_t *ub=NULL, *bp, *mp;
   int mbC, nbC;
   double avgK=0.0;
   if (!mb)
      return(NULL);
/*
 * We know ku=1 legal cleanup, use 3 most sig buts for U=4,3,2 resp,
 * 1 means an unroll of this size can be used for cleanup, 0 means not.
 * For an U to be legal for cleanup, it must be a multiple of all KB.
 */
   assert(32-MMF_MAXBIT >= 3); /* need 3 bits, ku=2,3,4 */
/*
 * K cleanup can be done by any kernel with the same storage for A/B/C, so it
 * doesn't have to be the exact same kernel implementation, assuming it can
 * handle all the required K values.  Therefore, start by finding all storage
 * types that require K-cleanup.
 */
   for (mp=mb; mp; mp = mp->next)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KCLN))
      {
         FindValidCleanKUs(mp, mb);
         if (MMKernCleansK(mp, mp))
            mp->ivar = 0;
         else if (!MMCompatKernPresent(ub, mp))/* no storg compat kern in ub */
         {
            ATL_mmnode_t *p;
            p = CloneMMNode(mp);
            p->next = ub;
            ub = p;
         }
      }
      else
         mp->ivar = 0;
   }
/*
 * We're going to use the same K-cleanup for all storage-compatible kernels,
 * which means it can be used for wildly varying dimensions.  We will use
 * the largest observed MB/NB, while using a small KB.  We'll set KB to roughly
 * half the size of the average KB it'll be cleaning.
 * For each K-clean target, examine all kernels it will clean for to find
 * dims to use.
 */
   for (bp=ub; bp; bp = bp->next)
   {
      int nk=0;
      mbC = nbC = 0;
      mp = mb;
      bp->flag |= 7<<29;
      while (mp = MMCompatKernPresent(mp, bp))
      {
         if (FLAG_IS_SET(mp->flag, MMF_KCLN))
         {
            const unsigned kb = mp->kbB;
            mbC = Mmax(mbC, mp->mbB);
            nbC = Mmax(mbC, mp->nbB);
            if (kb & 0x3)            /* if (mp->kbB % 4 != 0) */
               bp->flag &= ~(1<<31);
            if (kb % 3)
               bp->flag &= ~(1<<30);
            if (kb & 1)              /* if (mp->kbB % 2 != 0) */
               bp->flag &= ~(1<<29);
            avgK += mp->kbB;
            nk++;
         }
         mp = mp->next;
      }
      avgK /= (nk<<1);
      bp->mbB = mbC;
      bp->nbB = nbC;
      bp->kbB = avgK;
      assert(mbC && nbC && bp->kbB);
   }
   return(ub);
}

ATL_mmnode_t *BestCompatKern
(
   char pre,
   ATL_mmnode_t *dp, /* amm kerntype & dims desiring K cleanup */
   ATL_mmnode_t *mb  /* list of kernels to try */
)
/*
 * RETURNS: dp if it can do K cleanup on its own, else NULL if no kernel in mb
 *          can do cleanup, otherwise a clone of the best-performing node in mb
 *          with mflop set to cleanup performance.
 */
{
   ATL_mmnode_t *mp, *mpB=NULL;
   double mfB=0.0;
   int mbB, nbB, kbB, cnt, cntB=0;

   if (!dp)
      return(dp);
   if (MMKernCleansK(dp, dp))
   {
      printf("   ID=%d '%s' U=(%u,%u,%u) B=(%u,%u,%u), cleans itself!\n",
             dp->ID, dp->rout?dp->rout:"NULL", dp->mu, dp->nu, dp->ku,
             dp->mbB, dp->nbB, dp->kbB);
      return(dp);
   }
   if (!mb)
      return(NULL);
   mbB = dp->mbB;
   nbB = dp->nbB;
   kbB = dp->kbB;

   printf(
"   SEARCHING FOR BEST K-CLEANER FOR ID=%d, U=(%u,%u,%u) B=(%u,%u,%u) %c%d\n",
          dp->ID, dp->mu, dp->nu, dp->ku, dp->mbB, dp->nbB, dp->kbB,
          FLAG_IS_SET(dp->flag, MMF_KVEC)?'K':'M', dp->vlen);
   for (cnt=0,mp=mb; mp; mp = mp->next, cnt++)
   {
      if (MMKernCleansK(dp, mp))
      {
         double mf;
         mf = TimeMMKernel(0, 0, mp, pre, mbB, nbB, kbB, 0, 0, -1);
         printf("      ID=%d, IDX=%d, mf=%.0f\n", mp->ID, cnt, mf);
         if (mf > mfB)
         {
            mfB = mf;
            mpB = mp;
            cntB = cnt;
         }
      }
   }
   if (mpB)
   {
      printf("   DONE, BEST: ID=%d, IDX=%d, mf=%.0f.\n\n", mpB->ID, cntB, mfB);
      mpB = CloneMMNode(mpB);
      mpB->mbB = mbB;
      mpB->nbB = nbB;
      mpB->kbB = kbB;
      mpB->mflop[0] = mfB;
   }
   else
      printf("   DONE: NO KCLEAN FOR THIS KERNEL!\n\n");
   return(mpB);
}

ATL_mmnode_t *FindAllBestCompatKern
(
   char pre,
   ATL_mmnode_t *db, /* amm kerntype & dims desiring K cleanup */
   ATL_mmnode_t *mb  /* list of kernels to try */
)
{
   ATL_mmnode_t *kb=NULL, *dp;

   for (dp=db; dp; dp = dp->next)
   {
      ATL_mmnode_t *mpB;
      mpB = BestCompatKern(pre, dp, mb);
      if (mpB == dp)
         dp->ivar = -1;
      else
      {
         if (!mpB) /* no cleaner so far found */
         {
            mpB = CloneMMNode(dp);  /* keep placeholder */
            mpB->mflop[0] = -2.0;   /* and indicate no cleaner available yet */
         }
         mpB->next = kb;
         kb = mpB;
      }
   }
   return(ReverseMMQ(kb));  /* now in same order as db */
}

ATL_mmnode_t *FindK1Clean
(
   char pre,
   char *fout,       /* output filename */
   ATL_mmnode_t *db, /* storage formats & NBs desiring K cleanup */
   ATL_mmnode_t *mb  /* already-existing kernels to be scope for cleanup */
)
{
   ATL_mmnode_t *dp, *kb=NULL, *gb=NULL, *nb, *mp;
   unsigned int L1ELTS;
   char upr;

   if (pre == 'z')
      upr = 'd';
   else
      upr = (pre == 'c') ? 's' : pre;
   printf("FINDING K-CLEANERS AMONGST EXISTING KERNELS\n");
   kb = FindAllBestCompatKern(pre, db, mb);
/*
 * Now we consider adding kernels to do K-cleanup.  kb is best found so far
 * (-2 mflop indicates no cleanup so far, so forced new kern).
 * We have two options for K-clean: any working user-contributed kernel, or
 * a generated kernel.
 */
/*
 * Create generated cleaners for each kernel, using best-found flags from
 * prior searches.  Recall that gAMMRES first three cases are best cases
 * for (# is order in file): 0: A,B,C fit in L1, 1: B fits in L1, 2: mu*K panel
 * fits in L1 (L2 blocked).  Since M,N,K will vary by kernel, we can't be sure
 * of which of these to use, but since we can control only prefetch and nobcast,
 * its no big deal if we use the wrong one.
 */
   nb = ReadMMFileWithPath(pre, "res", "gAMMRES.sum"); /* tuned params */
   assert(nb); /* all 3 fit L1, ->B fits; ->->pan fits */
   L1ELTS = GetL1CacheElts(upr);
   for (dp=db; dp; dp = dp->next)
   {
      ATL_mmnode_t *gp, *op;
      const int M=dp->mbB, N=dp->nbB, K=dp->kbB;
      int flag;

      if (N*M + K*(N+M) <= L1ELTS)
         op = nb;
      else
         op = (K*N <= L1ELTS) ? nb->next : nb->next->next;
      flag = op->flag & MMF_ALLPF;
      if (!FLAG_IS_SET(dp->flag, MMF_KVEC))
         flag |= op->flag & (1<<MMF_NOBCAST);
      gp = MMGetGenCaseKClean(dp);
      gp->flag |= flag;
      gp->rout = DupString("ATL_tmp.c");
      gp->genstr = MMGetGenString(pre, gp);
      gp->next = gb;
      gb = gp;
   }
   KillAllMMNodes(nb);
/*
 * Get all working user-contributed kernels, and grab the ones that could
 * clean for selected kerns
 */
   nb = ReadMMFileWithPath(pre, "res", "WORKING.sum");
   if (nb)
   {
      mp = MMGetKCleanQ(db, nb);
      KillAllMMNodes(nb);
      if (mp)
      {
         nb = ATL_LastMMNode(mp);
         nb->next = gb;
         gb = mp;
      }
   }
/*
 * Now search for best new kernel to use amongst hand-tuned & generated
 */
   printf("FINDING K-CLEANERS AMONGST NEW KERNELS\n");
   nb = FindAllBestCompatKern(pre, db, gb);
   KillAllMMNodes(gb);
/*
 * nb (best new) & kb (best existing) are in same order, so we can compare
 * them directly.  Require a new kernel to be 4% faster than an existing
 * kernel before adding it, since expanding the list of supported kernels
 * is last thing we need to do.
 */
   gb = NULL;
   for (dp=kb, mp=nb; dp; dp=dp->next, mp=mp->next)
   {
      ATL_mmnode_t *newp;

      if (dp->mflop[0]*1.04 < mp->mflop[0])
      {
         newp = CloneMMNode(mp);
         newp->next = gb;
         newp->flag &= ~(1<<MMF_KCLN);
         newp->ivar = 0;
         gb = newp;
      }
   }

   KillAllMMNodes(kb);
   KillAllMMNodes(nb);
   if (gb)
   {
      printf("\nKERNELS ADDED TO MASTER LIST FOR KCLEAN:\n");
      PrintMMNodes(stdout, gb);
   }
   else
      printf("\nNO KERNELS ADDED TO MASTER LIST FOR KCLEAN.\n");
   return(gb);
}

int GetIdxOfKClean
(
   ATL_mmnode_t *dp, /* kern desiring K1 cleaner */
   ATL_mmnode_t *kb, /* queue containing K1 cleaner kerns */
   ATL_mmnode_t *ml  /* master list showing all kbs used by dp */
)
/*
 * dp is kern desiring K1 cleaner, kb is queue containing such cleaners
 */
{
   ATL_mmnode_t *kp;
   int cnt;
/*
 * Prioritize self-cleaning over other-cleaning so that we don't do unnec
 * instruction load for kb0 case
 */
   FindValidCleanKUs(dp, ml);
   if (MMKernCleansK(dp, dp))
   {
      for (cnt=0,kp=kb; kp; kp = kp->next, cnt++)
         if (kp == dp)
            return(cnt);
      assert(0);  /* should not be reached! */
   }
   for (cnt=0,kp=kb; kp; kp = kp->next, cnt++)
   {
      if (MMKernCleansK(dp, kp))
         return(cnt);
   }
   return(-1);
}

void NumberKClean
(
   ATL_mmnode_t *ml, /* non-unique list of all kernel usages */
   ATL_mmnode_t *mb  /* unique list of things requiring cleaning */
)
{
   ATL_mmnode_t *mp;
   int cnt;

   for (cnt=0,mp=mb; mp; mp = mp->next, cnt++)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KCLN))
      {
          int idx;
          idx = GetIdxOfKClean(mp, mb, ml);
          if (idx < 0)
          {
             fprintf(stderr, "WTF: cnt=%d\n", cnt);
             PrintMMLine(stderr, mp);
          }
          assert(idx >= 0);
          mp->ivar = idx + 1;
      }
   }
}

/*
 * This routine takes as input all amm kerns to be used in building ATLAS.
 * Any kernels with MMF_KCLN needs a K-cleaner.  We will search for cleaners
 * for any such kernel.  A K-cleaning kernel must have MMF_KRUNTIME &
 * ku <= X, where X is max of (mu-1,nu-1,4,KVEC?VLEN:1).
 * The best existing kernel will be compared to a generated kernel on perf.
 * If the genned kernel provides higher performance, then it will be added
 * to the list of output kernels (with KCLN=0, obviously).
 * For a given kernel match, only one K-cleaner is generated, which means that
 * multiple kernels may share.
 * In the output file, those kernels not requiring K-cleanup, or who can serve
 * as their own K-cleanup, have ivar=0, KCLN=0.  For those needing Kclean,
 * KCLN=1, and ivar=#, where # is the position in the file of the cleaning kern.
 *
 * This should be run as last step before generation, and the output file
 * will be all the output kernels to generate (more may have been added by
 * search in order to do cleanup).  It will use ivar to express how cleanup
 * is being done.  ivar=0 indicates that no cleanup is required (KCLN not
 * set or the function can handle all K-cleanup itself).
 * ivar > 0 means that kern@ ivar-1 in the given kernel order provided kcleanup.
 */
int main(int nargs, char **args)
{
   ATL_mmnode_t *mb, *db, *kb, *ml;
   unsigned int msk;
   char pre;
   char *fout;

   ml = GetFlags(nargs, args, &pre, &fout);

   db = GetDesiredK1(ml);  /* get list of K cleaners to make */
   mb = AddUniqueMMKernCompList(NULL, ml); /* remove repeated comps frm lst */
   printf("KERNELS REQUIRING NON-SELF KCLEANERS:\n");
   PrintMMNodes(stdout, db);
   printf("\n\n");
/*
 * Find K1 cleaners and add them to the list of final kernels to build
 */
   kb = FindK1Clean(pre, fout, db, mb);
   KillAllMMNodes(db);
   if (kb)
   {
      mb = AddUniqueMMKernCompList(mb, kb);
      KillAllMMNodes(kb);
   }
/*
 * Number ivar with location of Kcleaners, and use this to set str to
 * the actual node address, which reverts ivar to -1 for KCLN nodes; we
 * will renumber after reordering nodes to try to place Kcleaners right
 * next to node they clean up for.
 */
   NumberKClean(ml, mb);
   KillAllMMNodes(ml);
   WriteMMFile("tmp.sum", mb);
   MMIvar2str(mb);
/*
 * Start new queue in kb with all nodes needing non-self K cleaning, and have
 * their K-cleaners placed right after them when possible (an earlier kernel
 * using same K-cleaner will force us to use non-contiguous locations).
 */
   kb = NULL;
   KEEP_LOOKING:
   {
      ATL_mmnode_t *mp;
      mp = MMFindKCleanStr(mb); /* 1st node wth non-NULL,mp address */
      if (mp)
      {
         ATL_mmnode_t *p;
         mb = RemoveMMNodeFromQ(mb, mp);
         p = (void*)mp->str;
         mb = RemoveMMNodeFromQ(mb, p);
         mp->next = p;
         p->next = kb;
         kb = mp;
/*
 *       Now, remove all nodes requiring cleaner mp->str, and put them close
 */
         while ((p=MMFindStrMatch(mb, mp->str)))
         {
            mb = RemoveMMNodeFromQ(mb, p);
            p->next = kb;
            kb = p;
         }
         goto KEEP_LOOKING;
      }
   }
/*
 * Join K-clean and their cleaners together with non-cleaning nodes, and then
 * translate str back to ivar (str does not show up in file, ivar does).
 * Write out final output file that will be blindly used to generate master
 * kernel list for ATLAS.
 */
   if (mb)
      mb = ATL_JoinMMQs(kb, mb);
   MMstr2Ivar(mb);
/*
 * Get rid of temporary flag bits
 */
   msk = (1<<MMF_MAXBIT)-1;
   for (kb=mb; kb; kb=kb->next)
      kb->flag &= msk;
   printf("DONE K CLEANUP SEARCH, OUTPUT IN '%s'\n", fout);
   WriteMMFile(fout, mb);
   free(fout);
   KillAllMMNodes(mb);
   return(0);
}
