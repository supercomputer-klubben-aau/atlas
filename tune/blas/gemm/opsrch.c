#include "atlas_cache.h"
#include "atlas_genparse.h"
#include "atlas_mmtesttime.h"
#include <math.h>

static unsigned int FULLSRCH=0;
void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set type/precision prefix (d) \n");
   fprintf(stderr, "   -F 0/1: don't/do do full search (default 0) \n");
}

char GetFlags(int nargs, char **args)
{
   char pre = 'd';
   int i;

   for (i=1; i < nargs; i++)
   {
      int wch, *ip, **ipp, TST=0;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'F':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         FULLSRCH = atoi(args[i]);
         break;
      case 'p':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        pre = tolower(args[i][0]);
        assert(pre == 's' || pre == 'd' || pre == 'z' || pre == 'c');
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   return(pre);
}

void do_bfi(char pre, int K0, int KN, int iK, char *fin, char *fout,
            double runBon)
{
   int L, i;
   char *ln;
   L = strlen(fout) + 45;
   L += strlen(fin);
   L += NumDecDigits(K0);
   L += NumDecDigits(KN);
   L += NumDecDigits(iK);
   if (runBon > 0.0)
      L += 16;
   ln = malloc(L*sizeof(char));
   assert(ln);
   i = sprintf(ln, "./xbfisrch -p %c -s p O R -i %s -o res/%c%s",
               pre, fin, pre, fout);
   assert(i < L);
   if (K0 && KN && iK)
      i += sprintf(ln+i, " -Bk %u %u %u", K0, KN, iK);
   if (runBon > 0.0)
      i += sprintf(ln+i, " -r %e", runBon);
   assert(i < L);
   i = system(ln);
   if (i)
   {
      fprintf(stderr, "'%s' RET: %d!\n", ln, i);
      assert(0);
   }
   free(ln);
}
int opWorkSetOK(unsigned long csz, ATL_mmnode_t *mp, int mb, int nb, int kb)
/*
 * RETURNS: 0 if working set for usual outer product > csz, else 1
 */
{
   unsigned long sz;
   sz = mb*nb + kb*(mp->mu+mb+nb);
   return(sz <= csz);
}

void do_opdek(int verb, char pre, int NEK)
{
   ATL_mmnode_t *mp, *mb=NULL;
   int tflag = 0;
   unsigned int i, U, N;
   const char ch=NEK?'N':'M';
   char *fn=NEK?"opnek.sum":"opmek.sum";

   mb = TimeMMFileWithPath(pre, "res", fn, 0, verb|1, 0, 0, 0, -1);
   if (mb)
   {
      KillAllMMNodes(mb);
      return;
   }
   mp = ReadMMFileWithPath(pre, "res", "opASYMP.sum");
   if (NEK)
   {
      mp->flag = (mp->flag&(~MMF_MVSET)) |
                 ((1<<MMF_MVC)|(1<<MMF_CPC)|(1<<MMF_MVA));
      tflag |= 1<<31;
   }
   else
      mp->flag = (mp->flag&(~MMF_MVSET)) |
                 ((1<<MMF_MVC)|(1<<MMF_CPC)|(1<<MMF_MVB));
   mp->mflop[0] = 0;
   assert(mp);
   printf("FINDING BEST CASE FOR %c=K OUTER PRODUCT:\n", ch);
   MMExpandNK(verb, pre, tflag, 1, mp);
   printf("BEST %c=K: ID=%u, B=(%u,%u,%u), mf=%.2f\n\n", ch, mp->ID,
          mp->mbB, mp->nbB, mp->kbB, mp->mflop[0]);
   U = mp->ku;
   if (NEK)
   {
      U = ATL_iLCM(U, mp->nu);
      N = mp->nbB;
   }
   else
   {
      U = ATL_iLCM(U, mp->mu);
      N = mp->mbB;
   }
   U = (U == 1 || U == 2) ? 4 : U;
   for (i=U; i < N; i += U)
   {
      ATL_mmnode_t *p;
      unsigned int b;
      printf("   FINDING BEST CASE FOR %c=K=%u:\n", ch, i);
      p = CloneMMNode(mp);
      if (NEK)
         p->nbB = p->kbB = i;
      else
         p->mbB = p->kbB = i;
      MMExpandMorN(verb, pre, tflag, 1, p, !NEK, 512);
      printf("   BEST CASE FOR %c=K=%u: B=(%u,%u,%u), mf=%.2f\n", ch, i,
             p->mbB, p->nbB, p->kbB, p->mflop[0]);
      p->next = mb;
      mb = p;
   }
   mp->next = mb;
   mb = mp;
   mb = ReverseMMQ(mb);
   MMPruneMflopTol(mb, 0, 1.0);
   WriteMMFileWithPath(pre, "res", fn, mb);
   KillAllMMNodes(mb);
}

ATL_mmnode_t *opsrch_rng
   (int verb, char pre, int flag, double runBon, int kb0, int kbN, int maxB,
    unsigned long sz, int (*wrkSetOK)(unsigned long, ATL_mmnode_t*,int,int,int),
    ATL_mmnode_t *mpB, /* best kernel for top of this range */
    ATL_mmnode_t *mmb) /* list of other kerns to try */
/*
 * Does (relatively) quick search for outer product in range K=[kb0,kbN].
 * It will first choose a kernel for kbN.  This kernel will then be blindly
 * used for all K%ku==0 cases.  If ku <= 4, we will consider using it for
 * *all* cases.  We will try K=kBn-(ku-1), and if it wins, this one kernel
 * will be used for call cases, otherwise the winning kernel is chosen for
 * all % = ku-1 cases, and now we check ku-2, and so on.
 * RETURNS: list of kerns with specialized MB,NB for each [kb0,kbN], sorted
 *          in increasing KB order.
 */
{
   unsigned int KU=mpB->ku, kbB=mpB->kbB, k;
   ATL_mmnode_t *mp, *mkb=NULL;
   ATL_mmnode_t **bKs;

   bKs = calloc(KU, sizeof(ATL_mmnode_t*));
   assert(bKs);
   *bKs = mpB;
   printf("SEARCHING FOR DEFAULT CASES FOR ALL KU=%u REMAINDER CASES\n", KU);
   for (k=1; k < KU; k++) /* remainders */
   {
      int kb = kbB - k;
      printf("   FINDING BEST AMM KERNEL FOR KB=%u, REMAINDER=%u\n", kb, k);
      bKs[k] = bestNearMNK(pre, verb, mmb, mpB->mbB, mpB->nbB, kb, 4, 0, 1.015);
      mp = bKs[k];
      printf("   BEST ID=%d, '%s', B=(%d,%d,%d) mf=%.2f\n", mp->ID,
             GetMMLabelName(pre, mp), mp->mbB, mp->nbB, mp->kbB,mp->mflop[0]);
   }
   printf("DONE SEARCH OF ALL KU=%u CASES\n\n", KU);

   printf("SEARCHING FOR TUNED M&N FOR ALL %u <= K <= %u:\n", kb0, kbN);
   for (k=kbN; k >= kb0; k--)
   {
      mp = CloneMMNode(bKs[k%KU]);
      printf("   FINDING BEST M/N-B FOR KB=%u\n", k);
      mp->kbB = k;
      MMExpandMN(verb, pre, 4, 0, maxB, mp);
      printf("   BEST ID=%d, '%s', B=(%d,%d,%d) mf=%.2f\n", mp->ID,
             GetMMLabelName(pre, mp), mp->mbB, mp->nbB, mp->kbB,mp->mflop[0]);
      mp->next = mkb;
      mkb = mp;
   }
   printf("DONE TUNING FOR ALL %u <= K <= %u.\n\n", kb0, kbN);
   for (k=1; k < KU; k++)
      KillMMNode(bKs[k]);
   free(bKs);
   return(mkb);
}

/*
 * This search forms the backbone of ATLAS outer-product support.  It will
 * find all the kernels needed to support 3 categories of degenerate dims:
 *    opk.sum : K is degenerate, M&N can be varied for performance
 *    opmk.sum: M&K are degenerate, N can be varied for performance
 *    opnk.sum: N&K are degenerate, M can be varied for performance
 *    opsq.sum: all dims degenerate, force true square, does not cover all K
 *
 * For mn/nk, increasing the other dim will not improve performance, because
 * any intra-kernel reuse is already being achieved via inter-kernel reuse.
 * Therefore, we will increase the free dim only for very small probs where
 * it may help with function call overhead.
 *
 * For sq, all dims must be same so that we can use same storage for any
 * array.  This is the only search where some K cannot be handled, because
 * it should only be used by algorithms that can choose K
 *
 * For true rank-K, where only K is degenerate, increasing M will just
 * transfer inter-kernel reuse of B to intra-kernel reuse, and so it will
 * always be set to near-square with K.  Increasing N, however, will allow
 * more reuse of A, and so will be fully searched.
 */
/*
 * This search's job is to find all the kernels we are going to use in
 * building the rank-K update where K is the only degenerate dimension.
 * It searches using near-square block factors, where we split B into diff
 * regions.  The first region is tiny B where performance is discontinuous.
 * For this small region, we try all available kernels at all proposed sizes.
 *
 * For regions with larger blocks, we pick a midrange problem, and then
 * compare all available kernels at that rough size, and pick the best
 * performing.  Call this kernel RCK (Regional Champion Kernel).  The RCK
 * must have the ability to vary KB for the entire region (though it can have
 * any ku <= 32), either as runtime or compile-time KB.  Kernels that can
 * only handle a subset of the range  will be moved to a different queue,
 * and then later compared to the RCK in the subset of ranges where they work,
 * and they can displace the RCK's use for any place where they are faster.
 * Call these kernels RRK (Restricted Range Kernels).  All other kernels
 * are GPK (General Population Kernels).  These files are eventually written
 * out as op[rck,rrk,gpk].sum, though not until later refinements are finished.
 *
 * Of particular interest is the RCK of the largest region, which is used
 * to find the maximum KB that increases performance.
 *
 * If the RCK of two adjacent regions is the same kernel, then those regions
 * will be combined together.
 *
 * We now attempt to bound the regions more precisely, by comparing the RCK's
 * of adjacent regions, and finding when one takes over from the other as
 * precisely as possible.
 *
 * We now consider the RRK.  For each such kernel, we pick a particular
 * size that it handles, and compare it head-to-head vs the RCK at a mutually
 * compatible KB, if any.  If they have no compatible KB, it will be compared
 * to clost RKB of smaller size.  Any RRK that loses to the RCK is put back
 * into GPK.  We now tourney the RRK at all KB it handles against the RCK:
 * any KB for which it wins gets a new node in the oprrk.sum, which is
 * written out long with opgpk.sum at the end of this step.
 *
 * We now have good kernels for entire range, but only believe we've got the
 * best we can do when KB%ku, or we're at a KB where a RRK was used.  We now
 * discover all KUs that are not within 4 (or VLEN if the kern is vectorized
 * along the K dim) of the RRK+RCK kernels.  For such KBs, we discover the KUs
 * that could serve them, and then time all the GPK at a large KB, and
 * eliminate all GPK that don't win one of those KU contests.  We now have
 * the best kerns found for each KU, with possible kerns from RCK+RRK+GPK,
 * but only the kerns known to win somewhere.
 *
 * We now have kernels that can cover the entire range, and we now time all
 * surviving kernels at any unhandled KB and take the best.  This result is
 * written out as oprkk.sum.
 *
 * NOTE: stoped working here, because I suspect what we really want for rank-K
 *       is MB=MU, and NB=KB in order to maximize reuse of A while filling the
 *       cache with B, which is reused across amm calls.  Need to write a BFI
 *       search to prove/disprove.
 *
 * is written out as oprkk.sum, sorted by KB.  We do not presently have
 * all KBs in oprkk.sum, but merely all the kernels we expect to use to handle
 * the range.  We won't finalize this kernel list until we have sear
 *
 * The output file of this search will be a list of kernels sorted on KB.
 * Other dims are chosen near KB, but may not match KB, and so may not be
 * strictly increasing.  The first region is the only one with all KB filled
 * in.  Other regions have a single mention of the champion kernel when the
 * range begins, and that kernel is then used until the next region begins.
 * The final range will have a 2nd mention of the champion kernel at the
 * largest block factor that improved performance.  If all regions except the
 * first use the same kernel, then there will be only two regions according
 * to the ouput file.
 *
 * With this search complete, we know the largest square block factor that
 * We first do all rank-K searches using this regions idea, and all of the
 * used kernels will comprise our rank-K amm compilations.  This subset of
 * kernels will then be used to complete the kernels
 * greatest to least by KB (other dims are chosen near KB, but do not have
 * to be exactly KB; this means
 */
void do_opk(int verb, char pre, char *fnam, int maxK)
{
   long i, nelt_llpc, nelt_llc, nelt=L1C_SZ>>3;
   unsigned int b0=60;   /* largest kb where we brute-force search */
   unsigned int bP;      /* largest square block fitting in LLPC */
   unsigned int bL;      /* largest square block fitting in LLC */
   unsigned int MB, NB, KB, maxB;
   unsigned int eltsh=0;
   double mf;
   const char upr = (pre == 'z' || pre == 'd') ? 'd' : 's';
   ATL_mmnode_t *mb, *mp, *tb, *mpB;

   tb = TimeMMFileWithPath(pre, "res", "opFNL.sum", 0, verb|1, 0, 0, 0, -1);
   if (tb)
   {
      KillAllMMNodes(tb);
      return;
   }
   nelt_llc = (pre == 'd' || pre == 'z') ? (LLC_SZ>>3) : (LLC_SZ>>2);
   eltsh = (upr == 'd') ? 3 : 2;
   tb = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
   assert(tb);
/*
 * Get best-performing kernel, and add it and a kruntime-enabled copy of it
 * to list
 */
   mp = FindMaxMflopMMQ(tb, 0);
   mb = CloneMMNode(mp);
   if (!FLAG_IS_SET(mb->flag, MMF_KRUNTIME))
   {
      mp = CloneMMNode(mb);
      mp->flag |= 1<<MMF_KRUNTIME;
      mp->next = mb;
      mb = mp;
   }
   mb = AddUniqueMMKernCompList(mb, tb);
   KillAllMMNodes(tb);
   mp = ATL_LastMMNode(mb);
   mp->next = GetWorkingUserCases(verb, upr);
   for (mp=mb; mp; mp = mp->next)
      mp->flag = (mp->flag & ~MMF_MVSET) |
                 ((1<<MMF_MVA)|(1<<MMF_MVC)|(1<<MMF_CPC));
   WriteMMFile("res/tmp.sum", mb);
/*
 * If the L1 is the last level of private cache, fitting all blocks necessary
 * to generate no traffic to non-private busses beyond the forced movement
 * is of interest.  This fits both wC/C blks + A,B and next A in the cache:
 *   2MN + K*(MU+M+N) -> 2b^2 +2b^2 + MU*b -> 4b^2 + MU*b - L1elts = 0
 * Applying quadratic equations give: b = (-MU + sqrt(MU*MU+4*4*L1elts))/8
 * --> b = (-MU + sqrt(MU^2 + 16*L1elts))/8 ~= -1 + sqrt(1+16*L1elts))/8
 * --> b ~= sqrt(16*L1elts)/8 = 4/8 sqrt(L1elts) = sqrt(L1elts)/2
 * The above approximations will give us too large a block factor to fit
 * everything in theory, but in practice non-LRU caches will retain large
 * amounts of data well beyond this region, and this is just a search bound,
 * not a performance prediction, so it is OK.
 *
 * For the low end of this region, performance varies unpredictably, and
 * can occasionally do so even towards the end of the region.  Therefore,
 * since this is the shared-traffice minimizing region for this arch, brute-
 * force search the entire region.
 */
/*
 * Find best case for LLC.  We will first solve KB for largest square case
 * that fits, and use this size to find the best-performing kernel with a
 * fixed-K-unroll (to avoid i-cache thrashing for huge block sizes).
 */
   if (pre == 'z' || pre == 'c')
   {
      nelt=L1C_SZ >> eltsh;
      for (i=12; i*i <= nelt; i += 4);
      b0 = (i-4)>>1;
      #if !defined(LLPC_LVL) || LLPC_LVL == 1
         bP = b0;
      #else
         nelt = LLPC_SZ >> eltsh;
         for (i=b0; i*i <= nelt; i += 4);
         bP = (i-4)>>1;
      #endif
   }
   else
   {
   #if !defined(LLPC_LVL) || LLPC_LVL == 1
      nelt=L1C_SZ >> eltsh;
      for (i=12; i*i <= nelt; i += 4);
      bP = b0 = (i-4)>>1;
/*
 * If there is another private cache above L1, then this region is mainly
 * of interest for tiny problems, and is less critical overall, so just brute
 * force the search in the beginning region where perf is so unpredictable.
 */
   #else
      b0 = 32;
      nelt = LLPC_SZ >> eltsh;
      for (i=b0; i*i <= nelt; i += 4);
      bP = (i-4)>>1;
   #endif
   }
   do_bfi(pre, 3, b0, 1, "res/tmp.sum", "opL1.sum", -1.0);

/*
 * For large problems, fully unrolling K may result in i-cache thrashing,
 * so we want to explore very large problems using loops with KU < KB.
 * tb will contain all fully-unrolled kernels, while mb will be those with
 * fixed KU.  We will use mb to explore large blocking factors, but will
 * then need to time fully-unrolled to find best-performing case once rough
 * blocking is established.
 */
   nelt = LLC_SZ >> eltsh;
   for (i=bP; i*i <= nelt; i += 4);
   bL = (i-4)>>1;
   tb = NULL;
   MMSplitByFlagAny((1<<MMF_KUISKB), &mb, &tb); /* tb gets fully unrolled K */
   WriteMMFile("res/tmp.sum", mb);
/*
 * Find best with KB near max for square problem
 */
   printf("   FINDING BEST ROLLED KERNEL FOR BOUND SEARCH, B=%d\n", bL);
   mp = bestNearSquare(pre, verb, mb, bL, 4, 0, 1.015);
   mf = mp->mflop[0];
   printf("   BEST ID=%d, '%s', B=(%d,%d,%d) mf=%.2f\n", mp->ID,
          GetMMLabelName(pre, mp), mp->mbB, mp->nbB, mp->kbB, mp->mflop[0]);
/*
 * Now find largest blocking in all dim that improves performance
 */
   printf("   FINDING BEST MAX BLOCK FACTORS, current=(%d,%d,%d):\n",
          mp->mbB, mp->nbB, mp->kbB);
   MMExpandMNK(verb, pre, 4, 0, nelt_llc, opWorkSetOK, 512, mp);
   MB = mp->mbB; NB = mp->nbB; KB = mp->kbB;
   printf("   BEST B=(%d,%d,%d), mf=%.2f, speedup=%.2f\n", MB, NB, KB,
          mp->mflop[0], mp->mflop[0] / mf);
/*
 * Now find best performance at around this size using any kernel
 */
   mb = ATL_JoinMMQs(mb, tb);
   printf("   FINDING BEST KERNEL FOR ROUGH BLOCK=(%d,%d,%d):\n", MB, NB, KB);
   mpB = bestNearMNK(pre, verb, mb, MB, NB, KB, 4, 0, 1.015);
   printf("   BEST ID=%d, '%s', B=(%d,%d,%d) mf=%.2f\n", mpB->ID,
          GetMMLabelName(pre, mpB), mpB->mbB, mpB->nbB, mpB->kbB,mpB->mflop[0]);
/*
 * If we changed kernel, re-search best block factors
 */
   if (!MMKernsSame(mpB, mp))
   {
      KillMMNode(mp);
      mp = mpB;
      printf("   FINDING BEST MAX BLOCK FACTORS, current=(%d,%d,%d):\n",
             mp->mbB, mp->nbB, mp->kbB);
      MMExpandMNK(verb, pre, 4, 0, nelt_llc, opWorkSetOK, 512, mp);
      MB = mp->mbB; NB = mp->nbB; KB = mp->kbB;
      printf("   BEST B=(%d,%d,%d), mf=%.2f, speedup=%.2f\n", MB, NB, KB,
             mp->mflop[0], mp->mflop[0] / mf);
      mpB = bestNearMNK(pre, verb, mb, MB, NB, KB, 4, 0, 1.015);
      KillMMNode(mpB);
   }
   else
      KillMMNode(mpB);

   mpB = mp;
   WriteMMFileWithPath(pre, "res", "opASYMP.sum", mpB);
   maxB = mp->mbB;
   maxB = Mmax(maxB, mp->nbB);
   maxB = Mmax(maxB, mp->kbB);
/*
 * Get temporary queue, eliminate kerns not providing unique KU speedup,
 * then do search for all problems between [bP,bL]
 */
   if (bP >= KB)
      bP = b0;
   if (bP < KB)
   {
      tb = CloneMMQueue(mb);
      tb = MMWinnowByKU(tb, bP, 1.05);
      if (FULLSRCH)
      {
         WriteMMFile("res/tmp.sum", mb);
         printf("FINDING BEST RANK-K KERNELS FOR LLC:\n");
         do_bfi(pre, bP+1, KB, 1, "res/tmp.sum", "opLLC.sum", 1.015);
         printf("DONE\n");
      }
      else
      {
         ATL_mmnode_t *rb;
         rb = opsrch_rng(verb, pre, 4, 1.015, bP+1, KB, maxB, nelt_llc,
                         opWorkSetOK, mpB, tb);
         WriteMMFileWithPath(pre, "res", "opLLC.sum", rb);
         KillAllMMNodes(rb);
      }
      KillAllMMNodes(tb);
   }
   else
      assert(0);
   KillMMNode(mpB);
/*
 * Finally, find best-performing kernels in LLPC region.  Do our winnowing
 * search on problem well-within LLPC
 */
   #if defined(LLPC_LVL) && LLPC_LVL > 1 && LLPC_LVL < LLC_LVL
   if (bP > b0)
   {
      KB = bP - b0;
      if (KB > 16)
         KB = bP-8;
      else
         KB = b0 + (KB>>1);
      printf("   FINDING BEST KERNEL FOR LLPC BOUND SEARCH, B=%d\n", KB);
      mp = bestNearSquare(pre, verb, mb, KB, 4, 0, 1.015);
      mf = mp->mflop[0];
      printf("   BEST ID=%d, '%s', B=(%d,%d,%d) mf=%.2f\n", mp->ID,
             GetMMLabelName(pre, mp), mp->mbB, mp->nbB, mp->kbB, mp->mflop[0]);
/*
 *    Now find largest blocking in all dim that improves performance
 */
      printf("   FINDING BEST MAX BLOCK FACTORS, current=(%d,%d,%d):\n",
             mp->mbB, mp->nbB, mp->kbB);
      MMExpandMN(verb, pre, 4, 0, maxB, mp);
      MB = mp->mbB; NB = mp->nbB; KB = mp->kbB;
      printf("   BEST B=(%d,%d,%d), mf=%.2f, speedup=%.2f\n", MB, NB, KB,
             mp->mflop[0], mp->mflop[0] / mf);
/*
 *    Find best-performing kernel around this size
 */
      printf("   FINDING BEST KERNEL FOR ROUGH BLOCK=(%d,%d,%d):\n", MB,NB,KB);
      mpB = bestNearMNK(pre, verb, mb, MB, NB, KB, 4, 0, 1.015);
      printf("   BEST ID=%d, '%s', B=(%d,%d,%d) mf=%.2f\n", mpB->ID,
             GetMMLabelName(pre, mpB),mpB->mbB,mpB->nbB,mpB->kbB,mpB->mflop[0]);
/*
 *     If we changed kernel, re-search best block factors
 */
      if (!MMKernsSame(mpB, mp))
      {
         KillMMNode(mp);
         mp = mpB;
         printf("   FINDING BEST MAX BLOCK FACTORS, current=(%d,%d,%d):\n",
                mp->mbB, mp->nbB, mp->kbB);
         MMExpandMN(verb, pre, 4, 0, maxB, mp);
         MB = mp->mbB; NB = mp->nbB; KB = mp->kbB;
         printf("   BEST B=(%d,%d,%d), mf=%.2f, speedup=%.2f\n", MB, NB, KB,
                mp->mflop[0], mp->mflop[0] / mf);
         mpB = bestNearMNK(pre, verb, mb, MB, NB, KB, 4, 0, 1.015);
         KillMMNode(mpB);
      }
      else
         KillMMNode(mpB);
      mpB = mp;
      mb = MMWinnowByKU(mb, b0, 1.10);
      if (FULLSRCH)
      {
         WriteMMFile("res/tmp.sum", mb);
         printf("FINDING BEST RANK-K KERNELS FOR LLPC:\n");
         do_bfi(pre, b0+1, bP, 1, "res/tmp.sum", "opLLPC.sum", 1.015);
         printf("DONE\n");
      }
      else
      {
         ATL_mmnode_t *rb;
         tb = CloneMMQueue(mb);
         tb = MMWinnowByKU(tb, bP, 1.05);
         nelt_llpc = (pre == 'd' || pre == 'z') ? (LLPC_SZ>>3) : (LLPC_SZ>>2);
         rb = opsrch_rng(verb, pre, 4, 1.015, b0+1, bP, maxB, nelt_llpc,
                         opWorkSetOK, mpB, tb);
         KillAllMMNodes(tb);
         WriteMMFileWithPath(pre, "res", "opLLPC.sum", rb);
         KillAllMMNodes(rb);
      }
      KillMMNode(mpB);
   }
   #endif
   KillAllMMNodes(mb);
   GetMMLabelName(pre, NULL);
   mb = ReadMMFileWithPath(pre, "res", "opL1.sum");
   assert(mb);
   mp = ATL_LastMMNode(mb);
   tb = ReadMMFileWithPath(pre, "res", "opLLPC.sum");
   if (tb)
   {
      mp->next = tb;
      mp = ATL_LastMMNode(tb);
   }
   tb = ReadMMFileWithPath(pre, "res", "opLLC.sum");
   mp->next = tb;
   WriteMMFileWithPath(pre, "res", "opFNL.sum", mb);
   KillAllMMNodes(mb);
}

/*
 * For small sizes, want to just try all kernels because there can be
 * unpredictable effects due to cleanup and related.  For larger problems,
 * we'd like to restrict the number of kernels searched a bit, based on the
 * idea that once FLOPS dominate, only raw performance should matter, and
 * so we can take a rough size in the middle of a "region" and only retain
 * kernels that are competitive at that size for all timings within the
 * region.  We have several regions of interest (greatest-to-least):
 * (1) Outer product working set fits in LLC, b <= (sqrt(1+LLC_elts/4) - 1)/2
 * (2) B fits in LLC: b <= sqrt(LLPC_elts) -> N <= LLPC_elts/K
 * (2) Outer produce wrk set fits in LLPC, b <= (sqrt(1+LLPC_elts/4) - 1)/2
 * we can eliminate kernels wt inadequate performance on a certain szi
 */
int main(int nargs, char **args)
{
   unsigned long llc_elts;
   unsigned eltsz = 8;
   int maxK, verb=0;
   char pre, upr;
   pre = GetFlags(nargs, args);
/*
 * For outer product, K is fixed, and you can vary at most 2 dimensions.
 * in practice, wt K fixed, you mainly want to increase NB to maximize
 * A reuse within a amm call.  Changing MB has smaller affects on performance,
 * because if we reduce MB, we reduce B reuse within one kernel, but increase
 * its's reuse over the M/MB amm calls that compute a colpan of C.  These
 * timings not needed for for degenerate M, which forces you to compute a
 * row-panel of C instead (and thus potentially have TLB problems).
 *
 * For this search, we break the possible block factors into several regions.
 * Minimally, consider a small region where operands fit in the L1, and
 * one where operands fit in the last-level of cache.  Another important
 * region last level private cache (cache that cores don't share).  Caches
 * that are shared between cores may be much more strongly bandwidth-limited
 * than private caches when all cores are running at maixmal rates, and so
 * this region may provide the best performance on some systems.
 *
 * In practice, since caches are often high associativity and random replacement
 * the boundaries between regions will not be distinct, and best performance
 * may be found when operands exceed the cache.
 */
/*
 * Biggest outer-produce case that makes sense is: 2MN + K*(MU+M+N) == LLC_elts
 * In this way, we should be able to resuse both Cb (C in block form) and
 * Bb from cache, despite moving both C & A ptrs, and copying Cb->Cc.  In
 * practice, the copy to col-major C will knock stuff out due to cache
 * conflicts, which may cause unpredictable affects.  If we set b=M=N=K, and
 * assume the best case of MU=1, the above becomes:
 *    4b^2 + b - (LLC_elts/4), or b = (sqrt(1+LLC_elts/4) - 1)/2 frm quad equ.
 */
   if (pre == 's' || pre == 's')
      eltsz = 4;
   maxK = (sqrt(1 + ((LLC_SZ / eltsz)>>2)) - 3.0)/2.0;
   maxK = Mmin(480, maxK);
   do_opk(verb, pre, NULL, maxK);
   do_opdek(verb, pre, 0);
   do_opdek(verb, pre, 1);
   return(0);
}
