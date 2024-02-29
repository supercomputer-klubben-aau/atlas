#include "atlas_cache.h"
#include "atlas_genparse.h"
#include "atlas_mmtesttime.h"
#include <math.h>
#include <stdio.h>

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
   i = sprintf(ln, "./xbfisrch -p %c -s p i R -i %s -o res/%c%s",
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
/*
 * This routine suggests the largest square block that will fit the working
 * set of the computation in the cache of CS (bytes).  It first translates
 * CS to NELT (bytes -> elements).  The working set of one iteration is
 * mu*nu + kb*(mu+nu).  Our typical kernel used loop order MNK, which means
 * both C and B will move in inner loop, so for safety double these sizes:
 *   2*mu*nu + kb*(mu+2*nu) <= NELT ==> KB = (NELT - 2*mu*nu)/(mu+2*nu)
 * We will solve this equation for all kernels, and return the smallest kb.
 */
int maxSquareNB(char pre, size_t CS, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mp;
   unsigned int eltsh = (pre == 'd') ? 3 : 2, kbB=8000;
   const size_t NELT = CS >> eltsh;
   assert(pre == 'd' || pre == 's');
   assert(mmb);
   for (mp=mmb; mp; mp = mp->next)
   {
      const unsigned int mu=mp->mu, nu=mp->nu;
      const unsigned kb = (NELT-((mu+mu)*nu))/(mu+nu+nu);
      if (kb < kbB)
         kbB = kb;
   }
   return(kbB);
}

int minSquareNB(char pre, size_t CS, ATL_mmnode_t *mmb)
/*
 * This routine returns largest NB that will be able to fully reuse all
 * operands between calls assuming LRU cache. We will therefore want to
 * contain all three blocks in the cache, plus the working set of next
 * call for safety.  This is: M*N + K*(M+N+MU+NU) <= NELT, now if D=M=N=K:
 *    D^2 + D^2 + D^2 + D(MU+NU) = NELT --->  3D^2 + D(MU+NU) - NELT = 0
 * so quadratic equation gives: D = -(MU+NU) + sqrt((MU+NU)^2 + 12*NELT)/6
 * Do not allow it > 256 for double, or 512 for single.
 */
{
   ATL_mmnode_t *mp;
   unsigned int eltsh = (pre == 'd') ? 3 : 2, mupnu, D;
   const size_t NELT = CS >> eltsh;
   assert(pre == 'd' || pre == 's');
   assert(mmb);
   mupnu = mmb->mu + mmb->nu;
   for (mp=mmb; mp; mp = mp->next)
   {
      unsigned int mn=mp->mu + mp->nu;
      mupnu = (mn > mupnu) ? mn : mupnu;
   }
   D = sqrt(mupnu*mupnu + 3.0*(NELT<<2)) / 6.0;
   D -= mupnu;
   if (pre == 'd' || pre == 'z')
      D = Mmin(256,D);
   else
      D = Mmin(512,D);
   return(D);
}

void do_ipdek(int verb, char pre, int NEK)
{
   ATL_mmnode_t *mp, *mb=NULL;
   int tflag = 0;
   unsigned int i, U, N, ADJKUKB=0;
   const char ch=NEK?'N':'M';
   char *fn=NEK?"ipnek.sum":"ipmek.sum";

   mb = TimeMMFileWithPath(pre, "res", fn, 0, verb|1, 0, 0, 0, -1);
   if (mb)
   {
      KillAllMMNodes(mb);
      return;
   }
   mp = ReadMMFileWithPath(pre, "res", "ipgen.sum");
   assert(mp);
   while (mp->next)
      mp = KillMMNode(mp);
   mp->flag = (mp->flag&(~MMF_MVSET)) | ((1<<MMF_MVA)|(1<<MMF_MVB));
   if (mp->flag & (1<<MMF_KUISKB))
   {
      if (!mp->ID)
         ADJKUKB = 1;
      else if (!mp->kbmax || mp->kbmax > mp->kbmin)
         ADJKUKB = 1;
      else
        assert(0);  /* need gen a file/retime user here */
   }
   if (NEK)
      tflag |= 1<<31;
   mp->mflop[0] = 0;
   assert(mp);
   printf("FINDING BEST CASE FOR %c=K INNER PRODUCT:\n", ch);
   MMExpandNK(verb, pre, tflag, 1, mp);
   printf("BEST %c=K: ID=%u, B=(%u,%u,%u), mf=%.2f\n\n", ch, mp->ID,
          mp->mbB, mp->nbB, mp->kbB, mp->mflop[0]);
   U = (NEK) ? mp->nu : mp->mu;
   if (ADJKUKB)
   {
      if (U < mp->kbmin)
         U = ((mp->kbmin+U-1)/U)*U;
      N = (NEK) ? mp->nbB : mp->mbB;
   }
   else if (NEK)
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
   fprintf(stderr, "U=%d mp->kbmax=%d", U, mp->kbmax);
   if (mp->kbmax)
      assert(U <= mp->kbmax);
   for (i=U; i < N; i += U)
   {
      ATL_mmnode_t *p;
      unsigned int b;
      if (mp->kbmin && i < mp->kbmin)
          continue;
      if (mp->kbmax && i > mp->kbmax)
         break;
      printf("   FINDING BEST CASE FOR %c=K=%u:\n", ch, i);
      p = CloneMMNode(mp);
      if (NEK)
         p->nbB = p->kbB = i;
      else
         p->mbB = p->kbB = i;
      if (ADJKUKB)
         p->ku = i;
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

void findBestIP(char pre, int verb, ATL_mmnode_t *mmb, char *fnout)
/*
 * This search simply finds best performing kernel, produces ip.sum
 * with best found kernel and blocking parameter
 */
{
   ATL_mmnode_t *tb, *mp, *mpB, *mpU, *mpD;
   double mfB, mf, t0;
   unsigned int mbB, nbB, kbB, maxD=240, D, iD, iD0;
   char upr = pre;
   if (pre == 'z')
      upr = 'd';
   else if (pre == 'c')
      upr = 's';

   mp = TimeMMFileWithPath(pre, "res", "ipgen.sum", 0, 1, 0, 1, 0, -1);
   if (mp)
   {
      KillAllMMNodes(mmb);
      KillAllMMNodes(mp);
      return;
   }
   nbB = minSquareNB(upr, LLC_SZ, mmb);
   printf("bLLC = %u\n", nbB);
   if (nbB >= 132)
      nbB -= 12;
/*
 * Get rid of all kernels that aren't within 60% of max performance
 */
   D = ATL_CountNumberOfMMNodes(mmb);
   printf("ELIMINATING UNCOMPETITIVE KERNELS out of %u:\n", D);
   mmb = ElimSlowKern(pre, verb, mmb, nbB, nbB, nbB, 0, 0, 0.60);
   iD = ATL_CountNumberOfMMNodes(mmb);
   printf("DONE: ELIMINATED %u UNCOMPETITIVE KERNELS, %u LEFT.\n",
          D-iD, iD);
   PrintMMNodes(stdout, mmb);
   printf("LOOKING FOR BEST ASYMPTOTIC INNER-PRODUCT KERNEL:\n");
   mpB = mp = FindMaxMflopMMQ(mmb, 0);
   mfB = mp->mflop[0];
   mbB = mp->mbB;
   nbB = mp->nbB;
   kbB = mp->kbB;
   printf("\nBEST KERNEL AT START:\n");
   PrintMMLine(stdout, mpB);
   for (mp=mmb; mp; mp = mp->next)
   {
      mf = mp->mflop[0];
      if (mf*1.2 >= mfB)  /* try all competitive kerns */
      {
         double mfN;
         tb = mp->next;
         printf("   TUNE BLK for %s, current=(%u,%u,%u), mf=%.2f\n",
                GetMMLabelName(pre, mp), mp->mbB, mp->nbB, mp->kbB,
                mp->mflop[0]);
         mfN = MMExpandMNK(verb, pre, 0, 0, 0, NULL, maxD, mp);
         if (mfN > mfB)
         {
            printf(
               "   BEST SO FAR: (%u,%u,%u), mf=%.2f, spup=%.2f, spupB=%.2f\n",
                   mp->mbB, mp->nbB, mp->kbB, mfN, mfN/mf, mfN/mfB);
            PrintMMLine(stdout, mp);
            mpB = mp;
            mfB = mp->mflop[0];
            mbB = mp->mbB;
            nbB = mp->nbB;
            kbB = mp->kbB;
         }
         else
            printf("   NOT BEST: (%u,%u,%u), mf=%.2f, spdup=%.2f\n",
                   mp->mbB, mp->nbB, mp->kbB, mfN, mfN/mfB);
      }
      else
         printf("Not timing ID=%d %s (%u,%u,%u): mf=%.2f, mfB=%.2f\n",
                mp->ID, GetMMLabelName(pre, mp), mp->mbB, mp->nbB, mp->kbB,
                mf, mfB);
   }
   printf("CHOSEN: ID=%d %s (%u,%u,%u), mf=%.2f\n\n", mpB->ID,
           GetMMLabelName(pre, mpB), mbB, nbB, kbB, mfB);
   PrintMMLine(stdout, mpB);
   mpD = CloneMMNode(mpB);
   KillAllMMNodes(mmb);
   mpD->mflop[0] = mfB;
   mpD->mbB = mbB;
   mpD->nbB = nbB;
   mpD->kbB = kbB;
   printf("CHOSEN2: ID=%d %s (%u,%u,%u), mf=%.2f\n\n", mpD->ID,
          GetMMLabelName(pre, mpD), mpD->mbB,mpD->nbB,mpD->kbB, mpD->mflop[0]);
   mpU = CloneMMNode(mpD);
   printf("FINDING LARGE BLOCK FACTORS FOR CHOSEN KERNEL:\n");
   mf = MMExpandMNK(verb, pre, 0, 0, 0, NULL, 480, mpU);
   if (mf < mfB)
   {
      KillMMNode(mpU);
      mpU = CloneMMNode(mpD);
   }
   printf("LRGBLK CHOSEN: ID=%d %s (%u,%u,%u), mf=%.2f\n\n", mpU->ID,
           GetMMLabelName(pre, mpU), mpU->mbB, mpU->nbB, mpU->kbB, mf);

   PrintMMLine(stdout, mpU);
   iD0 = ATL_iLCM(mpU->mu, mpU->nu);
   maxD = Mmax(mpU->mbB, mpU->nbB);
   iD = maxD / iD0;
   if (iD*iD0 == maxD)
      maxD -= iD0;
   else
      maxD = iD*iD0;
   for (iD=iD0; iD < 24; iD += iD);
   printf("LOOKING FOR BEST KERNS WITH RESTRICTED M&N <= %u, STEPS=%u:\n",
          maxD, iD);
   mmb = mpU;
   mmb->flag |= (1<<MMF_KCLN);
   for (D=maxD; D >= iD; D -= iD)
   {
      printf("   LOOKING FOR BEST KERNS WITH M&N < %d:\n", D);
      mp = CloneMMNode(mpD);
      mf = MMExpandMN_K(verb, pre, 0, 0, 0, NULL, D, mp);
      mp->flag |= (1<<MMF_KCLN);
      mp->next = mmb;
      mmb = mp;
      printf("   FOUND ID=%d : B=(%d,%d,%d) mf=%.2f\n",  mp->ID, mp->mbB,
             mp->nbB, mp->kbB, mp->mflop[0]);
      PrintMMLine(stdout, mp);
   }
   if (iD != iD0)
   {
      for (D=D+iD-iD0; D >= iD0; D -= iD0)
      {
         printf("   LOOKING FOR BEST KERNS WITH M&N < %d:\n", D);
         mp = CloneMMNode(mpD);
         mf = MMExpandMN_K(verb, pre, 0, 0, 0, NULL, D, mp);
         mp->flag |= (1<<MMF_KCLN);
         mp->next = mmb;
         mmb = mp;
         printf("   FOUND ID=%d : B=(%d,%d,%d) mf=%.2f\n",  mp->ID, mp->mbB,
                mp->nbB, mp->kbB, mp->mflop[0]);
         PrintMMLine(stdout, mp);
      }
      if (iD0 > 16) /* didn't get small enough with LCM! */
      {
         iD = Mmax(mpU->mu, mpU->nu);
         for (D=iD0-iD; D >= iD; D -= iD)
         {
            printf("   LOOKING FOR BEST KERNS WITH M&N < %d:\n", D);
            mp = CloneMMNode(mpD);
            mf = MMExpandMN_K(verb, pre, 0, 0, 0, NULL, D, mp);
            mp->flag |= (1<<MMF_KCLN);
            mp->next = mmb;
            mmb = mp;
            printf("   FOUND ID=%d : B=(%d,%d,%d) mf=%.2f\n",  mp->ID, mp->mbB,
                   mp->nbB, mp->kbB, mp->mflop[0]);
            PrintMMLine(stdout, mp);
         }
      }
   }
   printf("DONE LOOKING FOR RESTRICTED D KERNELS.\n");
   MMWinnowByOrder(mmb, 0, 1.0);
   WriteMMFileWithPath(pre, "res", "ipgen.sum", mmb);

   KillAllMMNodes(mmb);
}

ATL_mmnode_t *bestUM(char pre, int verb, ATL_mmnode_t *IP)
   /* IP: fastest inner product kernel with mu not mult of nu */
{
   ATL_mmnode_t *bp, *prv, *mp, *bestp=NULL;
   const unsigned int MB=IP->mbB, NB=IP->nbB, KB=IP->kbB;
   double mfB=0.0;
   mp = TimeMMFileWithPath(pre, "res", "ipbestUM.sum", 0, 1, 0, 1, 0, -1);
   if (mp)
      return(mp);
/*
 * Get all working user cases, and remove any where M/N unroll not multiple or
 * K isn't runtime <= 4;
 */
   bp = ReadMMFileWithPath(pre, "res", "gSYRKUM.sum");
   assert(bp);
   bp->next = GetWorkingUserCases(verb, pre);
   prv = bp;
   mp = bp->next;
   while(mp)
   {
      ATL_mmnode_t *nxt=mp->next;
      const unsigned int mu=mp->mu, nu=mp->nu,
         maxKU=(FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen : 4);
      if (mu%nu && nu%mu || mp->ku > maxKU ||
          !FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
      {
         prv->next = nxt;
         KillMMNode(mp);
      }
      else
         prv = mp;
      mp = nxt;
   }
/*
 * Now, find fastest asymptotic sized kernel
 */
   printf("FINDING FASTEST AMM WITH MU/NU THAT DIVIDE EVENLY:\n");
   for (mp=bp; mp; mp = mp->next)
   {
      double mf;
      const unsigned int mu=mp->mu, nu=mp->nu, ku=mp->ku;
      unsigned int mb=(MB/mu)*mu, nb=(NB/nu)*nu, kb=(KB/ku)/ku;
      mb = (mb) ? mb : mu;
      nb = (nb) ? nb : nu;
      kb = (kb) ? kb : ku;
      mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, kb, 1, -1, -1);
      printf("   ID=%u: B=(%u,%u,%u), U=(%u,%u,%u), mf=%.2f\n", mp->ID,
             mb, nb, kb, mu, nu, ku, mf);
      if (mf > mfB)
      {
         mp->mbB = mb;
         mp->nbB = nb;
         mp->nbB = kb;
         bestp = mp;
         mfB = mf;
      }
   }
   assert(bestp);
   bestp->mflop[0] = mfB;
   mp = bestp;
   printf("CHOSE: ID=%u: B=(%u,%u,%u), mf=%.2f\n", mp->ID, mp->mbB, mp->nbB,
          mp->kbB, mfB);
   bp = RemoveMMNodeFromQ(bp, bestp);
   KillAllMMRouts(bp);
   WriteMMFileWithPath(pre, "res", "ipbestUM.sum", bestp);
   return(bestp);
}

void srch_menUM(char pre, int verb)
/*
 * Finds kernels for cases with M=N at fixed values, while K is tuned.
 * Used in SYRK.
 */
{
   ATL_mmnode_t *bp, *mmb=NULL;
   unsigned int U, D, d, maxK, maxKU;

   bp = TimeMMFileWithPath(pre, "res", "ipmenUM.sum", 0, 1, 0, 1, 0, -1);
   if (bp)
   {
      mmb = bp;
      bp = TimeMMFileWithPath(pre, "res", "ipsyrkUM.sum", 0, 1, 0, 1, 0, -1);
      if (!bp)
         goto TIME_SYRKUM;
      KillAllMMNodes(bp);
      KillAllMMNodes(mmb);
      return;
   }
   bp = ReadMMFileWithPath(pre, "res", "ipmen.sum");
   assert(bp);
   while (bp->next)
      bp = KillMMNode(bp);
   if (bp->mu % bp->nu && bp->nu % bp->mu)
   {
      ATL_mmnode_t *tp=bp;
      bp = bestUM(pre, verb, bp);
      KillMMNode(tp);
   }
   U = ATL_iLCM(bp->mu, bp->nu);
   D = Mmax(bp->mbB, bp->nbB);
   maxK = Mmax(bp->kbB, D);
   if (bp->kbmax)
      maxK = Mmin(maxK, bp->kbmax);
   for (d=U; d <= D; d += U)
   {
      ATL_mmnode_t *mp;
      mp = CloneMMNode(bp);
      mp->next = mmb;
      mmb = mp;
      mp->mbB = mp->nbB = d;
      MMExpandK(verb, pre, 0, 1, mp, maxK);
      printf("\n");
   }
   KillMMNode(bp);
/*   PrintMMNodes(stdout, mmb); */
   while (mmb->next && mmb->mflop[0] < mmb->next->mflop[0])
      mmb = KillMMNode(mmb);
   mmb = ReverseMMQ(mmb);
   MMPruneMflopTol(mmb, 0, 1.0);
   WriteMMFileWithPath(pre, "res", "ipmenUM.sum", mmb);
/*
 * Now, create a SYRK kernel with its own time in mflop[0], and gemm's in
 * mflop[1], and write it it to "ipsyrkUM.sum"
 */
TIME_SYRKUM:
   while (mmb->next)                      /* get rid of all but */
      mmb = KillMMNode(mmb);              /* asymptotic case */
   maxKU = FLAG_IS_SET(mmb->flag, MMF_KVEC) ? mmb->vlen : 4;
/*
 * Try generating with straight params from search
 */
   bp = MMGetNodeGEN(pre, 0, 0, mmb->mu, mmb->nu, mmb->ku, mmb->vlen,
                     FLAG_IS_SET(mmb->flag, MMF_KVEC), mmb->pref, mmb->pfLS,
                     NULL);
   bp->flag = ((mmb->flag | (1<<MMF_KRUNTIME)) &
               (~((1<<MMF_X87)|(1<<MMF_KUISKB))));
   bp->mbB = mmb->mbB;
   bp->nbB = mmb->nbB;
   bp->kbB = mmb->kbB;
   if (bp->ku > maxKU)
   {
      if (FLAG_IS_SET(bp->flag, MMF_KVEC))
         bp->ku = maxKU;
      else
      {
         unsigned int ku=maxKU, kb=mmb->kbB;
         while((kb/ku)*ku != kb)
            ku--;
         bp->ku = ku;
      }
   }
   bp->mflop[1] = mmb->mflop[0];
   bp->blask = ATL_KSYRK;
   if (bp->genstr)
      free(bp->genstr);
   bp->genstr = NULL;
/*
 * Duplicate this node so we can try BREG1 & NOBCAST
 */
   KillMMNode(mmb);
   mmb = CloneMMNode(bp);
   mmb->flag |= (1<<MMF_BREG1);
   bp->mflop[0] = TimeMMKernel(verb, 0, bp, pre, bp->mbB, bp->nbB, bp->kbB,
                               1, 0, -1);
   printf("   SYRK/GEMM DEFAULT MFLOP=%.4f\n", bp->mflop[0]/bp->mflop[1]);
   mmb->mflop[0] = TimeMMKernel(verb, 0, mmb, pre, bp->mbB, bp->nbB, bp->kbB,
                                1, 0, -1);
   printf("   SYRK/GEMM BREG1   MFLOP=%.4f\n", mmb->mflop[0]/mmb->mflop[1]);
/*
 * bp is best found case so far
 */
   if (mmb->mflop[0] >= bp->mflop[0])
   {
      ATL_mmnode_t *tp=mmb;
      mmb = bp;
      bp = tp;
   }
/*
 * If nu mult of vlen, try NOBCAST
 */
   if (!bp->vlen || (bp->nu % bp->vlen) == 0)
   {
      mmb->flag = (mmb->flag | (1<<MMF_NOBCAST)) & (~(1<<MMF_BREG1));
      mmb->mflop[0] = TimeMMKernel(verb, 0, mmb, pre, bp->mbB, bp->nbB, bp->kbB,
                                   1, 0, -1);
      printf("   SYRK/GEMM NOBCAST MFLOP=%.4f\n", mmb->mflop[0]/mmb->mflop[1]);
      if (mmb->mflop[0] >= bp->mflop[0])
      {
         ATL_mmnode_t *tp=mmb;
         mmb = bp;
         bp = tp;
      }
   }
   KillMMNode(mmb);
   printf("SYRK/GEMM MFLOP=%.4f\n", bp->mflop[0]/bp->mflop[1]);
   if (pre == 'c' || pre == 'z')
      bp->flag |= 1<<MMF_COMPLEX;
   WriteMMFileWithPath(pre, "res", "ipsyrkUM.sum", bp);
   KillMMNode(bp);
}

void srch_men(char pre, int verb)
/*
 * Finds kernels for cases with M=N at fixed values, while K is tuned.
 * Used in SYRK.
 */
{
   ATL_mmnode_t *bp, *mmb=NULL;
   unsigned int U, D, d, maxK;

   bp = TimeMMFileWithPath(pre, "res", "ipmen.sum", 0, 1, 0, 1, 0, -1);
   if (bp)
   {
      KillAllMMNodes(bp);
      return;
   }
   bp = ReadMMFileWithPath(pre, "res", "ipgen.sum");
   assert(bp);
   while (bp->next)
      bp = KillMMNode(bp);
   U = ATL_iLCM(bp->mu, bp->nu);
   D = Mmax(bp->mbB, bp->nbB);
   maxK = Mmax(bp->kbB, D);
   for (d=U; d <= D; d += U)
   {
      ATL_mmnode_t *mp;
      mp = CloneMMNode(bp);
      mp->next = mmb;
      mmb = mp;
      mp->mbB = mp->nbB = d;
      MMExpandK(verb, pre, 0, 1, mp, maxK);
      printf("\n");
   }
   KillMMNode(bp);
/*   PrintMMNodes(stdout, mmb); */
   while (mmb->next && mmb->mflop[0] < mmb->next->mflop[0])
      mmb = KillMMNode(mmb);
   mmb = ReverseMMQ(mmb);
   MMPruneMflopTol(mmb, 0, 1.0);
   WriteMMFileWithPath(pre, "res", "ipmen.sum", mmb);
   KillAllMMNodes(mmb);
}

int main(int nargs, char **args)
{
   ATL_mmnode_t *mmb, *mp;
   int verb=0, nkern0, nkern, flag=0;
   unsigned int bLLC;
   char pre, upr;
   upr = pre = GetFlags(nargs, args);
   if (pre == 'z')
      upr = 'd';
   else if (pre == 'c')
      upr = 's';

/*
 * Get user & gen cases, and try compile- and run-time versions
 */
   mmb = GetWorkingUserCases(verb, pre);
   mmb = MMApplyRules(mmb, flag, 0, 0);
   mp = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
   mp = MMApplyRules(mp, flag|(1<<MMR_MKCOMP), 0, 0);
   assert(mp);
   mmb = AddUniqueMMKernCompList(mmb, mp);
   KillAllMMNodes(mp);
   flag = (verb > 3) ? 3:(verb&3);
   flag |= 3 << 2;

   for (mp=mmb; mp; mp = mp->next)
      mp->flag = ((mp->flag) & (~MMF_MVSET)) | ((1<<MMF_MVA) | (1<<MMF_MVB));

   findBestIP(pre, verb, mmb, NULL);
   srch_men(pre, verb);
   srch_menUM(pre, verb);
   do_ipdek(verb, pre, 0);
   do_ipdek(verb, pre, 1);

   return(0);
}
