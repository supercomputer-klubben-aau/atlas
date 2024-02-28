/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 */
#define ATL_WANT ILCM 1
#include "atlas_iopt.h"
#include "atlas_cache.h"
#include "atlas_genparse.h"
#include "atlas_mmtesttime.h"

#define TSM_RIGHT  0
#define TSM_UPPER  1
#define TSM_TRANSA 2
void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set type/precision prefix (d) \n");
   fprintf(stderr, "   -S [L/R] : search Left or Right TRSM\n");
   exit(ierr ? ierr : -1);
}

ATL_mmnode_t *weedOutBadTrsm(char pre, int flag, ATL_mmnode_t *mb)
/*
 * Deletes any mmkerns that can't be used for the TRSM case given in flag,
 * and returns ptr to base of modified queue (original queue changed!)
 */
{
   ATL_mmnode_t *mp = mb;
   while (mp)
   {
      ATL_mmnode_t *nxt = mp->next;
      int KILL;
      KILL = !(mp->flag & (1<<MMF_KRUNTIME));  /* kill compile-time K */
      KILL |= mp->flag & (1<<MMF_KUISKB);  /* kill fixed ku=KB kerns */
      KILL |= (mp->kbmin > mp->ku) | mp->kbmax; /* kill K-restricted kerns */
      if (flag&(1<<TSM_RIGHT))      /* Right case requires */
         KILL |= (mp->nu % mp->ku); /* nu a multiple of ku */
      else                          /* Left case requires */
         KILL |= (mp->mu % mp->ku); /* mu a multiple of ku */
      if (KILL)
         mb = KillMMNodeFromQ(mb, mp);
      else
      {
         mp->blask = ATL_KTRSM;
         mp->flag = (mp->flag&(~MMF_MVSET))|MMF_MVDEF;
      }
      mp = nxt;
   }
   return(mb);
}

ATL_mmnode_t *forceTrsmCases(char pre, int flag)
{
   ATL_mmnode_t *mb, *mp, *gb=NULL;
   char *fn;
/*
 * TRMM has the same restrictions as TRSM, and TRMM has been searched for
 * a working generated case, so read the TRMM kernel in, and change it back
 * to a GEMM
 */
   fn = (flag&(1<<TSM_RIGHT)) ? "trmmKRL.sum" : "trmmKLL.sum";

   gb = ReadMMFileWithPath(pre, "res", fn);
   if (!gb)
   {
      char ln[32];
      sprintf(ln, "make res/%c%s", pre, fn);
      assert(!system(ln));
      gb = ReadMMFileWithPath(pre, "res", fn);
      assert(gb);
   }
   assert(!gb->next); /* expect only 1 kern from trmm */
   gb->blask = ATL_KGEMM;
/*
 * Bring in all working gemm kernels since they may beat trmm for trsm usage
 */
   gb->next = ReadMMFileWithPath((pre=='s'||pre=='c')?'s':'d', "res",
                                 "WORKING.sum");
   mb = AddUniqueMMKernsToList(NULL, gb);
   KillAllMMNodes(gb);
   mb = weedOutBadTrsm(pre, flag, mb);
/*
 * Remove duplicate cases after forcing TRSM reqs, then make them TRSM kerns
 */
   return(mb);
}

ATL_mmnode_t *getAllCandKerns(char pre, int flag)
{
   ATL_mmnode_t *mb, *mp;
   mb = ReadMMFileWithPath(pre, "res", "ipmek.sum");
   mb = ATL_JoinMMQs(mb, ReadMMFileWithPath(pre, "res", "ipnek.sum"));
   mb = ATL_JoinMMQs(mb, ReadMMFileWithPath(pre, "res", "ipgen.sum"));
   mb = ATL_JoinMMQs(mb, ReadMMFileWithPath(pre, "res", "ipmen.sum"));
   mp = AddUniqueMMKernsToList(NULL, mb);
   KillAllMMNodes(mb);
   mb = weedOutBadTrsm(pre, flag, mp);
   if (!mb) /* no candidates! */
      mb = forceTrsmCases(pre, flag);
   assert(mb);
   return(mb);
}

char GetFlags(int nargs, char **args, int *FLG)
{
   char pre = 'd';
   int flg=0;
   int i;

   for (i=1; i < nargs; i++)
   {
      int wch, *ip, **ipp, TST=0;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'S':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        wch = args[i][0];
        if (wch == 'R' || wch == 'r')
           flg |= 1<<TSM_RIGHT;
        else
           flg &= ~(1<<TSM_RIGHT);
        break;
      case 'U':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        wch = args[i][0];
        if (wch == 'U' || wch == 'u')
           flg |= 1<<TSM_UPPER;
        else
           flg &= ~(1<<TSM_UPPER);
        break;
      case 'A':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        wch = args[i][0];
        if (wch == 'T' || wch == 't')
           flg |= 1<<TSM_TRANSA;
        else
           flg &= ~(1<<TSM_TRANSA);
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
   *FLG = flg;
   return(pre);
}

double timeTRSM(char pre, int flag, int B, int R, ATL_mmnode_t *mmp)
{
   double mfB = 0.0;
   int it;
   const int RIGHT=(flag&(1<<TSM_RIGHT));
   int tflag= ((RIGHT) ? 4:0);
   const char TA[2] = {'N', 'T'};

   if (flag&(1<<TSM_UPPER))
      tflag |= 16;
   if (flag&(1<<TSM_TRANSA))
      tflag |= 32;
   for (it=0; it < 2; it++)
   {
      const int tfl = tflag | ((it)?8:0);
      int DOIT, NR;
      const int TR = (RIGHT) ? !it : it;
      if (TR)
      {
         #if 0 /* code does not presently handle! */
         if (mmp->flag & (1<<MMF_KVEC))
            DOIT = 1;
         else
         #endif
            DOIT = (mmp->nu % mmp->ku) == 0;
         NR = (R > mmp->mu) ? (R/mmp->mu)*mmp->mu : R;
      }
      else
      {
         #if 0 /* code does not presently handle! */
         if (mmp->flag & (1<<MMF_KVEC))
            DOIT = 1;
         else
         #endif
            DOIT = (mmp->mu % mmp->ku) == 0;
         NR = (R > mmp->nu) ? (R/mmp->nu)*mmp->nu : R;
      }
      if (DOIT)
      {
         double mf;
         int blask = mmp->blask;
         int mvflags = mmp->flag;
         ATL_MMF_MVPUT(mmp->flag, 4);
         mmp->blask = ATL_KTRSM;
         mf = TimeMMKernel(0, tfl, mmp, pre, NR, NR, B, 1, 0, -1);
         mmp->blask = blask;
         mmp->flag = mvflags;
         printf("      ID=%d : B=%u, NRHS=%u U=(%2u,%2u,%2u) TRANS=%c, "
                "mf=%.2f\n", mmp->ID, B, R, mmp->mu, mmp->nu, mmp->ku,
                TA[it], mf);
         if (mf > mfB)
         {
            mmp->TA = mmp->TB = it ? AtlasTrans : AtlasNoTrans;
            mfB = mmp->mflop[0] = mf;
         }
      }
   }
   return(mfB);
}

ATL_mmnode_t *timeAllTRSM(char pre, int flag, int B, int R, ATL_mmnode_t *mb)
/*
 * Times all TRSM in mb with kb=B and NRHS roughly R, returns pointer to
 * cloned node of best performing
 */
{
   double mfB=0.0;
   ATL_mmnode_t *mp, *mmB=NULL;
   const int RIGHT=(flag&(1<<TSM_RIGHT));

      printf("   FINDING BEST TRSM KERNEL B=(%u,%u):\n", B, R);
   for (mp=mb; mp; mp = mp->next)
   {
      double mf;
      mf = timeTRSM(pre, flag, B, R, mp);
      if (mf > mfB)
      {
         mfB = mf;
         mmB = mp;
      }
   }
   if (!mmB)
   {
      printf("   NO CASE FOUND FOR B=(%u,%u)\n", B, R);
      return(NULL);
   }
   mmB = CloneMMNode(mmB);
   assert(mmB->mflop[0] == mfB);
   mmB->mbB = mmB->nbB = mmB->kbB = B;
   if (RIGHT)
      mmB->mbB = R;
   else
      mmB->nbB = R;

   printf(
      "   BEST CASE FOR B=(%u,%u): ID=%u, U=(%u,%u,%u), TRANS=%c, mf=%.2f\n",
          B, R, mmB->ID, mmB->mu, mmB->nu, mmB->ku,
          (mmB->TA == AtlasTrans)?'T':'N', mfB);
   return(mmB);
}

ATL_mmnode_t *findBestU(char pre, int flag, int B, int R, int U)
/*
 * Finds the best case with ku=1, and either mu OR nu = U
 */
{
   ATL_mmnode_t *mmb, *mmp, *mmB=NULL;
   int mb, nb, transB=0;
   double mfB=0.0;
   const int RIGHT=(flag&(1<<TSM_RIGHT));

   mmp = mmb = ReadMMFileWithPath(pre, "res", "WORKING.sum");
   while (mmp)
   {
      ATL_mmnode_t *nxt=mmp->next;
      if ((mmp->mu != U && mmp->nu != U) || (U%(mmp->ku)) ||
          !FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
         mmb = KillMMNodeFromQ(mmb, mmp);
      mmp = nxt;
   }
   printf("\nTIMING U=%u CASES FOR B=%u, NRHS=%u:\n", U, B, R);
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      double mf;
      int it;
      mmp->blask = ATL_KTRSM;
      mf = timeTRSM(pre, flag, B, R, mmp);
      if (mf > mfB)
      {
         mmB = mmp;
         mfB = mf;
         transB = (mmp->TA == AtlasTrans) ? 1:0;
      }
   }
   if (mmB)
   {
      mmb = RemoveMMNodeFromQ(mmb, mmB);
      mmB->mbB = mmB->kbB = mmB->nbB = B;
      if (RIGHT)
         mmB->mbB = R;
      else
         mmB->nbB = R;
      mmB->mflop[0] = mfB;
      mmB->TA = mmB->TB = transB;
      printf(
      "BEST U=%u CASE FOR B=%u,%u: ID=%u, U=(%u,%u,%u), TRANS=%c, mf=%.2f\n\n",
             U, B, R, mmB->ID, mmB->mu, mmB->nu, mmB->ku, transB?'T':'N', mfB);
   }
   else
      printf("NO USER CASES FOUND FOR U=%u.\n\n", U);
   KillAllMMNodes(mmb);
   return(mmB);
}

void findBestTRSM(char pre, int flag)
{
   ATL_mmnode_t *bb;   /* blocking base det by left or right setting */
   ATL_mmnode_t *tb=NULL; /* trsm best case queue, 1 for each entry in opmek */
   ATL_mmnode_t *kb;   /* list of candidate gemm kernels to try */
   ATL_mmnode_t *bp, *mp;
   int R=1, maxB=1, BN, BNN;
   double mf;
   char fnout[16];

   sprintf(fnout, "trsm%c_%c%c.sum", (flag&(1<<TSM_RIGHT))?'R':'L',
           (flag&(1<<TSM_UPPER))?'U':'L',(flag&(1<<TSM_TRANSA))?'T':'N');
   tb = TimeMMFileWithPath(pre, "res", fnout, 0, 1, 0, 1, 0, -1);
   if (tb)
   {
      KillAllMMNodes(tb);
      return;
   }

   kb = getAllCandKerns(pre, flag);
   PrintMMNodes(stdout, kb);
/*
 * Get NB that we will tune TRSM for
 */
   if (flag & (1<<TSM_RIGHT))
      bb = ReadMMFileWithPath(pre, "res", "opnek.sum");
   else
      bb = ReadMMFileWithPath(pre, "res", "opmek.sum");
/*
 * Compute RHS (R) to use, and find maxNB we'll see
 */
   for (mp=kb; mp; mp = mp->next)
   {
      R = ATL_iLCM(R, mp->nu);
      R = ATL_iLCM(R, mp->mu);
   }
   for (mp=bb; mp; mp = mp->next)
   {
      maxB = Mmax(maxB, mp->kbB);
      R = ATL_iLCM(R, mp->nu);
      R = ATL_iLCM(R, mp->mu);
   }
   if (maxB > R)
      R = (maxB/R)*R;  /* NRHS to time */
   else
      R = maxB;
   BN = bb->kbB;
/*
 * This is commented out, because it is bringing in very slow kerns
 * and I don't want to hassle with new kerns showing up due to TRSM.
 * Once install is working, we'll extend so it can cause new amm to be
 * added, then we need to add a version of this that adds kerns only
 * when the perform well on small problems!
 */
   #if 0
/*
 * Find best user U=[1-3] cases, so we have kerns that can do very small
 * block factors w/o huge performance loss
 */
   mp = findBestU(pre, flag, bb->kbB, R, 1);
   if (mp)
   {
      mp->next = kb;
      kb = mp;
      mf = mp->mflop[0];
   }
   else
      mf = 0.0;
   if (bb->next)
   {
      BN = bb->next->kbB;
      if (bb->next->next)
         BNN = bb->next->next->kbB;
      else
         BNN = BN;
   }
   else
      BNN = BN;
   mp = findBestU(pre, flag, BN, R, 2);
   if (mp)
   {
      if (mp->mflop[0] > mf)
      {
         mf = mp->mflop[0];
         mp->next = kb;
         kb = mp;
      }
      else
          KillMMNode(mp);
   }
   mp = findBestU(pre, flag, BNN, R, 3);
   if (mp)
   {
      if (mp->mflop[0] > mf)
      {
         mp->next = kb;
         kb = mp;
      }
      else
          KillMMNode(mp);
   }
   #endif
/*
 * Now look thru block sizes, and pick best-performing
 */
   printf("FINDING BEST TRSM KERNEL FOR %u BLOCK FACTORS:\n",
          ATL_CountNumberOfMMNodes(bb));
   for (bp=bb; bp; bp = bp->next)
   {
      BN = bp->kbB;
      double mfB = 0.0;
      ATL_mmnode_t *mmB=NULL;
      mmB = timeAllTRSM(pre, flag, BN, R, kb);
      assert(mmB);
      mmB->next = tb;
      tb = mmB;
   }
   tb = ReverseMMQ(tb);
   printf("\nBEST CASES:\n");
   PrintMMNodes(stdout, tb);
   WriteMMFileWithPath(pre, "res", fnout, tb);
   KillAllMMNodes(tb);
   KillAllMMNodes(bb);
   KillAllMMNodes(kb);
}

int main(int nargs, char **args)
{
   char pre;
   int flag;

   pre = GetFlags(nargs, args, &flag);
   findBestTRSM(pre, flag);
   return(0);
}
