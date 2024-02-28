#define ATL_GETFLAGS
#include "atlas_genparse.h"
#include "atlas_mmtesttime.h"

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set type/precision prefix (d) \n");
   fprintf(stderr, "   -o <resfile>: (stdout) fastest kernels found\n");
   fprintf(stderr,
   "   -O : same as -o, but deletes any existing file with <resfile> name\n");
   fprintf(stderr, "   -i <input file>: list of kernel to consider\n");
   fprintf(stderr, "   -r <runBonus> : (1.0) timing bonus to give runtime-K\n");
   fprintf(stderr, "   -b[m,n,k] # b1 ... b#: specify desired blocking for"
                   " given dim as list\n");
   fprintf(stderr, "   -B[m,n,k] B0 BN incB: specify desired blocking for\n"
                   "                         given dim as inclusive range\n");
   fprintf(stderr, "   -s [f,p] [O,o,i] [C,R]: specify search type:\n");
   fprintf(stderr, "      p: precise search, performance may be reduced"
                   " by dims\n");
   fprintf(stderr, "      f: fuzzy search, storage may exceed dims\n");
   fprintf(stderr, "      O: outer product including C copy\n");
   fprintf(stderr, "      o: outer product gemm-time only\n");
   fprintf(stderr, "      i: inner product gemm-time only\n");
   fprintf(stderr, "      C: try timing runtime-K C files as compile-time K\n");
   fprintf(stderr, "      R: only time runtime-K as runtime\n");
   fprintf(stderr, "   -MNK <idx>: take kern idx frm file fnin &"
                   " tune all dims\n");
   fprintf(stderr, "   -MN <idx> <KB>: take kern idx frm fnin &"
                   " tune MB&NB\n");
   fprintf(stderr, "   -N <idx> <KB> <MB> <maxN>\n");
   fprintf(stderr, "   -e [exclusion clause]: exclude kernels from search:\n");
   fprintf(stderr, "      V=[K,M] : exclude K- or M-vectorized kernels\n");
   fprintf(stderr, "      V=# : exclude all kerns whose vlen != #\n");
   fprintf(stderr, "      K=C : exclude all compile-time K kernels\n");
   fprintf(stderr, "      K=F : exclude fully unrolled K-loop kernels\n");
   fprintf(stderr, "      K=# : exclude all kernels with KU >= #\n");
   fprintf(stderr, "      repeat -e to exclude multiple kern types\n");
   fprintf(stderr, "   -v <verb> : set verbosity (1)\n");
   exit(ierr ? ierr : -1);
}

#define FLG_FUZZY   2  /* 1st two bit pos are for verb */
#define FLG_OPCPY   3  /* outer-product search including C put time */
#define FLG_OP      4  /* if these two not set, use inner-prod search */
#define FLG_NOKVEC  5  /* don't use k-vectorized kernels */
#define FLG_NOMVEC  6  /* don't use m-vectorized kernels */
#define FLG_NOKCOMP 7  /* don't use kernels needing compile-time K */
#define FLG_NOKFULL 8  /* don't use fully-unrolled K kernels */
#define FLG_NOBIGKU 9  /* exclude kernels with KU > III */
#define FLG_MKCOMP  10 /* try compiling runtime-K as compile-time */
#define ALLSRCH ( (1<<FLG_EXPMN)|(1<<FLG_EXPMNK) )
#define FLG_EXPMN  30  /* do ExpandMN search on 1 kern */
#define FLG_EXPMNK 31  /* do ExpandMNK search on 1 kern */
char *GetFlags(int nargs, char **args, char *PRE, int *FLAG, int *MAXKU,
               int *REQVL, char **FNIN, int **MBS, int **NBS, int **KBS,
               double *runBon)
/*
 * RETURNS: name of output file as malloced string
 */
{
   char *fout=NULL;
   int *mbs=NULL, *nbs=NULL, *kbs=NULL;
   int i, flag, verb=0, FUZZY=1, OPCPY=1, OP=0, MKCOMP=0;
   int NOKVEC=0, NOMVEC=0, NOKCOMP=0, NOKFULL=0;
   char pre='d';

   *REQVL = *MAXKU = 0;
   *FNIN = NULL;
   *runBon = 1.0;
   for (i=1; i < nargs; i++)
   {
      int wch, *ip, **ipp, TST=0;
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
      case 'e': /* -e V=[K,M,#] or  K=[C,F,#] */
         if (++i >= nargs)
             PrintUsage(args[0], i-1, NULL);
         wch = args[i][0];
         if (wch != 'V' && wch != 'K')
            PrintUsage(args[0], i-1, "expecting 'V' or 'K' after -e");
         if (args[i][1] != '=')
            PrintUsage(args[0], i-1, "expecting '=' after -e V/K");
         if (wch == 'V')
         {
            wch = args[i][2];
            if (wch == 'K')
               NOKVEC = 1;
            else if (wch == 'M')
               NOMVEC = 1;
            else if (isdigit(wch))
               *REQVL = atoi(args[i]+2);
            else
               PrintUsage(args[0], i-1, "expecting M,K,# after -e V=");
         }
         else /* wch == 'K' K=[C,F,#] */
         {
            wch = args[i][2];
            if (wch == 'F')
               NOKFULL = 1;
            else if (wch == 'C')
               NOKCOMP = 1;
            else if (isdigit(wch))
               *MAXKU = atoi(args[i]+2);
            else
               PrintUsage(args[0], i-1, "expecting F,C,# after -e K=");
         }
         break;
      case 'i':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *FNIN = args[i];
        break;
      case 'v':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        verb = atoi(args[i]);
        verb = (verb > 3) ? 3 : verb;
        break;
      case 'O':
        TST=1;
      case 'o':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        fout = args[i];
        if (TST)
           remove(fout);
        break;
      case 's': /* -s [f,p] [O,o,i] [C,R] */
         if (++i >= nargs)
             PrintUsage(args[0], i-1, NULL);
         wch = args[i][0];
         if (wch != 'f' && wch != 'F' && wch != 'p' && wch != 'P')
             PrintUsage(args[0], i-1, "bad fuzzy/precice to -s");
         if (wch == 'p' || wch == 'P')
            FUZZY = 0;
         if (++i >= nargs)
             PrintUsage(args[0], i-1, NULL);
         wch = args[i][0];
         if (wch != 'O' && wch != 'o' && wch != 'i' && wch != 'I')
             PrintUsage(args[0], i-1, "bad inner/outer to -s");
         OPCPY = (wch == 'O');
         OP = (wch == 'o');
         if (++i >= nargs)
             PrintUsage(args[0], i-1, NULL);
         wch = args[i][0];
         if (wch != 'C' && wch != 'c' && wch != 'r' && wch != 'R')
             PrintUsage(args[0], i-1, "bad runcomp to -s");
         MKCOMP = (wch == 'C' || wch == 'c');
         break;
      case 'r':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *runBon = atof(args[i]);
        break;
      case 'M':  /* -MN idx KB or -MNK idx */
         if (args[i][2] != 'N')
            PrintUsage(args[0], i, args[i]);
         if (args[i][3] == 'K')
            TST=1;
         else if (args[i][3] == '\0')
            TST=2;
         else
            PrintUsage(args[0], i, args[i]);
         nbs = malloc(TST*sizeof(int));
         assert(nbs);
         for (wch=0; wch < TST; wch++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            nbs[wch] = atoi(args[i]);
         }
         if (TST == 1) /* for full search, put info in MBs not NBs */
         {
            mbs = nbs;
            nbs = NULL;
         }
         break;
      case 'N':
         nbs = malloc(4*sizeof(int));
         assert(nbs);
         for (wch=0; wch < 4; wch++)
         {
            if (++i >= nargs)
                PrintUsage(args[0], i-1, NULL);
            nbs[wch] = atoi(args[i]);
         }
         break;
      case 'B':                                 /* -B[m,n,k] B0, BN, incB */
         TST=1;
      case 'b':                                 /* -b[m,n,k] # b1 ... b# */
         wch = tolower(args[i][2]);
         if (wch == 'm')
            ipp = &mbs;
         else if (wch == 'n')
            ipp = &nbs;
         else if (wch == 'k')
            ipp = &kbs;
         else
         {
            fprintf(stderr, "Expecting dim of m n or k!");
            PrintUsage(args[0], i, args[i]);
         }
         if (TST)
         {
            int B0, BN, incB;
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            B0 = atoi(args[i]);
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            BN = atoi(args[i]);
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            incB = atoi(args[i]);
            ip = GF_IntRange2IntList(B0, BN, incB);
         }
         else
         {
            ip = GF_GetIntList(nargs, args, i, 1);
            i += *ip + 1;
         }
         if (*ipp)
            free(*ipp);
         *ipp = ip;
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   flag = 0;
   flag |= MKCOMP<<FLG_MKCOMP;
   flag |= OP<<FLG_OP;
   flag |= OPCPY<<FLG_OPCPY;
   flag |= FUZZY<<FLG_FUZZY;
   flag |= NOKFULL<<FLG_NOKFULL;
   flag |= NOKCOMP<<FLG_NOKCOMP;
   flag |= NOMVEC<<FLG_NOMVEC;
   flag |= NOKVEC<<FLG_NOKVEC;
   *FLAG = flag;
   if (fout)
      fout = DupString(fout);
   *PRE = pre;
   *MBS = mbs;
   *NBS = nbs;
   *KBS = kbs;
   return(fout);
}

ATL_mmnode_t *ApplyRules(ATL_mmnode_t *ob, int flag, int maxKU, int reqVL)
{
   ATL_mmnode_t *nb=NULL;
   while(ob)
   {
      int KILL=0, REGEN=0, KILLFIX=0;
      if (reqVL)
         KILL = (ob->vlen != reqVL);

      if (!KILL && maxKU && ob->ku > maxKU)
      {
         REGEN = KILL = 1;
         if (!ob->ID) /* for genned codes, see if we can reset KU */
         {
            if (FLAG_IS_SET(ob->flag, MMF_KVEC))
            {
               if (maxKU >= ob->vlen)
               {
                  KILL = 0;
                  ob->ku = ob->vlen;
               }
            }
            else
            {
               KILL = 0;
               ob->ku = 1;
            }
         }
      }
      if (!KILL&&(flag&(1<<FLG_NOKFULL))&&FLAG_IS_SET(ob->flag, MMF_KUISKB))
      {
         REGEN = KILL = 1;
         if (!ob->ID) /* genned codes poss changed to runtime */
         {
            KILL = 0;
            if (FLAG_IS_SET(ob->flag, MMF_KVEC))
            {
               ob->ku = ob->vlen;
               ob->flag ^= 1<<MMF_KUISKB;
            }
            else
            {
               ob->ku = 1;
               ob->flag ^= 1<<MMF_KUISKB;
            }
         }
      }
      if (!KILL&&(flag&(1<<FLG_NOKCOMP))&&!FLAG_IS_SET(ob->flag, MMF_KRUNTIME))
      {
         REGEN = KILL = 1;
         if (!ob->ID) /* genned codes poss changed to runtime */
         {
            KILL = 0;
            if (FLAG_IS_SET(ob->flag, MMF_KVEC))
            {
               if (ob->ku < ob->vlen)
                  ob->ku = ob->vlen;
               else
                  while (ob->ku % ob->vlen)
                     ob->ku--;
            }
            ob->flag |= 1<<MMF_KRUNTIME;
         }
      }
      if ((flag&(1<<FLG_NOKVEC)) && FLAG_IS_SET(ob->flag, MMF_KVEC))
         KILL = 1;
      else if ((flag&(1<<FLG_NOMVEC)) && !FLAG_IS_SET(ob->flag, MMF_KVEC))
         KILL = 1;
      if (KILL)
         ob = KillMMNode(ob);
      else /* add to queue, and set it to use selected search type */
      {
         ATL_mmnode_t *nxt=ob->next;
         ob->flag &= ~MMF_MVSET;
         if (flag&((1<<FLG_OPCPY)|(1<<FLG_OP)))
            ob->flag |= (1<<MMF_MVA) | (1<<MMF_MVC) | (1<<MMF_CPC);
         else
            ob->flag |= (1<<MMF_MVA) | (1<<MMF_MVB);
         if (!ob->ID)
         {
            if (ob->genstr)
            {
               free(ob->genstr);
               ob->genstr = NULL;
            }
            if (ob->rout)
               free(ob->rout);
            ob->rout = DupString("ATL_tmp.c");
            if (flag&(1<<FLG_MKCOMP)) /* want to try run- & compile-time K */
            {
               ATL_mmnode_t *np=NULL;
               if (FLAG_IS_SET(ob->flag, MMF_KVEC))
               {
                  if (ob->ku == ob->vlen)
                     np = CloneMMNode(ob);
               }
               else if (ob->ku == 1)
                  np = CloneMMNode(ob);
               if (np)
               {
                  if (FLAG_IS_SET(ob->flag, MMF_KRUNTIME))
                     np->flag ^= 1<<MMF_KRUNTIME;
                  else
                     np->flag |= 1<<MMF_KRUNTIME;
                  np->next = nb;
                  nb = np;
               }
            }
         }     /* end if this is generated kern */
         ob->next = nb;
         nb = ob;
         ob = nxt;
      }     /* end else for keeper case */
   }        /* end while(ob); */
   return(nb);
}

void prepOP(ATL_mmnode_t *mb)
/*
 * Modifies MOV bits for rank-K timings
 */
{
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
      mp->flag = (mp->flag&(~MMF_MVSET)) |
                 ((1<<MMF_MVA)|(1<<MMF_MVC)|(1<<MMF_CPC));
}
/*
 * This search takes idx kern from fnin, and then finds the best performing
 * dimensions.  If KB==0, all three dims are tuned, else KB is fixed, and
 * only MB&NB are varied.
 */
void DoSrchD(char pre, int flag, int KB, char *fnin, int idx, char *fnout)
{
   ATL_mmnode_t *cb, *mp;
   double mf;
   int i, nu, NB;
   const int verb = flag&3;
   const unsigned int flg = ((flag>>FLG_OPCPY)&1)<<2;
/*
 * Get the idx kernel from the input file, and use it for all timings
 */
   cb = ReadMMFile(fnin);
   assert(cb);
   for (i=0, mp=cb; mp && i < idx; i++, mp = mp->next);
   assert(mp);
   cb = RemoveMMNodeFromQ(cb, mp);
   KillAllMMNodes(cb);

   prepOP(mp);

   if (KB)
      mp->kbB = KB;
   printf("   TUNING B FOR %s, PRESENT B=(%u,%u,%u), mf=%.2f:\n",
          GetMMLabelName(pre, mp), mp->mbB, mp->nbB, mp->kbB,
          mp->mflop ? mp->mflop[0]:0.0);
   GetMMLabelName(pre, NULL);
   if (KB)
      mf = MMExpandMN(verb, pre, flg, 0, 512, mp);
   else
      mf = MMExpandMNK(verb, pre, flg, 0, 0, NULL, 512, mp);
   printf("   BEST B=(%u,%u,%u), MFLOPS=%.2f\n", mp->mbB, mp->nbB, mp->kbB, mf);
   PrintMMLine(stdout, mp);
   if (fnout)
      WriteMMFile(fnout, mp);
   printf("\nB=(%d,%d,%d), MFLOPS=%.2f\n", mp->mbB, mp->nbB, mp->kbB, mf);
   KillMMNode(mp);
}
/*
 * This search takes a kernel mp, sets MB=MB, KB=KB, and finds the best
 * performing NB, which is restricted only by maxN
 */
void DoSrchN(char pre, int flag, int MB, int KB, int maxN, char *fnin, int idx)
{
   ATL_mmnode_t *cb, *mp;
   double mf;
   int i, nu, NB;
   const int verb = flag&3;
   const unsigned int flg = ((flag>>FLG_OPCPY)&1)<<2;
/*
 * Get the idx kernel from the input file, and use it for all timings
 */
   cb = ReadMMFile(fnin);
   assert(cb);
   for (i=0, mp=cb; mp && i < idx; i++, mp = mp->next);
   assert(mp);
   cb = RemoveMMNodeFromQ(cb, mp);
   KillAllMMNodes(cb);

   assert(MB % mp->mu == 0);
   nu = mp->nu;
   i = maxN % nu;
   if (i)
      maxN += nu-i;
   mp->mbB = MB;
   mp->kbB = KB;
   prepOP(mp);

   printf("   TUNING NB FOR BEST KERN, PRESENT NB=%d:\n", mp->nbB);
   NB = MMVaryDim(verb, pre, flg, 0, maxN, 0, mp);
   mf = mp->mflop[1];
   printf("   BEST NB=%u, MFLOPS=%.2f\n", NB, mf);
   mp->nbB = NB;
   mp->mflop[0] = mf;
   mp->mflop[1] = 0.0;
   PrintMMLine(stdout, mp);
   printf("\nB=(%d,%d,%d), MFLOPS=%.2f\n", MB, KB, NB, mf);
}

/*
 * This search is used for outer produce with unrestricted N & M.  In this
 * case, we get our reuse of B between amm calls, leaving A reuse the only
 * reuse that is important intra-kernel.  Therefore, this search takes a
 * list of kernels, and a specific KB.  It then sorts the kern by ascending MU.
 * For each MU, it then sets N=MAX(maxN,maxNU), and compares all kerns wt
 * target MU.  The winning kern is then tuned to find the best NB given
 * M=MU & K=KB, and then placed in the bb (base ptr to Best list).
 * As we add larger MUs we refuse to add those that don't beat any smaller
 * MU that evenly divide the new MU.  First do hand-written case above.
 */

/*
 * This search used for outer-product.  A kernel is found for each listed KB.
 * MB can be be forced to roughly match KB, or freely chosen, in which case
 * we will usually make it roughly match KB unless KB is small, and then
 * we may enlarge it to help with loop overhead.
 */
void opsrch(char pre, int flag, char *outnm, int *KBs, double runBon)
{
   ATL_mmnode_t *tb=NULL, *bb=NULL;  /* try and best base ptrs */
   ATL_mmnode_t *mp;
   double mf, mfB;
   const int verb = flag&3;
   int i;
   const unsigned int nKB=(*KBs);
   unsigned int maxKB, maxD, flg, MB, NB;  /* flag for TimeMMKernel */

   flg = ((flag>>FLG_OPCPY)&1)<<2;
   if (outnm)  /* if outnm exists, just time negative mflop and return */
   {
      ATL_mmnode_t *mb;
      mb = TimeMMFile(pre, outnm, 0, 1, flg, 0, 0, -1);
      if (mb)
      {
         KillAllMMNodes(mb);
         return;
      }
   }
/*
 * Sort KBs greatest-to-least, then peel largest block factor in order to
 * find good initial MB/NB
 */
   ATL_isortG2L(nKB, ++KBs);
   maxKB = *KBs;
   maxD = Mmax(maxKB, 120);

   tb = GetWorkingUserCases(verb, pre);
/*
 * Our normal list is all user kerns + best generated kerns.  We may add
 * special generated kerns later.
 */
   mp = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
   tb = ATL_JoinMMQs(tb, mp);
   #if 0
   for (mp=tb; mp; mp = mp->next)
      mp->flag = (mp->flag & ~MMF_MVSET)|((1<<MMF_MVA)|(1<<MMF_MVC));
   #else
      tb = ApplyRules(tb, flag, 0, 0);
   #endif
   printf("SEARCHING FOR BEST KERNEL AMONGST %u KERNELS AND %u KBs:\n",
          ATL_CountNumberOfMMNodes(tb), nKB);

   printf("   SEARCHINGS KERNS FOR MAXKB=%u\n", maxKB);
   mp = bestNearSquare(pre, verb, tb, maxKB, flg, 0, runBon);
   assert(mp);  /* for now, later gen a case */
   mf = mp->mflop[0];
   printf("   BEST KERN KB=%u: ID=%d '%s' mf=%.2f\n", maxKB, mp->ID,
          GetMMLabelName(pre, mp), mf);
   printf("   TUNING M/NB FOR BEST KERN, PRESENT B=(%u,%u):\n",
          mp->mbB, mp->nbB);
   mfB = MMExpandMN(verb, pre, flg, 0, 512, mp);
   printf("   B=(%u,%u,%u) GIVES SPEEDUP: %.2f\n", mp->mbB, mp->nbB, mp->kbB,
          mfB/mf);
   KillMMNode(mp);
   MB = mp->mbB;
   NB = mp->nbB;
/*
 * We purposely re-search maxKB, because enlarged MB/NB often picks higher
 * performing kernel
 */
   for (i=0; i < nKB; i++)
   {
      const int kb=KBs[i];

      printf("   SEARCHINGS KERNS FOR KB=%u\n", kb);
      mp = bestNearMNK(pre, verb, tb, MB, NB, kb, flg, 0, runBon);
      assert(mp);  /* for now, later gen a case */
      mf = mp->mflop[0];
      printf("   BEST KERN KB=%u: ID=%d '%s' mf=%.2f\n", kb, mp->ID,
             GetMMLabelName(pre, mp), mf);
      printf("   TUNING M/NB FOR BEST KERN, PRESENT B=(%u,%u):\n",
             mp->mbB, mp->nbB);
      mfB = MMExpandMN(verb, pre, flg, 0, 512, mp);
      printf("   B=(%u,%u,%u) GIVES SPEEDUP: %.2f\n", mp->mbB, mp->nbB, mp->kbB,
             mfB/mf);
      mp->next = bb;
      bb = mp;
   }
   printf("DONE.\n\n");
   KillAllMMNodes(tb);

   GetMMLabelName(pre, NULL);
   WriteMMFile(outnm, bb);
   KillAllMMNodes(bb);
}
/*
 * Brute-fource and ignorance search just tries every kernel for every
 * candidate size.  It is primarily used for small sizes where things like
 * cleanup or doing strange things like targeting multiple operands to the L1
 * cache can have unpredictable performance affects.
 *
 * BFI supports two types of searches:
 *  fuzzy: finds kernel that performs best for given size, with non-matching
 *         unrolling flops not counted in the flop rate.  It often requires
 *         over-allocating for A/B/C workspace.
 *  exact: finds kernel performs best for precise size, but bad sizes may
 *         have terrible performance (eg., prime sizes eliminating almost
 *         all non-generated kernels, and forcing generated kernels to use
 *         terrible register blockings).
 */
int main(int nargs, char **args)
{
   char *fnin, *fnout;
   int *MBs, *NBs, *KBs;
   double runBon;
   int flag, maxKU, reqVL;
   char pre;

   fnout = GetFlags(nargs, args, &pre, &flag, &maxKU, &reqVL, &fnin,
                    &MBs, &NBs, &KBs, &runBon);

   if (KBs)
   {
      assert(!MBs && !NBs);          /* currently, unsupported */
      assert(KBs);
      opsrch(pre, flag, fnout, KBs, runBon);
   }
   else if (NBs)
      DoSrchD(pre, 1<<FLG_OPCPY, NBs[1], fnin, NBs[0], fnout);
   else if (MBs)
      DoSrchD(pre, 1<<FLG_OPCPY, 0, fnin, MBs[0], fnout);

   if (fnout)
      free(fnout);
   if (MBs)
      free(MBs);
   if (NBs)
      free(NBs);
   if (KBs)
      free(KBs);
   return(0);
}
