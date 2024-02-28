#include "atlas_misc.h"
#include "atlas_mmtesttime.h"
#define CON_NOKVEC 0
#define CON_NOMVEC 1
#define CON_NOCOMPK 2


void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -T #  : set pruning tolerance for slower kernels:\n");
   fprintf(stderr,
   "      # <= 0.0 : keep all legal kernel sizes regardless of performance\n");
   fprintf(stderr,
   "      # >=0.0 : delete all kerns wt perf*# < maxSmallerPerf\n");
   fprintf(stderr,
   "                maxSmallerPerf= max perf found in smaller-sized kern\n");
   fprintf(stderr, "   -p [s,d]: set precision prefix (d) \n");
   fprintf(stderr, "   -b # nb1 ... nb# : square NBs to force\n");
   fprintf(stderr, "   -B # mb1 nb1 kb1 ... mb# nb# kb#: dims to force\n");
   fprintf(stderr, "   -S B0 Bn Binc : do series [B0,Bn] wt Binc stride\n");
   fprintf(stderr, "   -N maxB : try all blks up to & including maxB\n");
   fprintf(stderr,
      "   -F # a/b/c1 ... a/b/c#: which matblks should be cache flushed\n");
   fprintf(stderr, "   -v <verb> : set verbosity (1)\n");
   fprintf(stderr, "   -C [kmc] : Constrain kernel choice:\n");
   fprintf(stderr, "      k : don't allow K-vectorized storage\n");
   fprintf(stderr, "      m : don't allow M-vectorized storage\n");
   fprintf(stderr, "      c : don't allow compile-time K kernels\n");
   fprintf(stderr, "   -K <1/0> : do/don't generate K cleanup\n");
   fprintf(stderr, "   -o <outfile> : [res/<pre>UAMMRES.sum]\n");

   exit(ierr ? ierr : -1);
}

char *GetFlags(int nargs, char **args, char *PRE, int *C, int *MVS,
               int *NN, int **MBS, int **NBS, int **KBS, float *TOL,
               int *KCLEAN)
{
   ATL_mmnode_t *mmb=NULL;
   int B0=0, BN, KI;
   int *mbs=NULL, *nbs=NULL, *kbs=NULL;
   int i, k, j, N=0, ALL=0, MV=3, b0, bn, incB;
   char *cs, *fout=NULL;
   static char fnam[32];

   *TOL = 0.0;
   *PRE = 'd';
   *C = 0;
   *KCLEAN = 1;
   for (i=1; i < nargs; i++)
   {
      int n;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'o':
         if (++i >= nargs)
             PrintUsage(args[0], i-1, NULL);
         fout = args[i];
         break;
      case 'p':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *PRE = tolower(args[i][0]);
         assert(*PRE == 's' || *PRE == 'd' || *PRE == 'z' || *PRE == 'c');
         break;
      case 'T':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *TOL = atof(args[i]);
         if (*TOL < 0.0)
            *TOL = 0.0;
         break;
      case 'F':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         n = atoi(args[i]);
         for (MV=k=0; k < n; k++)
         {
           if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            if (args[i][0] == 'a' || args[i][0] == 'A')
               MV |= 1;
            else if (args[i][0] == 'b' || args[i][0] == 'B')
               MV |= 2;
            else if (args[i][0] == 'c' || args[i][0] == 'C')
               MV |= 4;
            else
               PrintUsage(args[0], -i, "UNKNOWN MATRIX FOR -M");
         }
         break;
      case 'K':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *KCLEAN = atoi(args[i]);
         break;
      case 'C':  /* -C <constraint string */
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         cs = args[i];
         n = strlen(cs);
         for (k=0; k < n; k++)
         {
            switch(cs[k])
            {
            case 'm':
               *C |= 1<<CON_NOMVEC;
               break;
            case 'k':
               *C |= 1<<CON_NOKVEC;
               break;
            case 'c':
               *C |= 1<<CON_NOCOMPK;
               break;
            default:
               PrintUsage(args[0], -i, "UNKNOWN CONSTRAINT");
            }
         }
         break;
      case 'N': /* <maxNB> */
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         else
         {
            int b0, bn, incB, n, k;
            bn = atoi(args[i]);
            b0 = 4;
            incB = 1;
            n = 1 + (bn-b0)/incB;
            mbs = malloc(n*sizeof(int));
            nbs = malloc(n*sizeof(int));
            kbs = malloc(n*sizeof(int));
            assert(mbs && nbs && kbs);
            for (k=0,n=b0; n <= bn; k++, n += incB)
               mbs[k] = nbs[k] = kbs[k] = n;
            assert(k <= n);
            *NN = k;
         }
         break;
      case 'S':  /* <b0> <bN> <incB> */
         if (i+3 >= nargs)
            PrintUsage(args[0], i-1, NULL);
         else
         {
            int b0, bn, incB, n, k;

            b0 = atoi(args[i+1]);
            bn = atoi(args[i+2]);
            incB = atoi(args[i+3]);
            n = 1 + (bn-b0)/incB;
            mbs = malloc(n*sizeof(int));
            nbs = malloc(n*sizeof(int));
            kbs = malloc(n*sizeof(int));
            assert(mbs && nbs && kbs);
            for (k=0,n=b0; n <= bn; k++, n += incB)
               mbs[k] = nbs[k] = kbs[k] = n;
            assert(k <= n);
            i += 3;
            *NN = k;
         }
         break;
      case 'B':
         if (nbs)
         {
            free(mbs);
            free(nbs);
            free(kbs);
         }
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *NN = N = atoi(args[i]);
         assert(N > 0);
         nbs = malloc(N*sizeof(int));
         assert(nbs);
         mbs = malloc(N*sizeof(int));
         assert(mbs);
         kbs = malloc(N*sizeof(int));
         assert(kbs);
         for (k=0; k < N; k++)
         {
           if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
           mbs[k] = atoi(args[i]);
           if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
           nbs[k] = atoi(args[i]);
           if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
           kbs[k] = atoi(args[i]);
         }
         break;
      case 'b':
         if (nbs)
         {
            free(mbs);
            free(nbs);
            free(kbs);
         }
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         N = atoi(args[i]);
         if (N < 0)
         {
            ALL = (N == -2) ? 32 : -1;  /* -2: don't force MB */
            N = 36;
         }
         *NN = N;
         assert(N > 0);
         nbs = malloc(N*sizeof(int));
         assert(nbs);
         mbs = malloc(N*sizeof(int));
         assert(nbs);
         kbs = malloc(N*sizeof(int));
         assert(nbs);
         if (ALL)  /* special case to just try most possible NBs */
         {
            mbs[ 0]=nbs[ 0]=kbs[ 0]=4;
            mbs[ 1]=nbs[ 1]=kbs[ 1]=6;
            mbs[ 2]=nbs[ 2]=kbs[ 2]=8;
            mbs[ 3]=nbs[ 3]=kbs[ 3]=12;
            mbs[ 4]=nbs[ 4]=kbs[ 4]=14;
            mbs[ 5]=nbs[ 5]=kbs[ 5]=16;
            mbs[ 6]=nbs[ 6]=kbs[ 6]=20;
            mbs[ 7]=nbs[ 7]=kbs[ 7]=22;
            mbs[ 8]=nbs[ 8]=kbs[ 8]=24;
            mbs[ 9]=nbs[ 9]=kbs[ 9]=26;
            mbs[10]=nbs[10]=kbs[10]=28;
            mbs[11]=nbs[11]=kbs[11]=32;
            mbs[12]=nbs[12]=kbs[12]=36;
            mbs[13]=nbs[13]=kbs[13]=40;
            mbs[14]=nbs[14]=kbs[14]=44;
            mbs[15]=nbs[15]=kbs[15]=48;
            mbs[16]=nbs[16]=kbs[16]=52;
            mbs[17]=nbs[17]=kbs[17]=56;
            mbs[18]=nbs[18]=kbs[18]=60;
            mbs[19]=nbs[19]=kbs[19]=64;
            mbs[20]=nbs[20]=kbs[20]=72;
            mbs[21]=nbs[21]=kbs[21]=80;
            mbs[22]=nbs[22]=kbs[22]=84;
            mbs[23]=nbs[23]=kbs[23]=88;
            mbs[24]=nbs[24]=kbs[24]=96;
            mbs[25]=nbs[25]=kbs[25]=104;
            mbs[26]=nbs[26]=kbs[26]=112;
            mbs[27]=nbs[27]=kbs[27]=120;
            mbs[28]=nbs[28]=kbs[28]=128;
            mbs[29]=nbs[29]=kbs[29]=132;
            mbs[30]=nbs[30]=kbs[30]=144;
            mbs[31]=nbs[31]=kbs[31]=156;
            mbs[32]=nbs[32]=kbs[32]=168;
            mbs[33]=nbs[33]=kbs[33]=192;
            mbs[34]=nbs[34]=kbs[34]=216;
            mbs[35]=nbs[35]=kbs[35]=240;
            for (k=0; k < ALL; k++)     /* don't force any particular */
               mbs[k] = 0;              /* MB during search */
         }
         else
         {
            for (k=0; k < N; k++)
            {
              if (++i >= nargs)
                  PrintUsage(args[0], i-1, NULL);
               mbs[k] = kbs[k] = nbs[k] = atoi(args[i]);
            }
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!nbs)
      PrintUsage(args[0], -1, "Dimensional flag (-b, -B -N, or -S) required!");
   *MBS = mbs;
   *NBS = nbs;
   *KBS = kbs;
   *MVS = MV;
   if (!fout)
   {
      sprintf(fnam, "res/%cuAMMRES.sum", *PRE);
      fout = fnam;
   }
   return(fout);
}

ATL_mmnode_t *DoSrch(ATL_mmnode_t *mmb, char pre, int flag, int beta,
                     int cflush, int nB, int *mbs, int *nbs, int *kbs)
/* need args for: moves */
{
   ATL_mmnode_t *mmB=NULL;
   int ib;
   printf("\nFINDING BEST AMONGST %d KERNELS and %d BLOCK FACTORS:\n",
          ATL_CountNumberOfMMNodes(mmb), nB);
   for (ib=0; ib < nB; ib++)
   {
      ATL_mmnode_t *mmp;
      const int mb = mbs[ib], nb=nbs[ib], kb=kbs[ib];
      mmp = MMBestWithGenKB(1, 0, flag, mmb, pre, mb, nb, kb, beta, -1, cflush);
      if (mmp)
      {
         mmp->next = mmB;
         mmB = mmp;
         printf("   BEST FOR B=(%d,%d,%d): %d,%s, mf=%.2f\n",
                mb, nb, kb, mmp->ID, mmp->rout, mmp->mflop[0]);
      }
      else
         printf("   NO VALID KERNEL FOR B=(%d,%d,%d)!\n", mb, nb, kb);
   }
   return(ReverseMMQ(mmB));
}

int main(int nargs, char **args)
{
   ATL_mmnode_t *mmb, *mmB;
   int *mbs, *nbs, *kbs;
   char *fnout;
   int C, MVS, N, KCLEAN;
   float tol;
   char pre, upr;

   fnout = GetFlags(nargs, args, &pre, &C, &MVS, &N, &mbs, &nbs, &kbs, &tol,
                    &KCLEAN);
   if (pre == 'd' || pre == 's')
      upr = pre;
   else
      upr = (pre == 'z') ? 'd' : 's';
/*
 * Get all working user-supplied kerns & best found generated kerns
 * assumes default search already run to produce [WORKING,gAMMRES].sum
 */
   mmb = ReadMMFileWithPath(upr, "res", "WORKING.sum");
   ATL_LastMMNode(mmb)->next = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
   assert(mmb);
/*
 * Remove kernels that don't meet required constraints
 */
   if (C)
   {
      ATL_mmnode_t *mp=mmb;
      do
      {
         ATL_mmnode_t *nxt = mp->next;
         int kill;
         kill  = (C & (1<<CON_NOKVEC)) && FLAG_IS_SET(mp->flag, MMF_KVEC);
         kill |= (C & (1<<CON_NOMVEC)) && !FLAG_IS_SET(mp->flag, MMF_KVEC);
         kill |= (C & (1<<CON_NOCOMPK)) && !FLAG_IS_SET(mp->flag, MMF_KRUNTIME);
         if (kill)
         {
            mmb = RemoveMMNodeFromQ(mmb, mp);
            KillMMNode(mp);
         }
         mp = nxt;
      }
      while (mp);
      assert(mmb);
   }
   MMApplyMoves2Flags(mmb, MVS);  /* use requested flushing */
   mmB = DoSrch(mmb, pre, 0, 1, -1, N, mbs, nbs, kbs);
   KillAllMMNodes(mmb);
   free(mbs);
   free(nbs);
   free(kbs);
   assert(mmB);
   if (tol > 0.0)
      MMPruneMflopTol(mmB, 0, tol);
   WriteMMFile(fnout, mmB);
   KillAllMMNodes(mmB);
   return(0);
}
