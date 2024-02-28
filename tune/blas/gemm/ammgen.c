/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_mmgen.h"
#include "atlas_sys.h"

#define NMMLISTS 8
#define IGE   0
#define IGEK1 1
#define ISQ   2
#define ISQK1 3
#define IRKK  4
#define ITRSM 5
#define ISYRK 6
#define ISYRKT 7

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
void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -p [s,d]: set type/precision prefix (d) \n");
   fprintf(stderr, "      s/d will generate for complex (c/z) as well\n");
   fprintf(stderr, "   -d <outdir>: directory to dump files to\n");
   fprintf(stderr, "   -g [r,c] <rect kern file> <geKCleanFile>\n");
   fprintf(stderr, "   -r [r,c] <rank-K kernel file> \n");
   fprintf(stderr, "   -s [r,c] <square-case kernel file> <sqKCleanFile>\n");
   fprintf(stderr, "   -t [r,c] <trsm file> : should match -s\n");
   fprintf(stderr, "   -k syrk <syrk file> \n");
   exit(ierr ? ierr : -1);
}


/*
 * RETURNS: precision prefix [s,d,c,z]
 */
char *GetFlags(int nargs, char **args, char *PRE, ATL_mmnode_t **lists)
{
   int i, j=0, n, k;
   char pre='d', cpre, ch;
   char *outd=NULL;
   char *fn[8]={"geAMMRES.sum","geAMMKCLEAN.sum",
      "sqAMMRES.sum","sqAMMKCLEAN.sum", "rkAMMRES.sum",
      "tsAMMRES.sum","gAMSYRK.sum", "syrkK.sum"};

   for (i=0; i < 2*NMMLISTS; i++)
      lists[i] = NULL;
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
        if (pre == 'z')
           pre = 'd';
        else if (pre == 'c')
           pre = 's';
        assert(pre == 's' || pre == 'd');
        cpre = (pre == 'd') ? 'z' : 'c';
        break;
      case 's': /* [r,c] <amm file> <K1file> */
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        ch = args[i][0];
        if (ch == 'c')
           k = NMMLISTS + 2;
        else if (ch == 'r')
           k = 2;
        else
            PrintUsage(args[0], i, "1st arg to -s must be 'r' or 'c'!");
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        lists[k] = ReadMMFile(args[i]);
        assert(lists[k]);
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        lists[++k] = ReadMMFile(args[i]);
        assert(lists[k]);
        break;
      case 'g': /* [r,c] <amm file> <K1file> */
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        ch = args[i][0];
        if (ch == 'c')
           k = NMMLISTS + 0;
        else if (ch == 'r')
           k = 0;
        else
            PrintUsage(args[0], i, "1st arg to -g must be 'r' or 'c'!");
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        lists[k] = ReadMMFile(args[i]);
        assert(lists[k]);
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        lists[++k] = ReadMMFile(args[i]);
        assert(lists[k]);
        break;
      case 't': /* [r,c] <amm file> */
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        ch = args[i][0];
        if (ch == 'c')
           k = NMMLISTS + 5;
        else if (ch == 'r')
           k = 5;
        else
            PrintUsage(args[0], i, "1st arg to -t must be 'r' or 'c'!");
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        lists[k] = ReadMMFile(args[i]);
        assert(lists[k]);
      case 'r': /* [r,c] <amm file> */
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        ch = args[i][0];
        if (ch == 'c')
           k = NMMLISTS + 4;
        else if (ch == 'r')
           k = 4;
        else
            PrintUsage(args[0], i, "1st arg to -r must be 'r' or 'c'!");
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        lists[k] = ReadMMFile(args[i]);
        assert(lists[k]);
      case 'k':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        assert(!strcmp(args[i], "syrk"));
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        k = ISYRK;
        lists[k] = ReadMMFile(args[i]);
        break;
      case 'd':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        outd = DupString(args[i]);
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   *PRE = pre;
/*
 * Fill in any required list entry
 */
   for (i=0; i < 2*NMMLISTS; i++)
   {
      const char pr=(i<NMMLISTS) ? pre : cpre;
      const int k = (i<NMMLISTS) ? i : i-NMMLISTS;
      if (i == NMMLISTS+ITRSM) /* complex trsm not */
         continue;             /* presently used */
      if (!lists[i])
         lists[i] = ReadMMFileWithPath(pr, "res", fn[k]);
      if (!lists[i])
         fprintf(stderr, "CANNOT FIND FILE 'res/%c%s'!\n", pr, fn[k]);
      assert(lists[i]);
   }
   return(outd);
}

void GenSyrkH(char pre, char *outd, ATL_mmnode_t *mb, int RCSAME)
{
   ATL_cpnode_t *cb;
   FILE *fp;
   const int ntr = (pre == 'd' || pre == 's') ? 2 : 4;
   int ial, ibe, itr;
   const char be[3] = {'0', '1', 'n'};
   char bes[4] = {'1', 'n', 'X', '0'};
   char trs[4] = {'N', 'T', 'C', 'H'};
   char rcs[4] = {'c', 'r', 'c', 'r'};
   char cjs[4] = {' ', ' ', 'C', 'C'};


   if (!mb)
      return;
   assert(!mb->next);
   fp = OpenMMGenHeader(outd, 0, pre, NULL, "amm", "syrk", NULL);
   fprintf(fp, "#define SYRK_NB %d\n", mb->kbB);
   fprintf(fp, "#define ATL_SYRKK_VLEN %d\n", mb->vlen);
   fprintf(fp, "#define ATL_SYRKK_KVEC %d\n",
           FLAG_IS_SET(mb->flag, MMF_KVEC) ? mb->vlen:0);
   fprintf(fp, "#define ATL_SYRKK_NU %d\n", mb->nu);
   fprintf(fp, "#define ATL_SYRKK_KU %d\n", mb->ku);
/*
 * Prototype ATL_<pre>_[sq,um]syrkK kernel
 */
   if (RCSAME && (pre == 'c' || pre == 'z'))
   {
      const char upr = (pre == 'z') ? 'd' : 's';
      for (ibe=0; ibe < 3; ibe++)
         fprintf(fp, "#define ATL_%c%ssyrkK_b%c ATL_%c%ssyrkK_b%c\n",
                 pre, sh, be[ibe], upr, sh, be[ibe]);
   }
   for (ibe=0; ibe < 3; ibe++)
   {
      fprintf(fp,
"void ATL_%camsyrkK_b%c(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,\n",
              pre, be[ibe]);
      fprintf(fp,
      "                     TYPE*, const TYPE*, const TYPE*, const TYPE*);\n");
   }
/*
 * Prototype/rename all A cpy routs: ATL_<pre>cm2am_syrk<TA>
 */
#if 1
   cb = AddMMUniqueACopyFromMMNodes(pre, 'T', mb, NULL);
   fprintf(fp, "#define ATL_%ca2blk_syrkN %s_Na1\n",
           pre, cb->rout);
   fprintf(fp, "#define ATL_%ca2blk_syrkT %s_Ta1\n",
           pre, cb->rout);
   for (itr=0; itr < 2; itr++)
   {
      char tr = itr ? 'T' : 'N';
      fprintf(fp, "void ATL_%ca2blk_syrk%c", pre, tr);
      if (pre == 'd' || pre == 's')
         fprintf(fp,
         "(ATL_CSZT,ATL_CSZT,const TYPE,const TYPE*,ATL_CSZT,TYPE*);\n");
      else
         fprintf(fp,
         "(ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,ATL_CSZT,TYPE*,TYPE*);\n");
   }
   KillAllCopyNodes(cb);
#else
   for (itr=0; itr < ntr; itr++)
   {
      fprintf(fp, "#define ATL_%ccm2am_syrk%c ATL_%c%cm2am_a1_%dx%d%c\n",
              pre, trs[itr], pre, rcs[itr], mb->nu,
              FLAG_IS_SET(mb->flag, MMF_KVEC)?mb->vlen:0, cjs[itr]);
      fprintf(fp, "void ATL_%ccm2am_syrk%c", pre, trs[itr]);
      if (pre == 'd' || pre == 's')
         fprintf(fp,
         "(ATL_CSZT,ATL_CSZT,const TYPE,const TYPE*,ATL_CSZT,TYPE*);\n");
      else
         fprintf(fp,
         "(ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,ATL_CSZT,TYPE*,TYPE*);\n");
   }
#endif
/*
 * Prototype all C cpy routs: ATL_<pre>SyrkIntoC_a<alp>_b<bet>
 */
   for (ial=0; ial < 3; ial++)
   {
      for (ibe=0; ibe < 4; ibe++)
      {
         fprintf(fp, "void ATL_%cSyrkIntoC_a%c_b%c\n", pre,bes[ial],bes[ibe]);
         if (pre == 's' || pre == 'd')
            fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const SCALAR, TYPE *,ATL_CSZT);\n");
         else
            fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const TYPE*,const SCALAR, TYPE *,ATL_CSZT);\n");
      }
   }
   CloseGenHeader(fp);
}

double GenPerfH(char pre, char *outd, char *ip, char *cn, double mfMax,
                ATL_mmnode_t *mb, ATL_mmnode_t *k1)
{
   FILE *fp;
   ATL_mmnode_t *mp;
   #define NTHRSH 11
   int THRSH[NTHRSH] = {25, 33, 50, 66, 75, 80, 85, 90, 95, 98, 99};
   int idxT[NTHRSH];
   ATL_mmnode_t *mpT[NTHRSH];
   const int FNDMAX = (mfMax == 0.0);
   int i, n, idxMax=0;
   double fcnt;

   fp = OpenMMGenHeader(outd, 0, pre, ip, cn, "perf", mb);
/*
 * If not already set, compute the max perf & its index.  We'll use this
 * max as the denominator in our PERF array.
 */
   if (FNDMAX)
   {
      for (i=0; i < NTHRSH; i++)
         mpT[i] = NULL;
      for (n=0,mp=mb; mp; n++, mp = mp->next)
      {
         if (mp->mflop[0] > mfMax)
         {
            mfMax = mp->mflop[0];
            idxMax = n;
         }
      }
   }
/*
 * If our denom comes from other file, just count the number of entries
 */
   else
   {
      for (n=0,mp=mb; mp; n++, mp = mp->next);
      idxMax = -1;
   }
   fprintf(fp, "#define ATL_%sAMM_NCASES %d\n", cn, n);
   fprintf(fp, "#define ATL_%sAMM_DENOM %le /* (%.2f)*/ \n", cn,
           mfMax, mfMax);
   fprintf(fp, "#define ATL_%sAMM_MAXMFLOPIDX %d\n\n", cn, idxMax);
   if (FNDMAX)
   {
      for (n=0,mp=mb; mp; mp = mp->next, n++)
      {
         double mf = mp->mflop[0] / mfMax;
         for (i=0; i < NTHRSH; i++)
         {
            if (!mpT[i] && THRSH[i]*0.01*mfMax < mp->mflop[0])
            {
               mpT[i] = mp;
               idxT[i] = n;
            }
         }
      }
      for (i=0; i < NTHRSH; i++)
      {
         mp = mpT[i];
         fprintf(fp, "#define ATL_%sAMM_%dLCMU %d\n", cn, THRSH[i],
                 Mylcm(Mylcm(mp->mu,mp->nu),mp->ku));
         fprintf(fp, "#define ATL_%sAMM_%dLCMMN %d\n", cn, THRSH[i],
                 Mylcm(mp->mu,mp->nu));
         fprintf(fp, "#define ATL_%sAMM_%dMB %d\n", cn, THRSH[i], mp->mbB);
         fprintf(fp, "#define ATL_%sAMM_%dNB %d\n", cn, THRSH[i], mp->nbB);
         fprintf(fp, "#define ATL_%sAMM_%dKB %d\n", cn, THRSH[i], mp->kbB);
         fprintf(fp, "#define ATL_%sAMM_%dIDX %d\n", cn, THRSH[i], idxT[i]);
      }
      mp = ATL_FindLastNode(mb, GetOffset(&mb->next, mb));
      fprintf(fp, "#define ATL_%sAMM_LCMU %d\n", cn,
              Mylcm(Mylcm(mp->mu,mp->nu),mp->ku));
      fprintf(fp, "#define ATL_%sAMM_LCMMN %d\n", cn, Mylcm(mp->mu,mp->nu));
      fprintf(fp, "\n");
   }

   fprintf(fp, "#if !defined(NOPERF) && !defined(NOARRS)\n");
   fprintf(fp, "static const float ATL_%sAMM_PERF[%d] =\n", cn, n);
   fprintf(fp, "{  /* %% of performance of best kernel */\n");
   for (i=0,mp=mb; mp; i++,mp = mp->next)
      fprintf(fp, "   %f%c  /* IDX=%d, KB=%d, mf=%.0f */\n", mp->mflop[0]/mfMax,
              (mp->next)?',':' ', i, mp->kbB, mp->mflop[0]);
   fprintf(fp, "};\n#endif\n\n");

   if (k1)
   {
      ATL_mmnode_t *kp;
      fprintf(fp, "#if !defined(NOK1RATIO) && !defined(NOARRS)\n");
      fprintf(fp, "static const float ATL_%sAMM_K1RATIO[%d] =\n", cn, n);
      fprintf(fp, "{  /* ratio of %sAMM perf wt that of its K=1 K-cleaner */\n",
              cn);
      for (i=0,mp=mb,kp=k1; mp; i++,mp = mp->next,kp = kp->next)
         fprintf(fp, "   %f%c  /* IDX=%d, KB=%d IDs=[%d,%d]*/\n",
                 kp->mflop[0]/mp->mflop[0], (mp->next)?',':' ', i, mp->kbB,
                 mp->ID, kp->ID);
      fprintf(fp, "};\n#endif\n\n");
   }

   fprintf(fp, "#if !defined(NOTIME) && !defined(NOARRS)\n");
   fprintf(fp, "static const float ATL_%sAMM_TIME[%d] =\n", cn, n);
   fprintf(fp, "{  /* actual seconds to compute one block */\n");
   for (i=0,mp=mb; mp; i++,mp = mp->next)
   {
      if (mp->blask == 0)
         fcnt = (2.0*mp->mbB)*mp->nbB*mp->kbB;  /* gemm flop count */
      else
         fcnt = (1.0*mp->mbB)*mp->nbB*mp->kbB;  /* roughly right */
      fprintf(fp, "   %e%c  /* IDX=%d, B=(%d,%d,%d) */\n",
              fcnt / (mp->mflop[0]*1.0e6),
              (mp->next)?',':' ', i, mp->mbB, mp->nbB, mp->kbB);
   }
   fprintf(fp, "};\n#endif\n");
   CloseGenHeader(fp);
   return(mfMax);
}

void GenKernH
   (char pre, char *outd, char *ip, char *cn, ATL_mmnode_t *mb,
    ATL_mmnode_t *mkb, ATL_mmnode_t *ub)
{
   FILE *fp;
   int ib, inxt, iaut;
   char bes[3] = {'0', '1', 'n'};

   inxt = GetOffset(&mb->next, mb);
   iaut = GetOffset(&ub->auth, ub);
   fp = OpenMMGenHeader(outd, 0, pre, ip, cn, "kern", mb);
   for (ib=0; ib < 3; ib++)
      PrintMMProtos(fp, pre, "KERN", ub, iaut, bes[ib]);
   for (ib=0; ib < 3; ib++)
   {
      PrintStrArrAtOff(fp, pre, "KERN", mb, inxt, iaut, "ammkern_t",
                       0, 0, bes[ib]);
      if (mkb)
         PrintStrArrAtOff(fp, pre, "KERN_K1", mkb, inxt, iaut, "ammkern_t",
                          0, 0, bes[ib]);
   }
   CloseGenHeader(fp);
}

void PrintFlagArr(FILE *fp, ATL_mmnode_t *mb)
{
   ATL_mmnode_t *mp;
   int n;

   n = CountListEntries(mb, GetOffset(&mb->next, mb));
   fprintf(fp, "#ifndef NOKFLAG\n");
   fprintf(fp, "static const unsigned char ATL_AMM_KFLAG[%d] =\n{\n", n);
   for (n=0,mp=mb; mp; n++, mp = mp->next)
   {
      unsigned char flag=FLAG_IS_SET(mp->flag, MMF_KRUNTIME) ? 1 : 0;
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         flag |= 2;
      fprintf(fp, "%6d%c   /* IDX=%d */\n", flag, mp->next ? ',':' ', n);
   }
   fprintf(fp, "};\n");
   fprintf(fp, "   #define ATL_AMM_KRUNTIME(idx_) (ATL_AMM_KFLAG[idx_] & 1)\n");
   fprintf(fp, "   #define ATL_AMM_KMAJOR(idx_) (ATL_AMM_KFLAG[idx_] & 2)\n");
   fprintf(fp, "#endif\n");
}
void GenBlkH(char pre, char *outd, char *ip, char *cn, ATL_mmnode_t *mb)
{
   FILE *fp;
   ATL_mmnode_t *mp;
   int *KUs;
   int i, n, inxt;

   fp = OpenMMGenHeader(outd, 0, pre, ip, cn, "blk", mb);
   inxt = GetOffset(&mb->next, mb);
   mp = ATL_FindLastNode(mb, inxt);
   fprintf(fp, "#define ATL_%sAMM_LASTMB %d\n", cn, mp->mbB);
   fprintf(fp, "#define ATL_%sAMM_LASTNB %d\n", cn, mp->nbB);
   fprintf(fp, "#define ATL_%sAMM_LASTKB %d\n", cn, mp->kbB);
   fprintf(fp, "\n");
/*
 * Save original KUs, and overwrite compile-time-K KUs with kbB for printing
 */
   n = CountListEntries(mb, inxt);
   KUs = malloc(sizeof(int)*n*3);
   assert(KUs);
   for (i=0,mp=mb; mp; i++,mp = mp->next)
   {
      KUs[i] = mp->ku;
      KUs[i+n] = mp->kbmin;
      KUs[i+n+n] = mp->kbmax;
      #if 0
      if (FLAG_IS_SET(mp->flag,MMF_KUISKB))
         mp->kbmin = mp->kbmax = mp->ku = mp->kbB;
      else if (!FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
         mp->kbmin = mp->kbmax = mp->kbB;
      #endif
   }
   PrintIntArrAtOff(fp, pre, "MBs", mb, inxt, GetOffset(&mb->mbB, mb),0,0);
   PrintIntArrAtOff(fp, pre, "NBs", mb, inxt, GetOffset(&mb->nbB, mb),0,0);
   PrintIntArrAtOff(fp, pre, "KBs", mb, inxt, GetOffset(&mb->kbB, mb),0,0);
   PrintIntArrAtOff(fp, pre, "MUs", mb, inxt, GetOffset(&mb->mu, mb),0,0);
   PrintIntArrAtOff(fp, pre, "NUs", mb, inxt, GetOffset(&mb->nu, mb),0,0);
   PrintIntArrAtOff(fp, pre, "KUs", mb, inxt, GetOffset(&mb->ku, mb),0,0);
   PrintIntArrAtOff(fp, pre, "KBMINs", mb, inxt, GetOffset(&mb->kbmin, mb),0,0);
   PrintIntArrAtOff(fp, pre, "KBMAXs", mb, inxt, GetOffset(&mb->kbmax, mb),0,0);
   PrintIntArrAtOff(fp, pre, "VLENs", mb, inxt, GetOffset(&mb->vlen, mb),0,0);
   for (i=0,mp=mb; mp; i++,mp = mp->next)
   {
      mp->ku = KUs[i];
      mp->kbmin = KUs[i+n];
      mp->kbmax = KUs[i+n+n];
   }
   PrintFlagArr(fp, mb);
   free(KUs);
   CloseGenHeader(fp);
}


void PrintDiv(FILE *fp, char *exp, int v)
/*
 * Prints either exp/v or exp>>lg2(v), dep on v being power of 2 or not
 * caller must put parens where appropriate (if exp is not a simple variable
 * name, it must have parens to produce the right answer!).
 */
{
   int p;
   assert(v > 0);  /* 0 cannot be supported, no need for neg now */
   if (v == 1)
      fprintf(fp, "%s", exp);
   for (p=1; (1<<p) < v; p++);
   if (v == (1<<p))
      fprintf(fp, "%s>>%d", exp, p);
   else
      fprintf(fp, "%s/%d", exp, v);
}
void PrintCeil(FILE *fp, char *exp, int v)
/*
 * Prints either ((exp+v-1)/v)*v or same things wt shifts if power of 2.
 */
{
   int p;
   assert(v > 0);  /* 0 cannot be supported, no need for neg now */
   if (v == 1)
      fprintf(fp, "%s", exp);
   for (p=1; (1<<p) < v; p++);
   if (v == (1<<p))
      fprintf(fp, "((%s+%d)>>%d)", exp, v-1, p);
   else
      fprintf(fp, "((%s+%d)/%d)", exp, v-1, v);
}
void PrintCeilMul(FILE *fp, char *exp, int v)
{
   int p;
   assert(v > 0);  /* 0 cannot be supported, no need for neg now */
   if (v == 1)
      fprintf(fp, "%s", exp);
   for (p=1; (1<<p) < v; p++);
   if (v == (1<<p))
      fprintf(fp, "(((%s+%d)>>%d)<<%d)", exp, v-1, p, p);
   else
      fprintf(fp, "(((%s+%d)/%d)*%d)", exp, v-1, v, v);
}


double PrintSum(FILE *fp, char *cn, ATL_mmnode_t *b, double mfGE)
{
   int inxt, max, i, kvec;
   double mf;
   ATL_mmnode_t *mp;

   inxt = GetOffset(&b->next, b);
   fprintf(fp, "#define ATL_%sAMM_NCASES %d\n", cn, CountListEntries(b, inxt));
   mp = ATL_FindLastNode(b, inxt);
   fprintf(fp, "#define ATL_%sAMM_LASTMB %d\n", cn, mp->mbB);
   fprintf(fp, "#define ATL_%sAMM_LASTNB %d\n", cn, mp->nbB);
   fprintf(fp, "#define ATL_%sAMM_LASTKB %d\n", cn, mp->kbB);
   mf = mp->mflop[0];
   fprintf(fp, "#define ATL_%sAMM_LASTMU %d\n", cn, mp->mu);
   fprintf(fp, "#define ATL_%sAMM_LASTNU %d\n", cn, mp->nu);
   fprintf(fp, "#define ATL_%sAMM_LASTKU %d\n", cn, mp->nu);
   fprintf(fp, "#define ATL_%sAMM_LASTLCMMN %d\n\n", cn, Mylcm(mp->mu, mp->nu));
   fprintf(fp, "#define ATL_%sAMM_LASTLCMU %d\n\n", cn,
           Mylcm(Mylcm(mp->mu, mp->nu),mp->ku));
   #if 0
   fprintf(fp, "#define ATL_%sAMM_LASTMFLOP %e\n", cn, mf);
   if (mfGE > 0.0)
      fprintf(fp, "#define ATL_%sAMM_MFLOP_RATIO %e\n", cn, mf/mfGE);
   #endif
   for (kvec=0,mp=b; mp; mp = mp->next)
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         kvec = Mmax(kvec, mp->vlen);
   fprintf(fp, "#ifndef  ATL_%sAMM_MAXKVEC\n", cn);
   fprintf(fp, "   #define ATL_%sAMM_MAXKVEC %d\n", cn, kvec);
   fprintf(fp, "#endif\n");
   fprintf(fp, "\n");
   return(mf);
}

void GenSumH(char pre, char *outd, ATL_mmnode_t *gb, ATL_mmnode_t *sb,
             ATL_mmnode_t *rb)
/*
 * cn:[ge,sq,rk]
 * For [ge,sq,rk] report: ncases, max[nb,mb,kb]
 * -> for maxB, report: perf,time,ratio wt best ge
 */
{
   FILE *fp;
   double mf;
   int i;

   fp = OpenMMGenHeader(outd, 0, pre, NULL, "amm", "sum", gb);
   mf = PrintSum(fp, "ge", gb, 0.0);
   if (sb)
      PrintSum(fp, "sq", sb, mf);
   if (rb)
      PrintSum(fp, "rk", rb, mf);
   CloseGenHeader(fp);
}

void GenDegenH_IP(char pre, char *outd, char *cn, ATL_mmnode_t *mb)
/*
 * Used to generate timing/index file for when one or more dims are degenerate
 * for inner product.  mb should be ordered by the degen dim(s), ivar should
 * have the index of kern.h that corresponds to this kern.
 */
{
   FILE *fp;
   ATL_mmnode_t *mp;
   int i, n, idxMax=0, maxM=0, maxN=0, maxK=0;
   double fcnt, mfMax=0.0;

   fp = OpenMMGenHeader(outd, 0, pre, "ip", cn, "geIdx", mb);
/*
 * Compute max perf & M,N,K
 */
   for (n=0,mp=mb; mp; n++, mp = mp->next)
   {
      maxM = Mmax(maxM, mp->mbB);
      maxN = Mmax(maxN, mp->nbB);
      maxK = Mmax(maxK, mp->kbB);
      if (mp->mflop[0] > mfMax)
      {
         if (mp->ivar > 0)
            (mp->ivar)--;
         mfMax = mp->mflop[0];
         idxMax = n;
      }
   }
   fprintf(fp, "#ifdef NOARRS\n");
   fprintf(fp, "   #define NOPERF 1\n");
   fprintf(fp, "   #define NOTIME 1\n");
   fprintf(fp, "   #define NOKIDX 1\n");
   fprintf(fp, "   #define NONBs  1\n");
   fprintf(fp, "#endif\n");
   fprintf(fp, "#define ATL_%sAMM_NCASES %d\n", cn, n);
   fprintf(fp, "#define ATL_%sAMM_DENOM %le /* (%.2f)*/ \n", cn,
           mfMax, mfMax);
   fprintf(fp, "#define ATL_%sAMM_MAXMFLOPIDX %d\n", cn, idxMax);
   fprintf(fp, "#define ATL_%sAMM_MAXMB %d\n", cn, maxM);
   fprintf(fp, "#define ATL_%sAMM_MAXNB %d\n", cn, maxN);
   fprintf(fp, "#define ATL_%sAMM_MAXKB %d\n\n", cn, maxK);

   fprintf(fp, "#ifndef NOPERF\n");
   fprintf(fp, "static const float ATL_%sAMM_PERF[%d] =\n", cn, n);
   fprintf(fp, "{  /* %% of performance of best kernel */\n");
   for (i=0,mp=mb; mp; i++,mp = mp->next)
      fprintf(fp, "   %f%c  /* CNT=%d, KB=%d, mf=%.0f */\n", mp->mflop[0]/mfMax,
              (mp->next)?',':' ', i, mp->kbB, mp->mflop[0]);
   fprintf(fp, "};\n#endif\n\n");

   fprintf(fp, "#if !defined(NOTIME) && !defined(NOARRS)\n");
   fprintf(fp, "static const float ATL_%sAMM_TIME[%d] =\n", cn, n);
   fprintf(fp, "{  /* actual seconds to compute one block */\n");
   for (i=0,mp=mb; mp; i++,mp = mp->next)
   {
      if (mp->blask == 0)
         fcnt = (2.0*mp->mbB)*mp->nbB*mp->kbB;  /* gemm flop count */
      else
         fcnt = (1.0*mp->mbB)*mp->nbB*mp->kbB;  /* roughly right */
      fprintf(fp, "   %e%c  /* CNT=%d, B=(%d,%d,%d) */\n",
              fcnt / (mp->mflop[0]*1.0e6),
              (mp->next)?',':' ', i, mp->mbB, mp->nbB, mp->kbB);
   }
   fprintf(fp, "};\n#endif\n");

   fprintf(fp, "#define ATL_AMM_NBs ATL_pmnAMM_NBs\n");
   PrintIntArrAtOff(fp, pre, "NBs", mb, GetOffset(&mb->next, mb),
                    GetOffset(&mb->nbB, mb), 0, 0);
   fprintf(fp, "#undef ATL_AMM_NBs\n");
   fprintf(fp, "#define ATL_AMM_KBs ATL_pmnAMM_KBs\n");
   PrintIntArrAtOff(fp, pre, "KBs", mb, GetOffset(&mb->next, mb),
                    GetOffset(&mb->kbB, mb), 0, 0);
   fprintf(fp, "#undef ATL_AMM_KBs\n");
   PrintIntArrAtOff(fp, pre, "KIDX", mb, GetOffset(&mb->next, mb),
                    GetOffset(&mb->ivar, mb), 0, 0);
   fprintf(fp, "#ifdef NOARRS\n");
   fprintf(fp, "   #undef NOPERF\n");
   fprintf(fp, "   #undef NOTIME\n");
   fprintf(fp, "   #undef NOKIDX\n");
   fprintf(fp, "   #undef NONBs\n");
   fprintf(fp, "#endif\n");
   CloseGenHeader(fp);
}

void GenAllHeaders(char pre, char *outd, ATL_mmnode_t **lists)
{
   ATL_cpnode_t *cb;
   ATL_mmnode_t *mmb;
   char *nm="perf";
   char *cn[4]={"ge", "sq", "rk", "ts"};
   double mfMax;
   int i, k, SKSAME;
   char pr;
/*
 * Write degenerate case header files.  Kind of ugly hardwiring names here,
 * but I don't see point in taking these from command line.
 */
   pr=pre;
   DO_TYPE:
   if (pr == pre)
   {
      pr = (pre == 'd') ? 'z' : 'c';
      goto DO_TYPE;
   }
   SKSAME = MMKernsPerfSame(lists[ISYRK], lists[ISYRK+NMMLISTS]);
/*
 * Handle [perf,blk,flag,sum] which require only ordered lists
 */
   for (pr=pre,k=0; k < 2*NMMLISTS; k += NMMLISTS)
   {
      mfMax = GenPerfH(pr, outd, "ip", "ge", 0.0, lists[k+IGE], lists[k+IGEK1]);
      if (k < NMMLISTS)  /* ts only in real for now */
         GenPerfH(pr, outd, "tr","sm", mfMax, lists[k+ITRSM], NULL);
      KillAllMMNodes(lists[k+ITRSM]);  /* done with these! */
      lists[k+ITRSM] = NULL;
      GenPerfH(pr, outd, "sy","sk", mfMax, lists[k+ISYRKT], NULL);
      KillAllMMNodes(lists[k+ISYRKT]);  /* done with these! */
      lists[k+ISYRKT] = NULL;
      GenSyrkH(pr, outd, lists[k+ISYRK], SKSAME);

      GenPerfH(pr, outd, "op","sq", 0.0, lists[k+ISQ], lists[k+ISQK1]);
      GenPerfH(pr, outd, "op","rk", 0.0, lists[k+IRKK], NULL);

      GenBlkH(pr, outd, "ip", "ge", lists[k+IGE]);
      GenBlkH(pr, outd, "op", "sq", lists[k+ISQ]);
      GenBlkH(pr, outd, "op", "rk", lists[k+IRKK]);
      pr = (pre == 'd') ? 'z' : 'c';
   }
   GenSumH(pre, outd, lists[IGE], lists[ISQ], lists[IRKK]);
   GenSumH(pr, outd, lists[IGE+NMMLISTS], lists[ISQ+NMMLISTS],
           lists[IRKK+NMMLISTS]);
/*
 * Create lists of only the unique kernels for [ge,sq,rk] in real & cplx
 * for use in prototyping.  We combine kern & K1 lists for each.
 * We can now gen [kern].
 */
   for (pr=pre,k=0; k < 2; k++)
   {
      for (i=0; i < 3; i++)
      {
         ATL_mmnode_t *ulst;  /* unordered list of only unique kernels */
         const int h = i+i+k*NMMLISTS, u=k*3+i;

         ulst = AddUniqueMMKernCompList(NULL, lists[h]);
         if (i < 2)
         {
            char *ip = (i) ? "op":"ip";
            ulst = AddUniqueMMKernCompList(ulst, lists[h+1]);
            GenKernH(pr, outd, ip, cn[i], lists[h], lists[h+1], ulst);
         }
         else
            GenKernH(pr, outd, i==2?"op":"tr", cn[i], lists[h], NULL, ulst);
         KillAllMMNodes(ulst);  /* done with these! */
      }
      pr = (pre == 's') ? 'c' : 'z';
   }
}

static void AddAlpBetSuf(char *rt, int rl, int ial, int ibe)
/*
 * Add _aXbX suffix to rt, which is presently of length rl
 */
{
   rt[rl] = '_';
   rt[rl+1] = 'a';
   if (ial == -1)
      rt[rl+2] = 'n';
   else if (ial == 2)
      rt[rl+2] = 'X';
   else
      rt[rl+2] = ial+48;
   rt[rl+3] = 'b';
   if (ibe == -1)
      rt[rl+4] = 'n';
   else if (ibe > 1)
      rt[rl+4] = 'X';
   else
      rt[rl+4] = ibe+48;
}

void GenAllKerns(char pre, char *outd, ATL_mmnode_t *rb)
{
   ATL_mmnode_t *mp;
   char *sgen=NULL;
   int i, len, dlen, mID, mU;
   char pr=pre;
   dlen = strlen(outd);
   len = 0;
/*
 * Generate/copy all required kernels.  This list has files with compile-time
 * K repeated, but we'll check ID to avoid repeated copies of user-supplied,
 * and generated kernels need a new generation for each required KB.
 */
   printf("\nGENERATING AMM KERNS:\n");
   for (mp=rb; mp; mp = mp->next)
   {
      const int id=mp->ID;
      if (!id)
      {
         printf("   -> %s\n", mp->rout);
         assert(mp->genstr);
         if ( (i=system(mp->genstr)) )
         {
            fprintf(stderr, "GENSTR RETURNS %d:\n'%s'\n", i, mp->genstr);
            exit(i);
         }
      }
      else /* user-supplied kernel */
      {
         ATL_mmnode_t *p;
         for (p=rb; p != mp && p->ID != id; p = p->next);
         if (p == mp)  /* this is first mention of this ID */
         {
            printf("   %s -> %s\n", mp->genstr, mp->rout);
            i = strlen(mp->genstr) + strlen(mp->rout) + dlen + 16;
            if (i > len)
            {
               if (sgen)
                  free(sgen);
               sgen = malloc(i*sizeof(char));
               assert(sgen);
               len = i;
            }
            sprintf(sgen, "cp AMMCASES/%s %s/%s", mp->genstr, outd, mp->rout);

            if ( (i=system(sgen)) )
            {
               fprintf(stderr, "FAILED CP='%s'\n", sgen);
               exit(i);
            }
         }
         else /* they better have same filename! */
         {
            if (strcmp(mp->rout, p->rout))
            {
               printf("rout=(%s,%s)!\n", mp->rout, p->rout);
               exit(1);
            }
         }
      }
   }
   printf("DONE GENERATING AMM KERNS.\n");
   free(sgen);
}

void GenMake(char pre, char *outd, ATL_mmnode_t *mb)
/*
 * mb files have already been made at least compile-time unique (same source
 * file might occur multiple times due to need to compile with -DKB)
 */
{
   FILE *fp;
   char *fn;
   ATL_mmnode_t *mp;
   char *sals[3] = {"1", "N1", "X"};
   char als[3] = {'1', 'n', 'X'};
   char *sbes[4] = {"0", "1", "N1", "X"};
   char bes[4] = {'0', '1', 'n', 'X'};  /* use 1st 3 for mmkerns */
   char ctas[4] = {'N', 'T', 'C', 'H'};
   char dcomp[8] = {'$', '(', 'D', 'M', 'C', ')', '\0'};
   char dflags[12] = {'$','(','D','M','C','F','L','A','G','S',')','\0'};


   fn = malloc(strlen(outd)+11);
   assert(fn);
   sprintf(fn, "%s/%cMake_amm", outd, pre);
   fp = fopen(fn, "w");
   assert(fp);
   free(fn);
   fprintf(fp, "include ../Make.inc\n");
   fprintf(fp, "CDEFS2=$(CDEFS)\n\n");
/*
 * Spew out kernels to be compiled
 */
   fprintf(fp, "objs =");
   for (mp=mb; mp; mp = mp->next)
   {
      int ib;
      for (ib=0; ib < 3; ib++)
         fprintf(fp, " \\\n       %s_b%c.o", mp->auth, bes[ib]);
   }
   fprintf(fp, "\n");
/*
 * library make targets
 */
   fprintf(fp, "\n\nlib : %clib.grd\nall : %clib.grd\n%clib : %clib.grd\n",
           pre, pre, pre, pre);
   fprintf(fp, "%clib.grd : $(objs)\n", pre);
   fprintf(fp, "\t$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(objs)\n");
   fprintf(fp, "\t $(RANLIB) $(ATLASlib)\n");
   fprintf(fp, "\t touch %clib.grd\n", pre);
   fprintf(fp, "clean : %cclean\n", pre);
   fprintf(fp, "%cclean:\n\t- rm -f $(objs)\n", pre);
   fprintf(fp, "killall : %ckillall\n", pre);
   fprintf(fp, "%ckillall : %cclean\n", pre, pre);
   fprintf(fp, "\t- $(ARCHIVER) d $(ATLASlib) $(objs)\n");
   fprintf(fp, "\t $(RANLIB) $(ATLASlib)\n");
   fprintf(fp, "\t- rm -f ATL_%c*.[S,c]\n\n", pre);
/*
 * Make targets for amm kerns
 */
   fprintf(fp, "#\n#  AMM kernel rules\n#\n");
   fn = (pre == 'd' || pre == 'z') ?  "-DDREAL=1" : "-DSREAL";
   for (mp=mb; mp; mp = mp->next)
   {
      int ib;
      char *comp, *flgs;

      comp = GetMMKernComp(mp, dcomp, dflags, &flgs);
      for (ib=0; ib < 3; ib++)
      {
         char *sp=" ";
         fprintf(fp, "%s_b%c.o : %s\n", mp->auth, bes[ib], mp->rout);
         fprintf(fp, "\t%s $(CDEFS2) %s -DBETA%s=1 \\\n", comp, fn, sbes[ib]);
         if (!FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
            fprintf(fp, "        -DMB=%d -DNB=%d -DKB=%d",
                    mp->mbB, mp->nbB, mp->kbB);
         else
            sp = "        ";
         if (FLAG_IS_SET(mp->flag, MMF_MVA))
         {
            fprintf(fp, "%s-DATL_MOVEA", sp);
            sp = " ";
         }
         if (FLAG_IS_SET(mp->flag, MMF_MVB))
         {
            fprintf(fp, "%s-DATL_MOVEB", sp);
            sp = " ";
         }
         if (FLAG_IS_SET(mp->flag, MMF_MVC))
         {
            fprintf(fp, "%s-DATL_MOVEC", sp);
            sp = " ";
         }
         fprintf(fp, " \\\n        %s \\\n", flgs);
         fprintf(fp,
                 "        -DATL_USERMM=%s_b%c \\\n        -c -o %s_b%c.o \\\n",
                 mp->auth, bes[ib], mp->auth, bes[ib]);
         fprintf(fp, "        %s\n", mp->rout);
      }
   }
   fclose(fp);
}

void GenAllFiles(char pre, char *outd, ATL_mmnode_t **lists)
{
   FILE *fp;
   int i, RCSAME;
   ATL_mmnode_t *rb=NULL, *ib=NULL, *mp;
   ATL_cpnode_t *cpA=NULL, *cpC=NULL, *ncp;
   const char cpr = (pre == 'd') ? 'z' : 'c';
/*
 * First, generate performance files, which require routs separated by
 * various lists
 */
   printf("\nGENERATING HEADER FILES\n");
   GenAllHeaders(pre, outd, lists);
   printf("DONE HEADER GENERATION\n");

/*
 * Generating kernels & Makefile do not care which list kernels came from,
 * or their order, so combine them all into one list for each data type
 * (2 lists total: real,complex), and remove any redundancies.
 * NOTE: rank-K should be distinguished by flag setting (MMF_MV[A,B,C]),
 * and so doesn't need protection in new naming system.
 * NOTE2: I think we should generate copy routs at this time and then
 *        can pass info to GenAllKerns & GenMake.
 */
   for (i=0; i < ITRSM; i++)
   {
      rb = AddUniqueMMKernCompList(rb, lists[i]);
      KillAllMMNodes(lists[i]);
      ib = AddUniqueMMKernCompList(ib, lists[NMMLISTS+i]);
      KillAllMMNodes(lists[NMMLISTS+i]);
   }
/*
 * Now, kerns are compiled the same for real and complex, so reduce real
 * and complex to one unique list
 */
   rb = AddUniqueMMKernCompList(rb, ib);
   KillAllMMNodes(ib);
/*
 * SYRK never same as GEMM, so just add syrk kerns to list, but don't dup if
 * real/cplx are same kernel other than blocking
 */
   ib = lists[ISYRK];
   if (ib)
   {
      ATL_mmnode_t *p = lists[NMMLISTS+ISYRK];
      assert(p);
      if (!MMKernsPerfSame(lists[ISYRK], p))
      {
         ib->next = lists[NMMLISTS+ISYRK];
         ib->next->next = rb;
      }
      else
      {
         ib->next = rb;
         KillAllMMNodes(p);
      }
      rb = ib;
   }
   GenAllKerns(pre, outd, rb);
   GenMake(pre, outd, rb);
   KillAllMMNodes(rb);
   KillAllCopyNodes(cpA);
   KillAllCopyNodes(cpC);
}

int main(int nargs, char **args)
{
   char *outd;
   ATL_mmnode_t *lists[2*NMMLISTS];
   int i;
   char pre, cpr;

   outd = GetFlags(nargs, args, &pre, lists);
   cpr = (pre == 'd') ? 'z' : 'c';
/*
 * Prep file for generation.  Free present values, and replace with:
 * ->auth  : kernel name without _b[1,n,0] suffix
 * ->genstr: for ID=0: genstr, else user kernel name (came in ->rout)
 * ->rout  : correct present filename
 * TRSM/SYRKT entries we just set all these to NULL, since not needed.
 * TRSM/SYRKT used to [ts,sk]amm_perf.h, needs perf,NB,IDX
 */
   for (i=0; i < 2*NMMLISTS; i++)
   {
      if (i == ITRSM || i == ITRSM+NMMLISTS ||
          i == ISYRKT || i == ISYRKT+NMMLISTS)
         KillAllMMStrings(lists[i]);
      else if (i == ISYRK || i == ISYRK+NMMLISTS)
         PrepMMForGen((i==ISYRK)?pre:cpr, outd, "syrk", lists[i]);
      else
         PrepMMForGen(pre, outd, "amm", lists[i]);
   }

   GenAllFiles(pre, outd, lists);

   free(outd);
   return(0);
}
