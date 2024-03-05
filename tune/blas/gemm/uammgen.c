/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_mmgen.h"
#include "atlas_sys.h"

#define NMMLISTS 2
#define IUSR   0
#define IUSRK1 1

   static int UID=0, UIL=1;
   static char unam[16];

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
   fprintf(stderr,
      "   -I <ID>: unique non-negative ID for header/kern files\n");
   fprintf(stderr, "   -ub [r,c] <user-case kernel file> \n");
   fprintf(stderr, "   -uk [r,c] <user-case KClean file (if needed)>\n");
   exit(ierr ? ierr : -1);
}


/*
 * RETURNS: precision prefix [s,d,c,z]
 */
char *GetFlags(int nargs, char **args, char *PRE, ATL_mmnode_t **lists)
{
   int i, j=0, n, k;
   int GotAnySumFile = 0;
   char pre='d', cpre, ch;
   char *outd=NULL;
   char *fn[2]={"uAMMRES.sum","uAMMKCLEAN.sum"};

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
      case 'I':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         UID = atoi(args[i]);
         for (k=10; k<=UID; k*=10)
            UIL++;
         break;
      case 'u':
         GotAnySumFile = 1;
         switch(args[i][2])
         {
         case 'k': /* [r,c] <file> */
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            ch = args[i][0];
            if (ch == 'c')
               k = NMMLISTS + 1;
            else if (ch == 'r')
               k = 1;
            else
               PrintUsage(args[0],i, "1st arg to -uk must be 'r' or 'c'!");
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            lists[k] = ReadMMFile(args[i]);
            assert(lists[k]);
            break;
         case 'b': /* [r,c] <file> */
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            ch = args[i][0];
            if (ch == 'c')
               k = NMMLISTS + 0;
            else if (ch == 'r')
               k = 0;
            else
               PrintUsage(args[0],i, "1st arg to -ub must be 'r' or 'c'!");
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            lists[k] = ReadMMFile(args[i]);
            assert(lists[k]);
            break;
         default:
            PrintUsage(args[0], i, args[i]);
         }
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
   sprintf(unam, "u%d", UID);
/*
 * Fill in any required list entry
 */
   for (i=0; i < 2*NMMLISTS; i++)
   {
      const char pr=(i<NMMLISTS) ? pre : cpre;
      const int k = (i<NMMLISTS) ? i : i-NMMLISTS;
      if (!lists[i] && (!GotAnySumFile || (i&1==0)))
         lists[i] = ReadMMFileWithPath(pr, "res", fn[k]);
      if (!lists[i])
         fprintf(stderr, "CANNOT FIND FILE 'res/%c%s'!\n", pr, fn[k]);
      if (i&1 == 0) assert(lists[i]);
   }
   return(outd);
}


double GenPerfH(char pre, char *outd, char *cn, double mfMax,
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

   fp = OpenMMGenHeader(outd, 0, pre, cn, 0, 0, "perf", mb);
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
   (char pre, char *outd, char *cn, ATL_mmnode_t *mb, ATL_mmnode_t *mkb,
    ATL_mmnode_t *ub)
{
   FILE *fp;
   int ib, inxt, iaut;
   char bes[3] = {'0', '1', 'n'};

   inxt = GetOffset(&mb->next, mb);
   iaut = GetOffset(&ub->auth, ub);
   fp = OpenMMGenHeader(outd, 0, pre, cn, 0, 0, "kern", mb);
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
void GenBlkH(char pre, char *outd, char *cn, ATL_mmnode_t *mb)
{
   FILE *fp;
   ATL_mmnode_t *mp;
   int *KUs;
   int i, n, inxt;

   fp = OpenMMGenHeader(outd, 0, pre, cn, 0, 0, "blk", mb);
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

int StrFindReplace(char **org, char *find, char *replace)
{
   char *pre, *post, *final;
   char *orig = *(org);
   int idx = (strstr(orig, find) - orig) / sizeof(char);
   int len = strlen(orig);
   int flen = strlen(find);
   int rlen = strlen(replace);
   if (idx < 0) return(len); /* not found */
   pre = malloc(sizeof(char)*idx + 1);
   post = malloc(sizeof(char)*(len - idx - flen) + 1);
   strncpy(pre, orig, idx);
   pre[idx] = '\0';
   strcpy(post, orig+idx+flen);
   len = len - flen + rlen;
   free(orig);
   orig = malloc(sizeof(char)*(len) + 1);
   strcpy(orig, pre);
   strcat(orig, replace);
   strcat(orig, post);
   *org = orig;
   return(len);
}

void AdjustCPNamesWithID(char pre, ATL_cpnode_t *cb)
{
   char find[24] = "";
   char replace[24] = "";
   ATL_cpnode_t *cp;
   sprintf(find, "ATL_%c", pre);
   sprintf(replace, "ATL_%cu%d", pre, UID);
   for (cp=cb; cp; cp=cp->next)
      cp->rtlen = StrFindReplace(&cp->rout, find, replace);
}

void GenCpyA(char pre, char *outd, char *cn, char dir, char alp,
             ATL_mmnode_t *mb)
/* presently leaving alp='[1,n,X]', cn=[ge,sq,rk], dir='[F,T]' to caller */
{
   FILE *fp;
   const int nt = (pre == 'd' || pre == 's') ? 2 : 4;
   int flag=0, i, j, it, im;
   char tas[4] = {'N', 'T', 'C', 'H'};
   char mats[2] = {'A', 'B'};
   char arrnm[12]={'x', 'x', '2', 'B', 'L', 'K', '_', 'a', alp, '\0'};
   char fn[12];
   int mati, transi;

   if (dir == 'T')
   {
      sprintf(fn, "cm2am_a%c", alp);
      mati = 0; transi = 1;
   }
   else
   {
      sprintf(fn, "am2cm_a%c", alp);
      sprintf(arrnm, "BLK2xx_a%c", alp);
      mati = 4; transi = 5;
   }
   fp = OpenMMGenHeader(outd, 0, pre, cn, 0, 0, fn, mb);
   fn[6] = 't'; fn[7] = '\0';  /* switch fn to type name */

   for (im=0; im < 2; im++)
   {
      arrnm[mati] = mats[im];
      for (it=0; it < nt; it++)
      {
         ATL_cpnode_t *cb, *ucb;
         int i;
         const char ta = tas[it];

         arrnm[transi] = ta;
         cb = GetMMCopyFromMMNodes(CopyEncode(pre, dir, mats[im], ta), mb);
         assert(cb);
         AdjustCPNamesWithID(pre, cb);
         ucb = AddUniqueCopyNode(NULL, cb);
         arrnm[6] = '\0';
         PrintMMCpProtosA(fp, arrnm, ucb, dir, alp);
         KillAllCopyNodes(ucb);
         PrintStrArrAtOff(fp, pre, arrnm, cb, GetOffset(&cb->next, cb),
                          GetOffset(&cb->rout, cb), fn, ta, alp, 0);
         KillAllCopyNodes(cb);
      }
   }
   CloseGenHeader(fp);
}

void GenCpyC(FILE *fp, char pre, char *outd, char *cn, char dir, char *type,
             char alp, char bet, ATL_mmnode_t *mb)
{
   char *arrnm = (dir == 'T') ? "C2BLK" : "BLK2C";
   ATL_cpnode_t *cb, *ucb;
   int flag=0, i, j;

   cb = GetMMCopyFromMMNodes(CopyEncode(pre, dir, 'c', 'n'), mb);
   assert(cb);
   AdjustCPNamesWithID(pre, cb);
   #ifdef Debug
      i = CountListEntries(cb, GetOffset(&cb->next, cb));
      j = CountListEntries(mb, GetOffset(&mb->next, mb));
      fprintf(stderr, "ncpy=%d, nmm=%d!\n", i, j);
   #endif
   ucb = AddUniqueCopyNode(NULL, cb);
   PrintMMCpProtosC(fp, arrnm, ucb, dir, alp, bet);
   KillAllCopyNodes(ucb);
   PrintStrArrAtOff(fp, pre, arrnm, cb, GetOffset(&cb->next, cb),
                    GetOffset(&cb->rout, cb), type, 0, alp, bet);
   KillAllCopyNodes(cb);
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

void GenFullCpyC(char pre, char *outd, char dir, char *cn, ATL_mmnode_t *mb)
{
   ATL_mmnode_t *mp;
   FILE *fp;
   int ib, NEEDSZF=0;
   char rt[16];
   char bets[4] = {'0', '1', 'n', 'X'};

/*
 * Create atlas_<pre><cn>szC.h if necessary
 */
   for (mp=mb; mp; mp = mp->next)
      if (FLAG_IS_SET(mp->flag, MMF_KVEC) && ib == -1 &&
          (mp->mu*mp->nu) % mp->vlen != 0)
         break;
   if (mp)
   {
      ATL_cpnode_t *cb, *ucb, *cp;
      int n;
      NEEDSZF=1;
      fp = OpenMMGenHeader(outd, 0, pre, cn, 0, 0, "szC", mb);
      cb = GetMMCopyFromMMNodes(CopyEncode(pre, dir, 'c', 'n'), mb);
      assert(cb);
      AdjustCPNamesWithID(pre, cb);
      ucb = AddUniqueCopyNode(NULL, cb);
      KillAllCopyNodes(cb);
/*
 *    Generate needed funcs for size query
 */
      for (cp=ucb; cp; cp = cp->next)
      {
         const unsigned kvec=cp->kvec;
         const unsigned int mu=cp->mu, nu=cp->nu, b = mu*nu;
         const unsigned int B = (kvec > 1) ? ((b+kvec-1)/kvec)*kvec : 0;
         fprintf(fp, "#define uint unsigned int\n");
         if (b != B)
         {
            assert(!cp->ID);
            fprintf(fp, "static uint ATL_szC_%dx%d(uint M, uint N, "
                        "uint mu, uint nu, uint vlen)", mu, nu);
            fprintf(fp, "\n{\n");
            fprintf(fp, "   return(");
            PrintCeilMul(fp, "M", mu);
            fprintf(fp, "*");
            PrintCeilMul(fp, "N", nu);
            fprintf(fp, "*");
            fprintf(fp, "%d);\n", B);
            fprintf(fp, "}\n");
         }
         fprintf(fp, "typedef uint (*szC_t)(uint,uint,uint,uint,uint);\n");
         fprintf(fp, "#undef uint\n");
         for (n=0,cp=cb; cp; n++, cp = cp->next);
         fprintf(fp, "static szC_t ATL_AMM_SZCW[%d] =\n{\n", n);
         for (cp=cb,n=0; cp; cp = cp->next,n++)
         {
            const unsigned kvec=cp->kvec;
            const unsigned int mu=cp->mu, nu=cp->nu, b = mu*nu;
            const unsigned int B = (kvec > 1) ? ((b+kvec-1)/kvec)*kvec : 0;
            fprintf(fp, "/* IDX=%d */ ", n);
            if (b != B)
               fprintf(fp, "ATL_szC_%dx%d", mu, nu);
            else
               fprintf(fp, "NULL");
            if (cp->next)
               fprintf(fp, ",\n");
            else
               fprintf(fp, "\n");
         }
      }
      KillAllCopyNodes(cb);
      CloseGenHeader(fp);
   }
#if 0
#endif

   if (dir == 'T')
      strcpy(rt, "cmat2ablk");
   else
      strcpy (rt, "ablk2cmat");

   fp = OpenMMGenHeader(outd, 0, pre, cn, 0, 0, rt, mb);
   if(NEEDSZF)
   {
      fprintf(fp, "#ifndef NOSZC\n");
      fprintf(fp, "   #include \"atlas_%c%sszC.h\"\n", pre, cn);
      fprintf(fp, "#endif\n");
   }
   rt[9] = '_'; rt[10] = 't'; rt[11] = '\0';
   for (ib=0; ib < 4; ib++)
   {
      int ia;
      char alps[3] = {'1', 'n', 'X'};
      const char bet = bets[ib];

      for (ia=0; ia < 3; ia++)
         GenCpyC(fp, pre, outd, cn, dir, rt, alps[ia], bet, mb);
   }
   CloseGenHeader(fp);
}

void GenAllCpyC(char pre, char *outd, char dir, ATL_mmnode_t *ub)
{
   if (ub)
      GenFullCpyC(pre, outd, dir, unam, ub);
}
void GenAllCpyA(char pre, char *outd, char dir, ATL_mmnode_t *ub)
/* presently leaving dir='[F,T]' to caller */
{
   int ia;
   char alps[3] = {'1', 'n', 'X'};

   for (ia=0; ia < 3; ia++)
   {
      if (ub)
         GenCpyA(pre, outd, unam, dir, alps[ia], ub);
   }
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

void GenSumH(char pre, char *outd, ATL_mmnode_t *ub)
/*
 * cn:[u%d]
 * For [u%d] report: ncases, max[nb,mb,kb]
 * -> for maxB, report: perf,time,ratio wt best user-case
 */
{
   FILE *fp;
   double mf;
   int i;

   fp = OpenMMGenHeader(outd, 0, pre, unam, 0, 0, "sum", ub);
   mf = PrintSum(fp, unam, ub, 0.0);
   CloseGenHeader(fp);
}

void GenAllHeaders(char pre, char *outd, ATL_mmnode_t **lists)
{
   ATL_cpnode_t *cb;
   ATL_mmnode_t *mb;
   char *nm="perf";
   double mfMax;
   int i, k, SKSAME;
   char pr;

/*
 * Handle [perf,blk,flag,sum] which require only ordered lists
 */
   for (pr=pre,k=0; k < 2*NMMLISTS; k += NMMLISTS)
   {
      mfMax = GenPerfH(pr, outd, unam, 0.0, lists[k+IUSR], lists[k+IUSRK1]);

      GenBlkH(pr, outd, unam, lists[k+IUSR]);
      pr = (pre == 'd') ? 'z' : 'c';
   }
   GenSumH(pre, outd, lists[IUSR]);
   GenSumH(pr, outd, lists[IUSR+NMMLISTS]);
/*
 * Create lists of only the unique kernels for [ub] in real & cplx
 * for use in prototyping.  We combine kern & K1 lists for each.
 * We can now gen [kern].
 */
   for (pr=pre,k=0; k < 2; k++)
   {
      for (i=0; i < 1; i++)
      {
         ATL_mmnode_t *ulst;  /* unordered list of only unique kernels */
         const int h = i+i+k*NMMLISTS, u=k+i;
         ulst = AddUniqueMMKernCompList(NULL, lists[h]);
         ulst = AddUniqueMMKernCompList(ulst, lists[h+1]);
         GenKernH(pr, outd, unam, lists[h], lists[h+1], ulst);
         KillAllMMNodes(ulst);  /* done with these! */
      }
      pr = (pre == 's') ? 'c' : 'z';
   }
/*
 * [cmat2ablk_a[1,n,X]_b[0,1,n,X],ablk2cmat_a?_b?,am2cm_a[1,n,X]]
 * -> have never genned, can we?: am2rm_a[1,n,X],
 */
   GenAllCpyA(pre, outd, 'T', lists[IUSR]);
   GenAllCpyA(pre, outd, 'F', lists[IUSR]);
   GenAllCpyA(pr, outd, 'T', lists[IUSR+NMMLISTS]);
   GenAllCpyA(pr, outd, 'F', lists[IUSR+NMMLISTS]); /* for now */

   GenAllCpyC(pre, outd, 'F', lists[IUSR]);
   GenAllCpyC(pre, outd, 'T', lists[IUSR]);
   GenAllCpyC(pr, outd, 'F', lists[IUSR+NMMLISTS]);
   GenAllCpyC(pr, outd, 'T', lists[IUSR+NMMLISTS]); /* for now */
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

void GenAllKerns(char pre, char *outd, ATL_mmnode_t *rb,
                 ATL_cpnode_t *cpA, ATL_cpnode_t *cpC)
{
   ATL_cpnode_t *last, *cp;
   ATL_mmnode_t *mp;
   char *sgen;
   int i, len, dlen, mID, mU;
   char pr=pre;
/*
 * join all copy routs into one list temporarily
 */
   last = FindLastCopyNode(cpA);
   last->next = cpC;
/*
 * Get a string of max length for all system calls
 */
   mID=mU=len=i=0;
   for (cp=(!i)?cpA:cpC; cp; cp = cp->next)
   {
      len = Mmax(len, cp->rtlen);
      mID = Mmax(mID, cp->ID);
      mU = Mmax(mU, cp->mu);
      mU = Mmax(mU, cp->nu);
      mU = Mmax(mU, cp->kvec);
   }
   dlen = strlen(outd);
   len += dlen + 1;  /* outd */
   len += 17;  /* " > /dev/null 2>&1" */
   assert(mID == 0);  /* just until we support user copies */
   len += 64 - 2*4;  /* C format string */
   len += 4*NumDecDigits(mU); /* mu,nu,kvec */
   len += 2;                  /* extra for syrk */
   sgen = malloc(len+1);
   assert(sgen);
/*
 * Generate all needed copy routines
 */
   printf("\nGENERATING MMCOPY FILES\n");
   for (cp=cpA; cp; cp = cp->next)
   {
      int flag = cp->flag;
      int k;
      char pr;
      pr = CopyGetPre(cp->flag);
      printf("   -> %s\n", cp->rout);
      if (flag & (1<<CPF_CBLK))
      {
         unsigned int kvec = cp->kvec;
         if (flag & (1<<CPF_SYRK))
         {
            k = sprintf(sgen,
"make genall_syblk2C pre=%c vec=NA vlen=%d mu=%d nu=%d cpvlen=1 rt=%s/%s.c",
                        pr, kvec, cp->mu, cp->nu, outd, cp->rout);
         }
         else if (flag & (1<<CPF_CBLK))
         {
            k = sprintf(sgen,
      "make genall_%s pre=%c vec=NA vlen=%d mu=%d nu=%d cpvlen=1 rt=%s/%s.c",
                        (flag&(1<<CPF_TOBLK)) ? "C2blk":"blk2C", pr,
                        kvec, cp->mu, cp->nu, outd, cp->rout);
         }
      }
      else /* generate A/B copy */
      {
         unsigned int kvec = cp->kvec;
         char *targ;
         if (flag & (1<<CPF_REAL))
            targ = (flag & (1<<CPF_TOBLK)) ? "A2blk":"blk2A";
         else
            targ = (flag & (1<<CPF_TOBLK)) ? "cA2blk":"cblk2A";
         k = sprintf(sgen,
             "make genall_%s pre=%c kmaj=%d vlen=%d UR=%d cpvlen=1 rt=%s/%s.c",
                     targ, pr, kvec, kvec ? kvec:1,
                     cp->nu, outd, cp->rout);
      }
      assert(k < len);
      k += sprintf(sgen+k, " ");
      assert(k <= len);
      if ((k=system(sgen)))
      {
         fprintf(stderr, "\n\ncpgenstr='%s' returns %d!\n\n", sgen, k);
         exit(k);
      }
   }
   printf("DONE GENERATING MMCOPY FILES\n");
   last->next = NULL;  /* disconnect A/C lists */
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

void GenMake(char pre, char *outd, ATL_mmnode_t *mb,
             ATL_cpnode_t *cpA, ATL_cpnode_t *cpC)
/*
 * mb files have already been made at least compile-time unique (same source
 * file might occur multiple times due to need to compile with -DKB)
 */
{
   FILE *fp;
   char *fn;
   ATL_cpnode_t *cp;
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
   for (cp=cpC; cp; cp = cp->next)  /* C copy routs */
   {
      int ib;
      for (ib=0; ib < 4; ib++)
      {
         int ia;
         for (ia=0; ia < 3; ia++)
            fprintf(fp, " \\\n       %s_a%cb%c.o", cp->rout, als[ia], bes[ib]);
      }
   }

   for (cp=cpA; cp; cp = cp->next)  /* A/B copy routs */
   {
      const int NTA = (cp->flag & (1<<CPF_REAL)) ? 2:4;
      int it;
      for (it=0; it < NTA; it++)
      {
         int ia;
         for (ia=0; ia < 3; ia++)
            fprintf(fp, " \\\n       %s_%ca%c.o", cp->rout, ctas[it], als[ia]);
      }
   }
   fprintf(fp, "\n");
/*
 * library make targets
 */
   fprintf(fp, "\n\nlib : %clib.grd\nall : %clib.grd\n%clib : %clib.grd\n",
           pre, pre, pre, pre);
   fprintf(fp, "%clib.grd : $(objs)\n", pre);
   fprintf(fp, "\t$(ARCHIVER) $(ARFLAGS) $(UAMMlib) $(objs)\n");
   fprintf(fp, "\t $(RANLIB) $(UAMMlib)\n");
   fprintf(fp, "\t touch %clib.grd\n", pre);
   fprintf(fp, "clean : %cclean\n", pre);
   fprintf(fp, "%cclean:\n\t- rm -f $(objs)\n", pre);
   fprintf(fp, "killall : %ckillall\n", pre);
   fprintf(fp, "%ckillall : %cclean\n", pre, pre);
   fprintf(fp, "\t- $(ARCHIVER) d $(UAMMlib) $(objs)\n");
   fprintf(fp, "\t $(RANLIB) $(UAMMlib)\n");
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
/*
 * Make targets for C copy routines
 */
   fprintf(fp, "#\n#  C copy rules\n#\n");
   dflags[3] = dcomp[3] = 'K';
   for (cp=cpC; cp; cp = cp->next)  /* C copy routs */
   {
      char *styp;
      int ib;
      char pr;

      styp = CopyGetCompType(cp->flag);
      for (ib=0; ib < 4; ib++)
      {
         int ia;
         for (ia=0; ia < 3; ia++)
         {
            fprintf(fp, "%s_a%cb%c.o : %s.c\n", cp->rout, als[ia], bes[ib],
                    cp->rout);
            fprintf(fp, "\t%s %s -c $(CDEFS) -D%s=1 -DBETA%s=1 -DALPHA%s=1",
                    dcomp, dflags, styp, sbes[ib], sals[ia]);
            fprintf(fp,
            " \\\n           -DATL_USERCPMM=%s_a%c_b%c -DATL_MU=%d -DATL_NU=%d",
                    cp->rout, als[ia], bes[ib], cp->mu, cp->nu);
            fprintf(fp, " \\\n           -o %s_a%cb%c.o %s.c\n",
                    cp->rout, als[ia], bes[ib], cp->rout);
         }
      }
   }
/*
 * Make targets for A/B copy routines
 */
   fprintf(fp, "#\n#  A/B copy rules\n#\n");
   for (cp=cpA; cp; cp = cp->next)
   {
      char *styp, *cnj[2] = {"", "-DConj_=1"};
      const int flag = cp->flag;
      const int NC = (flag&(1<<CPF_REAL)) ? 1 : 2;
      int ic;
      char pr;
      styp = CopyGetCompType(cp->flag);
      for (ic=0; ic < NC; ic++)
      {
         int it;
         for (it=0; it < 2; it++)
         {
            int ia;
            const int itc = ic*NC+it;
            for (ia=0; ia < 3; ia++)
            {
               fprintf(fp, "%s_%ca%c.o : %s.c\n", cp->rout, ctas[itc],
                       als[ia], cp->rout);
               fprintf(fp, "\t%s %s -c $(CDEFS) -D%s=1 -DALPHA%s=1",
                       dcomp, dflags, styp, sals[ia]);
               fprintf(fp, " \\\n           -DATL_NU=%d -DTRANS%c_=1 %s",
                       cp->nu, ctas[it],  cnj[ic]);
               fprintf(fp, " \\\n           -DATL_USERCPMM=%s_%ca%c", cp->rout,
                       ctas[itc], als[ia]);
               fprintf(fp, " \\\n           -o %s_%ca%c.o %s.c\n",
                       cp->rout, ctas[itc], als[ia], cp->rout);
            }
         }
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
 * NOTE2: I think we should generate copy routs at this time and then
 *        can pass info to GenAllKerns & GenMake.
 */
   for (i=0; i < NMMLISTS; i++)
   {
      rb = AddUniqueMMKernCompList(rb, lists[i]);
      KillAllMMNodes(lists[i]);
      ib = AddUniqueMMKernCompList(ib, lists[NMMLISTS+i]);
      KillAllMMNodes(lists[NMMLISTS+i]);
   }
/*
 * Now create lists of copy routines for both real & complex
 */
   GetMMAllUniqueCopyFromMMNodes(pre, 'F', 'T', rb, &cpA, &cpC);
   GetMMAllUniqueCopyFromMMNodes(pre, 'T', 'F', rb, &cpA, &cpC);
   GetMMAllUniqueCopyFromMMNodes(cpr, 'F', 'T', ib, &cpA, &cpC);
   GetMMAllUniqueCopyFromMMNodes(cpr, 'T', 'F', ib, &cpA, &cpC); /* for now */
/*
 * Now, kerns are compiled the same for real and complex, so reduce real
 * and complex to one unique list
 */
   rb = AddUniqueMMKernCompList(rb, ib);
   KillAllMMNodes(ib);
   AdjustCPNamesWithID(pre, cpA);
   AdjustCPNamesWithID(cpr, cpA);
   AdjustCPNamesWithID(pre, cpC);
   AdjustCPNamesWithID(cpr, cpC);
   GenAllKerns(pre, outd, rb, cpA, cpC);
   GenMake(pre, outd, rb, cpA, cpC);
   KillAllMMNodes(rb);
   KillAllCopyNodes(cpA);
   KillAllCopyNodes(cpC);
}

void AdjustMMNamesWithID(char pre, ATL_mmnode_t *mb)
{
   char find[24] = "";
   char replace[24] = "";
   ATL_mmnode_t *mp;
   sprintf(find, "ATL_%c", pre);
   sprintf(replace, "ATL_%cu%d", pre, UID);
   for (mp=mb; mp; mp=mp->next)
      StrFindReplace(&mp->auth, find, replace);
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
 */
   for (i=0; i < 2*NMMLISTS; i++)
   {
      PrepMMForGen(pre, outd, "amm", lists[i]);
   }
   for (i=0; i < 2*NMMLISTS; i++)
      AdjustMMNamesWithID(pre, lists[i]);

   GenAllFiles(pre, outd, lists);

   free(outd);
   return(0);
}
