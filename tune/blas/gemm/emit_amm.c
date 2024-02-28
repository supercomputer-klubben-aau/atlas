/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015, 2013, 2012 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_mmparse.h"
#include "atlas_mmtesttime.h"
#include "atlas_sys.h"

ATL_mmnode_t *RemoveNonUniqueKernels(ATL_mmnode_t *mmb);

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set type/precision prefix (d) \n");
   fprintf(stderr, "   -d <outdir>: directory to dump files to\n");
   fprintf(stderr, "   -g <rect kern file> <geKCleanFile>\n");
   fprintf(stderr, "   -r <rank-K kernel file> \n");
   fprintf(stderr, "   -s <square-case kernel file> <sqKCleanFile>\n");
   fprintf(stderr, "   -t <trsm file> : should match -s\n");
   exit(ierr ? ierr : -1);
}

/*
 * RETURNS: precision prefix [s,d,c,z]
 */
char GetFlags(int nargs, char **args, char **DOUT,
              ATL_mmnode_t **GEMM, ATL_mmnode_t **GEK1,
              ATL_mmnode_t **SQMM, ATL_mmnode_t **SQK1,
              ATL_mmnode_t **RKMM, ATL_mmnode_t **TRSM,
              ATL_mmnode_t **cGEMM, ATL_mmnode_t **cGEK1,
              ATL_mmnode_t **cSQMM, ATL_mmnode_t **cSQK1,
              ATL_mmnode_t **cRKMM, ATL_mmnode_t **cTRSM)
{
   int i, j=0, n, k;
   char pre='d', cpre;
   ATL_mmnode_t *gemm=NULL, *gek1=NULL, *sqmm=NULL, *sqk1=NULL, *rkmm=NULL;
   ATL_mmnode_t *trsm=NULL, *mmb=NULL, *mmp, *mp;

   *DOUT = NULL;
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
      case 'g':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        gemm = ReadMMFile(args[i]);
        assert(gemm);
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        gek1 = ReadMMFile(args[i]);
        assert(gek1);
        break;
      case 't':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        trsm = ReadMMFile(args[i]);
        assert(trsm);
      case 's':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        sqmm = ReadMMFile(args[i]);
        assert(sqmm);
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        sqk1 = ReadMMFile(args[i]);
        assert(sqk1);
        break;
      case 'r':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        rkmm = ReadMMFile(args[i]);
        break;
      case 'd':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *DOUT = DupString(args[i]);
        break;
      case 'i':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        mmp = ReadMMFile(args[i]);
        if (mmb)
        {
           ATL_mmnode_t *mp;
           for (mp=mmb; mp->next; mp = mp->next);
           mp->next = mmp;
        }
        else
           mmb = mmp;
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!gemm)
   {
      assert(!gek1);
      gemm = ReadMMFileWithPath(pre, "res", "geAMMRES.sum");
      assert(gemm);
      gek1 = ReadMMFileWithPath(pre, "res", "geAMMKCLEAN.sum");
      assert(gek1);
   }
   if (!sqmm)
   {
      assert(!sqk1);
      sqmm = ReadMMFileWithPath(pre, "res", "sqAMMRES.sum");
      assert(sqmm);
      sqk1 = ReadMMFileWithPath(pre, "res", "sqAMMKCLEAN.sum");
      assert(sqk1);
   }
   if (!trsm)
   {
      trsm = ReadMMFileWithPath(pre, "res", "tsAMMRES.sum");
      assert(trsm);
   }
   if (!rkmm)
   {
      rkmm = ReadMMFileWithPath(pre, "res", "rkAMMRES.sum");
      assert(rkmm);
   }
   if (!(*DOUT))
      *DOUT = DupString("KERNEL");
   cpre = (pre == 'd') ? 'z' : 'c';
   *cGEMM = ReadMMFileWithPath(cpre, "res", "geAMMRES.sum");
   *cGEK1 = ReadMMFileWithPath(cpre, "res", "geAMMKCLEAN.sum");
   *cSQMM = ReadMMFileWithPath(cpre, "res", "sqAMMRES.sum");
   *cSQK1 = ReadMMFileWithPath(cpre, "res", "sqAMMKCLEAN.sum");
   *cRKMM = ReadMMFileWithPath(cpre, "res", "rkAMMRES.sum");
   *cTRSM = ReadMMFileWithPath(cpre, "res", "tsAMMRES.sum");
   *GEMM = gemm;
   *GEK1 = gek1;
   *SQMM = sqmm;
   *SQK1 = sqk1;
   *RKMM = rkmm;
   *TRSM = trsm;
   return(pre);
}

char *GetVecStr(char pre, int vlen)
{
   if (vlen == 1)
      return("scalar");
   #ifdef ATL_AVX
      if (pre == 'd' || pre == 'z')
      {
         if (vlen == 4)
            return("avx");
         else if (vlen == 2)
            return("sse");
      }
      else if (pre == 's' || pre == 'c')
      {
         if (vlen == 8)
            return("avx");
         else if (vlen == 4)
            return("sse");
      }
   #elif defined(ATL_SSE1)
      #ifdef ATL_SSE2
         if ((pre == 'd' || pre == 'z') && vlen == 2)
               return("sse");
      #endif
      if ((pre == 's' || pre == 'c') && vlen == 4)
         return("sse");
   #endif
/*
 * Any vector length > 1 that isn't one of our known cases uses gnuvec
 */
   return("gvec");
}

void PrintBegBlock(char pre, ATL_mmnode_t *mmb, char *pfx, char *nam, FILE *fp)
{
   ATL_mmnode_t *mp;
   char PRE = toupper(pre);
   int i;

   if (nam)
   {
      if (pfx)
         fprintf(fp,
         "#ifndef ATLAS_%c%sAMM_%s_H\n   #define ATLAS_%c%sAMM_%s_H\n\n",
                 PRE, pfx, nam, PRE, pfx, nam);
      else
         fprintf(fp,"#ifndef ATLAS_%cAMM_%s_H\n   #define ATLAS_%cAMM_%s_H\n\n",
                 PRE, nam, PRE, nam);
      fprintf(fp, "#include \"atlas_amm.h\"\n");
   }
   else
      fprintf(fp, "#ifndef ATLAS_%cAMM_H\n   #define ATLAS_%cAMM_H\n\n",
              PRE, PRE);
/*
 * Count mmb, and print def of NCASES
 */
   for (mp=mmb,i=0; mp; i++, mp = mp->next);
   if (pfx)
      fprintf(fp, "#define ATL_%c%sAMM_NCASES %d\n", pre, pfx, i);
   else if (!nam || strstr(nam, "RANKK") == NULL)
   {
      fprintf(fp, "#ifdef ATL_AMM_NCASES\n");
      fprintf(fp, "   #if ATL_AMM_NCASES != %d\n", i);
      fprintf(fp, "      #error \"NCASES MISMATCH!\"\n");
      fprintf(fp, "   #endif\n");
      fprintf(fp, "#else\n");
      fprintf(fp, "   #define ATL_AMM_NCASES %d\n", i);
      fprintf(fp, "#endif\n");
   }
}

char *GetHName(char pre, char *outd, char *pfx, char *bnam)
{
   int i, NOBASE=0;
   char *fnam;
   if (!bnam)
   {
      NOBASE = 1;
      bnam = "";
   }
   i = strlen(outd) + strlen(bnam) + strlen(pfx) + 16;

   fnam = malloc(i*sizeof(char));
   assert(fnam);
   if (NOBASE)
      sprintf(fnam, "%s/atlas_%c%samm.h", outd, pre, pfx);
   else
      sprintf(fnam, "%s/atlas_%c%samm_%s.h", outd, pre, pfx, bnam);
   return(fnam);
}

FILE *StandHStart(char pre, ATL_mmnode_t *mmb, char *pfx,
                  char *outd,  char *bnam)
{
   char *fnam;
   FILE *fp;
   int i;

   assert(outd);
   fnam = GetHName(pre, outd, pfx, bnam);
   fp = fopen(fnam, "w");
   assert(fp);
   if (bnam)
   {
      for (i=0; bnam[i]; i++)
         fnam[i] = toupper(bnam[i]);
      fnam[i] = '\0';
      if (strstr(bnam, "perf"))
         PrintBegBlock(pre, mmb, pfx, fnam, fp);
      else
         PrintBegBlock(pre, mmb, NULL, fnam, fp);
   }
   else
      PrintBegBlock(pre, mmb, NULL, NULL, fp);
   free(fnam);
   return(fp);
}

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

void GenAmmSum(char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *rkb,
               char *pfx, char *outd)
{
   ATL_mmnode_t *mp, *p66;
   char *fnam;
   FILE *fp;
   char *type = "unsigned short";
   int i, n, maxb, maxNB, maxMB, maxKB, maxkmaj;
   char PRE = toupper(pre), bc[3] = {'M', 'N', 'K'};
   double mfB;

   fp = StandHStart(pre, mmb, pfx, outd, "sum");
   maxkmaj = maxNB = maxKB = maxMB = 0;
   mfB = 0.0;
   for (n=0,mp=mmb; mp; n++, mp = mp->next)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         maxkmaj = Mmax(maxkmaj, mp->vlen);
      maxMB = Mmax(maxMB, mp->mbB);
      maxNB = Mmax(maxNB, mp->nbB);
      maxKB = Mmax(maxKB, mp->kbB);
      mfB = Mmax(mfB, mp->mflop[0]);
   }
   if (maxkmaj == 1)
      maxkmaj = 0;
   maxb = Mmax(maxMB, maxNB);
   maxb = Mmax(maxb, maxKB);
   fprintf(fp, "\n#define ATL_AMM_MAXMB %d\n", maxMB);
   fprintf(fp, "#define ATL_AMM_MAXNB %d\n", maxNB);
   fprintf(fp, "#define ATL_AMM_MAXKB %d\n", maxKB);
   fprintf(fp, "#ifndef  ATL_geAMM_MAXKVEC\n", cn);
   fprintf(fp, "   #define ATL_geAMM_MAXKVEC %d\n\n", maxkmaj);
   fprintf(fp, "#endif\n");

   for (mp=mmb; mp && mp->next; mp = mp->next);
   assert(mp);
   fprintf(fp, "#define ATL_AMM_LMU %d\n", mp->mu);
   fprintf(fp, "#define ATL_AMM_LNU %d\n", mp->nu);
   fprintf(fp, "#define ATL_AMM_LKU %d\n", mp->ku);
   fprintf(fp, "#define ATL_AMM_LLCMMN %d\n\n", Mylcm(mp->mu, mp->nu));
   fprintf(fp, "#define ATL_AMM_LLCMU %d\n\n",
           Mylcm(Mylcm(mp->mu, mp->nu),mp->ku));
/*
 * Find smallest case achieving 2/3 of maximal performance
 */
   for (i=0,mp=mmb; mp && mp->mflop[0]*1.5 < mfB; i++, mp = mp->next);
   assert(mp);
   fprintf(fp, "#define ATL_AMM_66IDX %d\n", i);
   fprintf(fp, "#define ATL_AMM_66MB %d\n", mp->mbB);
   fprintf(fp, "#define ATL_AMM_66NB %d\n", mp->nbB);
   fprintf(fp, "#define ATL_AMM_66KB %d\n", mp->kbB);
   fprintf(fp, "#define ATL_AMM_66LCMMN %d\n\n", Mylcm(mp->mu, mp->nu));
   fprintf(fp, "#define ATL_AMM_66LCMU %d\n\n",
           Mylcm(Mylcm(mp->mu, mp->nu),mp->ku));
   fprintf(fp, "#define ATL_AMM_66RATIO %1.4lf\n\n", mp->mflop[0]/mfB);
/*
 * Find smallest case achieving 98% of maximal performance
 */
   for (i=0,mp=mmb; mp && mp->mflop[0] < 0.98*mfB; i++, mp = mp->next);
   assert(mp);
   fprintf(fp, "#define ATL_AMM_98IDX %d\n", i);
   fprintf(fp, "#define ATL_AMM_98MB %d\n", mp->mbB);
   fprintf(fp, "#define ATL_AMM_98NB %d\n", mp->nbB);
   fprintf(fp, "#define ATL_AMM_98KB %d\n", mp->kbB);
   fprintf(fp, "#define ATL_AMM_98LCMMN %d\n\n", Mylcm(mp->mu, mp->nu));
   fprintf(fp, "#define ATL_AMM_98LCMU %d\n\n",
           Mylcm(Mylcm(mp->mu, mp->nu),mp->ku));
   fprintf(fp, "#define ATL_AMM_98RATIO %1.4lf\n\n", mp->mflop[0]/mfB);
   if (rkb)
   {
      maxkmaj = maxNB = maxKB = maxMB = 0;
      mfB = 0.0;
      for (n=0,mp=rkb; mp; n++, mp = mp->next)
      {
         if (FLAG_IS_SET(mp->flag, MMF_KVEC))
            maxkmaj = Mmax(maxkmaj, mp->vlen);
         maxMB = Mmax(maxMB, mp->mbB);
         maxNB = Mmax(maxNB, mp->nbB);
         maxKB = Mmax(maxKB, mp->kbB);
         mfB = Mmax(mfB, mp->mflop[0]);
      }
      if (maxkmaj == 1)
         maxkmaj = 0;
      fprintf(fp, "#define ATL_rkAMM_LASTMB %d\n", maxMB);
      fprintf(fp, "#define ATL_rkAMM_LASTNB %d\n", maxNB);
      fprintf(fp, "#define ATL_rkAMM_LASTKB %d\n", maxKB);
      fprintf(fp, "#define ATL_MAXKMAJ_RKK %d\n\n", maxkmaj);
/*
 *    Find smallest case achieving 2/3 of maximal performance
 */
      for (i=0,mp=rkb; mp && mp->mflop[0]*1.5 < mfB; i++, mp = mp->next);
      assert(mp);
      fprintf(fp, "#define ATL_66IDX_RKK %d\n", i);
      fprintf(fp, "#define ATL_66NB_RKK %d\n", mp->mbB);
      fprintf(fp, "#define ATL_66MB_RKK %d\n", mp->nbB);
      fprintf(fp, "#define ATL_66KB_RKK %d\n", mp->kbB);
      fprintf(fp, "#define ATL_66RATIO_RKK %1.4lf\n\n", mp->mflop[0]/mfB);
   }
   fprintf(fp, "#define ATL_AMMFLG_KRUNTIME(flg_) ((flg_) & 1)\n");
   fprintf(fp, "#define ATL_AMMFLG_KMAJOR(flg_) ((flg_) & 2)\n");

   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenPerfFile(char pre, ATL_mmnode_t *mmb, char *pfx, char *outd,
                 double mfMax)
{
   ATL_mmnode_t *mp;
   #define NTHRSH 11
   int THRSH[NTHRSH] = {25, 33, 50, 66, 75, 80, 85, 90, 95, 98, 99};
   int idxT[NTHRSH];
   ATL_mmnode_t *mpT[NTHRSH];
   char *fnam;
   FILE *fp;
   char *type = "float";
   int i, j, n, maxb, maxNB, maxMB, maxKB, idxMax=0;
   const int FNDMAX = (mfMax == 0.0);
   char PRE = toupper(pre), bc[3] = {'M', 'N', 'K'};

   fp = StandHStart(pre, mmb, pfx, outd, "perf");

   if (FNDMAX)
   {
      for (i=0; i < NTHRSH; i++)
         mpT[i] = NULL;
      for (n=0,mp=mmb; mp; n++, mp = mp->next)
      {
         if (mp->mflop[0] > mfMax)
         {
            mfMax = mp->mflop[0];
            idxMax = n;
         }
      }
   }
   else
   {
      for (n=0,mp=mmb; mp; n++, mp = mp->next);
      idxMax = -1;
   }
   fprintf(fp, "#define ATL_%sAMM_DENOM %le /* (%.2f)*/ \n", pfx,
           mfMax, mfMax);
   fprintf(fp, "#define ATL_%sAMM_MAXMFLOPIDX %d\n\n", pfx, idxMax);
   if (FNDMAX)
   {
      for (n=0,mp=mmb; mp; mp = mp->next, n++)
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
         fprintf(fp, "#define ATL_AMM_%dLCMU %d\n", THRSH[i],
                 Mylcm(Mylcm(mp->mu,mp->nu),mp->ku));
         fprintf(fp, "#define ATL_AMM_%dLCMMN %d\n", THRSH[i],
                 Mylcm(mp->mu,mp->nu));
         fprintf(fp, "#define ATL_AMM_%dKB %d\n", THRSH[i],
                 mp->kbB);
         fprintf(fp, "#define ATL_AMM_%dNB %d\n", THRSH[i],
                 mp->nbB);
         fprintf(fp, "#define ATL_AMM_%dMB %d\n", THRSH[i],
                 mp->mbB);
         fprintf(fp, "#define ATL_AMM_%dIDX %d\n", THRSH[i],
                 idxT[i]);
      }
      fprintf(fp, "\n");
   }

   fprintf(fp, "static const float ATL_%sAMM_PERF[%d] =", pfx, n);
   fprintf(fp, "   /* %% of performance of best kernel */\n{\n");
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
      fprintf(fp, "   %f%c  /* IDX=%d, KB=%d */\n", mp->mflop[0]/mfMax,
              (mp->next)?',':' ', j, mp->kbB);
   fprintf(fp, "};\n\n");

   fprintf(fp, "static const float ATL_%sAMM_TIME[%d] =", pfx, n);
   fprintf(fp, "   /* %% of performance of best kernel */\n{\n");
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
      fprintf(fp, "   %e%c  /* IDX=%d, KB=%d */\n",
              2.0*mp->mbB*mp->nbB*mp->kbB / (mp->mflop[0]*1.0e6),
              (mp->next)?',':' ', j, mp->kbB);
   fprintf(fp, "};\n\n");

   #if 0   /* found no use for this yet, so don't waste mem wt it */
   fprintf(fp, "static const float ATL_AMM_SPDUPNXT[%d] =", n);
   fprintf(fp, "   /* speedup of next higher NB */\n{\n");
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      double mf = (mp->next) ? mp->next->mflop[0] : mp->mflop[0];
      fprintf(fp, "   %f%c  /* IDX=%d, KB=%d vs. %d */\n", mf/mp->mflop[0],
              (mp->next)?',':' ', j, mp->kbB, mp->next?mp->next->kbB:mp->kbB);
   }
   fprintf(fp, "};\n\n");
   #endif

   fprintf(fp, "#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenBlockingFile(char pre, ATL_mmnode_t *mmb, char *pfx, char *outd)
{
   ATL_mmnode_t *mp;
   char *fnam;
   FILE *fp;
   char *type = "unsigned short";
   int i, n, maxb, maxNB, maxMB, maxKB, maxkmaj;
   char PRE = toupper(pre), bc[3] = {'M', 'N', 'K'};

   fp = StandHStart(pre, mmb, pfx, outd, "blk");
   maxkmaj = maxNB = maxKB = maxMB = 0;
   for (n=0,mp=mmb; mp; n++, mp = mp->next)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         maxkmaj = Mmax(maxkmaj, mp->vlen);
      maxMB = Mmax(maxMB, mp->mbB);
      maxNB = Mmax(maxNB, mp->nbB);
      maxKB = Mmax(maxKB, mp->kbB);
   }
   if (maxkmaj == 1)
      maxkmaj = 0;
   maxb = Mmax(maxMB, maxNB);
   maxb = Mmax(maxb, maxKB);
   fprintf(fp, "#define ATL_AMM_MAXMB %d\n", maxMB);
   fprintf(fp, "#define ATL_AMM_MAXNB %d\n", maxNB);
   fprintf(fp, "#define ATL_AMM_MAXKB %d\n", maxKB);
   fprintf(fp, "#ifndef  ATL_AMM_MAXKVEC\n", cn);
   fprintf(fp, "   #define ATL_AMM_MAXKVEC %d\n", maxkmaj);
   fprintf(fp, "#endif\n");
   fprintf(fp, "\n");

   if (maxb <= 255)
      type = "unsigned char";
   for (i=0; i < 3; i++)
   {
      int j;
      fprintf(fp, "static const %s ATL_AMM_%cBs[%d] =\n{\n",
              type, bc[i], n);
      for (j=0,mp=mmb; mp; j++,mp = mp->next)
      {
         int b;
         if (bc[i] == 'M')
            b = mp->mbB;
         else if (bc[i] == 'N')
            b = mp->nbB;
         else if (bc[i] == 'K')
            b = mp->kbB;
         if (mp->next)
            fprintf(fp, "%8d,  /* index %d */\n", b, j);
         else
            fprintf(fp, "%8d   /* index %d */\n", b, j);
      }
      fprintf(fp, "};\n\n");
   }
   for (i=0; i < 3; i++)
   {
      int j;
      fprintf(fp, "static const %s ATL_AMM_%cUs[%d] =\n{\n",
              type, bc[i], n);
      for (j=0,mp=mmb; mp; j++,mp = mp->next)
      {
         int b;
         if (bc[i] == 'M')
            b = mp->mu;
         else if (bc[i] == 'N')
            b = mp->nu;
         else if (bc[i] == 'K')
         {
            if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
               b = mp->ku;
            else
               b = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->ku : mp->kbB;
         }
         if (mp->next)
            fprintf(fp, "%8d,  /* index %d */\n", b, j);
         else
            fprintf(fp, "%8d   /* index %d */\n", b, j);
      }
      fprintf(fp, "};\n\n");
   }
   fprintf(fp, "static const %s ATL_AMM_KBMINs[%d] =\n{\n", type, n);
   for (i=0,mp=mmb; mp; i++,mp = mp->next)
   {
      if (mp->next)
         fprintf(fp, "%8d,  /* index %d */\n", mp->kbmin, i);
      else
         fprintf(fp, "%8d   /* index %d */\n", mp->kbmin, i);
   }
   fprintf(fp, "};\n\n");
   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenFlagH(char pre, ATL_mmnode_t *mmb, char *pfx, char *outd)
{
   FILE *fp;
   int j, n;
   ATL_mmnode_t *mp;

   fp = StandHStart(pre, mmb, pfx, outd, "flag");

   for (n=0,mp=mmb; mp; n++,mp = mp->next);

   fprintf(fp, "static const unsigned char ATL_AMM_KFLAG[%d] =\n{\n", n);
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      unsigned char flag=FLAG_IS_SET(mp->flag, MMF_KRUNTIME) ? 1 : 0;
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         flag |= 2;
      if (mp->next)
         fprintf(fp, "%6d,  /* index %d */\n", flag, j);
      else
         fprintf(fp, "%6d   /* index %d */\n", flag, j);
   }
   fprintf(fp, "};\n\n");
   fprintf(fp, "#define ATL_AMM_KRUNTIME(idx_) (ATL_AMM_KFLAG[idx_] & 1)\n");
   fprintf(fp, "#define ATL_AMM_KMAJOR(idx_) (ATL_AMM_KFLAG[idx_] & 2)\n");
   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void SpewForthC2MProto(char pre, FILE *fp0, FILE *fp1, int mu, int nu)
{
   char ac[3] = {'1', 'n', 'X'};
   char bc[4] = {'0', '1', 'n', 'X'};
   int ia, ib;
   for (ia=0; ia < 3; ia ++)
   {
      for (ib=0; ib < 4; ib++)
      {
         fprintf(fp0, "void ATL_%cablk2cmat_%dx%d_a%c_b%c\n",
                 pre, mu, nu, ac[ia], bc[ib]);
         if (pre == 'z' || pre == 'c')
            fprintf(fp0, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const TYPE*,const SCALAR,TYPE *,ATL_CSZT);\n");
         else
            fprintf(fp0, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const SCALAR,TYPE *,ATL_CSZT);\n");
         fprintf(fp1, "void ATL_%ccmat2ablk_%dx%d_a%c_b%c\n",
                 pre, mu, nu, ac[ia], bc[ib]);
         if (pre == 'z' || pre == 'c')
            fprintf(fp1, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,TYPE*);\n");
         else
            fprintf(fp1, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,const SCALAR,TYPE*);\n");
      }
   }
}

void SpewForthC2BDecl(char pre, ATL_mmnode_t *mmb, FILE *fp, char *rt,
                      char alp, char bet)
{
   ATL_mmnode_t *mp;
   int j;

   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      fprintf(fp, "   ATL_%c%s_%dx%d_a%c_b%c",
              pre, rt, mp->mu,mp->nu, alp, bet);
      if (mp->next)
         fprintf(fp, ",  /* index %d */\n", j);
      else
         fprintf(fp, "   /* index %d */\n", j);
      }
      fprintf(fp, "};\n\n");
}

void GenC2BLK(char pre, ATL_mmnode_t *mmb, char *pfx, char *outd)
{
   FILE *fp0, *fp1;
   ATL_mmnode_t *mp;
   int ia, ib;
   char ac[3] = {'1', 'n', 'X'};
   char bc[4] = {'0', '1', 'n', 'X'};
   char *fnam;

   fp0 = StandHStart(pre, mmb, pfx, outd, "ablk2cmat");
   fp1 = StandHStart(pre, mmb, pfx, outd, "cmat2ablk");
   fprintf(fp0, "\n");
   fprintf(fp1, "\n");
/*
 * Crank out prototypes
 */
   SpewForthC2MProto(pre, fp0, fp1, mmb->mu, mmb->nu);
   for (mp=mmb->next; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      const int mu=mp->mu, nu=mp->nu;
      for (p=mmb; p->mu != mu || p->nu != nu; p = p->next);
      if (p == mp)  /* first occurance of this mu,nu pair */
         SpewForthC2MProto(pre, fp0, fp1, mp->mu, mp->nu);
   }
   fprintf(fp0, "\n");
   fprintf(fp1, "\n");
/*
 * Now, crank out funcptr arrays
 */
   for (ia=0; ia < 3; ia ++)
   {
      for (ib=0; ib < 4; ib++)
      {
         fprintf(fp0,
            "static const ablk2cmat_t ATL_AMM_BLK2C_a%c_b%c[ATL_AMM_NCASES] =\n{\n",
                 ac[ia], bc[ib]);
         SpewForthC2BDecl(pre, mmb, fp0, "ablk2cmat", ac[ia], bc[ib]);
         fprintf(fp1,
            "static const cmat2ablk_t ATL_AMM_C2BLK_a%c_b%c[ATL_AMM_NCASES] =\n{\n",
                 ac[ia], bc[ib]);
         SpewForthC2BDecl(pre, mmb, fp1, "cmat2ablk", ac[ia], bc[ib]);
      }
   }
   fprintf(fp0, "\n#endif  /* end include file guard */\n");
   fclose(fp0);
   fprintf(fp1, "\n#endif  /* end include file guard */\n");
   fclose(fp1);
}

void SpewForthRevCpProto(char pre, FILE *fp, char alp, int u, int kmaj)
{
   const int G = (pre == 'c' || pre == 'z') ? 2 : 1;
   const char *cst[2] = {"", "C"};
   int g;

   for (g=0; g < G; g++)
   {
      if (kmaj > 1)
         fprintf(fp, "void ATL_%cam2cm_a%c_%dx%d%s\n",
                 pre, alp, kmaj, u, cst[g]);
      else
         fprintf(fp, "void ATL_%cam2cm_a%c_%d%s\n",pre, alp, u, cst[g]);
      if (pre == 'z' || pre == 'c')
         fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*,const TYPE*);\n");
      else
         fprintf(fp,
         "   (ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*);\n");
      if (kmaj > 1)
         fprintf(fp, "void ATL_%cam2rm_a%c_%dx%d%s\n",
                 pre, alp, kmaj, u, cst[g]);
      else
         fprintf(fp, "void ATL_%cam2rm_a%c_%d%s\n",pre, alp, u, cst[g]);
      if (pre == 'z' || pre == 'c')
         fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*,const TYPE*);\n");
      else
         fprintf(fp,
         "   (ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*);\n");
   }
}

void SpewForthCpProto(char pre, FILE *fp, char alp, int u, int kmaj)
{
   const int G = (pre == 'c' || pre == 'z') ? 2 : 1;
   const char *cst[2] = {"", "C"};
   int g;

   for (g=0; g < G; g++)
   {
      if (kmaj > 1)
         fprintf(fp, "void ATL_%ccm2am_a%c_%dx%d%s\n",
                 pre, alp, kmaj, u, cst[g]);
      else
         fprintf(fp, "void ATL_%ccm2am_a%c_%d%s\n", pre, alp, u, cst[g]);
      if (pre == 'z' || pre == 'c')
         fprintf(fp,
    "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*,TYPE*);\n");
      else
         fprintf(fp,
         "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*);\n");
      if (kmaj > 1)
         fprintf(fp, "void ATL_%crm2am_a%c_%dx%d%s\n",
                 pre, alp, kmaj, u, cst[g]);
      else
         fprintf(fp, "void ATL_%crm2am_a%c_%d%s\n", pre, alp, u, cst[g]);
      if (pre == 'z' || pre == 'c')
         fprintf(fp,
     "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*,TYPE*);\n");
      else
         fprintf(fp,
         "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*);\n");
   }
}

void SpewForthCpConjDecl(char pre, int REVERSE, ATL_mmnode_t *mmb, FILE *fp,
                         char *arr, char *rt, char alp, int u)
{
   ATL_mmnode_t *mp;
   int j;

   fprintf(fp, "static const %s_t %s_a%c[%d] =\n{\n",
           REVERSE?"am2cm":"cm2am", arr, alp,  ATL_CountNumberOfMMNodes(mmb));
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      const int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
      if (kmaj > 1)
         fprintf(fp, "   ATL_%c%s_a%c_%dx%dC", pre, rt, alp, kmaj,
                 u?mp->mu:mp->nu);
      else
         fprintf(fp, "   ATL_%c%s_a%c_%dC", pre, rt, alp, u?mp->mu:mp->nu);
      if (mp->next)
         fprintf(fp, ",");
      else
         fprintf(fp, " ");
      fprintf(fp, "  /* index %d */\n", j);
   }
   fprintf(fp, "};\n\n");
}

void SpewForthCpDecl(char pre, int REVERSE, ATL_mmnode_t *mmb, FILE *fp,
                     char *arr, char *rt, char alp, int u)
{
   ATL_mmnode_t *mp;
   int j;

   fprintf(fp, "static const %s_t %s_a%c[%d] =\n{\n",
           REVERSE?"am2cm":"cm2am", arr, alp, ATL_CountNumberOfMMNodes(mmb));
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      const int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
      if (kmaj > 1)
         fprintf(fp, "   ATL_%c%s_a%c_%dx%d", pre, rt, alp, kmaj,
                 u?mp->mu:mp->nu);
      else
         fprintf(fp, "   ATL_%c%s_a%c_%d", pre, rt, alp, u?mp->mu:mp->nu);
      if (mp->next)
         fprintf(fp, ",");
      else
         fprintf(fp, " ");
      fprintf(fp, "  /* index %d */\n", j);
   }
   fprintf(fp, "};\n\n");
}


void GenAMAJ2CMAJ(char pre, ATL_mmnode_t *mmb, char *suff, char *outd)
/*
 * 3. atlas_<pre>amm_am2cm_a[1,X,n]:
 *    defines: ATL_AMM_NCASES
 *    prototypes all am2rm & am2cm routines
 *    1 indexible array giving which to use for each block factor
 */
{
   char ac[3] = {'1', 'n', 'X'};
   int ia, j;
   char *fnam, *sp, *np;
   ATL_mmnode_t *mp;

   if (!suff)
      suff = "";
   ia = strlen(outd) + strlen(suff) + 23;
   fnam = malloc(ia*sizeof(char));
   assert(fnam);
   sprintf(fnam, "%s/atlas_%c%samm_am2cm_a1.h", outd, pre, suff);
   np = fnam+ia-23+12;
   assert(*np == 'a' && np[1] == 'm');
   sp = fnam+ia-4;
   assert(*sp == '1');

   for (ia=0; ia < 3; ia++)
   {
      char *rt[2] = {"am2cm", "am2rm"};
      FILE *fp;
      int kmaj = FLAG_IS_SET(mmb->flag, MMF_KVEC) ? mmb->vlen:0;

      if (kmaj == 1)
         kmaj = 0;
      *sp = ac[ia];
      fp = fopen(fnam, "w");
      assert(fp);
      sp[1] = '\0';
      PrintBegBlock(pre, mmb, NULL, np, fp);
      sp[1] = '.';
      fprintf(fp, "/*\n * mat2blk prototypes\n */\n");
      SpewForthRevCpProto(pre, fp, ac[ia], mmb->mu, kmaj);
      if (mmb->nu != mmb->mu || kmaj > 1)
         SpewForthRevCpProto(pre, fp, ac[ia], mmb->nu, kmaj);
      for (mp=mmb->next; mp; mp = mp->next)
      {
         ATL_mmnode_t *p;
         int mu = mp->mu, nu = mp->nu;
         int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
         for (p=mmb; p != mp; p = p->next)
         {
            int km = FLAG_IS_SET(p->flag, MMF_KVEC) ? p->vlen:0;
            if (mu == p->mu && km == kmaj)
               break;
         }
         if (p == mp) /* haven't seen before */
            SpewForthRevCpProto(pre, fp, ac[ia], mu, kmaj);
         for (p=mmb; p != mp; p = p->next)
         {
            int km = FLAG_IS_SET(p->flag, MMF_KVEC) ? p->vlen:0;
            if ((p->nu == nu && km == kmaj) ||
                (km == kmaj && p->mu == nu))
               break;
         }
         if (p == mp) /* haven't seen before */
            SpewForthRevCpProto(pre, fp, ac[ia], nu, kmaj);
      }
      fprintf(fp, "\n");
      SpewForthCpDecl(pre,1,mmb, fp, "ATL_AMM_BLK2A", "am2cm", ac[ia], 1);
      SpewForthCpDecl(pre,1,mmb, fp, "ATL_AMM_BLK2AT", "am2rm", ac[ia], 1);
      SpewForthCpDecl(pre,1,mmb, fp, "ATL_AMM_BLK2B", "am2cm", ac[ia], 0);
      SpewForthCpDecl(pre,1,mmb, fp, "ATL_AMM_BLK2BT", "am2rm", ac[ia], 0);
      if (pre == 'c' || pre == 'z')
      {
         SpewForthCpConjDecl(pre, 1, mmb, fp, "ATL_AMM_BLKC2A",
                             "am2cm", ac[ia], 1);
         SpewForthCpConjDecl(pre, 1, mmb, fp, "ATL_AMM_BLKH2A",
                             "am2rm", ac[ia], 1);
         SpewForthCpConjDecl(pre, 1, mmb, fp, "ATL_AMM_BLKC2B",
                             "am2cm", ac[ia], 0);
         SpewForthCpConjDecl(pre, 1, mmb, fp, "ATL_AMM_BLKH2B",
                             "am2rm", ac[ia], 0);
      }
      fprintf(fp, "\n#endif  /* end include file guard */\n");
      fclose(fp);
   }
   free(fnam);
}
void GenCMAJ2AMAJ(char pre, ATL_mmnode_t *mmb, char *pfx, char *outd)
/*
 * 3. atlas_<pre>amm_cm2am_a[1,X,n]:
 *    defines: ATL_AMM_NCASES
 *    prototypes all rm2am & cm2am routines
 *    1 indexible array giving which to use for each block factor
 */
{
   char ac[3] = {'1', 'n', 'X'};
   int ia, j;
   char *fnam, *sp, *np;
   ATL_mmnode_t *mp;

   ia = strlen(outd) + strlen(pfx) + 23;
   fnam = malloc(ia*sizeof(char));
   assert(fnam);
   sprintf(fnam, "%s/atlas_%c%samm_cm2am_a1.h", outd, pre, pfx);
   np = fnam+ia-23+12;
   assert(*np == 'c' && np[1] == 'm');
   sp = fnam+ia-4;
   assert(*sp == '1');

   for (ia=0; ia < 3; ia++)
   {
      char *rt[2] = {"cm2am", "rm2am"};
      FILE *fp;
      int kmaj = FLAG_IS_SET(mmb->flag, MMF_KVEC) ? mmb->vlen:0;

      if (kmaj == 1)
         kmaj = 0;
      *sp = ac[ia];
      fp = fopen(fnam, "w");
      assert(fp);
      sp[1] = '\0';
      PrintBegBlock(pre, mmb, NULL, np, fp);
      sp[1] = '.';
      fprintf(fp, "/*\n * mat2blk prototypes\n */\n");
      SpewForthCpProto(pre, fp, ac[ia], mmb->mu, kmaj);
      if (mmb->nu != mmb->mu || kmaj > 1)
         SpewForthCpProto(pre, fp, ac[ia], mmb->nu, kmaj);
      for (mp=mmb->next; mp; mp = mp->next)
      {
         ATL_mmnode_t *p;
         int mu = mp->mu, nu = mp->nu;
         int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
         for (p=mmb; p != mp; p = p->next)
         {
            int km = FLAG_IS_SET(p->flag, MMF_KVEC) ? p->vlen:0;
            if (mu == p->mu && km == kmaj)
               break;
         }
         if (p == mp) /* haven't seen before */
            SpewForthCpProto(pre, fp, ac[ia], mu, kmaj);
         for (p=mmb; p != mp; p = p->next)
         {
            int km = FLAG_IS_SET(p->flag, MMF_KVEC) ? p->vlen:0;
            if ((p->nu == nu && km == kmaj) ||
                (km == kmaj && p->mu == nu))
               break;
         }
         if (p == mp) /* haven't seen before */
            SpewForthCpProto(pre, fp, ac[ia], nu, kmaj);
      }
      fprintf(fp, "\n");
      SpewForthCpDecl(pre,0,mmb, fp, "ATL_AMM_AN2BLK", "cm2am", ac[ia], 1);
      SpewForthCpDecl(pre,0,mmb, fp, "ATL_AMM_AT2BLK", "rm2am", ac[ia], 1);
      SpewForthCpDecl(pre,0,mmb, fp, "ATL_AMM_BN2BLK", "cm2am", ac[ia], 0);
      SpewForthCpDecl(pre,0,mmb, fp, "ATL_AMM_BT2BLK", "rm2am", ac[ia], 0);
      if (pre == 'c' || pre == 'z')
      {
         SpewForthCpConjDecl(pre, 0, mmb, fp, "ATL_AMM_AC2BLK",
                             "cm2am", ac[ia], 1);
         SpewForthCpConjDecl(pre, 0, mmb, fp, "ATL_AMM_AH2BLK",
                             "rm2am", ac[ia], 1);
         SpewForthCpConjDecl(pre, 0, mmb, fp, "ATL_AMM_BC2BLK",
                             "cm2am", ac[ia], 0);
         SpewForthCpConjDecl(pre, 0, mmb, fp, "ATL_AMM_BH2BLK",
                             "rm2am", ac[ia], 0);
      }
      fprintf(fp, "\n#endif  /* end include file guard */\n");
      fclose(fp);
   }
   free(fnam);
}

int KernelIsExactSame(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if kernels are the same including KB, 0 otherwise
 */
{
/*
 * Kernels aren't the same if one is being compiled with specific KB,
 * and the other has runtime
 */
   if (FLAG_IS_SET(p0->flag, MMF_KRUNTIME) !=
       FLAG_IS_SET(p1->flag, MMF_KRUNTIME))
      return(0);
   if (!FLAG_IS_SET(p0->flag, MMF_KRUNTIME) && (p0->kbB != p1->kbB))
      return(0);
/*
 * Kernels aren't same if their vectorization strategy is different
 */
   if (FLAG_IS_SET(p0->flag, MMF_KVEC) != FLAG_IS_SET(p1->flag, MMF_KVEC))
      return(0);
   if (p0->vlen != p1->vlen)
      return(0);
/*
 * Kernels aren't same if the -DATL_MOVE bits don't match
 */
   if (ATL_MMF_MVGET(p0->flag) != ATL_MMF_MVGET(p1->flag))
      return(0);
/*
 * Two generated kernels are the same if mu,nu,kmaj,ku,VLEN,flag are the same.
 * NOTE: if we make generator handle muladd, etc, MUST UPDATE HERE!!!
 */
   if (p0->ID == 0 && p1->ID == 0)
      return(p0->mu == p1->mu && p0->nu == p1->nu && p0->ku == p1->ku &&
             p0->vlen == p1->vlen && p0->flag == p1->flag);
/*
 * If both are user kernels, then they may be repeats.  For user kernels,
 * they are the same if both ID and flag match, else they are not.
 */
   else if (p0->ID > 0 && p1->ID > 0)
      return(p0->ID == p1->ID && p0->flag == p1->flag);
   return(0);  /* Can't be the same if above criteria fails */
}

/*
 * RETURNS: flags that necessitate recompilation, not including KRUNTIME,
 * which is encoded in kb
 */
int GetCompTimeFlags(ATL_mmnode_t *mp)
{
   int iflg;
   iflg = ATL_MMF_MVGET(mp->flag);  /* MVbits change kern at comp time */
   iflg |=  (((mp->flag) & 1)<<3);  /* LDTOP/BOT could be compile-time dec */
   if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      iflg |= 1<<4;
   return(iflg);
}

void FillInGenStrings
(
   char pre,
   ATL_mmnode_t *mmb,  /* queue to look through */
   char *dir           /* output directory to generate into */
)
/*
 * Creates GenString for any ID=0 in mmb
 */
{
   ATL_mmnode_t *mp;

   if (pre == 'z')
      pre = 'd';
   else if (pre == 'c')
      pre = 's';

   for (mp=mmb; mp; mp = mp->next)
   {
      if (mp->ID == 0)  /* is generated file */
      {
         int lr, ld;
         char *sp, *sp0=mp->rout;

         assert(sp0);  /* should have been filled in by search */
         lr = strlen(sp0);
         ld = strlen(dir);
         sp = malloc(ld+ld+2);
         assert(sp);
         strncpy(sp, dir, ld);
         sp[ld] = '/';
         strcpy(sp+ld+1, sp0);
         mp->rout = sp;
         mp->genstr = MMGetGenString(pre, mp);
         mp->rout = sp0;
         free(sp);
      }
   }
}

int ExactKernelInList(ATL_mmnode_t *mmb, ATL_mmnode_t *p)
/*
 * RETURNS: 1 if p is duplicated in mmb, else 0
 */
{
   ATL_mmnode_t *mp;
   if (!p || !mmb)
      return(0);
   for (mp=mmb; mp; mp = mp->next)
      if (KernelIsExactSame(mp, p))
         return(1);
    return(0);
}

void SpewForthKernProto(FILE *fp, char pre, char *nm, ATL_mmnode_t *p, char bc)
{
   fprintf(fp, "void ATL_%c%s_%d_%d_%x_%dx%dx%d_b%c\n", pre, nm, p->ID,
           FLAG_IS_SET(p->flag, MMF_KRUNTIME)?0:p->kbB, GetCompTimeFlags(p),
           p->mu, p->nu, p->ku,bc);
   fprintf(fp,
      "   (ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,TYPE*,\n");
   fprintf(fp,
      "    const TYPE*,const TYPE*,const TYPE*);\n");
}

void SpewForthKernProtos(FILE *fp, char pre, char *nm, ATL_mmnode_t *mmb,
                         int nbet)
{
   ATL_mmnode_t *mp;
   for (mp=mmb; mp; mp = mp->next)
   {
      char bc[3] = {'0', '1', 'n'};  /* 0 must come first */
      int ib;
      for (ib=0; ib < nbet; ib++)
         SpewForthKernProto(fp, pre, nm, mp, bc[ib]);
   }
}

void SpewForthKernArray(FILE *fp, char pre, ATL_mmnode_t *mmb,
                        char *vnam, char cbet)
{
   ATL_mmnode_t *mp;
   int n;
   char *nm = (vnam[5] == 'R' && vnam[6] == 'K') ? "AMRK":"AMMM";

   for (n=0,mp=mmb; mp; n++, mp = mp->next);
   fprintf(fp, "static const ammkern_t ATL_AMM_%s[%d] =\n", vnam, n);
   fprintf(fp, "{\n");
   for (mp=mmb; mp; mp = mp->next)
   {
      fprintf(fp, "   ATL_%c%s_%d_%d_%x_%dx%dx%d_b%c", pre, nm, mp->ID,
              FLAG_IS_SET(mp->flag, MMF_KRUNTIME)?0:mp->kbB,
              GetCompTimeFlags(mp), mp->mu, mp->nu, mp->ku, cbet);
      if (mp->next)
         fprintf(fp, ",\n");
   }
   fprintf(fp, "\n};\n\n");
}

/*
 * RETURNS: new list of K-cleanup kernels for mmb taken from mmk;
 * mmb & mmk not changed
 */
ATL_mmnode_t *GetUniqueK1Kerns
(
   char pre,            /* z,c,d,s */
   ATL_mmnode_t *mmb,   /* list of kernels requiring K-cleanup */
   ATL_mmnode_t *mmk,   /* list of K-cleanup kernels */
   char *outd
)
{
   ATL_mmnode_t *kcb=NULL, *prev, *mp;
   for (mp=mmb; mp; mp = mp->next)
   {
      ATL_mmnode_t *kc;
      int kmaj = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
      for (kc=mmk; kc; kc = kc->next) /* search for matching cleanup routine */
      {
         int km = FLAG_IS_SET(kc->flag, MMF_KVEC) ? kc->vlen:0;
         if (kc->mu == mp->mu && kc->nu == mp->nu && km == kmaj)
            break;
      }
      if (!kc)  /* cleanup routine must be in original list! */
      {
         if (mp->ku == 1 && FLAG_IS_SET(mp->flag, MMF_KRUNTIME) &&
             (!mp->kbmin || mp->kbmin == 1))
            kc = mp;
         else
         {
            for (kc=mmb; kc; kc = kc->next)
            {
               int km = FLAG_IS_SET(kc->flag, MMF_KVEC) ? kc->vlen:0;
               if (kc->mu == mp->mu && kc->nu == mp->nu && km == kmaj
                   && kc->ku == 1 && FLAG_IS_SET(kc->flag, MMF_KRUNTIME) &&
                   (!kc->kbmin || kc->kbmin == 1) &&
                   (!kc->kbmax || kc->kbmax >= mp->kbB))
                  break;
            }
            assert(kc);
         }
      }
      assert(kc);
      kc = CloneMMNode(kc);
      kc->mbB = mp->mbB;
      kc->nbB = mp->nbB;
      kc->kbB = mp->kbB;
      if (kcb)
      {
         prev->next = kc;
         prev = prev->next;
      }
      else
         prev = kcb = kc;

   }
   FillInGenStrings(pre, kcb, outd);
   return(kcb);
}

/*
 * RETURNS: possibly updated list of all unique mu/nu comboes
 */
typedef struct mnur mnur_t;
struct mnur {int mu; int nu; int kmaj; mnur_t *next;};
mnur_t *GetUniqueMNUnrolls(ATL_mmnode_t *mmb, mnur_t *urb)
{
   ATL_mmnode_t *mp;
   mnur_t *up;
/*
 * For each node in mmb, add to urb if mu/nu combo not already there
 * kmaj only affects A/B copy, and this is for C put, so ignore kmaj
 */
   for (mp = mmb; mp; mp = mp->next)
   {
      for (up=urb; up; up = up->next)
         if (mp->mu == up->mu && mp->nu == up->nu)
            break;
      if (!up)
      {
         up = malloc(sizeof(mnur_t));
         assert(up);
         up->mu = mp->mu;
         up->nu = mp->nu;
         up->kmaj = 0;
         up->next = urb;
         urb = up;
      }
   }
   return(urb);
}

/*
 * RETURNS: list of just unique MUs from mnb not already in mub
 */
mnur_t *GetUniqueMUnrolls(ATL_mmnode_t *mnb, mnur_t *mub)
{
   mnur_t *mup;
   ATL_mmnode_t *mnp;
   if (!mnb)
      return(mub);
   for (mnp=mnb; mnp; mnp = mnp->next)
   {
      const int kmaj = FLAG_IS_SET(mnp->flag, MMF_KVEC) ? mnp->vlen:0;
      for (mup=mub; mup; mup = mup->next)
         if (mup->mu == mnp->mu && mup->kmaj == kmaj)
            break;
      if (!mup)  /* a new mu */
      {
         mup = malloc(sizeof(mnur_t));
         assert(mup);
         mup->nu = mup->mu = mnp->mu;
         mup->next = mub;
         mup->kmaj = kmaj;
         mub = mup;
      }
   }
   return(mub);
}
/*
 * RETURNS: list of just unique NUs from mnb not already in mub
 */
mnur_t *GetUniqueNUnrolls(ATL_mmnode_t *mnb, mnur_t *mub)
{
   mnur_t *mup;
   ATL_mmnode_t *mnp;
   if (!mnb)
      return(mub);
   for (mnp=mnb; mnp; mnp = mnp->next)
   {
      const int kmaj = FLAG_IS_SET(mnp->flag, MMF_KVEC) ? mnp->vlen:0;
      for (mup=mub; mup; mup = mup->next)
         if (mup->mu == mnp->nu && kmaj == mup->kmaj)
            break;
      if (!mup)  /* a new nu */
      {
         mup = malloc(sizeof(mnur_t));
         assert(mup);
         mup->nu = mup->mu = mnp->nu;
         mup->next = mub;
         mup->kmaj = kmaj;
         mub = mup;
      }
   }
   return(mub);
}

void KillUnrollList(mnur_t *b)
{
   mnur_t *p;
   while (b)
   {
      p = b->next;
      free(b);
      b = p;
   }
}

int MMKernsSame1(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if kernels are the same except for blocking, 0 otherwise
 */
{
/*
 * Kernels aren't the same if one is being compiled with specific KB,
 * and the other has runtime
 */
   if (FLAG_IS_SET(p0->flag, MMF_KRUNTIME) !=
       FLAG_IS_SET(p1->flag, MMF_KRUNTIME))
      return(0);
/*
 * Two generated kernels are the same if mu,nu,ku,VLEN,flag are the same.
 * NOTE: if we make generator handle muladd, etc, MUST UPDATE HERE!!!
 */
   if (p0->ID == 0 && p1->ID == 0)
      return(p0->mu == p1->mu && p0->nu == p1->nu && p0->ku == p1->ku &&
             p0->vlen == p1->vlen && p0->flag == p1->flag);
/*
 * If both are user kernels, then they may be repeats.  For user kernels,
 * they are the same if both ID and flag match, else they are not.
 */
   else if (p0->ID > 0 && p1->ID > 0)
      return(p0->ID == p1->ID && p0->flag == p1->flag);
   return(0);  /* Can't be the same if above criteria fails */
}


void GenRankKH
(
   char pre,
   ATL_mmnode_t *geb,  /* baseptr  of gemm kernels */
   ATL_mmnode_t *rkb,  /* rank-K kernels, one for each supported K */
   char *outd
)
{
   FILE *fp;
   mnur_t *putb, *cpyb, *up;
   ATL_mmnode_t *mp;
   int m, n, k;
   int ia, ib;
   char PRE = pre;
   char ac[3] =  {'1', 'n', 'X'};
   char bc[4] = {'1', 'n', '0', 'X'};
   if (pre == 'c')
      pre = 's';
   else  if (pre == 'z')
      pre = 'd';

   assert(geb && rkb);
   fp = StandHStart(PRE, rkb, "rk", outd, "kern");
   fprintf(fp, "\n");

   for (mp=geb; mp->next; mp = mp->next);
   fprintf(fp, "#define ATL_rkAMM_LASTMB %d\n", mp->mbB);
   fprintf(fp, "#define ATL_rkAMM_LASTNB %d\n", mp->nbB);
   fprintf(fp, "#define ATL_rkAMM_LASTKB %d\n\n", mp->kbB);
/*
 * Prototype needed copy routines
 */
   putb = GetUniqueMNUnrolls(rkb, NULL);
   fprintf(fp, "/*\n * cblk2mat put function prototypes\n */\n");
   for (up=putb; up; up = up->next)
   {
      for (ib=0; ib < 4; ib++)
      {
         fprintf(fp, "void ATL_%cablk2cmat_%dx%d_a1_b%c\n",
                 PRE, up->mu, up->nu, bc[ib]);
         if (PRE == 'c' || PRE == 'z')
            fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const TYPE*,const SCALAR,TYPE *,ATL_CSZT);\n");
         else
            fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const SCALAR,TYPE *,ATL_CSZT);\n");
      }
   }
   KillUnrollList(putb);
   fprintf(fp,
      "/*\n * Column-major to access-major copy function prototypes\n */\n");
   cpyb = GetUniqueMUnrolls(rkb, NULL);
   cpyb = GetUniqueNUnrolls(rkb, cpyb);
   for (up=cpyb; up; up = up->next)
   {
      for (ia=0; ia < 3; ia++)
         SpewForthCpProto(PRE, fp, ac[ia], up->mu, up->kmaj);
   }
   KillUnrollList(cpyb);
/*
 * Prototype the rank-K functions
 */
   fprintf(fp, "/*\n * rank-K AMMM kernel prototypes\n */\n");
   for (mp=rkb; mp; mp = mp->next)
   {
/*
 *    For runtime-K kernels, only prototype 1st time they are seen in list
 */
      if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
      {
         ATL_mmnode_t *p;
         if (mp->ID > 0)
         {
            for (p=rkb; p != mp; p = p->next)
               if (p->ID == mp->ID)
                  break;
         }
         else
         {
            for (p=rkb; p != mp; p = p->next)
               if (MMKernsSame(mp, p))
                  break;
         }
         if (mp == p)
         {
            SpewForthKernProto(fp, pre, "AMRK", mp, '0');
            SpewForthKernProto(fp, pre, "AMRK", mp, '1');
            SpewForthKernProto(fp, pre, "AMRK", mp, 'n');
         }
      }
/*
 *    compile-time-K kernels get prototyped for each invocation
 */
      else
      {
         SpewForthKernProto(fp, pre, "AMRK", mp, '0');
         SpewForthKernProto(fp, pre, "AMRK", mp, '1');
         SpewForthKernProto(fp, pre, "AMRK", mp, 'n');
      }
   }
/*
 * Now, crank out funcptr arrays
 */
   for (ib=0; ib < 4; ib++)
   {
      fprintf(fp,
         "\nstatic const ablk2cmat_t ATL_AMM_BLK2C_a1_b%c[%d] =\n{\n",
              bc[ib], ATL_CountNumberOfMMNodes(rkb));
      SpewForthC2BDecl(PRE, rkb, fp, "ablk2cmat", '1', bc[ib]);
   }
   for (ia=0; ia < 3; ia++)
   {
      fprintf(fp, "\n");
      SpewForthCpDecl(PRE, 0, rkb, fp, "ATL_AMM_AN2BLK", "cm2am", ac[ia], 1);
      SpewForthCpDecl(PRE, 0, rkb, fp, "ATL_AMM_AT2BLK", "rm2am", ac[ia], 1);
      SpewForthCpDecl(PRE, 0, rkb, fp, "ATL_AMM_BN2BLK", "cm2am", ac[ia], 0);
      SpewForthCpDecl(PRE, 0, rkb, fp, "ATL_AMM_BT2BLK", "rm2am", ac[ia], 0);
      if (PRE == 'z' || PRE == 'c')
      {
         SpewForthCpConjDecl(PRE,0,rkb, fp, "ATL_AMM_AC2BLK", "cm2am",ac[ia],1);
         SpewForthCpConjDecl(PRE,0,rkb, fp, "ATL_AMM_AH2BLK", "rm2am",ac[ia],1);
         SpewForthCpConjDecl(PRE,0,rkb, fp, "ATL_AMM_BC2BLK", "cm2am",ac[ia],0);
         SpewForthCpConjDecl(PRE,0,rkb, fp, "ATL_AMM_BH2BLK", "rm2am",ac[ia],0);
      }
   }
   SpewForthKernArray(fp, pre, rkb, "KERN_RKK", '0');
   SpewForthKernArray(fp, pre, rkb, "KERN_RKK_b1", '1');
   SpewForthKernArray(fp, pre, rkb, "KERN_RKK_bn", 'n');
   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenKernH
(
   char pre,
   char *outd,
   char *pfx,
   ATL_mmnode_t *mmb,  /* full kerns */
   ATL_mmnode_t *k1b   /* K cleanup kerns */
)
/*
 * 5. atlas_<pre>_<pfx>kern.h
 *    defines: ATL_AMM_NCASES
 *    prototypes all kernels, including K-cleanup
 *    1 indexible array gives kernel to use for each case as func ptr
 *    1 indexible array gives K-clean kernel
 */
{
   FILE *fp;
   ATL_mmnode_t *ukb;
   const char upr = (pre == 'z' || pre == 'd') ? 'd' : 's';

   fp = StandHStart(pre, mmb, pfx, outd, "kern");
   fprintf(fp, "\n");
/*
 * For prototyping, create a list of all full and cleanup kerns, and get
 * rid of duplicates
 */
   ukb = AddUniqueMMKernCompList(NULL, mmb);
   ukb = AddUniqueMMKernCompList(ukb, k1b);
   SpewForthKernProtos(fp, upr, "AMMM", ukb, 3);
   fprintf(fp, "\n");
   KillAllMMNodes(ukb);
/*
 * Dump out kernel ptr arrays
 */
   SpewForthKernArray(fp, upr, mmb, "KERN_b0", '0');
   SpewForthKernArray(fp, upr, mmb, "KERN_b1", '1');
   SpewForthKernArray(fp, upr, mmb, "KERN_bn", 'n');
   if (k1b)
   {
      SpewForthKernArray(fp, upr, k1b, "KERN_K1", '0');
      SpewForthKernArray(fp, upr, k1b, "KERN_K1_b1", '1');
      SpewForthKernArray(fp, upr, k1b, "KERN_K1_bn", 'n');
   }

   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenCmplxHeaders(char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *rkb,
                     char *pfx, char *outd)
{
   GenCMAJ2AMAJ(pre, mmb, pfx, outd);
   GenC2BLK(pre, mmb, pfx, outd);
   if (rkb)
   {
      GenCMAJ2AMAJ(pre, rkb, "rk", outd);
      GenC2BLK(pre, rkb, "rk", outd);
      GenRankKH(pre, mmb, rkb, outd);
   }
}

void GenHeaderFiles(char pre, char *outd,
                    ATL_mmnode_t *geb, ATL_mmnode_t *gek1b,
                    ATL_mmnode_t *sqb, ATL_mmnode_t *sqk1b,
                    ATL_mmnode_t *rkb, ATL_mmnode_t *trsmb)
/*
 * Header files required to build full gemm (no timing):
 *X1. atlas_<pre>amm_blk.h :
 *X   defines: ATL_AMM_NCASES, ATL_AMM_MAX[M,N,K]B
 *X   3 arrays indexed by case give blocking
 *X2  atlas_<pre>amm_flag.h
 *X   defines: ATL_AMM_NCASES
 *X   1 indexible array giving KRUNTIME for now
 *X3. atlas_<pre>amm_cm2am_a[1,X,n]:
 *X   defines: ATL_AMM_NCASES
 *X   prototypes all rm2am & cm2am routines
 *X   1 indexible array giving which to use for each block factor
 *X4. atlas_<pre>amm_ablk2cmat.h
 *X   defines: ATL_AMM_NCASES
 *X   prototypes all ablk2cmat routines
 *X   1 indexible array for each alpha,beta combination
 *X   -> 3*4 = 12 indexible arrays total
 *X5. atlas_<pre>amm_kerns.h
 *X   defines: ATL_AMM_NCASES
 *X   prototypes all kernels, including K-cleanup
 *X   1 indexible array gives kernel to use for each case as func ptr
 *X   1 indexible array gives K-clean kernel
 *X6. atlas_<pre>amm_cmat2ablk.h (I don't need, Rakib does)
 *X   defines: ATL_AMM_NCASES
 *X   prototypes all cmat2ablk routines
 *X   1 indexible array for each alpha,beta combination
 *X   -> 3*4 = 12 indexible arrays total
 */
{
   ATL_mmnode_t *mp;
   char *pfx;

   pfx = "ge";
   GenAmmSum(pre, geb, rkb, pfx, outd);
   GenBlockingFile(pre, geb, pfx, outd);
   GenPerfFile(pre, geb, pfx, outd, 0.0);
   GenFlagH(pre, geb, pfx, outd);
   GenCMAJ2AMAJ(pre, geb, pfx, outd);
   GenC2BLK(pre, geb, pfx, outd);
   GenKernH(pre, outd, pfx, geb, gek1b);

   pfx = "sq";
   GenAmmSum(pre, sqb, rkb, pfx, outd);
   GenBlockingFile(pre, sqb, pfx, outd);
   GenPerfFile(pre, sqb, pfx, outd, 0.0);
   GenFlagH(pre, sqb, pfx, outd);
   GenAMAJ2CMAJ(pre, sqb, pfx, outd);
   GenCMAJ2AMAJ(pre, sqb, pfx, outd);
   GenC2BLK(pre, sqb, pfx, outd);
   GenKernH(pre, outd, pfx, sqb, sqk1b);

   if (trsmb)
   {
      double maxMF=0.0;
      for (mp=sqb; mp; mp = mp->next)
         if (mp->mflop[0] > maxMF)
             maxMF = mp->mflop[0];
      GenPerfFile(pre, trsmb, "ts", outd, maxMF);
   }
   if (rkb)
   {
      pfx = "rk";
      GenBlockingFile(pre, rkb, pfx, outd);
      GenPerfFile(pre, rkb, pfx, outd, 0.0);
      GenFlagH(pre, rkb, pfx, outd);
      GenCMAJ2AMAJ(pre, rkb, pfx, outd);
      GenC2BLK(pre, rkb, pfx, outd);
      GenRankKH(pre, geb, rkb, outd);
   }
}

/*
 * Splits rkb into two lists: (1) Routines with runtime K (RUNB),
 * (2) Routines with fixed-KB (returned)
 * NOTE: original list rkb is unchanged
 */
ATL_mmnode_t *SplitRankK(ATL_mmnode_t *rkb, ATL_mmnode_t **RUNB)
{
   ATL_mmnode_t *runb=NULL, *fixb=NULL, *p, *np;
   for (p=rkb; p; p = p->next)
   {
      np = CloneMMNode(p);
      if (FLAG_IS_SET(np->flag, MMF_KRUNTIME))
      {
         np->next = runb;
         runb = np;
      }
      else
      {
         np->next = fixb;
         fixb = np;
      }
   }
   *RUNB = runb;
   return(fixb);
}

char *GetKernComp(ATL_mmnode_t *mmp, char *dcomp, char *dflags, char **flgs)
{
   char *comp = dcomp;
   if (mmp->comp)
   {
      comp = (mmp->comp[0] == 'g' && mmp->comp[1] == 'c' &&
              mmp->comp[2] == 'c' &&
             (mmp->comp[3] == '\0' || mmp->comp[3] == ' '))
             ? "$(GOODGCC)" : mmp->comp;
      *flgs = mmp->cflags;
   }
   else
      *flgs = dflags;
   return(comp);
}
void PrintKernComp
(
   FILE *fp,            /* file to print to */
   char pre,
   ATL_mmnode_t *mmp,   /* kernel compile rule is for */
   int UID,             /* user-ID for user-determined kerns */
   char *nam,           /* base name ("AMMM" or "AMRK") */
   char *comp,          /* compiler to use */
   char *cflags,        /* compiler flags to use */
   char *styp,          /* string defining type (eg. "-DSREAL") */
   char cbet,           /* character with beta name ('1', '0', 'n') */
   char *sbet           /* string wt full beta name ("1", "N1", "0") */
)
{
   const int kb = FLAG_IS_SET(mmp->flag, MMF_KRUNTIME)?0:mmp->kbB;
   const int flg = GetCompTimeFlags(mmp);
   fprintf(fp, "ATL_%c%s_%d_%d_%x_%dx%dx%d_b%c.o : %s\n", pre, nam,
           mmp->ID, kb, flg, mmp->mu, mmp->nu, mmp->ku, cbet, mmp->rout);
   fprintf(fp, "\t%s $(CDEFS2) %s -DBETA%s=1", comp, styp, sbet);
   if (!FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
      fprintf(fp, " -DMB=%d -DNB=%d, -DKB=%d", mmp->mbB, mmp->nbB, mmp->kbB);
   if (FLAG_IS_SET(mmp->flag, MMF_MVA))
      fprintf(fp, " -DATL_MOVEA");
   if (FLAG_IS_SET(mmp->flag, MMF_MVB))
      fprintf(fp, " -DATL_MOVEB");
   if (FLAG_IS_SET(mmp->flag, MMF_MVC))
      fprintf(fp, " -DATL_MOVEC");
         fprintf(fp, " -DATL_USERMM=ATL_%c%s_%d_%d_%x_%dx%dx%d_b%c", pre, nam,
                 mmp->ID, kb, flg, mmp->mu, mmp->nu, mmp->ku, cbet);
         fprintf(fp, " %s -o ATL_%c%s_%d_%d_%x_%dx%dx%d_b%c.o -c %s\n",
                 cflags, pre, nam, mmp->ID, kb, flg, mmp->mu, mmp->nu, mmp->ku,
                 cbet, mmp->rout);
}
void GenMakefile
(
   char pre,            /* type/precision prefix : s,d,c,z */
   char *outd,          /* directory to generate Makefile in */
   ATL_mmnode_t *geb,   /* main kernels for GEMM */
   ATL_mmnode_t *gek1b, /* K-cleanup for gemm kerns */
   ATL_mmnode_t *sqb,   /* kernels where MB=NB=KB */
   ATL_mmnode_t *sqk1b, /* K-cleanup for square kerns */
   ATL_mmnode_t *rkb    /* list of kernels to doing rank-K update */
)
{
   ATL_mmnode_t *mmb, *mmp, *p, *fixb, *runb;
   mnur_t *mnurb=NULL, *allub, *up;
   FILE *fp;
   char *comp, *cflags;
   char *ln;
   int i;
   char pres[2];
   char be[3] = {'1', 'n', '0'};
   char *bes[3] = {"1", "N1", "0"};
   char al[3] = {'1', 'n', 'X'};
   char dcomp[8] = {'$', '(', 'D', 'M', 'C', ')', '\0'};
   char dflags[12] = {'$', '(', 'D', 'M', 'C', 'F', 'L', 'A', 'G', 'S',
                     ')', '\0'};
   char *styps[2] = {"-DDREAL", "-DDCPLX"};
   char *styp = (pre == 'd' || pre == 'z') ? "-DDREAL" : "-DSREAL";
/*
 * Square & rect kernels handled same way: going to be compiled once for
 * each compile-time K, and once if runtime, so combine them into one big
 * queue for Make generation, w/o repeating same compile-time kernel
 */
   mmb = AddUniqueMMKernCompList(NULL, geb);
   mmb = AddUniqueMMKernCompList(mmb, gek1b);
   mmb = AddUniqueMMKernCompList(mmb, sqk1b);
   mmb = AddUniqueMMKernCompList(mmb, sqb);

   pres[0] = pre;
   if (pre == 's' || pre == 'c')
   {
      styps[0] = "-DSREAL";
      styps[1] = "-DSCPLX";
      pres[1] = 'c';
   }
   else
      pres[1] = 'z';

   ln = malloc((strlen(outd)+11)*sizeof(char));
   assert(ln);
   sprintf(ln, "%s/%cMake_amm", outd, pre);
   fp = fopen(ln, "w");
   assert(fp);
   free(ln);
   fprintf(fp, "include ../Make.inc\n");
   fprintf(fp, "CDEFS2=$(CDEFS)\n\n");
   if (pre == 'c')
   {
      fprintf(fp, "CMC=$(SMC)\n");
      fprintf(fp, "CKCFLAGS=$(SKCFLAGS)\n");
      fprintf(fp, "CMCFLAGS=$(SMCFLAGS)\n");
   }
   else if (pre == 'z')
   {
      fprintf(fp, "ZMC=$(DMC)\n");
      fprintf(fp, "ZKCFLAGS=$(DKCFLAGS)\n");
      fprintf(fp, "ZMCFLAGS=$(DMCFLAGS)\n");
   }
/*
 * Build list of all unique MU/NU combos for copy routines
 */
   mnurb = GetUniqueMNUnrolls(mmb, NULL);
   mnurb = GetUniqueMNUnrolls(rkb, mnurb);
   allub = GetUniqueMUnrolls(mmb, NULL);
   allub = GetUniqueMUnrolls(rkb, allub);
   allub = GetUniqueNUnrolls(mmb, allub);
   allub = GetUniqueNUnrolls(rkb, allub);
/*
 * Spew out all filenames that must be compiled
 */
   fprintf(fp, "objs =");
/*
 * Routines to copy from MU/NU-major to column major output array
 */
   for (up=mnurb; up; up = up->next)
   {
      int j;
      const int mu=up->mu, nu=up->nu;

      for (j=0; j < 3; j++)
      {
         int k;
         char *rtn[2] = {"ablk2cmat", "cmat2ablk"};
         for (k=0; k < 2; k++)
         {
            int h;
            for (h=0; h < 2; h++)
            {
               char pre = pres[h];
               fprintf(fp, " \\\n       ATL_%c%s_%dx%d_a%c_b1.o",
                       pre, rtn[k], mu, nu, al[j]);
               fprintf(fp, " ATL_%c%s_%dx%d_a%c_bX.o",
                       pre, rtn[k], mu, nu, al[j]);
               fprintf(fp, " \\\n       ATL_%c%s_%dx%d_a%c_b0.o",
                       pre, rtn[k], mu, nu, al[j]);
               fprintf(fp, " ATL_%c%s_%dx%d_a%c_bn.o",
                       pre, rtn[k], mu, nu, al[j]);
            }
         }
      }
   }
/*
 * Routines to copy back and forth from A and B
 */
   for (up=allub; up; up = up->next)
   {
      int h;
      for (h=0; h < 2; h++)
      {
         char pre=pres[h];
         int j;
         const int u = up->mu;
         for (j=0; j < 3; j++)
         {
            if (up->kmaj > 1)
            {
               fprintf(fp,
 " \\\n       ATL_%crm2am_a%c_%dx%d.o ATL_%ccm2am_a%c_%dx%d.o",
                       pre, al[j], up->kmaj, u, pre, al[j], up->kmaj, u);
               fprintf(fp,
 " \\\n       ATL_%cam2rm_a%c_%dx%d.o ATL_%cam2cm_a%c_%dx%d.o",
                       pre, al[j], up->kmaj, u, pre, al[j], up->kmaj, u);
            }
            else
            {
               fprintf(fp,
                       " \\\n       ATL_%crm2am_a%c_%d.o ATL_%ccm2am_a%c_%d.o",
                       pre, al[j], u, pre, al[j], u);
               fprintf(fp,
                       " \\\n       ATL_%cam2rm_a%c_%d.o ATL_%cam2cm_a%c_%d.o",
                       pre, al[j], u, pre, al[j], u);
            }
            if (pre == 'c' || pre == 'z')
            {
               if (up->kmaj > 1)
               {
                  fprintf(fp,
 " \\\n       ATL_%crm2am_a%c_%dx%dC.o ATL_%ccm2am_a%c_%dx%dC.o",
                          pre, al[j], up->kmaj, u, pre, al[j], up->kmaj, u);
                  fprintf(fp,
 " \\\n       ATL_%cam2rm_a%c_%dx%dC.o ATL_%cam2cm_a%c_%dx%dC.o",
                          pre, al[j], up->kmaj, u, pre, al[j], up->kmaj, u);
               }
               else
               {
                  fprintf(fp,
                     " \\\n       ATL_%crm2am_a%c_%dC.o ATL_%ccm2am_a%c_%dC.o",
                          pre, al[j], u, pre, al[j], u);
                  fprintf(fp,
                     " \\\n       ATL_%cam2rm_a%c_%dC.o ATL_%cam2cm_a%c_%dC.o",
                          pre, al[j], u, pre, al[j], u);
               }
            }
         }
      }
   }
/*
 * AMM kernel routines
 */
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      int kb = mmp->kbB;
/*
 *    Kernels that take runtime K are only compiled once, so don't repeat them
 *    for every KB.  Only generate a statement if this is the first one.
 */
      if (FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
      {
         const int id = mmp->ID;
         for (p=mmb; p != mmp; p = p->next)
            if (MMKernsSame(mmp, p))
               break;
         if (p != mmp)
            continue;
         kb = 0;
      }
/*
 *    ATL_<pre>AMMM_<ID>_<kb>_<flg>_<mu>x<nu>x<ku>_b<X>
 */
      for (i=0; i < 3; i++)
         fprintf(fp, " \\\n       ATL_%cAMMM_%d_%d_%x_%dx%dx%d_b%c.o",
                 pre, mmp->ID, kb, GetCompTimeFlags(mmp), mmp->mu, mmp->nu,
                 mmp->ku, be[i]);
   }
/*
 * rank-K AMM routines where KB is runtime variable
 */
   fixb = SplitRankK(rkb, &runb);
   runb = RemoveNonUniqueKernels(runb);
   for (mmp=runb; mmp; mmp = mmp->next)
      for (i=0; i < 3; i++)
         fprintf(fp, " \\\n       ATL_%cAMRK_%d_0_%x_%dx%dx%d_b%c.o", pre,
                 mmp->ID, GetCompTimeFlags(mmp),mmp->mu,mmp->nu,mmp->ku,be[i]);
/*
 * rank-K AMM routines where KB is compile-time, need obj for each unique kbB
 */
   for (mmp=fixb; mmp; mmp = mmp->next)
   {
      for (i=0; i < 3; i++)
         fprintf(fp, " \\\n       ATL_%cAMRK_%d_%d_%x_%dx%dx%d_b%c.o",
                 pre, mmp->ID, mmp->kbB, GetCompTimeFlags(mmp),
                 mmp->mu, mmp->nu, mmp->ku, be[i]);
   }
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
   fprintf(fp, "\t- rm -f ATL_%c*.[S,c]\n", pre);

/*
 * Print out the individual rules for each needed copy function
 */
   dcomp[2] = dflags[2] = toupper(pre);
   dflags[3] = dcomp[3] = 'K';
   fprintf(fp, "\ntsth.o : tsth.c\n");
   fprintf(fp, "\t%s %s $(CDEFS) %s -c tsth.c\n\n", dcomp, dflags, styp);
   fprintf(fp, "#\n# Data copy rules\n#\n");
/*
 * Print out 2-D ablk2Cmat, cmat2ablk and ammswp targets
 */
   for (up=mnurb; up; up = up->next)
   {
      const int mu=up->mu, nu = up->nu;
      char cbe[4] = {'0', '1', 'n', 'X'};
      int ibe[4] =  {0,    1,  -1,  2};
      int i, j;

      for (i=0; i < 4; i++)
      {
         int h;
         for (h=0; h < 2; h++)
         {
            char pre=pres[h];
            char *styp=styps[h];
            for (j=0; j < 3; j++)
            {
               int k;
               char *rtn[2] = {"ablk2cmat", "cmat2ablk"};
               for (k=0; k < 2; k++)
               {
                  char rn[64];
                  sprintf(rn, "ATL_%c%s_%dx%d_a%c_b%c",
                          pre, rtn[k], mu, nu, al[j], cbe[i]);
                  fprintf(fp, "%s.o : %s.c\n", rn, rn);
                  fprintf(fp, "\t%s %s $(CDEFS) %s -c %s.c\n",
                          dcomp, dflags, styp, rn);
               }
            }
         }
      }
   }
   KillUnrollList(mnurb);
/*
 * Print out 1-D copy-in routine rules
 */
   for (up=allub; up; up = up->next)
   {
      const int u = up->mu, kmaj=up->kmaj;
      int j;
      for (j=0; j < 3; j++)
      {
         int h;
         for (h=0; h < 2; h++)
         {
            char pre=pres[h];
            char *styp=styps[h];
            char *cst[2] = {"", "C"};
            char *cdef[2] = {"", "-DConj_=1 "};
            int g;
            const int G = (pre == 'c' || pre == 'z') ? 2:1;
            for (g=0; g < G; g++)
            {
               int kk, r;
               const int nmI = 5;
               char rout[32];
               char *routs[4] = {"rm2am", "cm2am", "am2rm", "am2cm"};
               for (r=0; r < 4; r++)
               {
                  if (kmaj > 1)
                     sprintf(rout, "ATL_%c%s_a%c_%dx%d",
                             pre, routs[r], al[j], kmaj, u);
                  else
                     sprintf(rout, "ATL_%c%s_a%c_%d", pre, routs[r], al[j], u);
                  fprintf(fp, "%s%s.o : %s.c\n", rout, cst[g], rout);
                  fprintf(fp, "\t%s %s $(CDEFS) %s%s -o %s%s.o -c %s.c\n",
                          dcomp, dflags, cdef[g], styp, rout, cst[g], rout);
               }
            }
         }
      }
   }
   KillUnrollList(allub);
/*
 * Print out the individual rules for each kernel compile
 */
   dflags[3] = dcomp[3] = 'M';
   fprintf(fp, "#\n#  AMM kernel rules\n#\n");
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
/*
 *    Kernels that take runtime K are only compiled once, so print rules on
 *    only the first encounter of that ID
 */
      if (FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
      {
         const int id = mmp->ID;
         for (p=mmb; p != mmp; p = p->next)
            if (MMKernsSame(mmp, p))
               break;
         if (p != mmp)
            continue;
      }
      comp = GetKernComp(mmp, dcomp, dflags, &cflags);
/*
 *    ATL_<pre>AMMM_<ID>_<kb>_<flg>_<mu>x<nu>x<ku>_b<X>,
 *    kerns in all 3 beta cases
 */
      for (i=0; i < 3; i++)
         PrintKernComp(fp, pre, mmp, -1, "AMMM", comp, cflags, styp,
                       be[i], bes[i]);
   }
/*
 * runtime-KB rank-K needs only one kernel compile rule
 */
   if (runb)
      fprintf(fp, "#\n#  rank-K kernels with run-time K\n#\n");
   for (mmp=runb; mmp; mmp = mmp->next)
   {
      comp = GetKernComp(mmp, dcomp, dflags, &cflags);
      for (i=0; i < 3; i++)
         PrintKernComp(fp, pre, mmp, -1, "AMRK", comp, cflags, styp,
                       be[i], bes[i]);
   }

/*
 * Compile-time rank-K must be recompiled for each kbB
 */
   if (fixb)
      fprintf(fp, "#\n#  rank-K kernels with compile-time KB\n#\n");
   for (mmp=fixb; mmp; mmp = mmp->next)
   {
      const int kmaj = FLAG_IS_SET(mmp->flag, MMF_KVEC) ? mmp->vlen:0;
      comp = GetKernComp(mmp, dcomp, dflags, &cflags);
      for (i=0; i < 3; i++)
      {
         int kbG = mmp->kbB;
         if (kmaj)
            kbG = ((mmp->kbB + kmaj-1)/kmaj)*kmaj;
         fprintf(fp, "ATL_%cAMRK_%d_%d_%x_%dx%dx%d_b%c.o : %s\n", pre, mmp->ID,
                 mmp->kbB, GetCompTimeFlags(mmp), mmp->mu, mmp->nu, mmp->ku,
                 be[i], mmp->rout);
         fprintf(fp, "\t%s $(CDEFS2) %s -DBETA%s=1", comp, styp, bes[i]);
         fprintf(fp, " -DKB=%d -DATL_USERMM=ATL_%cAMRK_%d_%d_%x_%dx%dx%d_b%c",
                 kbG, pre, mmp->ID, mmp->kbB, GetCompTimeFlags(mmp), mmp->mu,
                 mmp->nu, mmp->ku, be[i]);
         fprintf(fp, " -DATL_MOVEA -DATL_MOVEC");
         fprintf(fp, " %s -o ATL_%cAMRK_%d_%d_%x_%dx%dx%d_b%c.o -c %s\n",
                 cflags, pre, mmp->ID, mmp->kbB, GetCompTimeFlags(mmp),
                 mmp->mu, mmp->nu, mmp->ku, be[i], mmp->rout);
      }
   }
   KillAllMMNodes(fixb);
   KillAllMMNodes(runb);
   KillAllMMNodes(mmb);
   fclose(fp);
}

void GenKerns(char pre, char *outd, ATL_mmnode_t *geb, ATL_mmnode_t *gek1b,
              ATL_mmnode_t *sqb, ATL_mmnode_t *sqk1b, ATL_mmnode_t *rkb)
/*
 * Creates/copies all required matmul kernels into outd using specified names
 */
{
   ATL_mmnode_t *mmb, *mmp, *p;
   mnur_t *mnurb=NULL, *allub, *up;
   char *ln=NULL;
   int lnlen=0, dlen;
   char al[3] = {'1', 'n', 'X'};
   int ial[3] = {1,   -1,   2};
   char pres[2];
   pres[0] = pre;
   pres[1] = (pre == 's') ? 'c' : 'z';
/*
 * All non-rkK kerns are generated the same way, so combine them all into one
 * gigantic queue.  We then eliminate all the guys who aren't different
 * once we have compiled them, which still leaves duplicates that differ
 * only in their compilation options.  We'll go ahead and just generate
 * them multiple times (overwriting the file with the same file), so that
 * we don't have to write & support a stricter RemoveNonUnique
 */
   mmb = CloneMMQueue(geb);
   ATL_LastMMNode(mmb)->next = CloneMMQueue(sqb);
   ATL_LastMMNode(mmb)->next = CloneMMQueue(gek1b);
   ATL_LastMMNode(mmb)->next = CloneMMQueue(sqk1b);
   ATL_LastMMNode(mmb)->next = CloneMMQueue(rkb);
   mmb = RemoveNonUniqueKernels(mmb);
/*
 * Build list of all unique MU/NU combos for copy routines
 */
   mnurb = GetUniqueMNUnrolls(mmb, NULL);
   allub = GetUniqueMUnrolls(mmb, NULL);
   allub = GetUniqueNUnrolls(mmb, allub);
   dlen = strlen(outd);
/*
 * Extract every unique block-copy routine
 */
   for (up=mnurb; up; up = up->next)
   {
      const int mu=up->mu, nu=up->nu;
      char cbe[4] = {'0', '1', 'n', 'X'};
      int ibe[4] =  {0,    1,  -1,  2};
      int i, j;
      j = 64+8 + strlen(outd);
      j = (j > 128) ? j : 128;
      if (lnlen < j)
      {
         free(ln);
         lnlen = j;
         ln = malloc(j*sizeof(char));
         assert(ln);
      }
      for (i=0; i < 4; i++)
      {
         for (j=0; j < 3; j++)
         {
            char rn[64];
            int ierr;
            int k;
            for (k=0; k < 2; k++)
            {
               int h;
               for (h=0; h < 2; h++)
               {
                  char pre = pres[h];
                  if (!k)
                     sprintf(rn, "ATL_%cablk2cmat_%dx%d_a%c_b%c.c",
                             pre, mu, nu, al[j], cbe[i]);
                  else
                     sprintf(rn, "ATL_%ccmat2ablk_%dx%d_a%c_b%c.c",
                             pre, mu, nu, al[j], cbe[i]);
                  sprintf(ln,
                     "make %s pre=%c mu=%d nu=%d al=%c be=%c alpha=%d beta=%d",
                          rn, pre, mu, nu, al[j], cbe[i], ial[j], ibe[i]);
                  ierr = system(ln);
                  if (ierr)
                  {
                     fprintf(stderr, "FAILED CMND='%s'\n", ln);
                     exit(ierr);
                  }
                  sprintf(ln, "mv %s %s/.", rn, outd);
                  ierr = system(ln);
                  if (ierr)
                  {
                     fprintf(stderr, "FAILED CMND='%s'\n", ln);
                     exit(ierr);
                  }
               }
            }
         }
      }
   }
   KillUnrollList(mnurb);

/*
 * Routines to copy back and forth from A and B
 */
   for (up=allub; up; up = up->next)
   {
      const int u = up->mu, kmaj=up->kmaj;
      int j;
      j = 16 * 40 + (strlen(outd)<<1);
      j = (j > 90) ? j : 90;
      if (lnlen < j)
      {
         free(ln);
         lnlen = j;
         ln = malloc(j*sizeof(char));
         assert(ln);
      }
      for (j=0; j < 3; j++)
      {
         int h;
         for (h=0; h < 2; h++)
         {
            int ierr;
            char pre = pres[h];
            sprintf(ln, "make ATL_%crm2am_a%c_%d.c ATL_%ccm2am_a%c_%d.c "
                    "pre=%c UR=%d alpha=%d al=%c kmaj=%d",
                    pre, al[j], u, pre, al[j], u, pre, u, ial[j], al[j], kmaj);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            sprintf(ln, "make ATL_%cam2rm_a%c_%d.c ATL_%cam2cm_a%c_%d.c "
                    "pre=%c UR=%d alpha=%d al=%c kmaj=%d",
                    pre, al[j], u, pre, al[j], u, pre, u, ial[j], al[j], kmaj);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            if (kmaj > 1)
               sprintf(ln,
                       "mv ATL_%cam2rm_a%c_%d.c %s/ATL_%cam2rm_a%c_%dx%d.c",
                       pre, al[j], u, outd, pre, al[j], kmaj, u);
            else
               sprintf(ln,
                       "mv ATL_%cam2rm_a%c_%d.c %s/ATL_%cam2rm_a%c_%d.c",
                       pre, al[j], u, outd, pre, al[j], u);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            if (kmaj > 1)
               sprintf(ln,
                       "mv ATL_%cam2cm_a%c_%d.c %s/ATL_%cam2cm_a%c_%dx%d.c",
                       pre, al[j], u, outd, pre, al[j], kmaj, u);
            else
               sprintf(ln,
                       "mv ATL_%cam2cm_a%c_%d.c %s/ATL_%cam2cm_a%c_%d.c",
                       pre, al[j], u, outd, pre, al[j], u);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            if (kmaj > 1)
               sprintf(ln,
                       "mv ATL_%crm2am_a%c_%d.c %s/ATL_%crm2am_a%c_%dx%d.c",
                       pre, al[j], u, outd, pre, al[j], kmaj, u);
            else
               sprintf(ln,
                       "mv ATL_%crm2am_a%c_%d.c %s/ATL_%crm2am_a%c_%d.c",
                       pre, al[j], u, outd, pre, al[j], u);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            if (kmaj > 1)
               sprintf(ln,
                       "mv ATL_%ccm2am_a%c_%d.c %s/ATL_%ccm2am_a%c_%dx%d.c",
                       pre, al[j], u, outd, pre, al[j], kmaj, u);
            else
               sprintf(ln,
                       "mv ATL_%ccm2am_a%c_%d.c %s/ATL_%ccm2am_a%c_%d.c",
                       pre, al[j], u, outd, pre, al[j], u);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
         }
      }
   }
   KillUnrollList(allub);
/*
 * Copy/generate every unique file, but only once
 */
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      const int id=mmp->ID;
/*
 *    For generated files, just overwrite dups with same file, won't hurt
 */
      if (id == 0)
      {
         assert(mmp->genstr);
         assert(!system(mmp->genstr));
      }
      else  /* user-supplied files copied from AMMCASES directory */
      {
/*
 *       If this is the first time we've seen this ID, it must be copied
 */
         for (p=mmb; p != mmp && p->ID != id; p = p->next);
         if (p == mmp)
         {
            int i, ierr;
            i = strlen(mmp->rout) + dlen + 16;
            if (i > lnlen)
            {
               if (ln)
                  free(ln);
               ln = malloc(i*sizeof(char));
               assert(ln);
               lnlen = i;
            }
            sprintf(ln, "cp AMMCASES/%s %s/.", mmp->rout, outd);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
         }
      }
   }
/*
 * Copy/generate rank-K kerns
 */
   for (mmp=rkb; mmp; mmp = mmp->next)
   {
      const int id=mmp->ID;
/*
 *    For generated files, just generate it using genstr
 */
      if (id == 0)
      {
         assert(mmp->genstr);
         assert(!system(mmp->genstr));
      }
      else  /* user-supplied files copied from AMMCASES directory */
      {
         int i, ierr;
         i = strlen(mmp->rout) + dlen + 16;
         if (i > lnlen)
         {
            if (ln)
               free(ln);
            ln = malloc(i*sizeof(char));
            assert(ln);
            lnlen = i;
         }
         sprintf(ln, "cp AMMCASES/%s %s/.", mmp->rout, outd);
         ierr = system(ln);
         if (ierr)
         {
            fprintf(stderr, "FAILED CMND='%s'\n", ln);
            exit(ierr);
         }
      }
   }
   KillAllMMNodes(mmb);
   if (ln)
      free(ln);
}

int KernelInList(ATL_mmnode_t *mmb, ATL_mmnode_t *p)
/*
 * RETURNS: 1 if p is duplicated in mmb, else 0
 */
{
   ATL_mmnode_t *mp;
   if (!p || !mmb)
      return(0);
   for (mp=mmb; mp; mp = mp->next)
      if (MMKernsSame1(mp, p))
         return(1);
    return(0);
}

ATL_mmnode_t *StripNonUniqueKs(ATL_mmnode_t *ukb, ATL_mmnode_t *mmb)
/*
 * Deletes any ukb node that also appears in mmb,
 * RETURNS: possibly shortened ukb
 */
{
   ATL_mmnode_t *mp, *prev;
/*
 * Delete any repetitive nodes starting unique queue
 */
   while(ukb && KernelInList(mmb, ukb))
      ukb = KillMMNode(ukb);
   if (!ukb)
      return(ukb);
/*
 * Now, delete any internal non-unique K-cleanup
 */
   prev = ukb;
   mp = ukb->next;
   while (mp)
   {
      if (KernelInList(mmb, mp))
         mp = prev->next = KillMMNode(mp);
      else
      {
         prev = mp;
         mp = mp->next;
      }
   }
   return(ukb);
}
ATL_mmnode_t *StripExactMatchKs(ATL_mmnode_t *ukb, ATL_mmnode_t *mmb)
/*
 * Deletes any ukb node that also appears in mmb with same K-value,
 * RETURNS: possibly shortened ukb
 */
{
   ATL_mmnode_t *mp, *prev;
/*
 * Delete any repetitive nodes starting unique queue
 */
   while(ukb && ExactKernelInList(mmb, ukb))
      ukb = KillMMNode(ukb);
   if (!ukb)
      return(ukb);
/*
 * Now, delete any internal non-unique K-cleanup
 */
   prev = ukb;
   mp = ukb->next;
   while (mp)
   {
      if (ExactKernelInList(mmb, mp))
         mp = prev->next = KillMMNode(mp);
      else
      {
         prev = mp;
         mp = mp->next;
      }
   }
   return(ukb);
}
/*
 * RETURNS: mmb, with all repeated kernels removed
 */
ATL_mmnode_t *RemoveNonUniqueKernels(ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mp;
   if (!mmb || !mmb->next)
      return(mmb);
   for (mp=mmb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p=mp->next, *prev=mp;
      while (p)
      {
         if (MMKernsSame1(mp, p))
            prev->next = p = KillMMNode(p);
         else
         {
            prev = p;
            p = p->next;
         }
      }
   }
   return(mmb);
}

void GenAllFiles(char pre, char *outd, ATL_mmnode_t *geb, ATL_mmnode_t *gek1b,
                 ATL_mmnode_t *sqb, ATL_mmnode_t *sqk1b, ATL_mmnode_t *rkb,
                 ATL_mmnode_t *trsmb, ATL_mmnode_t *cgeb, ATL_mmnode_t *cgek1b,
                 ATL_mmnode_t *csqb, ATL_mmnode_t *csqk1b, ATL_mmnode_t *crkb,
                 ATL_mmnode_t *ctrsmb)
{
   ATL_mmnode_t *uge, *usq, *urk;  /* lists w/o object repeats */
   ATL_mmnode_t *mp;
   const char cpr = (pre == 'd') ? 'z' : 'c';

   if (!ctrsmb)
      ctrsmb = trsmb;
   GenHeaderFiles(pre, outd, geb, gek1b, sqb, sqk1b, rkb, trsmb);
   GenHeaderFiles(cpr, outd, cgeb, cgek1b, csqb, csqk1b, crkb, ctrsmb);
/*
 * Link complex and real together: real/cplx don't matter for generation
 * or compiling
 */
   ATL_LastMMNode(geb)->next = cgeb;
   ATL_LastMMNode(gek1b)->next = cgek1b;
   ATL_LastMMNode(sqb)->next = csqb;
   ATL_LastMMNode(sqk1b)->next = csqk1b;
   mp = AddUniqueMMKernCompList(NULL, rkb);
   KillAllMMNodes(rkb);
   rkb = AddUniqueMMKernCompList(mp, crkb);
   KillAllMMNodes(crkb);

   GenMakefile(pre, outd, geb, gek1b, sqb, sqk1b, rkb);
   GenKerns(pre, outd, geb, gek1b, sqb, sqk1b, rkb);
   KillAllMMNodes(geb);
   KillAllMMNodes(gek1b);
   KillAllMMNodes(sqb);
   KillAllMMNodes(sqk1b);
   KillAllMMNodes(rkb);
   KillAllMMNodes(trsmb);
}

void Cplx2RealMM(char pre, ATL_mmnode_t *mmb, char *outd)
{
   ATL_mmnode_t *mp;
   const int mask = ~(1<<MMF_COMPLEX);
   char pfx[8];
   const char cpr = (pre == 'd') ? 'z' : 'c';
   pfx[0] = 'A'; pfx[1] = 'T'; pfx[2] = 'L'; pfx[3] = '_';
   pfx[4] = cpr; pfx[5] = 0;
   for (mp=mmb; mp; mp = mp->next)
   {
      char *sp;
      mp->flag &= mask;
      if (mp->rout)
      {
         sp = strstr(mp->rout, pfx);
         if (sp)
            sp[4] = pre;
      }
   }
   FillInGenStrings(pre, mmb, outd);
}

int main(int nargs, char **args)
{
   char pre='d';
   int verb=1;
   int *nbs;
   char *outd, *ukin, *kcin, *rkin, *sqin;
   ATL_mmnode_t *geb, *gek1b, *sqb, *sqk1b, *rkb, *trsm;
   ATL_mmnode_t *cgeb, *cgek1b, *csqb, *csqk1b, *crkb, *ctrsm;
/*
 * We expect to get a list of rectangular blocking factor kernels for gemm,
 * K-cleanup kernels for gemm, and square kernels for building amm-based L3BLAS,
 * and their K-cleanup kernels, and finally a list of kernels for rank-K
 */
   pre = GetFlags(nargs, args, &outd, &geb, &gek1b, &sqb, &sqk1b, &rkb, &trsm,
                  &cgeb, &cgek1b, &csqb, &csqk1b, &crkb, &ctrsm);
   assert(geb && gek1b);
   assert(sqb && sqk1b);
   assert(rkb);
/*
 * Complex actually handled by real
 */
   if (pre == 'c')
      pre = 's';
   else if (pre == 'z')
      pre = 'd';
/*
 * Fill in correct generator strings for all lists
 */
   FillInGenStrings(pre, geb, outd);
   FillInGenStrings(pre, gek1b, outd);
   FillInGenStrings(pre, sqb, outd);
   FillInGenStrings(pre, sqk1b, outd);
   FillInGenStrings(pre, rkb, outd);
   Cplx2RealMM(pre, cgeb, outd);
   Cplx2RealMM(pre, cgek1b, outd);
   Cplx2RealMM(pre, csqb, outd);
   Cplx2RealMM(pre, csqk1b, outd);
   Cplx2RealMM(pre, crkb, outd);

   GenAllFiles(pre, outd, geb, gek1b, sqb, sqk1b, rkb, trsm,
               cgeb, cgek1b, csqb, csqk1b, crkb, ctrsm);

   return(0);
}
