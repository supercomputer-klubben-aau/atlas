#include "atlas_threads.h"
#include "atlas_misc.h"
#include "atlas_atmctr.h"
#include "atlas_gatmctr.h"
#define ATL_GETFLAGS 1
#include "atlas_genparse.h"
#ifdef ATL_USE_PERF_H
   #define NOPERF    1
   #define NOK1RATIO 1
   #include Mstr(Mjoin(ATLAS_PRE,geamm_perf.h))
#endif

void *ATL_atmctr_new_mut(long N);
void ATL_atmctr_free_mut(void *ac);
long ATL_atmctr_dec_mut(void *ac);


static double DELAY, timearr[ATL_NTHREADS];
static unsigned long cntPerRank[ATL_NTHREADS];
static void *ACNT;

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -P <nthreads> (max=%d)\n", ATL_NTHREADS);
   fprintf(stderr,
   "   -m <base> <multiplier> <max> (base in microsec, max in sec)\n");
   fprintf(stderr, "   -t <sec>: aim to spend sec seconds per delay line\n");
   exit (ierr ? ierr : -1);
}

FILE *GetFlags(int nargs, char **args, double *TPL, int *P, double *mula)
{
   int p=ATL_NTHREADS, i, verb=0;
   FILE *fpout=stdout;
   double tpl=8.0;

   mula[0] = 1e-6; mula[1] = 2.0; mula[2] = .01;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'm' : /* -m <base> <mul> <max> */
         if (i+3 >= nargs)
            PrintUsage(args[0], i, NULL);
         mula[0] = atof(args[i+1]) * 1e-6;
         mula[1] = atof(args[i+2]);
         mula[2] = atof(args[i+3]);
         assert(mula[0] > 0.0 && mula[2] >= mula[0] && mula[1] >= 0.0);
         i += 3;
         break;
      case 't':   /* -t <tpl> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         tpl = atof(args[i]);
         break;
      case 'P':   /* -P <p> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         p = atoi(args[i]);
         break;
      case 'o':  /* -o <outfile> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         if (!strcmp(args[i], "stdout"))
            fpout = stdout;
         else if (!strcmp(args[i], "stderr"))
            fpout = stderr;
         else
         {
            fpout = fopen(args[i], "w");
            assert(fpout);
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   *TPL = tpl;
   *P = p;
   return(fpout);
}

static void DoWork_glb(void *vpp, int rank, int vrank)
{
   double t0, t1, t2;
   unsigned long cnt=0;
/*
 * We don't trust walltime() to have enough accuracy to time one call to
 * DecAtom, so add all them together, and then subtract out cnt*DELAY;
 * This isn't perfect either, but should be better wt low-res clocks
 */
   t0 = ATL_walltime();
   while (ATL_gatmctr_dec(ACNT,vrank))
   {
      t1 = ATL_walltime();
      cnt++;
      do
         t2 = ATL_walltime();
      while (t2-t1 < DELAY);
      t0 += t2-t1 - DELAY;    /* correct any extra delay */
   }
   t1 = ATL_walltime();
   cntPerRank[rank] = cnt;
   timearr[rank] = t1 - t0 - DELAY*cnt;
}
static void DoWork_loc(void *vpp, int rank, int vrank)
{
   double t0, t1, t2;
   unsigned long cnt=0;
/*
 * We don't trust walltime() to have enough accuracy to time one call to
 * DecAtom, so add all them together, and then subtract out cnt*DELAY;
 * This isn't perfect either, but should be better wt low-res clocks
 */
   t0 = ATL_walltime();
   while (ATL_atmctr_dec(ACNT))
   {
      t1 = ATL_walltime();
      cnt++;
      do
         t2 = ATL_walltime();
      while (t2-t1 < DELAY);
      t0 += t2-t1 - DELAY;    /* correct any extra delay */
   }
   t1 = ATL_walltime();
   cntPerRank[rank] = cnt;
   timearr[rank] = t1 - t0 - DELAY*cnt;
}
static void DoWork_mut(void *vpp, int rank, int vrank)
{
   double t0, t1, t2;
   unsigned long cnt=0;
/*
 * We don't trust walltime() to have enough accuracy to time one call to
 * DecAtom, so add all them together, and then subtract out cnt*DELAY;
 * This isn't perfect either, but should be better wt low-res clocks
 */
   t0 = ATL_walltime();
   while (ATL_atmctr_dec_mut(ACNT))
   {
      t1 = ATL_walltime();
      cnt++;
      do
         t2 = ATL_walltime();
      while (t2-t1 < DELAY);
      t0 += t2-t1 - DELAY;    /* correct any extra delay */
   }
   t1 = ATL_walltime();
   cntPerRank[rank] = cnt;
   timearr[rank] = t1 - t0 - DELAY*cnt;
}

double PerCallSec(int P)
{
   double cost=0.0;
   if (P)
   {
      unsigned long cnt=0;
      int i;
      for (i=0; i < P; i++)
      {
         cnt += cntPerRank[i];
         cost += timearr[i];
      }
      cost /= cnt;
   }
   return(cost);
}

#define PRVDEF ATL_GAC_PRV
#define PUBDEF ATL_GAC_PUB
#define STLDEF ATL_GAC_MIX
int main(int nargs, char **args)
{
   FILE *fpout;
   int P;
   unsigned long cnt;
   double ma[3], del, max, mul, tpl;
   #ifdef ATL_USE_PERF_H
      #define NFILE 5
      char fnam[13]={'r','e','s','/',0, 0,0,0,'.','t','i','m'};
      char *NMs[NFILE]={"mut", "lac", "pub", "mix", "prv"};
      FILE *FPs[NFILE];
      int I;
   #endif

   fpout=GetFlags(nargs, args, &tpl, &P, ma);
   #ifdef ATL_USE_PERF_H
      for (I=0; I < NFILE; I++)
      {
         #ifdef DCPLX
            fnam[4] = 'z'
         #elif defined(SCPLX)
            fnam[4] = 'c';
         #elif defined(SREAL)
            fnam[4] = 's';
         #else
            fnam[4] = 'd';
         #endif
         fnam[5] = NMs[I][0];
         fnam[6] = NMs[I][1];
         fnam[7] = NMs[I][2];
         FPs[I] = fopen(fnam, "w");
         assert(FPs[I]);
         if (!I)
            fprintf(FPs[0], "#define ATL_AMM_NSCHED %d\n", ATL_geAMM_NCASES);
         fprintf(FPs[I], "#if !defined(NoARRS) && !defined(NOSCHED%s)\n",
                 NMs[I]);
         fprintf(FPs[I], "static const float ATL_AMM_SCHED%s[%d] =\n{\n",
                 NMs[I], ATL_geAMM_NCASES);
      }
   #endif
   fprintf(fpout,
"   DELAY      COUNT    ATMCTR_mut   ATMCTR    g_PUB    g_STL    g_PRV\n");
   fprintf(fpout,
"========  =========  ============  =======  =======  =======  =======\n");
   del = 0.0; mul = ma[1]; max = ma[2];
   #ifdef ATL_USE_PERF_H
   for (I=0; I < ATL_geAMM_NCASES; I++)
   #else
   while (del <= max)
   #endif
   {
      double tmut, tloc=0, tprv, tpub, tstl;
      unsigned long cnt;

      #ifdef ATL_USE_PERF_H
         del = ATL_geAMM_TIME[I];
      #endif
      DELAY = del;
      if (del <= 0.0)
         cnt = tpl / (5.0e-6);
      else
         cnt = tpl / (5.0*del);
      cnt = Mmax(1000, cnt);

      ACNT = ATL_gatmctr_new(P, cnt, STLDEF);
      ATL_goParallel(P, DoWork_glb, NULL, NULL, NULL);
      ATL_gatmctr_free(ACNT);
      tstl = PerCallSec(P);

      ACNT = ATL_gatmctr_new(P, cnt, PUBDEF);
      ATL_goParallel(P, DoWork_glb, NULL, NULL, NULL);
      ATL_gatmctr_free(ACNT);
      tpub = PerCallSec(P);

      ACNT = ATL_gatmctr_new(P, cnt, PRVDEF);
      ATL_goParallel(P, DoWork_glb, NULL, NULL, NULL);
      ATL_gatmctr_free(ACNT);
      tprv = PerCallSec(P);

      ACNT = ATL_atmctr_new(cnt);
      ATL_goParallel(P, DoWork_loc, NULL, NULL, NULL);
      ATL_atmctr_free(ACNT);
      tloc = PerCallSec(P);

      ACNT = ATL_atmctr_new_mut(cnt);
      ATL_goParallel(P, DoWork_mut, NULL, NULL, NULL);
      ATL_atmctr_free_mut(ACNT);
      tmut = PerCallSec(P);

      fprintf(fpout, "%.2e %10d  %12e %8.2f %8.2f %8.2f %8.2f\n", del,cnt,tmut,
              tmut/tloc, tmut/tpub, tmut/tstl, tmut/tprv);
      #ifdef ATL_USE_PERF_H
         fprintf(FPs[0], "   %e,  /* CNT=%u */\n", tmut, I);
         fprintf(FPs[1], "   %e,  /* CNT=%u */\n", tloc, I);
         fprintf(FPs[2], "   %e,  /* CNT=%u */\n", tpub, I);
         fprintf(FPs[3], "   %e,  /* CNT=%u */\n", tstl, I);
         fprintf(FPs[4], "   %e,  /* CNT=%u */\n", tprv, I);
      #endif
      #ifndef ATL_USE_PERF_H
         if (del == 0.0)
            del = ma[0];
         else
            del *= mul;
      #endif
   }
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
   #ifdef ATL_USE_PERF_H
      for (I=0; I < NFILE; I++)
      {
         fprintf(FPs[I], "};\n#endif\n");
         fclose(FPs[I]);
      }
      fnam[5] = 'd'; fnam[6] = 'l'; fnam[7] = 'y';
      fpout = fopen(fnam, "w");
      assert(fpout);
      fprintf(fpout, "#define ATL_AMM_NSCHED %d\n", ATL_geAMM_NCASES);
      fprintf(fpout, "#if !defined(NOARRS) && !defined(NOTIMEdly)\n");
      fprintf(fpout, "static const float ATL_AMM_SCHEDdly[%d] =\n{\n",
              ATL_geAMM_NCASES);
      for (I=0; I < ATL_geAMM_NCASES; I++)
         fprintf(fpout, "   %e,  /* CNT=%u */\n", ATL_geAMM_TIME[I], I);
      fprintf(fpout, "};\n#endif\n");
      fclose(fpout);
   #endif
   return(0);
}
