#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "atlas_threads.h"
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_level1.h"
#include "atlas_gentesttime.h"
#define ATL_GETFLAGS 1
#include "atlas_genparse.h"

static double wbeg[ATL_NTHREADS];

static void DoWork(void *vpp, int rank, int vrank)
{
   wbeg[vrank] = ATL_walltime();
}


void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -# <nreps>\n");
   fprintf(stderr, "   -P <nthreads> (max=%d)\n", ATL_NTHREADS);
   fprintf(stderr, "   -v <verb> [1] do/don't print indiv reps\n");
   exit (ierr ? ierr : -1);
}

FILE *GetFlags(int nargs, char **args, unsigned long *REPS, int *P, int *VERB)
{
   int p=ATL_NTHREADS, i, verb=0;
   unsigned long reps=40000;
   FILE *fpout=stdout;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case '#':   /* -# <reps> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         reps = atol(args[i]);
         break;
      case 'P':   /* -P <p> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         p = atoi(args[i]);
         break;
      case 'v':   /* -v <verb> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         verb = atoi(args[i]);
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
   *P = p;
   *VERB = verb;
   *REPS = reps;
   return(fpout);
}

void DoLine(FILE *fpout, int P, double t0)
{
   int i;
   fprintf(fpout, "%e", 1e6*(wbeg[0]-t0));
   t0 = wbeg[0];
   for (i=1; i < P; i++)
      fprintf(fpout, " %.2e", (wbeg[i]-wbeg[i-1])*1e6);
   fprintf(fpout, "\n");
   fflush(fpout);
}
void DoLine_file(FILE *fpout, int P, double t0)
{
   int i;
   fprintf(fpout, "%le", 1e6*(wbeg[0]-t0));
   t0 = wbeg[0];
   for (i=1; i < P; i++)
      fprintf(fpout, " %le", (wbeg[i]-wbeg[i-1])*1e6);
   fprintf(fpout, "\n");
   fflush(fpout);
}

int main(int nargs, char **args)
{
   FILE *fpout;
   unsigned long reps, r;
   double t0, maxInit=0.0, minInit=1e10, avgInit=0.0;
   double maxGap=0.0, minGap=1e10, avgGap=0.0;
   int P, verb;

   fpout = GetFlags(nargs, args, &reps, &P, &verb);
   void (*DoLn)(FILE *fpout, int P, double t0);

   DoLn = (fpout == stdout || fpout== stderr) ? DoLine : DoLine_file;
   ATL_goParallel(P, DoWork, NULL, NULL, NULL); /* start thread pool */
   fprintf(fpout, "TIMING IN MICROSECONDS FOR P=%d, R=%lu:\n", P, reps);
   for (r=0; r < reps; r++)
   {
      double d;
      int p;
      t0 = ATL_walltime();
      ATL_goParallel(P, DoWork, NULL, NULL, NULL);
      if (r == 0)
         continue;
      SortDoubles(ATL_NTHREADS, wbeg);
      d = wbeg[0] - t0;
      minInit = Mmin(minInit, d);
      maxInit = Mmax(maxInit, d);
      avgInit += d;
      for (p=1; p < P; p++)
      {
         d = wbeg[p] - wbeg[p-1];
         minGap = Mmin(minGap, d);
         maxGap = Mmax(maxGap, d);
         avgGap += d;
      }
      if (verb)
         DoLn(fpout, P, t0);
   }
   fprintf(fpout, "DONE\n");
   avgInit /= reps;
   avgGap /= reps;
   if (reps > 999)
   {
      FILE *fpres;
      unsigned long orep=0;
      fpres = fopen("res/thrbeg.sum", "r");
      if (fpres)
      {
         if (fscanf(fpres, "%lu", &orep) != 1)
            orep = 0;
         fclose(fpres);
      }
      if (orep < reps)
      {
         fpres = fopen("res/thrbeg.sum", "w");
         if (fpres)
         {
            fprintf(fpres, "%lu\n%le\n%le\n", reps, avgInit, avgGap);
            fclose(fpres);
         }
         fpres = fopen("res/thrbeg.h", "w");
         if (fpres)
         {
            fprintf(fpres, "   #define ATL_tstart_sec %le\n", avgInit);
            fprintf(fpres, "   #define ATL_tstartgap_sec %le\n", avgGap);
            fclose(fpres);
         }
      }
   }
   minInit *= 1.0e6;
   minGap *= 1.0e6;
   maxInit *= 1.0e6;
   maxGap *= 1.0e6;
   avgInit *= 1.0e6;
   avgGap *= 1.0e6;
   fprintf(fpout, "\nInit summary: max=%e, min=%e, avg=%e\n",
           maxInit, minInit, avgInit);
   fprintf(fpout, "Gap  summary: max=%e, min=%e, avg=%e\n",
           maxGap, minGap, avgGap);
   if (fpout != stderr && fpout != stdout)
      fclose(fpout);
   return(0);
}
