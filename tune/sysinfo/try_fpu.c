#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define INTV 0.2
#define DEFMF 200.0
/*
 * This routine written in assembly to allow useless computations to be used
 * in order to acheive FPU peak.
 * if nrep is 0, then the function returns the number of flops done in a single
 * loop iteration.  If nrep > 0, it performs that many loop iterations and
 * returns with junk as return value.
 * The array d is input/output of length at least 5000, and can be used
 * to make compiler-oriented routine unsure if flops are useless.
 */
unsigned long fpuStress(unsigned long nrep, double *d);
void printUsage(char *name, int iarg, char *flag)
{
   if (iarg)
      fprintf(stderr, "Error around argument %d (%s)!\n", iarg, flag);
   fprintf(stderr, "USAGE: %s [flags]\n", name);
   fprintf(stderr, "   -m <mflop> : set the number of MFLOP to time");
   fprintf(stderr, "   -t <sec> : set floor on timing interval");
   fprintf(stderr, "   -o <outfile> : default=stdout\n");
   exit(iarg ? iarg : -1);
}

FILE *getFlags(int nargs, char **args, double *mflop, double *time)
{
   FILE *fpout=stdout;
   int i;

   *mflop = *time = 0.0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         printUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'm':
         if (++i >= nargs)
            printUsage(args[0], i, "Out of args");
         *mflop = atof(args[i]);
         if (*mflop <= 0.0)  /* <=0 says use default MFLOP */
            *mflop = DEFMF;
         *time = 0.0;
         break;
      case 't':
         if (++i >= nargs)
            printUsage(args[0], i, "Out of args");
         *time = atof(args[i]);
         *mflop = 0.0;
         break;
      case 'o':
         if (++i >= nargs)
            printUsage(args[0], i, "Out of args");
         if (fpout != stdout)
            fclose(fpout);
         fpout = fopen(args[i], "w");
         assert(fpout);
         break;
      default:
         printUsage(args[0], i, args[i]);
      }
   }
   return(fpout);
}
/*
 * This program is linked to various assembly language backends to attempt
 * to exercise the FPU at its peak rate.  Once the correct backend is linked
 * in, we can use this peak performance to auto-discover what thread IDs
 * share an FPU, and thus which IDs should be used
 */
#define NTRY 3
int main(int nargs, char **args)
{
   double *d;
   double time, mf;
   FILE *fpout;
   double ATL_walltime(void);

   d = malloc(sizeof(double)*100);
   assert(d);
   fpout = getFlags(nargs, args, &mf, &time);
   if (mf > 0.0)
   {
      unsigned long nrep, fpi, k;
      double t0, t1, avg=0.0;

      fpi = fpuStress(0, d);      /* get # of flops in 1 loop iteration */
      nrep = mf/(1e-6*fpi);       /* compute # of its to get desired MFLOPs */

      for (k=0; k < NTRY; k++)
      {
         t0 = ATL_walltime();
         fpuStress(nrep, d);
         t1 = ATL_walltime() - t0;
         if (fpout == stdout)
            fprintf(fpout, "   MFLOP=%.2f (%e sec)\n", nrep*fpi*1e-6/t1, t1);
         avg += t1;
      }
      avg /= (double)NTRY;
      if (fpout == stdout)
         fprintf(fpout, "AVG MFLOP=%.2f\n", nrep*fpi*1e-6/avg);
      else
         fprintf(fpout, "%le\n",  nrep*fpi*1e-6 / avg);
   }
   else
     assert(0);
   free(d);
   return(0);
}
