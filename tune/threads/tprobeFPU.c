#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atlas_fopen.h"
#define ATL_DEF_RUNTHR 1
#include "atlas_tprim.h"
/*
 * The idea here is to find cores that share significant FPU hardware
 * resources, and eliminate them as threading targets for ATLAS, which is
 * essentially FPU-bound.  This should eliminate using hyperthreading,
 * for instance.  We read in the fpu.S produced by probeFPU, which is
 * assembly code that usually drives the FPU at peak rate while accessing
 * almost no memory and requiring negligible ALU usage.
 *
 * We will first run this code with only 1 thread active on core=0 to get
 * the base case.  We will then add core 1, and retime: if 0's performance
 * drops by >= 8%, we declare that as a virtual thread and eliminate it
 * from the set of threads to use.  We keep adding core IDs until every
 * core has either been added to the non-aliased list or eliminated.
 */
static int nact=1, actRanks[ATL_AFF_NUMID]={0};
static unsigned long NREP=0;
static volatile unsigned chkin[ATL_AFF_NUMID]={0};
static double PEAK=27065.0, times[ATL_AFF_NUMID];
#define NTRY 3

/*
 * This routine written in assembly to allow useless computations to be used
 * in order to achieve FPU peak.
 * if nrep is 0, then the function returns the number of flops done in a single
 * loop iteration.  If nrep > 0, it performs that many loop iterations and
 * returns with junk as return value.
 * The array d is input/output of length at least 5000, and can be used
 * to make compiler-oriented routine unsure if flops are useless.
 */
unsigned long fpuStress(unsigned long nrep, double *d);
double ATL_walltime(void);

static void *DoWork(void *vp)
{
   ATL_thread_t *tp=vp;
   double avg=0.0, t0;
   unsigned long nrep;
   unsigned int i;
   const unsigned int rank = tp->rank;

/*   printf("%u,%u:\n", rank, vrank); */
   nrep = NREP;
/*
 * Barrier until all active & inactive threads arrive
 */
   if (!rank)
   {
      for (i=1; i < ATL_AFF_NUMID; i++)
         while(!chkin[i]);
      *chkin = 1;
   }
   else
   {
      chkin[rank] = 1;
      while (!(*chkin));
   }
/*
 * If my core isn't active in this test, just return;
 */
   for (i=0; rank != actRanks[i] && i < nact; i++);
   if (i == nact)
   {
/*      printf("RETURNING: nact=%u, rank=%u\n", nact, rank); */
      return(NULL);
   }
/*
 * Take average of NTRY timings
 */
   for (i=0; i < NTRY; i++)
   {
      t0 = ATL_walltime();
      fpuStress(nrep, NULL);
      avg += ATL_walltime() - t0;
   }
   avg /= NTRY;
/*printf("times[%d]=%e ntry=%u, nrep=%lu)\n", rank, avg, NTRY, nrep); */
   times[rank] = avg;

   if (nact > 1)
   {
      if (!rank)
      {
         for (i=1; i < nact; i++)
            while(!chkin[i]);
         *chkin = 1;
      }
      else
      {
         chkin[rank] = 1;
         while (!(*chkin));
      }
   }
   return(NULL);
}

void PrintUsage(char *name, int iarg, char *arg)
{
   fprintf(stderr, "\nERROR around arg %d (%s).\n", iarg,
           arg ? arg : "unknown");
   fprintf(stderr, "USAGE: %s [flags] where flags are:\n", name);
   fprintf(stderr, "   -o <affinity in/outfile>\n");
   exit(iarg? iarg:-1);
}
char *GetFlags(int nargs, char **args)
{
   int i;
   char *fnam=NULL;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'o':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         fnam = args[i];
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   return(fnam);
}

int main(int nargs, char **args)
{
   int i;
   unsigned long nf;
   double mf;
   char *fnam;
   fnam = GetFlags(nargs, args);
/*
 * Compute how many reps to run around .05 seconds at detected peak
 */
   nf = fpuStress(0, NULL);            /* find # of flops in 1 it of loop */
   NREP = (0.05*PEAK*1e6+nf-1)/nf;     /* .05 above most wallclock res */
   mf = NREP*1e-6*nf;
/*
 * First time with only ID=0 running
 */
   ATL_runThreads(ATL_AFF_NUMID, DoWork, NULL);
   printf("\nSingle-node timings = %e sec, mf=%.0f\n", *times, mf/(*times));
/*
 * Now try adding each affID in turn, and look for performance regressions
 */
   for (i=1; i < ATL_AFF_NUMID; i++)
   {
      double oldt[ATL_AFF_NUMID];
      int k, NO=0;
      for (k=0; k < nact; k++)
         oldt[k] = times[k];
      actRanks[nact++] = i;
      for (k=0; k < ATL_AFF_NUMID; k++)
         chkin[k] = 0;
      ATL_runThreads(ATL_AFF_NUMID, DoWork, NULL);
      for (k=0; k < nact-1; k++)
      {
/*
 *       If new timing more than 10% slower than old, declare these threads
 *       to share too much hardware to be used as indepedent thread
 */
         if (oldt[k]*1.1 < times[k])
         {
            NO++;
            printf("      affID %u seems %.0f%% aliased with affID %u "
                   "(old=%.2e,new=%.2e)!\n",
            #ifndef ATL_RANK_IS_PROCESSORID
                   ATL_affinityIDs[i], (oldt[k]/times[k])*100.0,
                   ATL_affinityIDs[k], oldt[k], times[k]);
            #else
                   i, (oldt[k]/times[k])*100.0, k, oldt[k], times[k]);
            #endif
         }
      }
      if (NO)
      {
         #ifndef ATL_RANK_IS_PROCESSORID
            printf("   Dependent thread ID %u rejected!\n", ATL_affinityIDs[i]);
         #else
            printf("   Dependent thread ID %u rejected!\n", i);
         #endif
         nact--;
      }
      else
         #ifndef ATL_RANK_IS_PROCESSORID
            printf("   Independent thread ID %u accepted.\n",
                   ATL_affinityIDs[i]);
         #else
            printf("   Independent thread ID %u accepted.\n", i);
         #endif

   }
   #ifndef ATL_RANK_IS_PROCESSORID
      printf("Independent IDs={%u", ATL_affinityIDs[actRanks[0]]);
      for (i=1; i < nact; i++)
         printf(",%u", ATL_affinityIDs[actRanks[i]]);
   #else
      printf("Independent IDs={0");
      for (i=1; i < nact; i++)
         printf(",%u", actRanks[i]);
   #endif
   printf("}\n");
/*
 * If anything changed, and we've been told to overwrite the affinity include
 * file, do so
 */
   if (nact < ATL_AFF_NUMID && fnam)
   {
      char ln[128];
      FILE *fpin, *fpout;
      int gap=0, rankIsID;
      assert(nact);
      for (i=0; i < nact; i++)
      #ifndef ATL_RANK_IS_PROCESSORID
         if (ATL_affinityIDs[actRanks[i]] != i)
      #else
         if (actRanks[i] != i)
      #endif
            break;
      rankIsID = (i == nact);
/*
 *    See if IDs start at zero, and then go with a constant stride: if so,
 *    we can use NUMID,IDSTRIDE,RANK_IS_PROCESSORID
 *    Otherwise, must use explicit affinityIDs list
 */
      #ifndef ATL_RANK_IS_PROCESSORID
         if (ATL_affinityIDs[actRanks[0]] == 0)
         {
            if (nact > 2)
            {
               gap = ATL_affinityIDs[actRanks[1]]-ATL_affinityIDs[actRanks[0]];
               for (i=2; i < nact; i++)
                  if ((ATL_affinityIDs[actRanks[i]]-
                      ATL_affinityIDs[actRanks[i-1]]) != gap)
                     break;
               gap = (i < nact) ? 0 : gap;
            }
         }
      #else
         if (actRanks[0] == 0)
         {
            if (nact > 2)
            {
               gap = actRanks[1] - actRanks[0];
               for (i=2; i < nact; i++)
                  if (actRanks[i]-actRanks[i-1] != gap)
                     break;
               gap = (i < nact) ? 0 : gap;
            }
         }

      #endif
      assert(!rename(fnam, "res/atlas_affinity_aliased.h"));
      fpin = fopen("res/atlas_affinity_aliased.h", "r");
      assert(fpin);
      fpout = fopen(fnam, "w");
      assert(fpout);
      while(fgets(ln, 128, fpin))
      {
/*
 *       Find line(s) of file that needs to change
 */
         if (strstr(ln, "ATL_AFF_NUMID"))
         {
            fprintf(fpout, "#define ATL_AFF_NUMID %u\n", nact);
            if (rankIsID)
               fprintf(fpout, "#define ATL_RANK_IS_PROCESSORID 1 "
                       "/* good IDs [0,%u] */\n", nact-1);
         }
         else if (strstr(ln, "ATL_AFF_IDSTRIDE"))
            fprintf(fpout, "#define ATL_AFF_IDSTRIDE %u\n", gap);
         else if (strstr(ln, "ATL_RANK_IS_PROCESSORID"))
            ;  /* line handled above, so don't output it */
         else if (strstr(ln, "static int ATL_affinityIDs"))
         {
            while (!strstr(ln, "};")) /* eat lns until end of array found */
            {
               assert(fgets(ln, 128, fpin));
            }
            if (!rankIsID)
            {
               fprintf(fpout, "static int ATL_affinityIDs[%u]\n", nact);
               #ifndef ATL_RANK_IS_PROCESSORID
                  fprintf(fpout, "   = {%u", ATL_affinityIDs[actRanks[0]]);
                  for (i=1; i < nact; i++)
                     fprintf(fpout, ",%u", ATL_affinityIDs[actRanks[i]]);
               #else
                  fprintf(fpout, "   = {%u", actRanks[0]);
                  for (i=1; i < nact; i++)
                     fprintf(fpout, ",%u", actRanks[i]);
               #endif
               fprintf(fpout, "};\n");
            }
         }
         else
            fputs(ln, fpout);
      }
      fclose(fpout);
      fclose(fpin);
/*
 *    Now patch NPROC in Make.inc, and re-install archdefs
 */
      sprintf(ln, "make patch_archdefs nproc=%u", nact);
      assert(system(ln) == 0);
   }
   return(0);
}
