#include "atlas_threads.h"
#include "atlas_misc.h"
#include <assert.h>
static int NTHR=0, NREP=100000;
volatile char *bchk;
static double *timearr=NULL;
#define ATL_membarrier()

#define TimePT 1
#ifdef TimePT
   #include <pthread.h>
   static pthread_barrier_t ptb;
#endif
void ATL_cbc_barrier
(
   ATL_CUINT P,    /* # of threads to barrier */
   ATL_CUINT iam   /* rank of calling thread in barrier */
)
{
   ATL_UINT d;
   const char newv = !bchk[iam];
/*
 * Any thread that did not participate in this job, simply updates
 * its entry to match what the active threads will change to, and
 * returns w/o waiting on barrier (they are not part of barrier!).
 * NOTE: all cores must call this routine to avoid having the boolean array
 *       get out of sync, or this array must be reset manually during serial
 *       execution.
 *
 */
#if 0
   if (iam > P)
   {
      bchk[iam] = newv;
      return;
   }
#endif
   for (d=1; iam+d < P; d <<= 1)
   {
      if ((iam>>(d-1))&1) /* partner on right is done */
         break;
      while (bchk[iam+d] != newv); /* await partner signal */
   }
   bchk[iam] = newv;
   if (iam)
      while (*bchk != newv);
}

void ATL_cbc_barrlinGap
(
   ATL_CUINT P,    /* # of threads to barrier */
   ATL_CUINT iam   /* rank of calling thread in barrier */
)
{
   ATL_CUINT II = iam << 7;
   const char newv = !bchk[II];

   if (iam)
   {
      bchk[II] = newv;
      while (*bchk != newv);
   }
   else
   {
      int i;
      for (i=1; i < P; i++)
      {
         ATL_CUINT d = i<<7;
         while (bchk[d] != newv);
      }
      *bchk = newv;
   }
}

void ATL_cbc_barrlin
(
   ATL_CUINT P,    /* # of threads to barrier */
   ATL_CUINT iam   /* rank of calling thread in barrier */
)
{
   const char newv = !bchk[iam];

   if (iam)
   {
      bchk[iam] = newv;
      while (*bchk != newv);
   }
   else
   {
      int i;
      for (i=1; i < P; i++)
         while (bchk[i] != newv);
      *bchk = newv;
   }
}

void DoPTB(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp = vpp;
   int i;
   double t0;
   t0 = ATL_walltime();
   for (i=0; i < NREP; i++)
   {
      pthread_barrier_wait(&ptb);
   }
   timearr[rank] = ATL_walltime() - t0;
}
void DoCBC(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp = vpp;
   int i;
   double t0;
   t0 = ATL_walltime();
   for (i=0; i < NREP; i++)
   {
      ATL_cbc_barrier(pp->nthr, vrank);
   }
   timearr[rank] = ATL_walltime() - t0;
}
void DoCBC_lin(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp = vpp;
   int i;
   double t0;
   t0 = ATL_walltime();
   for (i=0; i < NREP; i++)
   {
      ATL_cbc_barrlin(pp->nthr, vrank);
   }
   timearr[rank] = ATL_walltime() - t0;
}

void DoCBC_gap(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp = vpp;
   int i;
   double t0;
   t0 = ATL_walltime();
   for (i=0; i < NREP; i++)
   {
      ATL_cbc_barrlinGap(pp->nthr, vrank);
   }
   timearr[rank] = ATL_walltime() - t0;
}
void DoStart(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp = vpp;
   ATL_cbc_barrier(pp->nthr, vrank);
}

double PrintTimes()
{
   int i;
   double tmin, tmax, tsum;
   tmin = tmax = tsum = timearr[0];
   for (i=1; i < NTHR; i++)
   {
      const double t0 = timearr[i];
      tsum += t0;
      if (t0 > tmax)
         tmax = t0;
      if (t0 < tmin)
         tmin = t0;
   }
   printf("   TIMES: MAX=%e, MIN=%e, AVG=%e\n", tmax, tmin, tsum/NTHR);
   return(tsum);
}
int main(int nargs, char **args)
{
   int i, reps=100000;
   int nthr=ATL_NTHREADS;
   ATL_CUINT nthr128 = nthr<<7;
   double tlin, tlog, tpt=0, tgap;
   if (nargs > 1)
   {
      nthr = atoi(args[1]);
      if (nargs > 2)
         reps = atoi(args[2]);
   }
   NTHR = ATL_NTHREADS;
   NREP = reps;
   timearr = malloc(sizeof(double)*nthr);
   assert(timearr);
   bchk = calloc(nthr128, sizeof(char));
   assert(bchk);
/*
 * Start up ATLAS's threadpool by doing single barrier; we don't want spawn cost
 * in any of our timings.
 */
   ATL_goParallel(nthr, DoStart, NULL, NULL, NULL);
   printf("FINDING SPEED OF %s BARRIER USING %d THREADS AND %d REPS\n",
          "log2(CBC)", nthr, reps);
   ATL_goParallel(nthr, DoCBC, NULL, NULL, NULL);
   tlog = PrintTimes();

   printf("FINDING SPEED OF %s BARRIER USING %d THREADS AND %d REPS\n",
          "lin(CBC)", nthr, reps);
   ATL_goParallel(nthr, DoCBC_lin, NULL, NULL, NULL);
   tlin = PrintTimes();

   for (i=0; i < nthr; i++)
      bchk[i] = 0;
   printf("FINDING SPEED OF %s BARRIER USING %d THREADS AND %d REPS\n",
          "linCL(CBC)", nthr, reps);
   ATL_goParallel(nthr, DoCBC_gap, NULL, NULL, NULL);
   tgap = PrintTimes();

   #ifdef TimePT
      assert(pthread_barrier_init(&ptb, NULL, nthr) == 0);
      printf("FINDING SPEED OF %s BARRIER USING %d THREADS AND %d REPS\n",
             "PTHREAD", nthr, reps);
      ATL_goParallel(nthr, DoPTB, NULL, NULL, NULL);
      tpt = PrintTimes();
      assert(pthread_barrier_destroy(&ptb) == 0);
   #endif
   printf("LIN/LOG = %0.2f, LOG/LIN = %0.2f\n", tlin/tlog, tlog/tlin);
   printf("compressed/seperated = %0.2f\n", tlin/tgap);
   if (tpt != 0.0)
      printf("PTHREAD/LINCBC = %0.2f, LINCBC/PTHREAD = %0.2f\n",
             tpt/tlin, tlin/tpt);
   return(0);
}
