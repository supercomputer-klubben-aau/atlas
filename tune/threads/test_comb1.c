#include "atlas_misc.h"
#include <math.h>  /* for fabs */
#include "atlas_threads.h"
#include "atlas_cbc.h"
#include "atlas_tst.h"

TYPE *DATP, *MAXP, *MINP, *SUMP;
int *NERRS;
int NPROB, NREP, NTHR;

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -p <#> : # of problems to cycle through\n");
   fprintf(stderr, "   -r <#> : # of reps to test\n");
   fprintf(stderr, "   -T <#> : # of threads to use\n");
   exit(ierr ? ierr : -1);
}

/*
 * returns: nthreads to use
 */
int GetFlags(int nargs, char **args, int *NP, int *R)
{
   int i, nthr = ATL_NTHREADS;

   *R = 50000;
   *NP = 100;
   for (i=1; i < nargs; i++)
   {
      int ii;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'p':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *NP = atoi(args[i]);
         break;
      case 'r':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *R = atoi(args[i]);
         break;
      case 'T':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         nthr = atoi(args[i]);
         break;
      }
   }
   return(nthr);
}

ATL_igegen(int M, int BAD1, int *X, int BAD2, int seed)
{
   int i;
   srand(seed);
   for (i=0; i < M; i++)
      X[i] = rand();
}


void DoCombs(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp = vpp;
   const int nthr = NTHR, nprob=NPROB, nrep=NREP;
   int r, p, nerr=0;
   TYPE *dp, *datp=DATP, *maxp=MAXP, *minp=MINP, *sump=SUMP;
   if (rank > NTHR)
      return;
/*
 * Perform nrep calls, repeating problems as necessary
 */
   dp = datp;
   for (r=p=0; r < nrep; r++)
   {
      TYPE max, min, sum;
      #ifdef TCPLX
         TYPE cmax[2], cmin[2], csum[2];
         cmax[0] = cmin[0] = csum[0] = dp[rank+rank];
         cmax[1] = cmin[1] = csum[1] = dp[rank+rank+1];
         Mjoin(PATL,comb_min)(nthr, rank, cmin, NULL);
         Mjoin(PATL,comb_max)(nthr, rank, cmax, NULL);
         Mjoin(PATL,comb_sum)(nthr, rank, csum), NULL;
         if (cmin[0] != minp[p+p] || cmin[1] != minp[p+p+1])
            nerr++;
         if (cmax[0] != maxp[p+p] || cmax[1] != maxp[p+p+1])
            nerr++;
         if (fabs(csum[0]-sump[p+p]) > 2.0*nthr*EPS ||
             fabs(csum[1]-sump[p+p+1]) > 2.0*nthr*EPS)
            nerr++;
      #else
         min = Mjoin(PATL,comb_min)(nthr, rank, dp[rank], NULL);
         max = Mjoin(PATL,comb_max)(nthr, rank, dp[rank], NULL);
         sum = Mjoin(PATL,comb_sum)(nthr, rank, dp[rank], NULL);
         if (min != minp[p])
            nerr++;
         if (max != maxp[p])
            nerr++;
         #ifdef SINT
            if (sum != sump[p])
               nerr++;
         #else
            if (fabs(sum-sump[p]) > 2.0*nthr*EPS)
               nerr++;
         #endif
      #endif
      if (++p != nprob)
         dp += nthr SHIFT;
      else
      {
         p = 0;
         dp = datp;
      }
   }
   NERRS[rank] = nerr;
}

int DoTests(int nthr, int nprob, int nrep)
{
   ATL_CUINT setsz = nthr * nprob;
   TYPE *datp, *maxp, *minp, *sump, *dp;
   #ifdef TCPLX
      TYPE max, min, rsum, isum;
   #else
      TYPE max, min, sum;
   #endif
   int *idx, *ip, *nerrs;
   int p, r, i;

   NPROB = nprob;
   NREP = nrep;
   NTHR = nthr;
/*
 * allocate space to solve nprob problems and store that many answers
 */
   DATP = datp = malloc(ATL_MulBySize(setsz+nprob+nprob+nprob));
   NERRS = nerrs = malloc(sizeof(int)*nthr);
   ATL_assert(datp && nerrs);
   MAXP = maxp = datp + (setsz SHIFT);
   MINP = minp = maxp + (nprob SHIFT);
   SUMP = sump = minp + (nprob SHIFT);
   Mjoin(PATL,gegen)(setsz, 1, datp, setsz, setsz+nthr);
/*
 * Precompute correct answers to all problems
 */
   for (dp=datp,p=0; p < nprob; p++, dp += (nthr SHIFT))
   {
      #ifdef TCPLX
         rsum = *dp;
         isum = dp[1];
         maxp[p+p]   = minp[p+p]   = rsum;
         maxp[p+p+1] = minp[p+p+1] = isum;
         max = min = fabs(rsum) + fabs(isum);
         for (i=2; i < nthr+nthr; i += 2)
         {
            const TYPE rdv = dp[i], idv = dp[i+1];
            const TYPE dv = fabs(rdv) + fabs(idv);
            rsum += rdv;
            isum += idv;
            if (dv > max)
            {
               max = dv;
               maxp[p+p] = rdv;
               maxp[p+p+1] = idv;
            }
            if (dv < min)
            {
               min = dv;
               minp[p+p] = rdv;
               minp[p+p+1] = idv;
            }
         }
         sump[p+p]   = rsum;
         sump[p+p+1] = isum;
      #else
         max = min = sum = *dp;
         for (i=1; i < nthr; i++)
         {
            const TYPE dv = dp[i];
            sum += dv;
            min = Mmin(min, dv);
            max = Mmax(max, dv);
         }
         maxp[p] = max;
         minp[p] = min;
         sump[p] = sum;
      #endif
   }
   ATL_goParallel(nthr, DoCombs, NULL, NULL, NULL);
/* get mutex for weakly-ordered! */
   free(datp);
   for (p=i=0; i < nthr; i++)
      p += NERRS[i];
   free(NERRS);
   return(p);
}

int main(int nargs, char **args)
{
   int nerr, nthr, nprob, nrep;

   nthr = GetFlags(nargs, args, &nprob, &nrep);
   nerr = DoTests(nthr, nprob, nrep);
   if (!nerr)
      printf("ALL TESTS PASS\n");
   else
      printf("FAILED %d TESTS\n", nerr);
   return(nerr);
}
