#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "atlas_genparse.h"
#define ulong unsigned long
#define uint  unsigned int
#define KSH 10 /* kilo-byte (1024) shift */
#define MSH 20 /* mega-byte (104857) shift */

int *getPivAndSort(int N, long *szs)
/*
 * Finds a pivot that will put szs in order.  Each entry in ipiv says what row
 * should be swapped to mov the next smallest size into that position.
 * Sorts szs at same time using selection sort.
 */
{
   int *ipiv;
   int i;

   ipiv = malloc(sizeof(int)*N);
   for (i=0; i < N-1; i++)
   {
      long minsz = szs[i];
      int j, imin = i;
      for (j=i+1; j < N; j++)
      {
         long sz = szs[j];
         if (sz >= minsz)
            continue;
         imin = j;
         minsz = szs[j];
      }
      ipiv[i] = imin;
      if (imin != i)
      {
         szs[imin] = szs[i];
         szs[i] = minsz;
      }
   }
   ipiv[N-1] = N-1;
   return(ipiv);
}

void pivFloats(int N, int *ipiv, float *f)
{
   int i;
   N--;
   for (i=0; i < N; i++)
   {
      int k = ipiv[i];
      if (k != i)
      {
         float tmp = f[i];
         f[i] = f[k];
         f[k] = tmp;
      }
   }
}

void PrintAllTimes(unsigned int N, unsigned int P, long *szs, float *times,
                   float *var)
{
   int i;
   for (i=0; i < N; i++)
   {
      int j;
      float *f = times + i;;
      printf("%7lu", szs[i]>>KSH);
      for (j=0; j < P; j++, f += N)
         printf(" %.2e", *f);
      printf("\n");
   }
   if (var)
   {
      printf("varian:");
      for (i=0; i < P; i++)
         printf(" %.2e", var[i]);
      printf("\n");
   }
}

float smoothRes(unsigned int N, float *fa, float *raw)
/*
 * Logically, small-case performance must always be faster than large.  Results
 * showing greater time for smaller buffers must be an error due to either
 * timing variance or due to an associativity problem (which we don't want
 * to consider here).  Before this call, raw will be sorted by size, so this
 * routine goes thru raw, and any result that is less than a prior result gets
 * set to the prior result.
 * RETURNS: ratio of greatest error to its denominator.
 */
{
   unsigned int i;
   float errRat = 0.0, max;

   if (N < 2)
      return(0.0);
   *fa = max = *raw;
   for (i=1; i < N; i++)
   {
      const float rw = raw[i];
      if (rw < max)  /* new result less than last */
      {
         float er;
         fa[i] = max;
         er = (max - rw) / max;
         errRat = (errRat >= er) ? errRat : er;
      }
      else
      {
         fa[i] = raw[i];
         max = rw;
      }
   }
   return(errRat);
}

void printRow(unsigned int N, long *szs, float *tim)
{
   unsigned int j;
   if (N < 1)
      return;
   if (szs)
   {
      printf("%8lu", *szs);
      for (j=1; j < N; j++)
         printf(" %8lu", szs[j]);
      printf("\n");
   }
   printf("%.2e", *tim);
   for (j=1; j < N; j++)
      printf(" %.2e", tim[j]);
   printf("\n");
}

void printCol(unsigned int N, long *szs, float *tim)
{
   unsigned int j;
   if (N < 1)
      return;
   if (szs)
   {
      printf("%9lu %e\n", *szs, *tim);
      for (j=1; j < N; j++)
         printf("%9lu %e\n", szs[j], tim[j]);
   }
   else
   {
      for (j=0; j < N; j++)
         printf(" %e\n", tim[j]);
   }
}

void *printRes(unsigned int N, long *szs, float *tim)
{
   void *ret=szs;
   unsigned long len;
   len = (szs) ? 8+(N-1)*9 : N*9;
   if (len <= 80) /* can it fit in a vt320 window? */
   {
      printRow(N, szs, tim);
      ret=NULL;
   }
   else
      printCol(N, szs, tim);
   return(ret);
}

int nextLvl(unsigned int N, float tol, float *tim)
/*
 * Compares adjacent times, t0 (smaller size), t1 (larger).  If
 * t0*tol < t1, declares it to be a cache level.  We will usually set
 * tol to something like 1.15, to indicate a less than 15% drop in
 * performance isn't likely a cache level (associativity effects for instance).
 * RETURNS: -1 for no detected level, else last in-cache index.
 */
{
   int i;
   for (i=0; i < N-1; i++)
      if (tim[i]*tol < tim[i+1])
         return(i);
   return(-1);
}

ulong *readMemWalk(char *fnam, uint *N0, uint *P0, float **RES, float **VAR)
/*
 * Reads fnam as results of a search run.  N is the number of sizes, P # cores.
 * All calls fill in the N*P array RES, and if VAR is set, VAR will be
 * allocated and filled in with the observed variance in running the smallest
 * timed case.  NOTE: VAR must be non-null for power-of-two search!
 * RETURNS: array of sizes in bytes as ulong, NULL if file can't be read.
 */
{
   unsigned long *szs;
   float *var=NULL, *ra, *f;
   int *ipiv;
   FILE *fp;
   unsigned long varsz;
   unsigned int P, N, i, j;

   fp = fopen(fnam, "r");
   if (!fp)
      return(NULL);
   assert(fscanf(fp, "%u %u", &P, &N) == 2);
   assert(N > 1 && P > 0);
   if (VAR)
   {
      var = malloc(P*sizeof(float));
      assert(var);
      N--;
   }
   ra = malloc(N*P*sizeof(float));
   szs = malloc(N*sizeof(long));
   assert(ra && szs);

   if (VAR)
      assert(fscanf(fp, " %lu", &varsz) == 1);
   for (i=0; i < N; i++)
      assert(fscanf(fp, " %lu", szs+i) == 1);
   ipiv  = getPivAndSort(N, szs);
   if (VAR)
      assert(varsz == szs[0]);

   for (f=ra, j=0; j < P; j++, f += N)
   {
      if (VAR)
         assert(fscanf(fp, " %e", var+j) == 1);
      for (i=0; i < N; i++)
         assert(fscanf(fp, " %e", f+i) == 1);
      pivFloats(N, ipiv, f);
   }
   fclose(fp);
   free(ipiv);

   if (VAR)
      *VAR = var;
   *P0 = P;
   *N0 = N;
   *RES = ra;
   return(szs);
}

ulong *doMemWalk(int flg, uint nrep, uint N0, uint NN, uint incN, uint *P0,
                 uint *Nt, float **RES, float **VAR)
/*
 * Calls xtmemwalk, doing writes if WRT is set, in range [N0,NN].
 * Sets:  *P0: # of cores;  *Nt: # of times; *RES: Nt*P timings;
 *   VAR is NULL if non-power2 search, else variance for smallest size.
 * RETURNS: all sizes timed.
 */
{
   char *fnam, *cmnd;
   unsigned long *szs=NULL;
   unsigned int dlen=2, flen, clen, tlen, i;
   const int WRT=flg&1, FRC=(flg>>1)&1;

   if (incN)
      dlen = NumDecDigits(incN) + 1;
   dlen += NumDecDigits(N0);
   dlen += NumDecDigits(NN);
   dlen += NumDecDigits(nrep);
   flen = 16;
   clen = 29;
   tlen = flen + clen + (dlen<<1) + 1;
   cmnd = malloc(tlen*sizeof(char));
   assert(cmnd);
   i = sprintf(cmnd,
       "./xtmemwalk -w %u -m %u -S %u %u %u -o res/cachetim%u_%u_%u_%ux%u",
               WRT, nrep, N0, NN, incN, WRT, N0, NN, incN, nrep);
   assert(i < tlen);
   fnam = cmnd + clen + dlen;
   printf("fnam='%s' cmnd='%s\n", fnam, cmnd);
   if (FRC)
      i = system(cmnd);
   else
   {
      i = 0;
      szs = readMemWalk(fnam, Nt, P0, RES, VAR);
      if (!szs)
         i = system(cmnd);
   }
   if (!szs)
      szs = readMemWalk(fnam, Nt, P0, RES, VAR);
   if (i || !szs)
      fprintf(stderr, "\nERROR: cmnd='%s'\n\n", cmnd);
   assert(!i);
   free(cmnd);
   return(szs);
}

ulong findL1(unsigned int nrep)
/*
 * This routine assumes that the maximum L1 data cache size is 128KB.  Since
 * bandwidth issues should not affect missing the L1 (assuming at least an L2),
 * we use reads rather than writes so that write-through caches show strong
 * performance loss when the L1 is exceeded.  In theory, L1 detection could
 * be as hard as further levels, because it could have pseudo-random
 * replacement that resists perfection.  Traditionally, L1 caches are
 * low associativity, and do a good approximiation of LRU, which means we
 * can find them on most architectures.  Even on modern intel, with high
 * associativity, the L1 cache is highly perfectable, and so this usually
 * works even when the larger caches defy precise preciction.  We declare the
 * first substantial drop in performance to be the L1; in the worst case,
 * this might be only some part of the L1 due to random replacement, but if
 * this algorithm sees a drop, so will our gemm, so its still a good enough
 * guess as long as the drop is repeatable.
 *
 * RETURNS: size of detected L1 on all cores in bytes.  If cores had different
 *          detected L1 sizes, returns 0.
 */
{
   unsigned long *szs, *sz;
   float *tim, *var, *fa, *f;
   unsigned long l1sz=0;
   int P, Nt, j;

   sz = szs = doMemWalk(0, nrep, 2<<KSH, 256<<KSH, 0, &P, &Nt, &tim, &var);
   fa = malloc(Nt * sizeof(float));
   assert(fa);
   f = tim;
   for (j=0; j < P; j++, f += Nt)
   {
      int il1;
      float tol, t0, v0;
      v0 = 1.0 + smoothRes(Nt, fa, f);
      sz = printRes(Nt, sz, fa);
      t0 = tim[j*Nt];
      tol = (t0 + var[j])/t0;
      printf(" --> %u: tol=%.4f v0=%.4f", j, tol*100.0, v0*100.0);
      tol = (v0 > tol) ? v0 : tol;
      tol = (tol > 1.1) ? tol : 1.1;
      il1 = nextLvl(Nt, tol, fa);
      assert(il1 >= 0 && il1 < Nt);
      printf(" CacheSz=%lu\n", szs[il1]);
      if (!l1sz)
         l1sz = il1;
      else if (il1 != l1sz)
      {
         free(szs);
         free(fa);
         free(tim);
         free(var);
         return(0);
      }
   }
   l1sz = szs[l1sz];
   free(szs);
   free(fa);
   free(tim);
   free(var);
   return(l1sz);
}

int FindLargestDrop(int P, int Nt, ulong *szs, float *tim, uint ldt, float *var)
{
   float *fa, *f;
   int idx=(-1), MISMATCH=0, j;
   unsigned long *sz=szs;

   fa = malloc(Nt * sizeof(float));
   assert(fa);
   f = tim;
   for (j=0; j < P; j++, f += ldt)
   {
      unsigned int i, ix;
      float tol, t0, v0, maxgap;
      v0 = 1.0 + smoothRes(Nt, fa, f);
      sz = printRes(Nt, sz, fa);
      t0 = tim[j*Nt];
      tol = (t0 + var[j])/t0;
      printf(" --> %u: tol=%.4f v0=%.4f", j, tol*100.0, v0*100.0);
      tol = (v0 > tol) ? v0 : tol;
      tol = (tol > 1.1) ? tol : 1.1;

      ix = 0;
      maxgap = fa[1] - fa[0];
      for (i=1; i < Nt-1; i++)
      {
         const float gap = fa[i+1] - fa[i];
         if (gap > maxgap)
         {
            maxgap = gap;
            ix = i;
         }
      }
      if (fa[ix]*tol >= fa[ix+1]) /* Is gap meaningful? */
         ix = -2;                 /* no real cache edge found */
      printf(" idx=%u\n", ix);
      if (idx == -1)
         idx = ix;
      else if (ix != idx)  /* see if LLC is close enough */
         MISMATCH = 1;
   }
   free(fa);
   if (MISMATCH)  /* need to see if times are close enough */
      return(-1); /* for now, freak out */
   return(idx);
}

/*
 * Finding the edge of an LRU cache is easy and precise: you should get constant
 * performance until the size is reached, at which point you see the perf
 * of the next level of the heirarchy.  If you plot a series of runs with
 * chain size on X axis & time/elt on Y, you will see a sharp step, where
 * each step up corresponds to a memory heirarchy level.
 *
 * FIFO replacement will also cleanly miss all lines once it begins missing
 * in the cache, providing you have doubled the CS between timings, since
 * we have no reuse in the chain following (LRU & FIFO same w/o reuse).
 * If the actual cache size is not a power of two (eg, NWAYS=3), then FIFO
 * will behave differently from LRU, since there was some reuse in one,
 * and possibly both, size timings.
 *
 * In caches with multiple ways and random replacement, the situation is worse.
 * You begin to miss in the cache even before you reach the cache size,
 * because some of your lines randomly conflict even prior to hitting the
 * cache edge.  Once the size is exceeded, the fact that you throw lines
 * randomly away means that some lines missing in LRU will be retained.
 * This results in a steadily declining performance once CS / nways is exceeded.
 * In this scenario, you should see a series of slopes to new maximums, where
 * new maxes come from exceeding CS/NWAYS (bringing in more interference),
 * while the slope comes from an increasing number of conflicts (i.e. in 1st
 * CS/NWAYS, no conflicts possible, then 2 conflicts poss for first set,
 * then two for 2nd set, etc., then 3 conflicts possible, etc.).
 *
 * For truly random replacement, the cache should be *perfectable*.  The idea
 * is that you take some data set that does fit in the cache, and you repeatedly
 * access it.  On the second access, you have unnecessary misses due to random
 * replacement, and some lines go into the unused ways, but others go into
 * the used ways, causing some more misses for 3rd access.  If repeated
 * enough times, you should eventually get all ways used in the cache just
 * by random chance, and since you aren't accessing any data outside it,
 * you stop doing replacement, and now the line between in- and out-of-cache
 * is sharp.  In practice, this does not seem to work, particularly around
 * the boundary (trying to distinguish between aggregate cache size of
 * 2MB+256KB vs 2MB, for instance).  At a guess, OS and other processes
 * get swapped in, polluting both the instruction & data caches (further
 * levels are shared, remember), and this is enough to prevent perfecting in
 * practice.  However, increasing the reps does tend to improve the accuracy,
 * so tmemwalk takes the -m # arg, which says do the timing # times, and take
 * the minimum walk time as the result.
 *
 * It is therefore pretty much impossible to reliably detect the exact cache
 * size w/o knowing the number of ways and replacement policy.  We therefore
 * accept innaccuracy, and time a region between twice the L1 size up to
 * 64MB / core, which we assume should be completely out-of-cache.
 * We time all powers of 2, using writes during the chain in order to
 * maximize bandwidth requirements and thus magnify the gap when the LLC
 * (last level cache) is exceeded.  We will then take the largest drop in
 * performance in absolute terms as the effective LLC.  We will then examine
 * the region between the detected LLC and the L1 to see if there is a
 * LLPC (last level private cache).
 *
 * Since size refinement doesn't work reliably, we just accept the power-of-two
 * estimate, which may be correct for pure power-of-two sizes.  On modern
 * Intel, for instance, caches are exclusive, so it'll miss some stuff.
 * I.e. on a system with 2MB LLC, and 256 private L2, the real size is
 * roughly 2MB+246KB, but we'll declare it as 2MB.
 *
 * We won't blindly use these cache numbers, but rather tune using them as
 * bound estimates, so this should be close enough for our needs.
 */
ulong findLLC(uint nrep, ulong L1sz, ulong *LLPsz)
/*
 * *LLPsz is set to the last level of private cache found, or 0 if we can't find
 * RETURNS: last-level cache discovered on all cores, or 0 if they disagree
 */
{
   unsigned long *szs;
   unsigned long LLsz;
   float *tim, *var;
   int ilast, ipriv, P, Nt;

   printf("SEARCHING FOR LAST LEVEL OF CACHE\n");
   szs = doMemWalk(1, nrep, L1sz<<1, 64<<MSH, 0, &P, &Nt, &tim, &var);
   ilast = FindLargestDrop(P, Nt, szs, tim, Nt, var);
   assert(ilast < Nt-1);
   if (ilast >= 0)
      LLsz = szs[ilast];
   else
      LLsz = 0;
   if (LLsz)
      printf("LLC ESTIMATED AS %luKB (%luMB).\n", LLsz>>KSH, LLsz>>MSH);
   else
      printf("LLC CANNOT BE FOUND.\n");
   if (LLPsz)
   {
      if (ilast > 1 && ilast < Nt-1)
         ipriv = FindLargestDrop(P, ilast+1, szs, tim, Nt, var);
      else
         ipriv = -1;

      if (ipriv >= 0 && ipriv < ilast)
         *LLPsz = szs[ipriv];
      else
         *LLPsz = 0;
      if (ipriv < ilast && ipriv >= 0)
         printf("LAST LEVEL OF PRIVATE CACHE ESTIMATED AS: %luKB.\n",
                (*LLPsz)>>KSH);
      else
         printf("LAST LEVEL OF PRIVATE CACHE CANNOT BE FOUND\n");
   }
   free(tim);
   free(var);
   free(szs);
   return(LLsz);
}

void goToTown(void)
{
   unsigned long l1sz, LLsz, LLPsz, cs;
   unsigned long szs[3];
   int i, lvl;
   FILE *fp;

   printf("FINDING THE SIZE OF YOUR L1 DATACACHE:\n");
   l1sz = findL1(16);
   if (!l1sz)
      l1sz = findL1(64);
   printf("L1 SIZE DETECTED AS %lu for all cores.\n", l1sz);
   fflush(stdout);
   assert(l1sz);
   LLsz = findLLC(16, l1sz, &LLPsz);
   if (!LLsz)
      LLsz = findLLC(64, l1sz, &LLPsz);
   fflush(stdout);
   assert(LLsz);

   fp = fopen("res/atlas_cache.h", "w");
   assert(fp);
   fprintf(fp, "/* generated by ATLAS/tune/sysinfo/cachesrch.c */\n");
   fprintf(fp, "#ifndef ATLAS_CACHE_H\n   #define ATLAS_CACHE_H 1\n\n");

   szs[0] = l1sz; szs[1] = LLPsz; szs[2] = LLsz;
   for (lvl=0,i=0; i < 3; i++)
   {
      cs = szs[i];
      if (cs)
      {
         lvl++;
         fprintf(fp, "#define L%uC_SZ %lu\n", lvl, cs);
         fprintf(fp, "#ifdef SREAL\n   #define L%uC_ELTS %lu\n", lvl, cs>>2);
         fprintf(fp, "#elif defined(DREAL) || defined(SCPLX)\n   "
                     "#define L%uC_ELTS %lu\n", lvl, cs>>3);
         fprintf(fp, "#elif defined(DCPLX)\n   #define L%uC_ELTS %lu\n",
                 lvl, cs>>4);
         fprintf(fp, "#endif\n\n");
      }
   }
   i = (szs[1] == 0) ? 1 : 2;
   fprintf(fp, "#define LLPC_LVL %u\n", i);
   fprintf(fp, "#define LLPC_SZ  L%uC_SZ\n", i);
   fprintf(fp, "#define LLPC_ELTS  L%uC_ELTS\n", i);
   fprintf(fp, "#define LLC_LVL %u\n", ++i);
   fprintf(fp, "#define LLC_SZ  L%uC_SZ\n", i);
   fprintf(fp, "#define LLC_ELTS  L%uC_ELTS\n", i);
   fprintf(fp, "\n#endif /* end multiple inclusion guard */\n");
   fclose(fp);
}

int main(int nargs, char **args)
{
   FILE *fp;
/*
 * Make sure we don't already have result we need
 */
   fp = fopen("res/atlas_cache.h", "r");
   if (fp)
   {
      fclose(fp);  /* already have answer */
      printf("Cache already mapped.  "
             "Delete res/atlas_cache.h to force research.\n");
      return(0);   /* so we are done */
   }
/*
 * If we are on an x86, try using cpuinfo to map cache
 */
   #if defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)
   if (!system("make findCache_x86"))
      return(0);
   #endif
   goToTown();
   return(0);
}
