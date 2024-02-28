#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#define ATL_DEF_RUNTHR 1
#include "atlas_tprim.h"
#ifndef MINLS
   #define MINLS 16        /* smallest cache line size */
#endif
#ifndef MAXSZ
   #define MAXSZ (1<<30)   /* allow at most 2GB total allocation */
#endif
#ifndef VARCNT
   #define VARCNT 8
#endif
/*
 * large mmap calls never succeed for me, unless I allow 0 return, which
 * doesn't work, obviously.
 */
#define NOMAP 1
#ifndef NOMAP
   #include <sys/mman.h>
#endif
#define ulong unsigned long
#define KSH 10 /* kilo-byte (1024) shift */
#define MSH 20 /* mega-byte (104857) shift */
int DOWRITE=1;
volatile void *VP=NULL;
size_t NSZP=0;
ulong SZ0=0, SZN=0, SZINC=0;
float *SZTIMES=NULL;
unsigned int MINREP=16, NTHR=ATL_NTHREADS;
volatile unsigned char chkin[ATL_NTHREADS]={0};

double ATL_walltime(void);

static void my_barrier(int rnk)
{
   unsigned char newv = !chkin[rnk];
   if (rnk)
   {
      chkin[rnk] = newv;
      while (*chkin != newv);
   }
   else
   {
      unsigned int i;
      for (i=1; i < NTHR; i++)
         while (chkin[i] != newv);
      *chkin = newv;
   }
}

volatile void **findAnyNull(size_t N, volatile char *cb)
{
   size_t i;
   for (i=0; i < N; i += MINLS)
   {
      volatile void **vp = (volatile void **)(cb+i);
      if (*vp == NULL)
         return(vp);
   }
   return(NULL);
}
unsigned long countNulls(size_t N, volatile char *cb)
{
   size_t i;
   unsigned long cnt=0;
   for (i=0; i < N; i += MINLS)
   {
      volatile void **vp = (volatile void **)(cb+i);
      if (*vp == NULL)
         cnt++;
   }
   return(cnt);
}

volatile void **findNull(size_t N, volatile void **b, volatile void **pI)
/*
 * Tries to find a NULL ptr near location pI in buffer b
 */
{
   const size_t beg=(size_t)b, end = beg + N-MINLS, p=(size_t)pI;
   size_t d;

   for (d=MINLS; d < N; d += MINLS)
   {
      size_t ad = p + d;
      volatile void **vp = (volatile void**)ad;
      if (ad <= end)
      {
         if (*vp == NULL)
            return(vp);
      }
      else /* we are past end of array */
      {
         size_t ad;
         if (p <= d)  /* we are past valid memory on low end! */
            break;    /* and end of array, stop looking */
         ad = p - d;
         if (ad < beg)  /* we are past beginning & end of array */
            break;      /* so stop looking */
      }
      if (p > d)
      {
         ad = p - d;
         vp = (volatile void**) ad;
         if (ad >= beg)
            if (*vp == NULL)
               return(vp);
      }
   }
   return(NULL);
}

void initChain(size_t N, volatile char *b)
/*
 * Initializes b with a circular chain of pointers moving between random
 * locations, and zeros the character next to the pointer.  N is length of b,
 * and must evenly divide by MINLS.
 */
{
   volatile void **p0, **pn;
   volatile char *cp;
   const size_t n = N / MINLS;
   size_t i;
   unsigned int seed;

   assert(sizeof(char)+sizeof(void*) <= MINLS);
   assert(n*MINLS == N);
/*
 * NULL all chain ptrs, and zero char read/write var
 */
   cp = b + sizeof(char*);
   for (i=0; i < N; i += MINLS)
   {
      *((char**)(b+i)) = NULL;
      cp[i] = 0;
   }
/*
 * Now use random jumps to create circular chain that visits all MINLS locations
 */
   seed = ((size_t)cp) | N;
   p0 = (volatile void**) b;
   for (i=0; i < n-1; i++)
   {
      int k;
      k = rand_r(&seed) % n;
      *p0 = (volatile void*) 0x1;  /* don't let pn be this loc */
      pn = (volatile void**)(b + k*MINLS); /* srch pt for next link in chain */
      if (*pn != NULL) /* can't use this exact spot, find nearby NULL */
         pn = findNull(N, (volatile void**)b, pn); /* find next link in chain */
      assert(pn != NULL);
      assert(*pn == NULL);
      *p0 = pn;
      p0 = pn;
   }
   *p0 = (volatile void*)b;

   i = countNulls(N, b);
   if (i)
   {
      size_t idx;
      idx = (size_t) findAnyNull(N, b);
      idx -= (size_t) b;
      idx /= MINLS;
      fprintf(stderr, "FOUND %lu NULLS, first one at %lu out of %lu!\n",
              (unsigned long)i, (unsigned long)idx, (unsigned long)n);
      assert(0);
   }
}

void printChain(size_t N, char *b)
{
   size_t i;

   printf("INDEX:");
   for (i=0; i < N; i += MINLS)
      printf("%5lu,", i/MINLS);
   printf("\n");
   printf("CHAIN:");
   for (i=0; i < N; i += MINLS)
   {
      char **p = (char**) (b+i);
      printf("%5lu,", ((unsigned long)(*p - b))/MINLS);
   }
   printf("\n\n");
}
#define takeChain(N_, p_, r_) \
{  \
   size_t i; \
   void **pp=(void**)(p_); \
   for (i=0; i < N_; i += MINLS) \
      pp = *pp; \
   r_ = (void*)pp; \
}

#define addChain(N_, p_) \
{  \
   size_t i; \
   void **pp=(void**)(p_); \
   for (i=0; i < N_; i += MINLS) \
   {  char *cp=((char*)pp)+sizeof(void*);\
      *cp += *cp; \
      pp = *pp; \
   } \
}


double doTime(int rnk, int nmin, size_t N, volatile char *cp)
/*
 * NOTE: fact that we do writes may make it harder to detect write-thru L1,
 * but makes last-level & shared caches easier to see, and we care more about
 * these.  We can always do a read-only test for L1 if we want.
 */
{
   double mintime=10e15;
   size_t r,i;
   int m;
/*
 * Using walltime, so minimum result most accurate (no context switch, etc.)
 */
   for (m=0; m < nmin; m++)
   {
      double t0;
      void *vp;

      my_barrier(rnk);            /* agree we are all ready for next timing */
      if (DOWRITE)
      {
         t0 = ATL_walltime();
         addChain(N, cp);
         t0 = ATL_walltime() - t0;
      }
      else
      {
         t0 = ATL_walltime();
         takeChain(N, cp, vp);
         t0 = ATL_walltime() - t0;
         assert(vp);  /* use vp so takeChain not dead code */
      }
      mintime = (mintime <= t0) ? mintime : t0;
      assert (t0 > 0.0);
   }
   return(mintime);
}

void *doWork(void *vp)
{
   ATL_thread_t *tp = vp;
   const int rnk = tp->rank;
   float *mytim = SZTIMES + NSZP*rnk;
   volatile char *cp;
   const size_t N = ((64<<MSH) > SZN) ? (64<<MSH) : SZN;
   size_t mysz = SZN, nrep=1, cnt=0, nmin=MINREP;
   size_t i, k;
   double gapMax=0.0, last;

   cp = VP;
   cp += SZN * rnk;
   if (!SZINC)
   {
/*
 *    Loop over decreasing sizes
 */
      do
      {
/*
 *       Reassure auditor we haven't hung; flush i/o before barrier so it
 *       doesn't affect timings, also so auditor can see it real-time.
 */
         if (!rnk)
         {
            printf("   attempting to populate the ways . . .\n");
            fflush(stdout);
         }
/*
 *       In case this fits into some level of high-associativity cache, minimize
 *       random replacement conflicts by repeatedly reading cache before timing.
 *       Random replacement can only occur on a miss, and with enough reads you
 *       could therefore have no misses.  In practice, its a shared i/d cache,
 *       and other processes are running, so this can't work perfectly.
 */
         initChain(mysz, cp);
         for (k=0; k < VARCNT; k++)
            assert(!countNulls(mysz, cp));
         if (!rnk)
         {
            printf("   timing %lu KB (%lu MB) with %lu trials.\n",
                   (unsigned long)(mysz>>KSH), (unsigned long)(mysz>>MSH),
                   (unsigned long) nmin);
            fflush(stdout);
         }
         mytim[cnt++] = doTime(rnk, nmin, N, cp);
         mysz >>= 1;
         nrep <<= 1;
      }
      while (mysz >= SZ0);
/*
 *    Now, measure the variance on in-L1 by retiming the last size VARCNT times
 *    Ways should already be populated by prior calls.
 */
      mysz <<= 1;
      nrep >>= 1;
      for (k=0; k < VARCNT; k++)
      {
         double t0;
         t0 = doTime(rnk, nmin, N, cp);
         if (k)
         {
            double gap;
            gap = t0 - last;
            if (gap < 0.0)
               gap = -gap;
            if (gap > gapMax)
               gapMax = gap;
         }
         last = t0;
      }
      mytim[cnt] = gapMax;
   }
   else /* doing a search between [SZ0,SZN,SZINC] */
   {
      double t0;
      ulong i, sz=SZ0;

      for (i=0; i < NSZP; i++, sz += SZINC)
      {
         int m;
         if (i == NSZP-1)
            sz = SZN;
         nrep = (64<<MSH) / sz;
         initChain(sz, cp);
         for (k=0; k < VARCNT; k++)
            assert(!countNulls(sz, cp));
         t0 = doTime(rnk, nmin, N, cp);
         mytim[cnt++] = t0 / nrep;
      }
   }
}

ulong NumSizes(ulong sz0, ulong szN, ulong szInc)
{
   size_t cnt=0;
   if (!szInc)
   {
      do
      {
         szN >>= 1;
         cnt++;
      }
      while(szN >= sz0);
   }
   else
   {
      ulong i;
      for (i=sz0; i < szN; i += szInc)
         cnt++;
      cnt++;  /* always timing szN, even if it isn't sz0+X*szInc */
   }
   return(cnt);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr,"USAGE: %s [flags]:\n", name);
   fprintf(stderr,"   -o <resfile> : (stdout)\n");
   fprintf(stderr,"   -m <minrepeats> : (8)\n");
   fprintf(stderr,"   -w 0/1 : don't/do write during buffer traversal\n");
   fprintf(stderr,"   -p # : use only # threads (0 < # <= %u)\n", ATL_NTHREADS);
   fprintf(stderr,
      "   -S sz0 szN szINC set szINC=0 for power-of-2 search (default)\n");
   fprintf(stderr, "   -o <output file>\n");
   exit(ierr ? ierr : -1);
}

FILE *GetFlags(int nargs, char **args)
{
   unsigned long sz0=4<<KSH, szn=64<<MSH, szinc=0;
   int mrep=8;
   FILE *fpout=stdout;
   int i;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'o':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        fpout = fopen(args[i], "w");
        assert(fpout);
        break;
      case 'p':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        NTHR = atoi(args[i]);
        break;
      case 'w':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        DOWRITE = atoi(args[i]);
        break;
      case 'm':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        mrep = atoi(args[i]);
        break;
      case 'S':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        sz0 = atol(args[i]);
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        szn = atol(args[i]);
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        szinc = atol(args[i]);
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   MINREP = mrep;
   SZ0 = sz0;
   SZN = szn;
   SZINC = szinc;
   return(fpout);
}

void *GetMem
(
   ulong sz0,  /* low size per thread to try (power of 2) */
   ulong szN,  /* max size per thread.  Power of 2 in bytes */
   ulong szInc,/* 0: do power-of-two search, else increment between sizes */
   ulong pgsz, /* page size.  Core's mem must start on pg bound. */
   size_t *TSZ /* total size allocated */
)
{
   void *vp;
   size_t i, totsz;
   ulong nsz;

   nsz = NumSizes(sz0, szN, szInc);
   if (!szInc)  /* want large power of two search */
   {
      if (szN > (1<<MSH))
         printf(
            "Finding size per core I can allocate, starting at size %luMB\n",
                (unsigned long) (szN>>MSH));
      else
         printf(
            "Finding size per core I can allocate, starting at size %luKB\n",
                (unsigned long) (szN>>KSH));
      do
      {
         szN = ((szN+pgsz-1)/pgsz)*pgsz;
         totsz = pgsz + NTHR*(szN + nsz*sizeof(float));
         vp = NULL;
         if (totsz <= MAXSZ)
         {
            #ifdef NOMAP
               vp = malloc(totsz);
            #else
               vp = mmap(NULL, totsz, PROT_READ|PROT_WRITE,
                         MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_LOCKED,
                         -1, -1);
            #endif
         }
         szN >>= 1;
         nsz--;
      }
      #ifdef NOMAP
      while(!vp);
      #else
      while (vp == NULL || vp == MAP_FAILED);
      #endif
      szN <<= 1;
      nsz += 2;     /* put back 1, have 1 slot for saving variance timeing */
      NSZP = nsz;
      SZN = szN;
      if (szN > (1<<MSH))
         printf("Allocated %lu MB per core (%luMB total), NSIZES= %lu.\n",
                (unsigned long) (szN>>MSH), (unsigned long) totsz>>MSH,
                (unsigned long) nsz);
       else
         printf("Allocated %lu KB per core (%luKB total), NSIZES= %lu.\n",
                (unsigned long) (szN>>KSH), (unsigned long) totsz>>KSH,
                (unsigned long) nsz);
   }
   else /* doing simple range search */
   {    /* max size already known to work from pwr2 search */
      totsz = pgsz + NTHR*(szN + nsz*sizeof(float));
      vp = malloc(totsz);
      NSZP = nsz;
      assert(vp);
   }
   SZTIMES = vp;
   i = (size_t) vp;
   i += nsz * sizeof(float);
   i = ((i+pgsz-1)/pgsz)*pgsz;
   VP = (void*)i;
   *TSZ = totsz;
   return(vp);
}

int main(int nargs, char **args)
{
   void *vp;
   unsigned long pgsz;
   size_t totsz, nsz;
   FILE *fp;

   fp = GetFlags(nargs, args);
   pgsz = sysconf(_SC_PAGESIZE);
   if (pgsz == -1)
      pgsz = 4096;

   vp = GetMem(SZ0, SZN, SZINC, pgsz, &totsz);
   nsz = NSZP;
   ATL_runThreads(NTHR, doWork, NULL);
   if (SZINC)
   {
      ulong i, sz=SZ0;
      for (i=0; i < nsz; i++, sz += SZINC)
      {
         int j;
         sz = (i != nsz-1) ? sz : SZN;
         printf("%6lu ", (unsigned long)(sz>>KSH));
         for (j=0; j < NTHR; j++)
               printf(" %.2e", SZTIMES[i+j*nsz]/sz);
         printf("\n");
      }
   }
   else
   {
      size_t sz;
      int i;
      for (sz=SZN, i=0; i < nsz-1; sz>>=1, i++)
      {
         int j;
         printf("%6lu ", (unsigned long)(sz>>KSH));
         for (j=0; j < NTHR; j++)
               printf(" %.2e", SZTIMES[i+j*nsz]);
         printf("\n");
      }
      printf("varien:");
      for (i=0; i < NTHR; i++)
         printf(" %.2e", SZTIMES[nsz-1+i*nsz]);
   }
   printf("\n");
   if (fp != stderr && fp != stdout)
   {
      size_t sz;
      float *f = SZTIMES;

      if (!SZINC)
      {
         int j;
         fprintf(fp, "%u %lu\n", NTHR, nsz);
         fprintf(fp, "%lu\n", SZ0);
         for (sz=SZN; sz >= SZ0; sz >>= 1)
            fprintf(fp, "%lu\n", sz);
         fprintf(fp, "\n");
         for (j=0; j < NTHR; j++, f += nsz)
         {
            int i;
            fprintf(fp, "%e\n", f[nsz-1]);
            for (i=0; i < nsz-1; i++)
               fprintf(fp, "%e\n", f[i]);
            if (j != NTHR-1)
               fprintf(fp, "\n");
         }
      }
      else
      {
         int j;
         fprintf(fp, "%u %lu\n", NTHR, nsz);
         fprintf(fp, "%lu\n", SZ0);
         for (j=0; j < nsz; j++)
            fprintf(fp, "%lu\n", (j != nsz-1) ? (SZ0+j*SZINC) : SZN);
         fprintf(fp, "\n");
         for (j=0; j < NTHR; j++, f += nsz)
         {
            int i;
            for (i=0; i < nsz; i++)
               fprintf(fp, "%e\n", f[i]);
            if (j != NTHR-1)
               fprintf(fp, "\n");
         }
      }
      fclose(fp);
   }
   #ifdef NOMAP
      free(vp);
   #else
      munmap(vp, totsz);
   #endif
   return(0);
}
