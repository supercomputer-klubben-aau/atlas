#include "atlas_threads.h"
#include "assert.h"
#include "atlas_tbitvec.h"

size_t GNREP=0;
ATL_BV_t *GBV=NULL;
double timearr[ATL_NTHREADS];
void *GMUT=NULL;

void PrintUsage(char *exe)
{
   fprintf(stderr, "USAGE: %s [-r <reps>] [-n <nbits]\n", exe);
   exit(-1);
}

double GetMax(int P)
{
   double max=timearr[0];
   int i;
   for (i=1; i < P; i++)
      if (timearr[i] > max)
         max = timearr[i];

   return(max);
}
size_t GetFlags(int nargs, char **args, size_t *NREP)
{
   size_t nbits = 1000, nrep = 5555;
   int i;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'n':
         if (++i >= nargs)
            PrintUsage(args[0]);
         nbits = atoll(args[i]);
         break;
      case 'r':
         if (++i >= nargs)
            PrintUsage(args[0]);
         nrep = atoll(args[i]);
         break;
      #if 0
      case 'o':
         if (++i >= nargs)
            PrintUsage(args[0]);
         *outfile = args[i];
         break;
      #endif
      default:
         PrintUsage(args[0]);
      }
   }
   *NREP = nrep;
   return(nbits);
}

#ifdef SET_ATOM
int setFirst(void *vp, unsigned int skp)
{
   int was;
   int ATL_SetBitAtomic(void *bv, unsigned char pos);
   do
   {
      ATL_BV_t *bv=vp;
      skp = ATL_FindFirstUnsetBitBV(bv, skp);
      if (skp == -1)
         return(-1);
      bv += (skp>>shBV) + 1;
      was = ATL_SetBitAtomic(bv, skp&modmskBV);
   }
   while(was == 1);
   return(0);
}

int unsetFirst(void *vp, unsigned int skp)
{
   int was;
   int ATL_UnsetBitAtomic(void *bv, unsigned char pos);
   do
   {
      ATL_BV_t *bv=vp;
      skp = ATL_FindFirstSetBitBV(bv, skp);
      if (skp == -1)
         return(-1);
      bv += (skp>>shBV) + 1;
      was = ATL_UnsetBitAtomic(bv, skp&modmskBV);
   }
   while(was == 0);
   return(0);
}
#endif
void DoWork_ast(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   double t0;
   const size_t nrep=GNREP;
   ATL_BV_t *bv=GBV;
   size_t i, nbits, nbitsP=1, ibeg=0;
   unsigned int P;
   int (*fndChgBit)(void *, unsigned int);

   nbits = ATL_GetTotBitsBV(bv);
   P = (pp) ? pp->nworkers : 1;
   nbitsP = nbits / P;
   ibeg = nbitsP*rank;
   t0 = ATL_walltime();
   for (i=0; i < nrep; i++)
   {
      int b=ibeg;
      if (i&1)
      {
         fndChgBit = setFirst;
      }
      else
      {
         fndChgBit = unsetFirst;
      }
      do
      {
         b = fndChgBit(bv, b);
         if (b == -1)
            b = fndChgBit(bv, 0);
      }
      while (b != -1);
   }
   timearr[rank] = ATL_walltime() - t0;
}

void DoWork_mut(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   double t0;
   const size_t nrep=GNREP;
   ATL_BV_t *bv=GBV;
   size_t i, nbits, nbitsP=1, ibeg=0;
   unsigned int P;
   int (*changeBit)(ATL_BV_t *, unsigned int);
   int (*findBit)(ATL_BV_t *, unsigned int);

   nbits = ATL_GetTotBitsBV(bv);
   P = (pp) ? pp->nworkers : 1;
   nbitsP = nbits / P;
   ibeg = nbitsP*rank;
   t0 = ATL_walltime();
   for (i=0; i < nrep; i++)
   {
      int b=ibeg;
      if (i&1)
      {
         changeBit = ATL_UnsetBitBV;
         findBit = ATL_FindFirstSetBitBV;
      }
      else
      {
         changeBit = ATL_SetBitBV;
         findBit = ATL_FindFirstUnsetBitBV;
      }
      do
      {
         b = findBit(bv, b);
         if (b == -1)
            b=findBit(bv, 0);
         if (b != -1)
         {
            ATL_mutex_lock(GMUT);
            b = findBit(bv, b);
            if (b == -1)
               b=findBit(bv, 0);
            if (b != -1)
               changeBit(bv, b);
            ATL_mutex_unlock(GMUT);
         }
      }
      while (b != -1);
   }
   timearr[rank] = ATL_walltime() - t0;
}

void DoWork_thr(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   double t0;
   const size_t nrep=GNREP;
   ATL_BV_t *bv=GBV;
   size_t i, nbits, nbitsP=1, ibeg=0;
   unsigned int P;
   int (*fndChgBit)(void *, unsigned int);

   nbits = ATL_GetTotBitsBV(bv);
   P = (pp) ? pp->nworkers : 1;
   nbitsP = nbits / P;
   ibeg = nbitsP*rank;
   t0 = ATL_walltime();
   for (i=0; i < nrep; i++)
   {
      int b=ibeg;
      if (i&1)
      {
         fndChgBit = ATL_tSetFirstUnsetBV;
      }
      else
      {
         fndChgBit = ATL_tUnsetFirstSetBV;
      }
      do
      {
         b = fndChgBit(bv, b);
         if (b == -1)
            b = fndChgBit(bv, 0);
      }
      while (b != -1);
   }
   timearr[rank] = ATL_walltime() - t0;
}

void DoWork_ser(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   double t0;
   const size_t nrep=GNREP;
   ATL_BV_t *bv=GBV;
   size_t i, nbits, nbitsP=1, ibeg=0;
   unsigned int P;
   int (*changeBit)(ATL_BV_t *, unsigned int);
   int (*findBit)(ATL_BV_t *, unsigned int);

   nbits = ATL_GetTotBitsBV(bv);
   P = (pp) ? pp->nworkers : 1;
   nbitsP = nbits / P;
   ibeg = nbitsP*rank;
   t0 = ATL_walltime();
   for (i=0; i < nrep; i++)
   {
      int b=ibeg;
      if (i&1)
      {
         changeBit = ATL_UnsetBitBV;
         findBit = ATL_FindFirstSetBitBV;
      }
      else
      {
         changeBit = ATL_SetBitBV;
         findBit = ATL_FindFirstUnsetBitBV;
      }
      do
      {
         b = findBit(bv, b);
         if (b == -1)
            b=findBit(bv, 0);
         if (b != -1)
            changeBit(bv, b);
      }
      while (b != -1);
   }
   timearr[rank] = ATL_walltime() - t0;
}

int main(int nargs, char **args)
{
   size_t nbit;
   double tser_s, tthr_s, tmut_s, tast_s;
   double tser, tthr, tmut, tast;
   nbit = GetFlags(nargs, args, &GNREP);
   GBV = ATL_NewBV(nbit);
   assert(ATL_FindFirstSetBitBV(GBV, 0) == -1); /* bring into cache */
   GMUT = ATL_mutex_init();

   DoWork_ser(NULL, 0, 0);
   tser_s = timearr[0];
   printf("   serial  BV time = %e\n", tser_s);
   ATL_UnsetAllBitsBV(GBV);

   DoWork_thr(NULL, 0, 0);
   tthr_s = timearr[0];
   printf("   serial tBV time = %e\n", tthr_s);
   ATL_UnsetAllBitsBV(GBV);

   DoWork_mut(NULL, 0, 0);
   tmut_s = timearr[0];
   printf("   serial mut time = %e\n", tmut_s);
   ATL_UnsetAllBitsBV(GBV);

#ifndef SET_ATOM
   printf("SPEEDUP OVER MUTEX: tBV=%.2f, unsafeBV=%.2f\n\n",
          tmut_s/tthr_s, tmut_s/tser_s);
#else
   DoWork_ast(NULL, 0, 0);
   tast_s = timearr[0];
   printf("   serial ast time = %e\n", tast_s);
   printf("SPEEDUP OVER MUTEX: tBV=%.2f, ast=%.2f, unsafeBV=%.2f\n\n",
          tmut_s/tthr_s, tmut_s/tast_s, tmut_s/tser_s);
   ATL_UnsetAllBitsBV(GBV);
#endif


   ATL_goParallel(ATL_NTHREADS, DoWork_ser, NULL, NULL, NULL);
   tser = GetMax(ATL_NTHREADS);
   printf("   P=%4d  BV time = %e\n", ATL_NTHREADS, tser);
   ATL_UnsetAllBitsBV(GBV);

   ATL_goParallel(ATL_NTHREADS, DoWork_thr, NULL, NULL, NULL);
   tthr = GetMax(ATL_NTHREADS);
   printf("   P=%4d tBV time = %e\n", ATL_NTHREADS, tthr);
   ATL_UnsetAllBitsBV(GBV);

   ATL_goParallel(ATL_NTHREADS, DoWork_mut, NULL, NULL, NULL);
   tmut = GetMax(ATL_NTHREADS);
   printf("   P=%4d mut time = %e\n", ATL_NTHREADS, tmut);
#ifdef SET_ATOM
   ATL_UnsetAllBitsBV(GBV);
   ATL_goParallel(ATL_NTHREADS, DoWork_ast, NULL, NULL, NULL);
   tast = GetMax(ATL_NTHREADS);
   printf("   P=%4d ast time = %e\n", ATL_NTHREADS, tast);
   printf("SERIAL SPEEDUP: unsf=%.2f, iBV=%.2f, ast=%.2f, mut=%.2f\n",
          tser/tser_s, tthr/tthr_s, tast/tast_s, tmut/tmut_s);
   printf("SPEEDUP OVER MUTEX: tBV=%.2f, ast=%.2f, unsafeBV=%.2f\n\n",
          tmut/tthr, tmut/tast, tmut/tser);
#else
   printf("SERIAL SPEEDUP: unsf=%.2f, iBV=%.2f, mut=%.2f\n",
          tser/tser_s, tthr/tthr_s, tmut/tmut_s);
   printf("SPEEDUP OVER MUTEX: tBV=%.2f, unsafeBV=%.2f\n\n",
          tmut/tthr, tmut/tser);
#endif
   return(0);
}
