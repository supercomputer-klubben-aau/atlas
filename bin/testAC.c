#ifdef TEST_GLOBAL
   #include "atlas_gatmctr.h"
#else
   #include "atlas_atmctr.h"
   #include "atlas_threads.h"
   #define ATL_gatmctr_new(p_, cnt_, flg_) ATL_atmctr_new(cnt_)
   #define ATL_gatmctr_free(ac_) ATL_atmctr_free(ac_)
   #define ATL_gatmctr_print(fp_, ac_)
   #define ATL_gatmctr_dec(ac_, rnk_) ATL_atmctr_dec(ac_)
   #define ATL_GAC_MIX 0
#endif
#include "atlas_bitvec.h"
#define ATL_GETFLAGS 1
#include "atlas_genparse.h"

static long GN=500, NERR=0;
volatile static unsigned char *CNT;
static void *AC;
void DoWorkAC(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   const long N=GN;
   long cnt;
   while ((cnt = ATL_gatmctr_dec(AC, vrank)))
   {
      if (cnt < 1 || cnt > N)
      {
         fprintf(stderr, "   AtmCtr out of range: %lu\n", cnt);
         NERR++;
      }
      else
         CNT[N-cnt]++;
   }
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -v <verb> [1]\n");
   fprintf(stderr, "   -f <flag> [NOPOLL]\n");
   fprintf(stderr, "   -F # <flag1> ... <flagN>\n");
   fprintf(stderr, "   -p <nthreads> [ATL_NTHREAD]\n");
   fprintf(stderr, "   -P <nthr0> <nthrN> <nthrINC> [ATL_NTHREAD]\n");
   fprintf(stderr, "   -n <cnt> [8111]\n");
   fprintf(stderr, "   -N <cnt0> <cntN> <cntINC> [8111]\n");
   fprintf(stderr, "   -g <gap> : nctrs = P-gap [0]\n");
   exit(ierr ? ierr : -1);
}

int GetFlags(int nargs, char **args, int *N0, int *NN, int *NINC,
             int *P0, int *PN, int *PINC, int **FLAGs, int *GAP)
{
   int i, P, verb=1;

   *GAP = 0;
   *P0 = *PN = *PINC = ATL_NTHREADS;
   *N0 = *NN = *NINC = 8111;
   *FLAGs=NULL;
   for (i=1; i < nargs; i++)
   {
      int WH=0;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'p':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *P0 = *PN = *PINC = atoi(args[i]);
         if (*P0 == 0)
            *PN = -1;
         break;
      case 'g':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *GAP = atoi(args[i]);
         break;
      case 'n':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *N0 = *NN = *NINC = atoi(args[i]);
         break;
      case 'N':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *N0 = atoi(args[i]);
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *NN = atoi(args[i]);
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *NINC = atoi(args[i]);
         break;
      case 'P':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *P0 = atoi(args[i]);
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *PN = atoi(args[i]);
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *PINC = atoi(args[i]);
         break;
      case 'v':         /* -p <nthr> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         verb = atoi(args[i]);
         break;
      case 'F':
         *FLAGs = GF_GetIntList(nargs, args, i, 1);
         i += (*FLAGs)[0] + 1;
         break;
      case 'f':         /* -f <flag> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *FLAGs = GF_GetIntList1(atoi(args[i]));
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (*FLAGs == NULL)
      *FLAGs = GF_GetIntList1(ATL_GAC_MIX);
   return(verb);
}

int main(int nargs, char **args)
{
   void *ac;
   long i;
   int p, f, verb, N0, NN, incN, P0, PN, incP, *flgs, gap;
   unsigned long nerr=0;
   verb = GetFlags(nargs, args, &N0, &NN, &incN, &P0, &PN, &incP, &flgs, &gap);

/*
 * First just check that we can alloc/free (use valgrind to validate)
 */
   if (verb)
      fprintf(stderr, "%d of %s\n", __LINE__, __FILE__);
   ac = ATL_gatmctr_new(ATL_NTHREADS, ATL_NTHREADS*(ATL_NTHREADS+2), flgs[1]);
   ATL_gatmctr_free(ac);
   if (verb)
      fprintf(stderr, "%d of %s\n", __LINE__, __FILE__);
   ac = ATL_gatmctr_new(ATL_NTHREADS, ATL_NTHREADS, flgs[1]);
   ATL_gatmctr_free(ac);
   ac = ATL_gatmctr_new(ATL_NTHREADS, ATL_NTHREADS>>1, flgs[1]);
   if (verb)
      fprintf(stderr, "%d of %s\n", __LINE__, __FILE__);
   ATL_gatmctr_free(ac);
   if (verb)
      fprintf(stderr, "%d of %s\n", __LINE__, __FILE__);
   i = (PN > 1) ? PN : 2;
   ac = ATL_gatmctr_new(i, i-1, flgs[1]);
   ATL_gatmctr_free(ac);
   if (verb)
      fprintf(stderr, "%d of %s\n", __LINE__, __FILE__);
   ac = ATL_gatmctr_new(0, 0, flgs[1]);
   ATL_gatmctr_free(ac);
/*
 * Test basic functionality
 */
   for (f=0; f < flgs[0]; f++)
   {
      int flg = flgs[f+1];
      for (p=P0; p <= PN; p += incP)
      {
         int nctr = (p > gap) ? p-gap : 1;
         for (GN=N0; GN <= NN; GN += incN)
         {
            int err=0;
            CNT = calloc(GN, sizeof(char));
            AC = ac = ATL_gatmctr_new(nctr, GN, flg);
            fprintf(stderr, "\n\n");
            ATL_gatmctr_print(stderr, ac);
            fprintf(stderr, "TESTING P=%d, N=%d:\n", p, GN);
            ATL_goParallel(p, DoWorkAC, NULL, NULL, NULL);
            ATL_gatmctr_free(ac);
            for (i=0; i < GN; i++)
            {
               if (CNT[i] != 1)
               {
                  err++;
                  if (verb)
                     fprintf(stderr, "   CNT[%d] = %d, expected 1\n",i,CNT[i]);
               }
            }
            err += NERR;
            NERR = 0;
            fprintf(stderr, "DONE    P=%d, N=%d, NERR=%d\n", p, GN, err);
            nerr += err;
            free((void*)CNT);
         }
      }
   }
   fprintf(stderr, "\nDONE ALL TESTING, NERR=%d\n", nerr);
   free(flgs);
   return(nerr);
}
