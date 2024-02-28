#include "atlconf.h"

void GetPmake(int verb, enum OSTYPE OS, int ncpu, char *pmake)
{
   FILE *fp;
   if (ncpu == 0) ncpu = 1;
   if (!OSIsWin(OS))
   { /* using gnu make */
      char *res;
      res = atlsys_1L(NULL, "make DoNothing --version | fgrep GNU", verb, 1);
      if (res)
      {
         free(res);
         #if defined(ATL_GCCCLANG) || defined(ATL_GCC3P)
            sprintf(pmake, "$(MAKE) -j 1");
         #else
            if (ncpu == 1) sprintf(pmake, "$(MAKE) -j 2");
            else sprintf(pmake, "$(MAKE) -j %d", ncpu);
         #endif
      }
   }
   else
   {  /* AIX, HP-UX, SunOS make do not have parallel option */
      switch(OS)
      {
      case OSIRIX:
         #if defined(ATL_GCCCLANG) || defined(ATL_GCC3P)
            strcpy(pmake, "$(MAKE)");
         #else
            strcpy(pmake, "$(MAKE) -P");
         #endif
         break;
      case OSFreeBSD:
         #if defined(ATL_GCCCLANG) || defined(ATL_GCC3P)
            sprintf(pmake, "$(MAKE) -j 1");
         #else
            if (ncpu == 1) sprintf(pmake, "$(MAKE) -j 2");
            else sprintf(pmake, "$(MAKE) -j %d", ncpu);
         #endif
         break;
      default:
         strcpy(pmake, "$(MAKE)");
      }
   }
}

/*
 * Defined some macros for checking booleans passed by -Si
 */
#define SI_FKO 0
#define SI_BOZOL1 1
#define SI_ARCHDEF 2
#define SI_IEEE 3
#define SI_LATUNE 4
#define SI_NOF77  5
#define SI_NOCYGWIN 6
#define SI_OMP 6
#define SI_LAREF 7
#define SI_ADCLOSEP 8
#define SI_SKPTHRCHK 9
#define SI_INDTHR 10
#define SI_IS_TRUE(bv_, b_) (((bv_)>>(b_))&1)
#define SI_SET_BIT(bv_, b_) (bv_) |= 1 << (b_)
#define SI_UNSET_BIT(bv_, b_) (bv_) &= ~(1 << (b_))

void PrintUsage(char *name, int iarg, char *arg)
{
   fprintf(stderr, "\nERROR around arg %d (%s).\n", iarg,
           arg ? arg : "unknown");
   fprintf(stderr, "USAGE: %s [flags] where flags are:\n", name);
   fprintf(stderr, "   -v <verb> : verbosity level\n");
   fprintf(stderr, "   -O <enum OSTYPE #>  : set OS type\n");
   fprintf(stderr,
           "   -t <#> : set # of threads (-1: autodect; 0: no threading)\n");
   fprintf(stderr,
      "NOTE: enum #s can be found by : make xprint_enums ; ./xprint_enums\n");
   exit(iarg);
}

void GetFlags(int nargs,                /* nargs as passed into main */
              char **args,              /* args as passed into main */
              int *verb,                /* verbosity setting */
              enum OSTYPE *OS,          /* OS to assume */
              int *nthreads,           /* # of threads */
              char **targ             /* mach to ssh to*/
             )
{
   int i, k, k0, kn, DoInt;
   char *sp, *sp0;

   *verb = 0;
   *targ = NULL;

   *OS = 0;
   *verb = 0;
   *nthreads = -1;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 't':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *nthreads = atoi(args[i]);
         break;
      case 'O':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *OS = atoi(args[i]);
         break;
      case 'v':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *verb = atoi(args[i]);
         break;
      case 'T':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *targ = args[i];
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
}

int main(int nargs, char **args)
{
   int verb, OS, arch, ncpu;
   char *targ;
   char pmake[256];
   GetFlags(nargs, args, &verb, &OS, &ncpu, &targ);
   GetPmake(verb, OS, ncpu, pmake);
   printf("PMAKE='%s'\n", pmake);
   return(0);
}
