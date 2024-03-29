/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2006 R. Clint Whaley
 */
#include "atlconf.h"

int RunISAProbe(char *isaxnam, int verb, char *targ, char *opt)
{
   char *cmnd, *res, *frm;
   int i=1;
   if (targ)
   {
      frm = "make IRun_%s atlrun=atlas_runX targ=%s MYFLAGS='%s' | grep -F SUCCESS";
      i += strlen(targ);
   }
   else
      frm = "make IRun_%s MYFLAGS='%s' | grep -F SUCCESS";
   i += strlen(frm) + strlen(isaxnam) + strlen(opt);
   cmnd = malloc(sizeof(char)*i);
   assert(cmnd);
   if (targ)
      sprintf(cmnd, frm, isaxnam, targ, opt);
   else
      sprintf(cmnd, frm, isaxnam, opt);
   res = atlsys_1L(targ, cmnd, verb, 0);
   if(res)
   {
      if (strstr(res, "SUCCESS"))
      {
         if (verb)
            fprintf(stdout, "   %s: DETECTED!\n", isaxnam);
         free(res);
         return(1);
      }
   }
   if (verb > 1)
      fprintf(stdout, "   cmnd='%s' out='%s'\n", cmnd, res);
   free(cmnd);
   free(res);
   if (verb)
      fprintf(stdout, "   %s: NO.\n", isaxnam);
   return(0);
}

int GetAllISAExt(int verb, char *targ, enum OSTYPE OS, enum ASMDIA asmb)
{
   int i, iret=0;
   char ln[256];

   if (verb)
      fprintf(stdout, "\nProbing for supported ISA extensions:\n");

/*
 * For OS X, throw try throwing their random-ass annoyance flag
 */
   if (OS == OSOSX)
   {
      if (RunISAProbe(ISAXNAM[ISA_AV], verb, targ, ln))
         iret |= (1<<ISA_AV);
   }
   sprintf(ln, "-DATL_OS_%s -DATL_%s", osnam[OS], ASMNAM[asmb]);
   for (i=1; i < NISA; i++)
   {
      if (i >= ISA_NEON && i <= ISA_VFP3D16MAC)
      {
         if (asmb != gas_arm)
         {
            printf("   %s: SKIPPED (bad assembly)\n", ISAXNAM[i]);
            continue;
         }
         #if defined(ATL_NONIEEE) && ATL_NONIEEE != 0
            if (i >= ISA_NEON && i <= ISANEON16)
            {
               printf("   %s: SKIPPED (non_IEEE)\n", ISAXNAM[i]);
               continue;
            }
         #endif
      }
      else if (i >= ISA_AVXZ && i <= ISA_3DNow)
      {
         if (asmb != gas_x86_32 && asmb != gas_x86_64)
         {
            printf("   %s: SKIPPED (bad assembly)\n", ISAXNAM[i]);
            continue;
         }
      }
      if (RunISAProbe(ISAXNAM[i], verb, targ, ln))
         iret |= (1<<i);
   }
   return(iret);
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
   fprintf(stderr, "   -s <enum ASMDIA #>  : set assembly dialect\n");
   fprintf(stderr,
      "NOTE: enum #s can be found by : make xprint_enums ; ./xprint_enums\n");
   exit(iarg);
}

void GetFlags(int nargs,                /* nargs as passed into main */
              char **args,              /* args as passed into main */
              int *verb,                /* verbosity setting */
              enum OSTYPE *OS,          /* OS to assume */
              enum ASMDIA *asmb,        /* assembly dialect to assume */
              char **targ             /* mach to ssh to*/
             )
{
   int i, k, k0, kn, DoInt;
   char *sp, *sp0;

   *verb = 0;
   *targ = NULL;

   *asmb = 0;
   *OS = 0;
   *verb = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 's':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *asmb = atoi(args[i]);
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
   int verb, iret;
   enum OSTYPE OS;
   enum ASMDIA asmb;
   char *targ;

   GetFlags(nargs, args, &verb, &OS, &asmb, &targ);
   iret = GetAllISAExt(verb, targ, OS, asmb);
   printf("VECFLAG=%d\n", iret);
   return(0);
}
