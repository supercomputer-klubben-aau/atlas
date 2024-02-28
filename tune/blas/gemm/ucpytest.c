#include "atlas_cpparse.h"
#include "atlas_cptesttime.h"

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -i <copy.idx>: (stdin) list of user-supplied copies\n");
   fprintf(stderr, "   -o <cpyWORK.CPS>: (stdout) user copy kernels\n");
   exit(ierr ? ierr : -1);
}

ATL_cpnode_t *GetFlags(int nargs, char **args, char *PRE, char **FOUT)
/*
 * RETURNS: list of all kerns/blk factors to tune copies for
 */
{
   ATL_cpnode_t *cb=NULL, *cp;
   int i;
   int INDONE=0;
   char pre = 'd';

   *FOUT = NULL;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'p':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         pre = tolower(args[i][0]);
         assert(pre == 's' || pre == 'd' || pre == 'z' || pre == 'c');
         break;
      case 'i':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        cb = ReadCPFile(args[i]);
        INDONE = 1;
        break;
      case 'o':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *FOUT = args[i];
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!INDONE)
      cb = ReadCPFile(NULL);
   if (cb)
      PrefixStrAllNodes(cb, GetOffset(&cb->next, cb), GetOffset(&cb->rout, cb),
                       "CPYCASES/");
   *PRE = pre;
   return(cb);
}

int FailAnyClaimedCase(char pre, ATL_cpnode_t *cp)
{
   const unsigned int flg = cp->flag, TA=(flg>>CPF_TRANS)&1,
      CBLK=(flg>>CPF_CBLK)&1;
   int ia;
   const int ialps[3] = {1,-1,2};
   const int ibets[4] = {1,-1,2,0};
   fprintf(stderr, "   TESTING %d-'%s':\n", cp->ID, cp->rout);
   if (!CBLK)  /* For A/B copy, NU is M or N unroll, while MU is ku */
   {
      for (ia=CPF_AL1; ia <= CPF_ALX; ia++)
      {
         int ialp;
         if (flg & (1<<ia) == 0)
            continue;
         ialp = ialps[ia-CPF_AL1];
         if (CPKernelFailsTest(pre, TA, 367, 233, ialp, 1, cp))
         {
            fprintf(stderr, "   FAILED %d: alpha=%d!\n", cp->ID, ialp);
            return(1);
         }
      }
   }
   else
   {
      for (ia=CPF_AL1; ia <= CPF_ALX; ia++)
      {
         int ib, ialp;
         if (flg & (1<<ia) == 0)
            continue;
         ialp = ialps[ia-CPF_AL1];
         for (ib=CPF_BE1; ib <= CPF_BE0; ib++)
         {
            int ibet;
            if (flg & (1<<ib) == 0)
               continue;
            ibet=ibets[ib-CPF_BE1];
            if (CPKernelFailsTest(pre, TA, 367, 233, ialp, 1, cp))
            {
               fprintf(stderr, "   FAILED %d: alpha=%d, beta=%d!\n",
                       cp->ID, ialp, ibet);
               return(1);
            }
         }
      }
   }
   fprintf(stderr, "   PASS:   %d-'%s'\n", cp->ID, cp->rout);
   return(0);
}
int main(int nargs, char **args)
{
   char *fnout;
   char pre;
   ATL_cpnode_t *cb, *cp;

   cb = GetFlags(nargs, args, &pre, &fnout);
   while (cb && FailAnyClaimedCase(pre, cb))
      cb = KillCPNode(cb);
   if (cb)
   {
      ATL_cpnode_t *prev=cb;
      cp = cb->next;
      while(cp)
      {
         ATL_cpnode_t *next = cp->next;
         if (FailAnyClaimedCase(pre, cp))
         {
            cp =  KillCPNode(cp);
            prev->next = cp;
         }
         else
            prev = cp;
         cp = next;
      }
   }
   fprintf(stderr, "\nSURVIVING COPY KERNELS, PRE='%c':\n", pre);
   if (fnout)
      WriteCPFile("stderr", cb);
   WriteCPFile(fnout, cb);
   return(0);
}
