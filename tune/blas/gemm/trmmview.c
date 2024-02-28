/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#define ATL_WANT ILCM 1
#include "atlas_iopt.h"
#include "atlas_cache.h"
#include "atlas_genparse.h"
#include "atlas_mmtesttime.h"

#define BSR_RIGHT  0
#define BSR_UPPER  1
#define BSR_TRANSA 2
void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set type/precision prefix (d) \n");
   fprintf(stderr, "   -S [L/R] : search Left or Right TRMM\n");
   fprintf(stderr, "   -U [U/L] : search Upper or Lower TRMM\n");
   exit(ierr ? ierr : -1);
}

char GetFlags(int nargs, char **args, int *FLG)
{
   char pre = 'd';
   int flg=0;
   int i;

   for (i=1; i < nargs; i++)
   {
      int wch, *ip, **ipp, TST=0;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'S':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        wch = args[i][0];
        if (wch == 'R' || wch == 'r')
           flg |= 1<<BSR_RIGHT;
        else
           flg &= ~(1<<BSR_RIGHT);
        break;
      case 'U':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        wch = args[i][0];
        if (wch == 'U' || wch == 'u')
           flg |= 1<<BSR_UPPER;
        else
           flg &= ~(1<<BSR_UPPER);
        break;
      case 'A':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        wch = args[i][0];
        if (wch == 'T' || wch == 't')
           flg |= 1<<BSR_TRANSA;
        else
           flg &= ~(1<<BSR_TRANSA);
        break;
      case 'p':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        pre = tolower(args[i][0]);
        assert(pre == 's' || pre == 'd' || pre == 'z' || pre == 'c');
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   *FLG = flg;
   return(pre);
}

int main(int nargs, char **args)
{
   ATL_mmnode_t *mb, *tk, *mp, *tb, *tp;
   char *onm;
   int flag, mmflg=0;
   char pre;

   pre = GetFlags(nargs, args, &flag);
   if (flag&(1<<BSR_RIGHT))      /* Right case requires  N=K */
   {
      mmflg = 4;
      mb = ReadMMFileWithPath(pre, "res", "ipnek.sum");
      tk = ReadMMFileWithPath(pre, "res", "trmmKRL.sum");
      tk->flag |= (1<<MMF_RIGHT);
      if (flag&(1<<BSR_UPPER))
      {
         tk->flag |= 1<<MMF_UPPER;
         onm = "trmmRU.sum";
      }
      else
         onm = "trmmRL.sum";
   }
   else
   {
      mb = ReadMMFileWithPath(pre, "res", "ipmek.sum");
      tk = ReadMMFileWithPath(pre, "res", "trmmKLL.sum");
      if (flag&(1<<BSR_UPPER))
      {
         tk->flag |= 1<<MMF_UPPER;
         onm = "trmmLU.sum";
      }
      else
         onm = "trmmLL.sum";
   }
   tk->blask = ATL_KTRMM;
   tb = TimeMMFileWithPath(pre, "res", onm, 0, 0, 0, 1, 0, -1);
   if (tb)
   {
      WriteMMFileWithPath(pre,"res", onm, tb);
      KillAllMMNodes(tp);
      return(0);
   }
/*
 * Time all trmm at block factors forced by gemm to make views
 */
   printf("TRMMK TO BE TIMED FOR SIDE=%c:\n", (flag&(1<<BSR_RIGHT)) ? 'R':'L');
   PrintMMLine(stdout, tk);
   printf("\nTIMING TRMMK FOR ALL GEMM BLOCK FACTORS:\n");
   for (mp=mb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      double mf;
      p = CloneMMNode(tk);
      p->mbB = mp->mbB;
      p->nbB = mp->nbB;
      p->kbB = mp->kbB;
      mf = TimeMMKernel(0, mmflg, p, pre, p->mbB, p->nbB, p->kbB, 1, 0,-1);
      printf("   (%u,%u,%u) = %.2f\n", p->mbB, p->nbB, p->kbB, mf);
      p->mflop[0] = mf;
      if (tb)
         tp->next = p;
      else
         tb = p;
      tp = p;
   }
   KillAllMMNodes(mb);
   KillMMNode(tk);
   WriteMMFileWithPath(pre, "res", onm, tb);
   KillAllMMNodes(tb);
   return(0);
}
