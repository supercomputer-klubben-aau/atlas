/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_mmparse.h"
void PrintUsage(char *name, int ierr, char *flag)
{
   fprintf(stderr,
"This routine adds a demand for every kernel in input to get a\n"
"K-dimension cleanup routine.\n"
"NOTE: it is legal for both input and output file to have same name.\n");
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr,"USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -i : input sumfile (stdin); can be repeated\n");
   fprintf(stderr, "   -o : output sumfile (stdout)\n");
   exit(ierr ? ierr : -1);
}

ATL_mmnode_t *GetFlags(int nargs, char **args, char **FOUT)
{
   FILE *fpin=stdin;
   char *fout=NULL;
   ATL_mmnode_t *mb=NULL, *mp;
   int i;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'i':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        mp = ReadMMFile(args[i]);
        assert(mp);
        if (!mb)
           mb = mp;
        else
        {
           ATL_mmnode_t *p;
           p = ATL_LastMMNode(mb);
           p->next = mp;
        }
        break;
      case 'o':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        fout = args[i];
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!mb)
   {
      mb = ReadMMFile(NULL);
      assert(mb);
   }
   *FOUT = fout;
   return(mb);
}

int main(int nargs, char **args)
{
   ATL_mmnode_t *mb, *mp;
   char *fout;
   mb = GetFlags(nargs, args, &fout);

   for (mp=mb; mp; mp = mp->next)
      mp->flag |= (1<<MMF_KCLN);

   WriteMMFile(fout, mb);
   KillAllMMNodes(mb);
   return(0);
}
