#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "atlas_cpparse.h"
void PrintUsage(char *name)
{
   fprintf(stderr,
           "USAGE: %s <files> : negate mflops in standard CP files\n",
           name);
   exit(-1);
}
int main(int nargs, char **args)
{
   int i, k;
   double *mfs;
   ATL_cpnode_t *p, *pb;
   if (nargs < 2)
      PrintUsage(args[0]);
   for (i=1; i < nargs; i++)
   {
      pb = ReadCPFile(args[i]);
      assert(pb);
      for (p=pb; p; p = p->next)
      {
         mfs = p->mflop;
         for (k=0; k < 4; k++)
            mfs[k] = (mfs[k] > 0.0) ? -mfs[k] : mfs[k];
      }
      ResubGoodGccInCPNodes(pb);
      WriteCPFile(args[i], pb);
      KillAllCPNodes(pb);
   }
   exit(0);
}
