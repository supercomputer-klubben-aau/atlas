#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
void PrintPref(FILE *fp, int ipf)
{
   if (ipf & 2)
   {
      if (ipf & 64)
         fprintf(fp, " Ab to L1");
      else if (ipf & 512)
         fprintf(fp, " Ab to LLC");
      else
         fprintf(fp, " Ab to L2");
   }
   else
      fprintf(fp, " NO Ab PREF");
   if (ipf & 4)
   {
      if (ipf & 128)
         fprintf(fp, ", Bb to L1");
      else if (ipf & 1024)
         fprintf(fp, ", Bb to LLC");
      else
         fprintf(fp, ", Bb to L2");
   }
   else
      fprintf(fp, ", NO Bb PREF");
   if (ipf & 8)
   {
      if (ipf & 2048)
         fprintf(fp, ", Ak to L1");
      else if (ipf & 8192)
         fprintf(fp, ", Ak to LLC");
      else
         fprintf(fp, ", Ak to L2");
   }
   else
      fprintf(fp, ", NO Ak PREF");
   if (ipf & 16)
   {
      if (ipf & 4096)
         fprintf(fp, ", Bk to L1");
      else if (ipf & 16384)
         fprintf(fp, ", Bk to LLC");
      else
         fprintf(fp, ", Bk to L2");
   }
   else
      fprintf(fp, ", NO Bk PREF");
   if (ipf & 1)
   {
      if (ipf & 32)
         fprintf(fp, ", C to L1");
      else if (ipf &256)
         fprintf(fp, ", C to LLC");
      else
         fprintf(fp, ", C to L2");
   }
   else
      fprintf(fp, ", NO C PREF");
}
void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -b <ipf> : translate pref bitvec to string\n");
   fprintf(stderr, "   -p <Ab> <Bb> <C> <Ak> <Bk>\n");
   fprintf(stderr, "      Each arg to -p is the cache level to prefetch to.\n");
   fprintf(stderr, "      0 means do not prefetch, 3 means pfX\n");
   exit(ierr ? ierr : -1);
}

int GetFlags(int nargs, char **args, int ipfs[5])
{
   int i, k, ipf;
   if (nargs < 2)
   {
      fprintf(stderr, "\nMUST SPECIFY EITHER -p OR -b!\n");
      PrintUsage(args[0], 1, NULL);
   }
   ipf = 0;
   ipfs[0] = -1;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'b':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         ipf = atoi(args[i]);
         break;
      case 'p':
         for (k=0; k < 5; k++)
         {
            int j;
            if (++i >= nargs)
                PrintUsage(args[0], i-1, NULL);
            j = ipfs[k] = atoi(args[i]);
            if (j > 3 || j < 0)
               PrintUsage(args[0], i, "Cache level out of range");
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (ipf && ipfs[0] != -1)
   {
      printf("\nCANNOT SPECIFY BOTH -p AND -b!\n");
      PrintUsage(args[0], 1, NULL);
   }
   return(ipf);
}

int main(int nargs, char **args)
{
   int ipf, ipfs[5];

   ipf = GetFlags(nargs, args, ipfs);
   if (ipfs[0] == -1)  /* translating bitvec to string */
   {
      printf("IPF=%d:", ipf);
      PrintPref(stdout, ipf);
      printf("\n");
   }
   else  /* translate pref setting to bitvec */
   {     /* ipfs = 0:Ab 1:Bb 2:C 3:Ak 4:Bk */
      ipf = 0;
      if (ipfs[0])
      {
         ipf |= 2;
         if (ipfs[0] == 1)
            ipf |= 64;
         else if (ipfs[0] == 3)
            ipf |= 512;
      }
      if (ipfs[1])
      {
         ipf |= 4;
         if (ipfs[1] == 1)
            ipf |= 128;
         else if (ipfs[1] == 3)
            ipf |= 1024;
      }
      if (ipfs[2])
      {
         ipf |= 1;
         if (ipfs[2] == 1)
            ipf |= 32;
         else if (ipfs[2] == 3)
            ipf |= 256;
      }
      if (ipfs[3])
      {
         ipf |= 8;
         if (ipfs[3] == 1)
            ipf |= 2048;
         else if (ipfs[3] == 3)
            ipf |= 8192;
      }
      if (ipfs[4])
      {
         ipf |= 16;
         if (ipfs[4] == 1)
            ipf |= 4096;
         else if (ipfs[4] == 3)
            ipf |= 16384;
      }
      printf("Ab=%d, Bb=%d, C=%d, Ak=%d, Bk=%d: bitvec=%d\n",
             ipfs[0],ipfs[1],ipfs[2],ipfs[3],ipfs[4], ipf);
   }/* ipfs = 0:Ab 1:Bb 2:C 3:Ak 4:Bk */
   return(0);
}
