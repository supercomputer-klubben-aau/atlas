#include "atlas_tvec.h"

void PrintUsage(char *name, char *arg, int i)
{
   fprintf(stderr,
"This routine takes tvecs and spits them out in platicus prefab data format\n");
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' ON %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s <flags> ; flags include:\n", name);
   fprintf(stderr, "   -i <file> : (stdin) file with vecs to reduce\n");
   fprintf(stderr, "   -o <file>  : (stdout) platicus prefab data file out\n");
   fprintf(stderr, "   -C # <nam1> ... <nam#>: vectors (in order) to use\n");
   fprintf(stderr, "   -h 1/0: do (don't) print headers\n");
   exit (i ? i : -1);
}

char **GetFlags         /* RETURNS: array of names to put into output file */
(
   int nargs,
   char **args,
   int *Nv,             /* # of vecs to put in output file */
   int *head,           /* print headers or not */
   FILE **fpin,         /* input stream */
   FILE **fpout         /* output stream */
)
{
   char **vs=NULL, *sp;
   int i, j, n=0;
   FILE *fp;

   *fpin = stdin;
   *fpout = stdout;
   *head = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], "no '-' preceding flag!", i);
      switch(args[i][1])
      {
      case 'i':    /* -i <file> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -i ", i-1);
         *fpin = fopen(args[i], "r");
         assert(*fpin);
         break;
      case 'o':    /* -o <file> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -i ", i-1);
         fp = fopen(args[i], "w");
         if (!fp)
         {
            fprintf(stderr, "unable to open file %s!\n", args[i]);
            assert(fp);
         }
         *fpout = fp;
         break;
      case 'h':    /* -h 0/1 */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -h ", i-1);
         *head = atoi(args[i]);
         break;
      case 'C':    /* -C # <nam1> ... <nam#> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -C ", i-1);
         n = atoi(args[i]);
         assert(n > 0 && n < 2048);
         vs = malloc(sizeof(char*)*n);
         assert(vs);
         for (j=0; j < n; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -C ", i-1);
            vs[j] = args[i];
         }
         break;
      default :
         PrintUsage(args[0], args[i], i);
      }                                         /* end switch over flags */
   }                                            /* end for over flags */
   if (!vs)
   {
      n = 2;
      vs = malloc(2*sizeof(char*));
      assert(vs);
      vs[0] = "MFLOP";
      vs[1] = "N";
   }
   *Nv = n;
   return(vs);
}

int main(int nargs, char **args)
{
   FILE *fpin, *fpout;
   char **vnams, *cmnt;
   int N, i, j, RNGINC=0, head;
   ATL_tvec_t *tp, *np, *nb=NULL;

   vnams = GetFlags(nargs, args, &N, &head, &fpin, &fpout);

/*
 * Grab only the vectors to be output (in the order the user has specified)
 * from list, and free all unused vectors
 */

   tp = ATL_ReadTvecFile(fpin, &cmnt, &i);
   if (fpin != stdin)
      fclose(fpin);
   free(cmnt);
   nb = ATL_PullNamedVecsFromList(N, vnams, &tp);
   free(vnams);
   if (tp)
      ATL_KillAllTvecs(tp);
/*
 * Print them out, and we are done
 */
   ATL_PrintTvecsInRow(fpout, nb, "\t", head);
   ATL_KillAllTvecs(nb);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
   return(0);
}
