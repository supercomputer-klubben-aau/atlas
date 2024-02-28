/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "atlas_misc.h"
#include "atlas_lvl3.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),sysinfo.h))

#define dumb_seed(iseed_) srand(iseed_)
#ifndef RAND_MAX  /* rather dangerous non-ansi workaround */
   #define RAND_MAX ((unsigned long)(1<<30))
#endif
#define dumb_rand() ( 0.5 - ((double)rand())/((double)RAND_MAX) )

double time00(void);

#ifndef DENMAT
   #define DENMAT 200
#endif

#ifdef Right_
   char Side = 'R';
#else
   char Side = 'L';
#endif
#ifdef Upper_
   char Uplo = 'U';
#else
   char Uplo = 'L';
#endif
#ifdef Transpose_
   char Tran = 'T';
#else
   char Tran = 'N';
#endif
#ifdef UnitDiag_
   char Diag = 'U';
#else
   char Diag = 'N';
#endif


void PrintUsage(char *nam)
{
   fprintf(stderr,
           "usage: %s -m <M> -N <N0> <NN> <incN> -a <alpha> -f <filename>\n",
           nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, int *M, int *N0, int *NN, int *incN,
              TYPE *alpha,  char *file)
{
   int i;
   char *in, *out=file;

   file[0] = '\0';
   *M = ATL_opsq_66MB;
   *N0 = 100;
   *NN = 2000;
   *incN = 100;
   *alpha = 1.0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'f':
         in = args[++i];
         while (*file++ = *in++);
         break;
      case 'm':
         *M = atoi(args[++i]);
         break;
      case 'a':
         *alpha = atof(args[++i]);
      case 'N':
         *N0 = atoi(args[++i]);
         *NN = atoi(args[++i]);
         *incN = atoi(args[++i]);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
}

int main(int nargs, char *args[])
{
   char fnam[256];
   TYPE alpha;
   int M, N0, NN, incN, n, nn, k;
   FILE *fpout;

   GetFlags(nargs, args, &M, &N0, &NN, &incN, &alpha, fnam);
   if (fnam[0])
   {
      fpout = fopen(fnam, "a");
      assert(fpout);
      nn = 3*ATL_opsq_66MB;
      fprintf(fpout, "#define TRSM_%c%c%c%c_Xover %d\n",
              Side, Uplo, Tran, Diag, nn);
   }
   fprintf(stdout, "\n\nXover point at NB=%d, N=%d\n\n", ATL_opsq_66MB, nn);
   return(0);
}

