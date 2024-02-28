/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "atlas_type.h"
#include "atlas_fopen.h"
#include "atlas_sys.h"
int GetL1CacheSize(int MaxL1Size)
{
   FILE *L1f;
   char ln[80];
   int L1Size;

   if (!FileExists("res/L1CacheSize"))
   {
      sprintf(ln, "make RunL1 MaxL1=%d\n",MaxL1Size);
      assert(system(ln) == 0);
   }
   L1f = fopen("res/L1CacheSize", "r");
   assert(L1f != NULL);
   assert(fscanf(L1f, "%d", &L1Size) == 1);
   fclose(L1f);
   fprintf(stderr, "\n      Read in L1 Cache size as = %dKB.\n",L1Size);
   return(L1Size);
}

size_t GetPageSize(void)
{
   size_t pgsz=4096;
   char *sp;
   sp = atlsys_1L(NULL, "getconf PAGESIZE", 0, 1);
   if (!sp)
      sp = atlsys_1L(NULL, "getconf PAGE_SIZE", 0, 1);
   if (sp)
   {
      unsigned long sz;
      if (sscanf(sp, "%lu", &sz) == 1)
         pgsz = sz;
      free(sp);
   }
   return(pgsz);
}
void getfpinfo0(char pre, int *muladd, int *lat, int *lbnreg, double *mf)
{
   char fnam[128];
   FILE *fp;

   if (pre == 'z') pre = 'd';
   else if (pre == 'c') pre = 's';

   sprintf(fnam, "res/%cMULADD", pre);
   if (!FileExists(fnam))
   {
      sprintf(fnam, "make res/%cMULADD pre=%c\n", pre, pre);
      assert(system(fnam) == 0);
      sprintf(fnam, "res/%cMULADD", pre);
   }
   fp = fopen(fnam, "r");
   assert(fp);
   assert(fscanf(fp, " %d", muladd) == 1);
   assert(fscanf(fp, " %d", lat) == 1);
   assert(fscanf(fp, " %lf", mf) == 1);
   assert(fscanf(fp, " %d", lbnreg) == 1);
   fclose(fp);
}

void getfpinfo(char pre, int *muladd, int *lat, int *lbnreg, int *nkflop)
{
   double mf;
   char ln[32];

   getfpinfo0(pre, muladd, lat, lbnreg, &mf);
   if (mf <= 0.0)
   {
      if (pre == 'c') pre = 's';
      else if (pre == 'z') pre = 'd';
      sprintf(ln, "make RunMADef pre=%c\n", pre);
      assert(system(ln) == 0);
      getfpinfo0(pre, muladd, lat, lbnreg, &mf);
      assert(mf >= 0.0);
   }
   mf = mf * 0.75 * 750;
   #ifdef ATL_ARCH_ATHLON  /* kludge for athlon defaults */
      mf *= 3.0;
   #endif
   if (mf > (double)(1<<30)) *nkflop = ~(1<<31);
   else *nkflop = mf;
}

void CreateHeader(char pre, char *fnam, int L1Size, int muladd, int lat,
                  int lbnreg, int nkflop, int mmnreg, size_t pgsz)
{
   FILE *fpout;

   fpout = fopen(fnam, "w");
   assert(fpout);
   fprintf(fpout, "#ifndef ATL_%cSYSINFO_H\n   #define ATL_%cSYSINFO_H\n\n",
           toupper(pre), toupper(pre));
   if (muladd) fprintf(fpout, "#define ATL_MULADD\n");
   else fprintf(fpout, "#define ATL_NOMULADD\n");
   fprintf(fpout, "#define ATL_fplat  %d\n", lat);
   fprintf(fpout, "#define ATL_lbnreg %d\n", lbnreg);
   fprintf(fpout, "#define ATL_mmnreg %d\n", mmnreg);
   fprintf(fpout, "#define ATL_nkflop %d\n", nkflop);
   fprintf(fpout, "#define ATL_pgsz %lu\n", pgsz);
   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
}

int getmmnreg(char pre)
{
   char fnam[128];
   int mmnregs;
   FILE *fp;

   if (pre == 'z') pre = 'd';
   else if (pre == 'c') pre = 's';
   sprintf(fnam, "res/%cnreg", pre);
   if (!FileExists(fnam))
   {
      sprintf(fnam, "make res/%cnreg\n", pre);
      assert(system(fnam) == 0);
      sprintf(fnam, "res/%cnreg", pre);
   }
   fp = fopen(fnam, "r");
   assert(fp);
   assert(fscanf(fp, " %d", &mmnregs) == 1);
   fclose(fp);
   return(mmnregs);
}

int main(int nargs, char **args)
{
   size_t pgsz;
   int MaxL1Size;
   int muladd, lat, lbnreg, L1Size, mmnreg, nkflop;
   FILE *fpout;
   char pre;

   if (nargs != 3)
   {
      fprintf(stderr, "USAGE: %s <pre> <file>\n", args[0]);
      exit(-1);
   }
   pre = *args[1];

   getfpinfo(pre, &muladd, &lat, &lbnreg, &nkflop);
   pgsz = GetPageSize();
   CreateHeader(pre, args[2], L1Size, muladd, lat, lbnreg, nkflop, 0, pgsz);
   mmnreg = getmmnreg(pre);
   CreateHeader(pre, args[2], L1Size, muladd, lat, lbnreg, nkflop, mmnreg,pgsz);
   return(0);
}
