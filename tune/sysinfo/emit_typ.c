/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 1997, 2017 R. Clint Whaley
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

int GetPower2(int n)
{
   int pwr2, i;

   if (n == 1) return(0);
   for (pwr2=0, i=1; i < n; i <<= 1, pwr2++);
   if (i != n) pwr2 = 0;
   return(pwr2);
}

#define ShiftThresh 2
char *GetDiv(int N, char *inc)
{
   static char ln[256];
   int pwr2 = GetPower2(N);
   if (N == 1) sprintf(ln, "%s", inc);
   else if (pwr2) sprintf(ln, "((%s) >> %d)", inc, pwr2);
   else sprintf(ln, "((%s) / %d)", inc, N);
   return(ln);
}

char *GetInc(int N, char *inc)
{
   static char ln0[256];
   char ln[256];
   char *p=ln;
   int i, n=N, iPLUS=0;

   if (n == 0)
   {
      ln[0] = '0';
      ln[1] = '\0';
   }
   while(n > 1)
   {
      for (i=0; n >= (1<<i); i++);
      if ( (1 << i) > n) i--;
      if (iPLUS++) *p++ = '+';
      sprintf(p, "((%s) << %d)", inc, i);
      p += strlen(p);
      n -= (1 << i);
   }
   if (n == 1)
   {
      if (iPLUS++) *p++ = '+';
      sprintf(p, "%s", inc);
   }
   if (iPLUS > ShiftThresh) sprintf(ln0, "(%d*(%s))", N, inc);
   else if (iPLUS) sprintf(ln0, "(%s)", ln);
   else sprintf(ln0, "%s", ln);
   return(ln0);
}

void PrintTypeHead(FILE *fpout)
{
   int pwr2;

   fprintf(fpout, "#ifndef ATLAS_TYPE_H\n");
   fprintf(fpout, "#define ATLAS_TYPE_H\n\n");

   if (sizeof(int) == sizeof(void*))
      fprintf(fpout, "#define ATL_iptr_t unsigned int\n");
   else if (sizeof(long) == sizeof(void*))
      fprintf(fpout, "#define ATL_iptr_t unsigned long\n");
   else if (sizeof(long long) == sizeof(void*))
      fprintf(fpout, "#define ATL_iptr_t unsigned long long\n");
   else if (sizeof(short) == sizeof(void*))
      fprintf(fpout, "#define ATL_iptr_t unsigned short\n");
   else if (sizeof(char) == sizeof(void*))
      fprintf(fpout, "#define ATL_iptr_t unsigned char\n");
   else
   {
      fprintf(stderr, "No integral type same length as pointer!\n");
      assert(0);
   }
   if (sizeof(void*) >= sizeof(double))
      fprintf(fpout, "#define ATL_cparr_t ATL_iptr_t\n"
              "#define ATL_cpflt_t double\n");
   else if (sizeof(void*) >= sizeof(float))
      fprintf(fpout, "#define ATL_cparr_t ATL_iptr_t\n"
              "#define ATL_cpflt_t float\n");
   else if (sizeof(int) >= sizeof(float))
      fprintf(fpout, "#define ATL_cparr_t unsigned int\n"
              "#define ATL_cpflt_t float\n");
   else if (sizeof(long) >= sizeof(float))
      fprintf(fpout, "#define ATL_cparr_t unsigned long\n"
              "#define ATL_cpflt_t float\n");
   else if (sizeof(long long) >= sizeof(float))
      fprintf(fpout, "#define ATL_cparr_t unsigned long long\n"
              "#define ATL_cpflt_t float\n");
   else
   {
      fprintf(stderr, "\nNo integral type longer than float!\n");
      assert(0);
   }
   fprintf(fpout, "\n#define ATL_PSIZE %d\n", (int) sizeof(void*));
   fprintf(fpout, "#define ATL_LSIZE %d\n", (int) sizeof(long));
   fprintf(fpout, "#define ATL_ISIZE %d\n", (int) sizeof(int));
   fprintf(fpout, "#define ATL_SSIZE %d\n", (int) sizeof(float));
   fprintf(fpout, "#define ATL_DSIZE %d\n", (int) sizeof(double));
   fprintf(fpout, "#define ATL_CSIZE %d\n", (int) (2*sizeof(float)));
   fprintf(fpout, "#define ATL_ZSIZE %d\n", (int) (2*sizeof(double)));
   fprintf(fpout, "#define ATL_psize ((size_t)%d)\n", (int) sizeof(void*));
   fprintf(fpout, "#define ATL_lsize ((size_t)%d)\n", (int) sizeof(long));
   fprintf(fpout, "#define ATL_isize ((size_t)%d)\n", (int) sizeof(int));
   fprintf(fpout, "#define ATL_ssize ((size_t)%d)\n", (int) sizeof(float));
   fprintf(fpout, "#define ATL_dsize ((size_t)%d)\n", (int) sizeof(double));
   fprintf(fpout, "#define ATL_csize ((size_t)%d)\n", (int) (2*sizeof(float)));
   fprintf(fpout, "#define ATL_zsize ((size_t)%d)\n", (int) (2*sizeof(double)));
   fprintf(fpout, "#define ATL_%cMulBySize(N_) %s\n", 'p',
           GetInc(sizeof(void*), "((size_t)(N_))"));
   fprintf(fpout, "#define ATL_%cMulBySize(N_) %s\n", 'l',
           GetInc(sizeof(long), "((size_t)(N_))"));
   fprintf(fpout, "#define ATL_%cMulBySize(N_) %s\n", 'i',
           GetInc(sizeof(int), "((size_t)(N_))"));
   fprintf(fpout, "#define ATL_%cMulBySize(N_) %s\n", 's',
           GetInc(sizeof(float), "((size_t)(N_))"));
   fprintf(fpout, "#define ATL_%cMulBySize(N_) %s\n", 'd',
           GetInc(sizeof(double), "((size_t)(N_))"));
   fprintf(fpout, "#define ATL_%cMulBySize(N_) %s\n", 'c',
           GetInc(2*sizeof(float), "((size_t)(N_))"));
   fprintf(fpout, "#define ATL_%cMulBySize(N_) %s\n", 'z',
           GetInc(2*sizeof(double), "((size_t)(N_))"));


   pwr2 = GetPower2(sizeof(int));
   if (pwr2)
   {
      fprintf(fpout, "#define ATL_ishift %d\n", pwr2);
      fprintf(fpout, "#define ATL_iDivBySize(N_) ((N_) >> %d)\n", pwr2);
   }
   else
      fprintf(fpout, "#define ATL_iDivBySize(N_) ((N_) / %d)\n",
              sizeof(int));
   pwr2 = GetPower2(sizeof(long));
   if (pwr2)
   {
      fprintf(fpout, "#define ATL_lshift %d\n", pwr2);
      fprintf(fpout, "#define ATL_lDivBySize(N_) ((N_) >> %d)\n", pwr2);
   }
   else
      fprintf(fpout, "#define ATL_lDivBySize(N_) ((N_) / %d)\n",
              sizeof(long));
   pwr2 = GetPower2(sizeof(void*));
   if (pwr2)
   {
      fprintf(fpout, "#define ATL_pshift %d\n", pwr2);
      fprintf(fpout, "#define ATL_pDivBySize(N_) ((N_) >> %d)\n", pwr2);
   }
   else
      fprintf(fpout, "#define ATL_pDivBySize(N_) ((N_) / %d)\n",
              sizeof(void*));
   pwr2 = GetPower2(sizeof(float));
   if (pwr2)
   {
      fprintf(fpout, "#define ATL_sshift %d\n", pwr2);
      fprintf(fpout, "#define ATL_sDivBySize(N_) ((N_) >> %d)\n", pwr2);
      fprintf(fpout, "#define ATL_cshift %d\n", pwr2+1);
      fprintf(fpout, "#define ATL_cDivBySize(N_) ((N_) >> %d)\n", pwr2+1);
   }
   else
   {
      fprintf(fpout, "#define ATL_sDivBySize(N_) ((N_) / %d)\n",
              sizeof(float));
      fprintf(fpout, "#define ATL_cDivBySize(N_) ((N_) / %d)\n",
              2*sizeof(float));
   }
   pwr2 = GetPower2(sizeof(double));
   if (pwr2)
   {
      fprintf(fpout, "#define ATL_dshift %d\n", pwr2);
      fprintf(fpout, "#define ATL_dDivBySize(N_) ((N_) >> %d)\n", pwr2);
      fprintf(fpout, "#define ATL_zshift %d\n", pwr2+1);
      fprintf(fpout, "#define ATL_zDivBySize(N_) ((N_) >> %d)\n", pwr2+1);
   }
   else
   {
      fprintf(fpout, "#define ATL_dDivBySize(N_) ((N_) / %d)\n",
              sizeof(double));
      fprintf(fpout, "#define ATL_zDivBySize(N_) ((N_) / %d)\n",
              2*sizeof(double));
   }
   fprintf(fpout, "\n#endif\n");
}
int main(int nargs, char *args[])
{
   FILE *fpout=NULL;
   if (nargs == 1) fpout = stdout;
   else if (nargs != 2)
   {
      fprintf(stderr, "usage: %s <file out>\n", args[0]);
      exit(-1);
   }
   if (fpout == NULL)
   {
      fpout = fopen(args[1], "w");
      assert(fpout != NULL);
   }
   PrintTypeHead(fpout);
   if (fpout != stdout) fclose(fpout);
   return(0);
}
