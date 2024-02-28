#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(__MINGW32__) || defined(__MINGW64__)

int slashdrivesub(char *ln)
/*
 * replaces \\c\ with c:\, returns change in string length
 * this version required for older cygwins
 */
{
   char *sp, *lp=ln, ctmp;
   int nrep=0;
   do
   {
      sp = strstr(lp, "\\\\");
      if (sp && strlen(sp) > 3)
      {
         if (sp[2] == 'a' || sp[2] == 'b' || sp[2] == 'c' || sp[2] == 'd' ||
             sp[2] == 'e' || sp[2] == 'f' || sp[2] == 'g' || sp[2] == 'h')
         {
            if (sp[3] == '\\')
            {
               ctmp = sp[2];
               sp[0] = sp[2];
               sp[1] = ':';
               sp[2] = '\\';
               for (lp=sp+3; *lp = lp[1]; lp++);
               lp = sp + 3;
               nrep++;
            }
            else lp = sp + 2;
         }
         else lp = sp + 2;
      }
      else lp = sp + 2;
   }
   while (sp);
   return(-nrep);
}

int cygdrivesub(char *ln)
/*
 * replaces \cygdrive\c\ with c:\, returns change in string length
 * this version works cygnus version 1.1.0
 */
{
   char *sp;
   int i=0;

   while(sp = strstr(ln, "\\cygdrive\\"))
   {
      i++;
      sp[0] = sp[10];
      sp[1] = ':';
      sp[2] = '\\';
      sp += 3;
      while (*sp = sp[9]) sp++;
   }
   return( slashdrivesub(ln) - (i*9) );
}

void slashsub(char *ln)
/*
 * changes forward slash of unix to backslash of windoze
 */
{
   int i;
   for (i=0; ln[i]; i++) if (ln[i] == '/') ln[i] = '\\';
}

#endif
int sComputeRound(void)
/*
 * Blind translation of netlib LAPACK LAMCH's rounding computation
 * RETURNS: 1 if numbers are correctly rounded, 0 if they are truncated
 */
{
   volatile float a, b, c, f;
   int rnd=0;
   b = a = 1.0;
   do
   {
      a *= 2.0;
      c = a + b;
      c = c - a;
   }
   while (c == 1.0);
   b = 0.5*FLT_RADIX;
   c = .01*(-FLT_RADIX);
   f = b + c;
   c = f + a;
   rnd = (c == a);

   c = 0.01*(FLT_RADIX);
   f = b + c;
   c = f + a;
   if (rnd && c == a)
      rnd = 0;
   return(rnd);
}

float sComputeSafmin(void)
/*
 * BFI translation of netlib LAPACK LAMCH's safmin calc, adapted to use float.h
 * RETURNS: LAMCH's safmin
 */
{
   volatile float small, sfmin, eps;
   const float ONE = 1.0;
/*
 * NOTE:
 *   TINY(X) : smallest positive number storable in type -> [FLT,DBL]_MIN
 *   HUGE(X) : largest positive number storable in type  -> [FLT,DBL]_MAX
 */
      eps = 0.5;
      eps *= FLT_EPSILON;
      sfmin = FLT_MIN;
      small = ONE / FLT_MAX;
      if (small >= sfmin)
         sfmin = small*(ONE+FLT_EPSILON);
   return(sfmin);
}

void emit_slamch(char *path)
{
   FILE *fpout;
   char *name;
   volatile float f, under, over;
   int len = 16, bad;

   if (path)
   {
      len += strlen(path);
      name = malloc(len);
      assert(name);
      sprintf(name, "%s/atlas_slamch.h", path);
      fpout = fopen(name, "w");
      assert(fpout);
      free(name);
   }
   else
      fpout = stdout;
   fprintf(fpout, "/* generated by %s */\n\n", __FILE__);
   fprintf(fpout, "#ifndef ATLAS_SLAMCH_H\n");
   fprintf(fpout, "   #define ATLAS_SLAMCH_H\n\n");

   fprintf(fpout, "#define ATL_slaMANTDIG     %d\n", FLT_MANT_DIG);
   fprintf(fpout, "#define ATL_slaMINEXP      %d\n", FLT_MIN_EXP);
   fprintf(fpout, "#define ATL_slaMAXEXP      %d\n", FLT_MAX_EXP);
   fprintf(fpout, "#define ATL_slaBASE        %d\n", FLT_RADIX);
   f = 0.5;
   f *= FLT_EPSILON;
   fprintf(fpout, "#define ATL_slaEPSILON     %60.53e\n", f);
   f = 0.5 * FLT_RADIX;
   f *= FLT_EPSILON;
   fprintf(fpout, "#define ATL_slaPRECISION   %60.53e\n", f);
   fprintf(fpout, "#define ATL_slaUNDERTHRESH %60.53e\n", FLT_MIN);
   under = FLT_MIN;
   f = pow(FLT_RADIX, FLT_MAX_EXP-2)*(4.0-2.0*FLT_EPSILON);
   fprintf(fpout, "#define ATL_slaOVERTHRESH  %60.53e\n", f);
   over = f;
   fprintf(fpout, "#define ATL_slaSAFMIN      %60.53e\n",
           sComputeSafmin());
   fprintf(fpout, "#define ATL_slaROUND       %d\n", sComputeRound());
/*
 * Blind translation of LAPACK's _LABAD test
 */
   f = 2000.0;
   bad = (log10(over) > f);
   fprintf(fpout, "#define ATL_slaBAD %d\n", bad);
   if (bad)
   {
      fprintf(fpout, "#define ATL_slabadUNDERTHRESH %60.53e\n",
              sqrt(under));
      fprintf(fpout, "#define ATL_slabadOVERTHRESH  %60.53e\n",
              sqrt(over));
   }
   else
   {
      fprintf(fpout,
              "#define ATL_slabadUNDERTHRESH ATL_slaUNDERTHRESH\n");
      fprintf(fpout,
              "#define ATL_slabadOVERTHRESH ATL_slaOVERTHRESH\n");
   }

   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
}

int dComputeRound(void)
/*
 * Blind translation of netlib LAPACK LAMCH's rounding computation
 * RETURNS: 1 if numbers are correctly rounded, 0 if they are truncated
 */
{
   volatile double a, b, c, f;
   int rnd=0;
   b = a = 1.0;
   do
   {
      a *= 2.0;
      c = a + b;
      c = c - a;
   }
   while (c == 1.0);
   b = 0.5*FLT_RADIX;
   c = .01*(-FLT_RADIX);
   f = b + c;
   c = f + a;
   rnd = (c == a);

   c = 0.01*(FLT_RADIX);
   f = b + c;
   c = f + a;
   if (rnd && c == a)
      rnd = 0;
   return(rnd);
}

double dComputeSafmin(void)
/*
 * BFI translation of netlib LAPACK LAMCH's safmin calc, adapted to use float.h
 * RETURNS: LAMCH's safmin
 */
{
   volatile double small, sfmin, eps;
   const double ONE = 1.0;
/*
 * NOTE:
 *   TINY(X) : smallest positive number storable in type -> [FLT,DBL]_MIN
 *   HUGE(X) : largest positive number storable in type  -> [FLT,DBL]_MAX
 */
      eps = 0.5;
      eps *= DBL_EPSILON;
      sfmin = DBL_MIN;
      small = ONE / DBL_MAX;
      if (small >= sfmin)
         sfmin = small*(ONE+DBL_EPSILON);
   return(sfmin);
}

void emit_dlamch(char *path)
{
   FILE *fpout;
   char *name;
   volatile double f, under, over;
   int len = 16, bad;

   if (path)
   {
      len += strlen(path);
      name = malloc(len);
      assert(name);
      sprintf(name, "%s/atlas_dlamch.h", path);
      fpout = fopen(name, "w");
      assert(fpout);
      free(name);
   }
   else
      fpout = stdout;
   fprintf(fpout, "/* generated by %s */\n\n", __FILE__);
   fprintf(fpout, "#ifndef ATLAS_DLAMCH_H\n");
   fprintf(fpout, "   #define ATLAS_DLAMCH_H\n\n");

   fprintf(fpout, "#define ATL_dlaMANTDIG     %d\n", DBL_MANT_DIG);
   fprintf(fpout, "#define ATL_dlaMINEXP      %d\n", DBL_MIN_EXP);
   fprintf(fpout, "#define ATL_dlaMAXEXP      %d\n", DBL_MAX_EXP);
   fprintf(fpout, "#define ATL_dlaBASE        %d\n", FLT_RADIX);
   f = 0.5;
   f *= DBL_EPSILON;
   fprintf(fpout, "#define ATL_dlaEPSILON     %60.53e\n", f);
   f = 0.5 * FLT_RADIX;
   f *= DBL_EPSILON;
   fprintf(fpout, "#define ATL_dlaPRECISION   %60.53e\n", f);
   fprintf(fpout, "#define ATL_dlaUNDERTHRESH %60.53e\n", DBL_MIN);
   under = DBL_MIN;
   f = pow(FLT_RADIX, DBL_MAX_EXP-2)*(4.0-2.0*DBL_EPSILON);
   fprintf(fpout, "#define ATL_dlaOVERTHRESH  %60.53e\n", f);
   over = f;
   fprintf(fpout, "#define ATL_dlaSAFMIN      %60.53e\n",
           dComputeSafmin());
   fprintf(fpout, "#define ATL_dlaROUND       %d\n", dComputeRound());
/*
 * Blind translation of LAPACK's _LABAD test
 */
   f = 2000.0;
   bad = (log10(over) > f);
   fprintf(fpout, "#define ATL_dlaBAD %d\n", bad);
   if (bad)
   {
      fprintf(fpout, "#define ATL_dlabadUNDERTHRESH %60.53e\n",
              sqrt(under));
      fprintf(fpout, "#define ATL_dlabadOVERTHRESH  %60.53e\n",
              sqrt(over));
   }
   else
   {
      fprintf(fpout,
              "#define ATL_dlabadUNDERTHRESH ATL_dlaUNDERTHRESH\n");
      fprintf(fpout,
              "#define ATL_dlabadOVERTHRESH ATL_dlaOVERTHRESH\n");
   }

   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
}

int main (int nargs, char **args)
{
   char *path = "res/";
   if (nargs > 1)
      path = args[1];
   #if defined(__MINGW32__) || defined(__MINGW64__)
   {
      char *winpath;
      winpath = malloc(sizeof(char)*(strlen(path)+1));
      assert(winpath);
      strcpy(winpath, path);
      slashsub(winpath);
      cygdrivesub(winpath);
      emit_dlamch(winpath);
      emit_slamch(winpath);
      free(winpath);
   }
   #else
      emit_dlamch(path);
      emit_slamch(path);
   #endif
   return(0);
}
