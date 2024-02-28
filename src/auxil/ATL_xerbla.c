/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
void ATL_xerbla(int p, char *rout, char *form, ...)
{
   va_list argptr;

   va_start(argptr, form);
#ifdef GCCWIN
   if (p) printf("Parameter %d to routine %s was incorrect\n", p, rout);
   vprintf(form, argptr);
#else
   if (p)
      fprintf(stderr, "Parameter %d to routine %s was incorrect\n", p, rout);
   vfprintf(stderr, form, argptr);
#endif
   va_end(argptr);
   exit(-1);
}
