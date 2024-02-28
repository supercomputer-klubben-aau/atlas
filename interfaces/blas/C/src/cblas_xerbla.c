/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias3.h"
#endif
#include "atlas_level3.h"
#include "cblas.h"

#include <stdarg.h>
void cblas_xerbla(int p, const char *rout, const char *form, ...)
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
